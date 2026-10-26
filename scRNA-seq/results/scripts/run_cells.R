#!/usr/bin/env Rscript
# run_cells.R -- execute the tutorial's R cells in ONE persistent session, notebook-style.
# Behaviour: a failing cell is recorded and execution continues with the next cell
# (same as pressing Shift+Enter through a notebook). Nothing is auto-corrected.
#
# usage: Rscript run_cells.R <cells_file> <out_dir> <label>

args   <- commandArgs(trailingOnly = TRUE)
cells_file <- args[[1]]
out_dir    <- args[[2]]
label      <- if (length(args) >= 3) args[[3]] else "run"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
log_con <- file(file.path(out_dir, paste0(label, "_cell_log.txt")), open = "wt")
say <- function(...) { txt <- paste0(...); cat(txt, "\n", sep = "")
  writeLines(txt, log_con); flush(log_con) }

src <- readLines(cells_file, warn = FALSE)
starts <- grep("^#@@CELL [0-9]+$", src)
if (!length(starts)) stop("no '#@@CELL <id>' sentinels found")
ids <- as.integer(sub("^#@@CELL ", "", src[starts]))
ends <- c(starts[-1] - 1L, length(src))
cells <- lapply(seq_along(starts), function(k)
  list(id = ids[k], code = src[(starts[k] + 1L):ends[k]]))

say("### driver: run_cells.R | label=", label)
say("### R: ", R.version.string, " | Seurat: ", as.character(packageVersion("Seurat")),
    " | cells: ", length(cells))
say("### wd: ", getwd())
say("### start: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))

probe <- function() {
  out <- character(0)
  add <- function(s) out <<- c(out, s)
  for (nm in c("sc_data", "seurat_obj", "merged_seurat", "ctrl_seurat", "stim_seurat",
               "seurat_obj.markers", "naive_cd4_markers", "memory_vs_naive_markers")) {
    if (!exists(nm, envir = .GlobalEnv)) next
    obj <- get(nm, envir = .GlobalEnv)
    if (inherits(obj, "Seurat")) {
      add(sprintf("%s: %d genes x %d cells", nm, nrow(obj), ncol(obj)))
      add(sprintf("%s: reductions=[%s]", nm, paste(Reductions(obj), collapse = ",")))
      if ("seurat_clusters" %in% colnames(obj@meta.data))
        add(sprintf("%s: clusters=%d", nm,
                    length(unique(obj$seurat_clusters))))
      if (length(Idents(obj)) && !all(is.na(Idents(obj))))
        add(sprintf("%s: idents=%d [%s]", nm, nlevels(Idents(obj)),
                    paste(levels(Idents(obj)), collapse = " | ")))
    } else if (inherits(obj, c("dgCMatrix", "dgTMatrix", "matrix"))) {
      add(sprintf("%s: %d x %d matrix", nm, nrow(obj), ncol(obj)))
    } else if (is.data.frame(obj)) {
      add(sprintf("%s: data.frame %d rows x %d cols", nm, nrow(obj), ncol(obj)))
    } else if (is.character(obj)) {
      add(sprintf("%s: character(%d) = %s", nm, length(obj),
                  paste(head(obj, 4), collapse = ",")))
    } else {
      add(sprintf("%s: <%s>", nm, class(obj)[1]))
    }
  }
  if (length(dev.list())) add(sprintf("open_devices=%s (cur=%d)",
                                      paste(dev.list(), collapse = ","), dev.cur()))
  out
}

summary_rows <- list()
for (k in seq_along(cells)) {
  id <- cells[[k]]$id; code <- cells[[k]]$code
  say(""); say("===== CELL ", id, " =====")
  for (ln in code) say("  | ", ln)
  warns <- character(0); msgs <- character(0)
  t0 <- Sys.time()
  stdout_txt <- capture.output({
    res <- tryCatch(
      withCallingHandlers(
        eval(parse(text = code), envir = .GlobalEnv),
        warning = function(w) { warns <<- c(warns, conditionMessage(w))
                                invokeRestart("muffleWarning") },
        message = function(m) { msgs <<- c(msgs, conditionMessage(m))
                                invokeRestart("muffleMessage") }
      ),
      error = function(e) e
    )
  }, type = "output")
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  ok <- !inherits(res, "error")
  if (length(stdout_txt)) { say("  --stdout--"); for (l in stdout_txt) say("  > ", l) }
  if (length(msgs))  { say("  --messages(", length(msgs), ")--")
    for (m in unique(msgs)) say("  * ", m) }
  if (length(warns)) { say("  --warnings(", length(warns), ")--")
    for (w in unique(warns)) say("  ! ", w) }
  if (!ok) say("  --ERROR-- ", conditionMessage(res))
  st <- probe(); if (length(st)) { say("  --state--"); for (s in st) say("  # ", s) }
  say(sprintf("  --status-- %s in %.2fs", if (ok) "OK" else "ERROR", elapsed))
  summary_rows[[length(summary_rows) + 1]] <- data.frame(
    cell = id, status = if (ok) "OK" else "ERROR", seconds = round(elapsed, 2),
    n_warn = length(warns), n_msg = length(msgs),
    n_unique_warn = length(unique(warns)),
    error = if (ok) "" else substr(gsub("[\r\n]+", " ", conditionMessage(res)), 1, 300),
    warn_first = if (length(warns)) substr(gsub("[\r\n]+", " ", unique(warns)[1]), 1, 300) else "",
    state = paste(st, collapse = " ;; "),
    stringsAsFactors = FALSE)
}

summ <- do.call(rbind, summary_rows)
write.table(summ, file.path(out_dir, paste0(label, "_cell_summary.tsv")),
            sep = "\t", quote = FALSE, row.names = FALSE)
say(""); say("===== SUMMARY =====")
for (i in seq_len(nrow(summ))) say(sprintf("cell %-3s %-5s %7.2fs warn=%d msg=%d",
                                          summ$cell[i], summ$status[i], summ$seconds[i],
                                          summ$n_unique_warn[i], summ$n_msg[i]))
ei <- file(file.path(out_dir, paste0(label, "_sessionInfo.txt")), open = "wt")
writeLines(capture.output(sessionInfo()), ei); close(ei)
say("### end: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
close(log_con)
