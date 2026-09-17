#!/usr/bin/env Rscript
# ============================================================================
# HCR Probe Designer v41 — Electron sidecar bootstrap runner.
#
# This script is executed by the embedded/bundled R runtime that the Electron
# main process spawns. It:
#   1. Resolves the directory that holds HCR_probe_design_v41.R.
#   2. Points R's package library at the bundled `library/` tree (if present).
#   3. Sets HCR_REFERENCE_DIR to the (copied) reference-data folder.
#   4. Sources the app script (defines `ui` and `server` — no side effects).
#   5. Picks a free TCP port and writes `HCR_PORT=<n>` to stdout so Electron
#      can attach when the Shiny server is listening.
#
# Usage (from Electron):
#   R --no-save --no-restore --no-environ --no-site-file --no-init-file \
#     --slave -f runner.R --args <app_dir> [ <reference_dir> ]
# ============================================================================

args <- commandArgs(trailingOnly = TRUE)
app_dir <- if (length(args) >= 1 && nzchar(args[1])) args[1] else "."
if (!dir.exists(app_dir)) app_dir <- getwd()

# ---------------------------------------------------------------------------
# 1. Library path: prefer a bundled library next to the app payload.
#    Layout: <resources>/r/{R.framework|R-x.y}/...  and <resources>/r/library
# ---------------------------------------------------------------------------
bundled_lib <- file.path(dirname(app_dir), "r", "library")
if (dir.exists(bundled_lib)) {
  .libPaths(c(bundled_lib, .libPaths()))
}
# Report which library will actually serve `shiny` — good for diagnostics.
cat("HCR_LIBS=", paste(.libPaths(), collapse = ";"), "\n", sep = "")

# ---------------------------------------------------------------------------
# 2. Reference data directory (HCR_REFERENCE_DIR env var wins over arg).
# ---------------------------------------------------------------------------
ref_dir <- if (length(args) >= 2 && nzchar(args[2])) args[2] else ""
if (nzchar(ref_dir)) {
  Sys.setenv(HCR_REFERENCE_DIR = ref_dir)
}
cat("HCR_REFERENCE_DIR=", Sys.getenv("HCR_REFERENCE_DIR", unset = ""), "\n", sep = "")

# ---------------------------------------------------------------------------
# 3. Source the application script. It must define ui/server in the global env.
# ---------------------------------------------------------------------------
app_file <- file.path(app_dir, "HCR_probe_design_v41.R")
if (!file.exists(app_file)) {
  stop("HCR app script not found: ", app_file,
       "\nExpected HCR_probe_design_v41.R in: ", app_dir)
}
if (getRversion() < "4.0") {
  stop("HCR Probe Designer requires R >= 4.0 (found ", getRversion(), ")")
}
source(app_file, local = FALSE, chdir = TRUE)

if (!exists("ui", envir = globalenv()) || !exists("server", envir = globalenv())) {
  stop("App script did not define ui/server objects as expected")
}

# ---------------------------------------------------------------------------
# 4. Choose a free port and launch Shiny headlessly.
# ---------------------------------------------------------------------------
if (!requireNamespace("shiny", quietly = TRUE)) {
  stop("The 'shiny' package is not installed in the bundled library ",
       "(looked in ", paste(.libPaths(), collapse = ", "), "). ",
       "Re-run the bundling step: `npm run bundle:mac` / `npm run bundle:win`.")
}

port <- NULL
for (i in seq_len(25)) {
  p <- try(httpuv::randomPort(30000L, 60000L), silent = TRUE)
  if (inherits(p, "try-error") || is.null(p) || is.na(p) || length(p) != 1L) next
  port <- as.integer(p)
  break
}
if (is.null(port)) {
  stop("Could not find a free TCP port to host the Shiny server")
}

shiny_app <- shiny::shinyApp(ui = ui, server = server)

# Handshake line that Electron watches for
cat("HCR_PORT=", port, "\n", sep = "")
flush(stdout())

shiny::runApp(
  shiny_app,
  host = "127.0.0.1",
  port = port,
  launch.browser = FALSE
)