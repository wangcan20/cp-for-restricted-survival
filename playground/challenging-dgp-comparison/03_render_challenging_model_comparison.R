args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript 03_render_challenging_model_comparison.R <result_dir>")
}

result_dir <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]
script_dir <- normalizePath(dirname(sub("^--file=", "", script_arg)),
                            winslash = "/", mustWork = TRUE)
qmd <- file.path(script_dir, "02_analysis_challenging_model_comparison.qmd")
tag <- basename(result_dir)
html_dir <- file.path(script_dir, "reports", tag)
dir.create(html_dir, recursive = TRUE, showWarnings = FALSE)

Sys.setenv(RESULT_DIR = result_dir)

cmd <- c(
  "render",
  qmd,
  "--to", "html",
  "--output", "challenging-model-comparison-analysis.html",
  "--output-dir", html_dir
)

status <- system2("quarto", cmd)
if (!identical(status, 0L)) {
  stop("quarto render failed with exit status ", status)
}

src_html <- file.path(html_dir, "challenging-model-comparison-analysis.html")
dst_html <- file.path(result_dir, "challenging-model-comparison-analysis.html")
file.copy(src_html, dst_html, overwrite = TRUE)
