# 03_render_model_comparison.R
#
# Render the playground working-model comparison report.

suppressPackageStartupMessages({
  library(stringr)
})

parse_cli <- function(args) {
  cfg <- list(
    results_dir = "playground/model-comparison/results/pilot",
    output_file = "model-comparison-analysis.html"
  )
  for (a in args) {
    if (!startsWith(a, "--") || !str_detect(a, "=")) next
    kv <- str_split_fixed(sub("^--", "", a), "=", 2)
    key <- kv[1]
    val <- kv[2]
    if (key %in% names(cfg)) cfg[[key]] <- val
  }
  cfg
}

is_absolute_path <- function(path) {
  grepl("^(/|[A-Za-z]:[/\\\\])", path)
}

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  cfg <- parse_cli(args)

  script_dir <- normalizePath(dirname(sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])),
                              winslash = "/", mustWork = TRUE)
  project_root <- normalizePath(file.path(script_dir, "..", ".."), winslash = "/", mustWork = TRUE)

  qmd_path <- file.path(script_dir, "02_analysis_model_comparison.qmd")
  results_dir <- if (is_absolute_path(cfg$results_dir)) cfg$results_dir else file.path(project_root, cfg$results_dir)
  results_dir <- normalizePath(results_dir, winslash = "/", mustWork = TRUE)
  output_dir <- file.path(results_dir, "report")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  cmd <- sprintf(
    "quarto render %s --to html --output %s --output-dir %s -P results_dir:%s",
    shQuote(qmd_path),
    shQuote(cfg$output_file),
    shQuote(output_dir),
    shQuote(results_dir)
  )

  message("Rendering report from: ", results_dir)
  status <- system(cmd, intern = FALSE, ignore.stdout = FALSE, ignore.stderr = FALSE)
  if (!identical(status, 0L)) stop("Quarto render failed with status ", status)

  final_path <- file.path(output_dir, cfg$output_file)
  dep_dir_name <- paste0(tools::file_path_sans_ext(basename(qmd_path)), "_files")
  dep_dir_in_report <- file.path(output_dir, dep_dir_name)

  candidate_paths <- c(
    final_path,
    file.path(results_dir, cfg$output_file),
    file.path(script_dir, cfg$output_file),
    file.path(project_root, cfg$output_file),
    file.path(project_root, "playground", "model-comparison", "results", cfg$output_file)
  )
  hit <- candidate_paths[file.exists(candidate_paths)]
  if (length(hit) == 0L) {
    warning("Render succeeded but report path was not found in expected locations.")
  } else {
    src <- normalizePath(hit[1], winslash = "/", mustWork = TRUE)
    if (!identical(src, normalizePath(final_path, winslash = "/", mustWork = FALSE))) {
      ok_report <- file.copy(src, final_path, overwrite = TRUE)
      if (!ok_report) warning("Failed to copy report into the report/ directory.")
    }

    root_html <- file.path(results_dir, cfg$output_file)
    ok_html <- suppressWarnings(file.copy(final_path, root_html, overwrite = TRUE))

    if (dir.exists(dep_dir_in_report)) {
      root_dep <- file.path(results_dir, dep_dir_name)
      if (dir.exists(root_dep)) unlink(root_dep, recursive = TRUE, force = TRUE)
      ok_dep <- suppressWarnings(file.copy(dep_dir_in_report, results_dir, overwrite = TRUE, recursive = TRUE))
    } else {
      ok_dep <- TRUE
    }

    if (!isTRUE(ok_html) || any(!ok_dep)) {
      message("Root-level mirror may be incomplete; the report copy under report/ is still available.")
    }
    message("Report saved at: ", normalizePath(root_html, winslash = "/", mustWork = TRUE))
  }
}

if (sys.nframe() == 0L) {
  main()
}
