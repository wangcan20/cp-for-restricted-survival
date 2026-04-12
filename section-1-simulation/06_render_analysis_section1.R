# 06_render_analysis_section1.R
#
# Render the Section 1 analysis report from saved simulation outputs.
#
# Example:
# Rscript section-1-simulation/06_render_analysis_section1.R \
#   --results_dir=/Users/pql/Desktop/thesis/section-1-simulation/results/pilot

suppressPackageStartupMessages({
  library(stringr)
})

parse_cli <- function(args) {
  cfg <- list(
    results_dir = "section-1-simulation/results/pilot",
    output_file = "section1-analysis.html"
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
  project_root <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = TRUE)

  qmd_path <- file.path(script_dir, "05_analysis_section1.qmd")
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
    file.path(project_root, cfg$output_file)
  )
  hit <- candidate_paths[file.exists(candidate_paths)]
  if (length(hit) == 0L) {
    warning("Render succeeded but report path was not found in expected locations.")
  } else {
    src <- normalizePath(hit[1], winslash = "/", mustWork = TRUE)
    if (!identical(src, normalizePath(final_path, winslash = "/", mustWork = FALSE))) {
      ok <- file.copy(src, final_path, overwrite = TRUE)
      if (!ok) warning("Failed to copy report to expected report directory.")
    }
    if (file.exists(final_path)) {
      message("Report saved at: ", normalizePath(final_path, winslash = "/", mustWork = TRUE))
    } else {
      message("Report saved at: ", src)
    }

    # Convenience mirror: keep a readable copy at results_dir root too.
    root_html <- file.path(results_dir, cfg$output_file)
    if (!identical(normalizePath(root_html, winslash = "/", mustWork = FALSE),
                   normalizePath(final_path, winslash = "/", mustWork = FALSE))) {
      ok_html <- file.copy(final_path, root_html, overwrite = TRUE)
      if (dir.exists(dep_dir_in_report)) {
        root_dep <- file.path(results_dir, dep_dir_name)
        if (!dir.exists(root_dep)) dir.create(root_dep, recursive = TRUE, showWarnings = FALSE)
        ok_dep <- file.copy(list.files(dep_dir_in_report, full.names = TRUE), root_dep,
                            overwrite = TRUE, recursive = TRUE)
      } else {
        ok_dep <- TRUE
      }
      if (!all(ok_html, ok_dep)) warning("Some root-level mirror files were not copied successfully.")
    }
  }
}

if (sys.nframe() == 0L) {
  main()
}
