# 02_make_crossfit_summary_figure.R
#
# Build a compact figure and summary table data comparing the original
# CP + Cox method against the cross-fitted variant.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
})

get_script_dir <- function() {
  file_arg <- commandArgs(trailingOnly = FALSE)
  key <- "--file="
  hit <- grep(key, file_arg, value = TRUE)
  if (length(hit) == 0) return(getwd())
  normalizePath(dirname(sub(key, "", hit[1])), winslash = "/", mustWork = TRUE)
}

is_absolute_path <- function(path) {
  grepl("^(/|[A-Za-z]:[/\\\\])", path)
}

parse_cli <- function(args) {
  cfg <- list(
    results_dir = "playground/crossfit-support-fit-comparison/results/main-grid-crossfit-support-fit-v1",
    figure_name = "crossfit_vs_original_by_tau.pdf"
  )
  for (a in args) {
    if (!startsWith(a, "--") || !grepl("=", a, fixed = TRUE)) next
    kv <- strsplit(sub("^--", "", a), "=", fixed = TRUE)[[1]]
    key <- kv[1]
    val <- kv[2]
    if (key %in% names(cfg)) cfg[[key]] <- val
  }
  cfg
}

main <- function() {
  cfg <- parse_cli(commandArgs(trailingOnly = TRUE))
  script_dir <- get_script_dir()
  project_root <- normalizePath(file.path(script_dir, "..", ".."), winslash = "/", mustWork = TRUE)

  results_dir <- if (is_absolute_path(cfg$results_dir)) cfg$results_dir else file.path(project_root, cfg$results_dir)
  results_dir <- normalizePath(results_dir, winslash = "/", mustWork = TRUE)

  raw_dir <- file.path(results_dir, "raw")
  derived_dir <- file.path(results_dir, "derived")
  fig_dir <- file.path(results_dir, "figures")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

  rep_df <- read_csv(file.path(raw_dir, "rep_level_results.csv"), show_col_types = FALSE)

  tau_summary <- rep_df %>%
    group_by(method_label, tau_quantile) %>%
    summarise(
      coverage_mean = mean(coverage, na.rm = TRUE),
      length_mean = mean(mean_len, na.rm = TRUE),
      .groups = "drop"
    )

  overall_summary <- rep_df %>%
    group_by(method_label) %>%
    summarise(
      tau_quantile = NA_real_,
      coverage_mean = mean(coverage, na.rm = TRUE),
      length_mean = mean(mean_len, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(regime = "Overall")

  tau_summary_out <- tau_summary %>%
    mutate(regime = paste0("tau_q=", formatC(tau_quantile, digits = 2, format = "f")))

  summary_out <- bind_rows(tau_summary_out, overall_summary) %>%
    relocate(regime, .before = method_label)

  write_csv(summary_out, file.path(derived_dir, "tau_overall_summary.csv"))

  plot_df <- tau_summary %>%
    transmute(
      method_label,
      tau_quantile,
      Coverage = coverage_mean,
      `Mean Interval Length` = length_mean
    ) %>%
    pivot_longer(
      cols = c(Coverage, `Mean Interval Length`),
      names_to = "metric",
      values_to = "value"
    ) %>%
    mutate(
      method_label = factor(method_label, levels = c("CP + Cox", "CP + Cox (Cross-Fit)")),
      metric = factor(metric, levels = c("Coverage", "Mean Interval Length"))
    )

  p <- ggplot(plot_df, aes(x = tau_quantile, y = value, color = method_label, group = method_label)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2.2) +
    facet_wrap(~metric, ncol = 1, scales = "free_y") +
    scale_x_continuous(breaks = sort(unique(tau_summary$tau_quantile))) +
    scale_color_manual(values = c("CP + Cox" = "#386cb0", "CP + Cox (Cross-Fit)" = "#d95f02")) +
    labs(
      x = "tau quantile used in the DGP",
      y = NULL,
      color = "Method"
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "bottom",
      strip.background = element_rect(fill = "gray95")
    )

  p <- p +
    geom_hline(
      data = data.frame(metric = "Coverage", yint = 0.9),
      aes(yintercept = yint),
      linetype = "dashed",
      color = "gray35",
      linewidth = 0.5
    )

  ggsave(file.path(fig_dir, cfg$figure_name), p, width = 7.2, height = 7.8)
}

if (sys.nframe() == 0L) {
  main()
}
