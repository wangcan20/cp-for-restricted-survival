# 04_run_experiments_section1.R
#
# Batch runner for Section 1 conformal PI simulations.
# It saves traceable outputs (scenario manifest, per-rep metrics,
# sampled test-level intervals, and aggregated summaries) under:
#   section-1-simulation/results/<tag>/

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(stringr)
  library(purrr)
})

get_script_dir <- function() {
  file_arg <- commandArgs(trailingOnly = FALSE)
  key <- "--file="
  hit <- grep(key, file_arg, value = TRUE)
  if (length(hit) == 0) {
    return(getwd())
  }
  normalizePath(dirname(sub(key, "", hit[1])), winslash = "/", mustWork = TRUE)
}

is_absolute_path <- function(path) {
  grepl("^(/|[A-Za-z]:[/\\\\])", path)
}

parse_cli <- function(args) {
  cfg <- list(
    tag = "pilot",
    out_dir = NULL,
    n_rep = 8L,
    n = 800L,
    test_frac = 0.30,
    p = 3L,
    high_dim_p = 8L,
    alpha = 0.10,
    B = 500L,
    seed0 = 20260212L,
    sample_n = 200L
  )

  for (a in args) {
    if (!startsWith(a, "--") || !str_detect(a, "=")) next
    kv <- str_split_fixed(sub("^--", "", a), "=", 2)
    key <- kv[1]
    val <- kv[2]

    if (key %in% c("tag", "out_dir")) cfg[[key]] <- val
    if (key %in% c("n_rep", "n", "p", "high_dim_p", "B", "seed0", "sample_n")) cfg[[key]] <- as.integer(val)
    if (key %in% c("test_frac", "alpha")) cfg[[key]] <- as.numeric(val)
  }

  cfg
}

build_scenarios <- function(p_base = 3L, p_high_dim = 8L) {
  dgp_grid <- tibble(
    dgp_id = c("base", "strong_signal", "high_shape", "high_dim"),
    dgp_label = c(
      "Base Weibull PH",
      "Stronger covariate effect",
      "Higher Weibull shape",
      "Higher covariate dimension"
    ),
    gamma = c(1.2, 1.2, 1.8, 1.2),
    kappa = c(1.0, 1.0, 1.0, 1.0),
    beta_scale = c(0.30, 0.60, 0.30, 0.30),
    p = c(as.integer(p_base), as.integer(p_base), as.integer(p_base), as.integer(p_high_dim))
  )

  censor_grid <- tibble(cens_target = c(0.10, 0.20, 0.35))
  tau_grid <- tibble(tau_quantile = c(0.75, 0.90, 0.97))

  crossing(dgp_grid, censor_grid, tau_grid) %>%
    mutate(setting_id = sprintf("s%03d", row_number())) %>%
    relocate(setting_id)
}

summarise_with_ci <- function(rep_df, alpha) {
  target <- 1 - alpha

  rep_df %>%
    group_by(
      setting_id, dgp_id, dgp_label, p, beta_scale, gamma, kappa,
      cens_target, tau_quantile, n, test_frac, alpha, B
    ) %>%
    summarise(
      n_rep = n(),
      tau_mean = mean(tau, na.rm = TRUE),
      tau_sd = sd(tau, na.rm = TRUE),
      censor_rate_emp_mean = mean(censor_rate_emp, na.rm = TRUE),
      admin_mass_emp_mean = mean(admin_mass_emp, na.rm = TRUE),
      coverage_mean = mean(coverage, na.rm = TRUE),
      coverage_sd = sd(coverage, na.rm = TRUE),
      length_mean = mean(mean_len, na.rm = TRUE),
      length_sd = sd(mean_len, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      tcrit = if_else(n_rep > 1, qt(0.975, df = n_rep - 1), NA_real_),
      coverage_se = coverage_sd / sqrt(n_rep),
      length_se = length_sd / sqrt(n_rep),
      coverage_ci_low = coverage_mean - tcrit * coverage_se,
      coverage_ci_high = coverage_mean + tcrit * coverage_se,
      length_ci_low = length_mean - tcrit * length_se,
      length_ci_high = length_mean + tcrit * length_se,
      target_coverage = target,
      abs_cov_error = abs(coverage_mean - target)
    ) %>%
    arrange(dgp_id, cens_target, tau_quantile)
}

run_grid <- function(scenarios, cfg) {
  n_total <- nrow(scenarios) * cfg$n_rep
  rep_rows <- vector("list", n_total)
  sample_rows <- list()
  row_ptr <- 1L

  start_time <- Sys.time()

  for (i in seq_len(nrow(scenarios))) {
    s <- scenarios[i, ]
    beta <- rep(s$beta_scale, s$p)

    message(sprintf(
      "[%s] Running %s (%d/%d): dgp=%s, cens=%.2f, tau_q=%.2f",
      format(Sys.time(), "%H:%M:%S"),
      s$setting_id, i, nrow(scenarios),
      s$dgp_id, s$cens_target, s$tau_quantile
    ))

    for (r in seq_len(cfg$n_rep)) {
      seed <- cfg$seed0 + i * 10000L + r

      out <- run_one_rep(
        n = cfg$n,
        test_frac = cfg$test_frac,
        p = s$p,
        beta = beta,
        gamma = s$gamma,
        kappa = s$kappa,
        tau_quantile = s$tau_quantile,
        cens_target = s$cens_target,
        alpha = cfg$alpha,
        B = cfg$B,
        seed = seed
      )

      rep_rows[[row_ptr]] <- tibble(
        setting_id = s$setting_id,
        rep = as.integer(r),
        seed = as.integer(seed),
        dgp_id = s$dgp_id,
        dgp_label = s$dgp_label,
        p = s$p,
        beta_scale = s$beta_scale,
        gamma = s$gamma,
        kappa = s$kappa,
        cens_target = s$cens_target,
        tau_quantile = s$tau_quantile,
        n = cfg$n,
        test_frac = cfg$test_frac,
        alpha = cfg$alpha,
        B = cfg$B,
        tau = out$tau,
        censor_rate_emp = out$censor_rate_emp,
        admin_mass_emp = out$admin_mass_emp,
        coverage = out$coverage,
        mean_len = out$mean_len,
        sd_len = out$sd_len
      )

      if (r == 1L) {
        n_take <- min(cfg$sample_n, nrow(out$res))
        sample_rows[[length(sample_rows) + 1L]] <- out$res %>%
          slice_sample(n = n_take) %>%
          mutate(
            setting_id = s$setting_id,
            rep = as.integer(r),
            seed = as.integer(seed),
            dgp_id = s$dgp_id,
            dgp_label = s$dgp_label,
            beta_scale = s$beta_scale,
            gamma = s$gamma,
            kappa = s$kappa,
            cens_target = s$cens_target,
            tau_quantile = s$tau_quantile
          ) %>%
          select(
            setting_id, rep, seed, dgp_id, dgp_label, beta_scale, gamma, kappa,
            cens_target, tau_quantile, tau, y, delta, T_true, T_star,
            pi_lower, pi_upper, covered, length, everything()
          )
      }

      row_ptr <- row_ptr + 1L
    }
  }

  run_minutes <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))

  list(
    rep_df = bind_rows(rep_rows),
    sample_df = bind_rows(sample_rows),
    run_minutes = run_minutes
  )
}

write_outputs <- function(scenarios, results, cfg, out_dir) {
  raw_dir <- file.path(out_dir, "raw")
  derived_dir <- file.path(out_dir, "derived")
  dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(derived_dir, recursive = TRUE, showWarnings = FALSE)

  setting_summary <- summarise_with_ci(results$rep_df, alpha = cfg$alpha)

  overall_summary <- results$rep_df %>%
    summarise(
      total_settings = n_distinct(setting_id),
      total_reps = n(),
      target_coverage = 1 - first(alpha),
      coverage_mean = mean(coverage, na.rm = TRUE),
      coverage_sd = sd(coverage, na.rm = TRUE),
      length_mean = mean(mean_len, na.rm = TRUE),
      length_sd = sd(mean_len, na.rm = TRUE)
    )

  run_config <- tibble(
    tag = cfg$tag,
    out_dir = out_dir,
    method = "CP + Cox (Upper=tau)",
    n_rep = cfg$n_rep,
    n = cfg$n,
    test_frac = cfg$test_frac,
    p = cfg$p,
    high_dim_p = cfg$high_dim_p,
    alpha = cfg$alpha,
    B = cfg$B,
    seed0 = cfg$seed0,
    sample_n = cfg$sample_n,
    run_minutes = results$run_minutes,
    created_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  )

  write_csv(scenarios, file.path(raw_dir, "scenario_manifest.csv"))
  write_csv(results$rep_df, file.path(raw_dir, "rep_level_results.csv"))
  write_csv(results$sample_df, file.path(raw_dir, "test_level_samples_rep1.csv"))
  write_csv(run_config, file.path(raw_dir, "run_config.csv"))

  write_csv(setting_summary, file.path(derived_dir, "setting_summary.csv"))
  write_csv(overall_summary, file.path(derived_dir, "overall_summary.csv"))

  writeLines(capture.output(sessionInfo()), con = file.path(raw_dir, "sessionInfo.txt"))
}

main <- function() {
  script_dir <- get_script_dir()
  project_root <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = TRUE)

  source(file.path(script_dir, "00_utils.R"))
  source(file.path(script_dir, "01_dgp_weibull_ph.R"))
  source(file.path(script_dir, "02_conformal_cox_tau.R"))
  source(file.path(script_dir, "03_simulation_section1.R"))

  cfg <- parse_cli(commandArgs(trailingOnly = TRUE))
  scenarios <- build_scenarios(p_base = cfg$p, p_high_dim = cfg$high_dim_p)

  if (is.null(cfg$out_dir) || identical(cfg$out_dir, "")) {
    out_dir <- file.path(project_root, "section-1-simulation", "results", cfg$tag)
  } else if (is_absolute_path(cfg$out_dir)) {
    out_dir <- cfg$out_dir
  } else {
    out_dir <- file.path(project_root, cfg$out_dir)
  }

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  message("Output directory: ", out_dir)
  message("Total settings: ", nrow(scenarios), "; reps per setting: ", cfg$n_rep)

  results <- run_grid(scenarios, cfg)
  write_outputs(scenarios, results, cfg, out_dir)

  message("Saved outputs to: ", out_dir)
  message(sprintf("Elapsed time: %.2f minutes", results$run_minutes))
}

if (sys.nframe() == 0L) {
  main()
}
