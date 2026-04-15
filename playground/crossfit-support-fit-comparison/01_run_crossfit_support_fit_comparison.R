# 01_run_crossfit_support_fit_comparison.R
#
# Compare the final main CP + Cox upper-tau method against an exploratory
# cross-fitted variant that uses K training folds to decouple support estimation
# from working-model fitting during calibration.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(stringr)
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
    tag = "main-grid-crossfit-support-fit-v1",
    out_dir = NULL,
    n_rep = 30L,
    n = 2000L,
    test_frac = 0.30,
    p = 3L,
    high_dim_p = 9L,
    alpha = 0.10,
    B = 1500L,
    seed0 = 20260412L,
    sample_n = 220L,
    n_folds = 5L
  )

  for (a in args) {
    if (!startsWith(a, "--") || !str_detect(a, "=")) next
    kv <- str_split_fixed(sub("^--", "", a), "=", 2)
    key <- kv[1]
    val <- kv[2]

    if (key %in% c("tag", "out_dir")) cfg[[key]] <- val
    if (key %in% c("n_rep", "n", "p", "high_dim_p", "B", "seed0", "sample_n", "n_folds")) {
      cfg[[key]] <- as.integer(val)
    }
    if (key %in% c("test_frac", "alpha")) {
      cfg[[key]] <- as.numeric(val)
    }
  }

  cfg
}

build_scenarios <- function(p_base = 3L, p_high_dim = 9L) {
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

evaluate_interval_fit <- function(fit_obj, method_label, method_id) {
  xnames <- grep("^x\\d+$", names(fit_obj$intervals), value = TRUE)
  eta_hat <- as.numeric(as.matrix(fit_obj$intervals[, xnames, drop = FALSE]) %*% fit_obj$cox_fit$coefficients)

  res <- fit_obj$intervals %>%
    mutate(
      eta_hat = eta_hat,
      covered = as.integer(T_star >= pi_lower & T_star <= pi_upper),
      length = pmax(pi_upper - pi_lower, 0),
      method_id = method_id,
      method_label = method_label
    )

  list(
    coverage = mean(res$covered),
    mean_len = mean_na(res$length),
    sd_len = sd_na(res$length),
    res = res
  )
}

run_one_rep_both_methods <- function(n,
                                     test_frac,
                                     p,
                                     beta,
                                     gamma,
                                     kappa,
                                     tau_quantile,
                                     cens_target,
                                     alpha,
                                     B,
                                     n_folds,
                                     seed) {
  sim <- simulate_rcensored_finite_tau(
    n = n,
    p = p,
    beta = beta,
    gamma = gamma,
    kappa = kappa,
    tau_quantile = tau_quantile,
    cens_target = cens_target,
    seed = seed
  )
  dat <- sim$dat
  tau <- sim$tau

  set.seed(seed + 10000L)
  idx <- sample.int(nrow(dat))
  n_test <- floor(test_frac * nrow(dat))
  test_id <- idx[seq_len(n_test)]
  train_id <- idx[(n_test + 1):nrow(dat)]

  train_df <- dat[train_id, , drop = FALSE]
  test_df <- dat[test_id, , drop = FALSE]

  current_fit <- conformal_pi_cox_Tstar(
    train_df = train_df,
    test_df = test_df,
    tau = tau,
    alpha = alpha,
    B = B,
    seed = seed + 20000L
  )

  crossfit_fit <- conformal_pi_cox_Tstar_crossfit(
    train_df = train_df,
    test_df = test_df,
    tau = tau,
    alpha = alpha,
    B = B,
    seed = seed + 30000L,
    n_folds = n_folds,
    fold_seed = seed + 40000L
  )

  current_eval <- evaluate_interval_fit(
    fit_obj = current_fit,
    method_label = "CP + Cox",
    method_id = "cp_cox_full"
  )

  crossfit_eval <- evaluate_interval_fit(
    fit_obj = crossfit_fit,
    method_label = "CP + Cox (Cross-Fit)",
    method_id = "cp_cox_crossfit"
  )

  list(
    tau = tau,
    censor_rate_emp = mean(train_df$delta == 0),
    admin_mass_emp = mean(train_df$delta == 0 & is_tau(train_df$y, tau)),
    train_n = nrow(train_df),
    test_n = nrow(test_df),
    method_rows = bind_rows(
      tibble(
        method_id = "cp_cox_full",
        method_label = "CP + Cox",
        coverage = current_eval$coverage,
        mean_len = current_eval$mean_len,
        sd_len = current_eval$sd_len
      ),
      tibble(
        method_id = "cp_cox_crossfit",
        method_label = "CP + Cox (Cross-Fit)",
        coverage = crossfit_eval$coverage,
        mean_len = crossfit_eval$mean_len,
        sd_len = crossfit_eval$sd_len
      )
    ),
    sample_rows = bind_rows(current_eval$res, crossfit_eval$res)
  )
}

summarise_by_method <- function(rep_df, alpha) {
  target <- 1 - alpha

  rep_df %>%
    group_by(
      setting_id, dgp_id, dgp_label, p, beta_scale, gamma, kappa,
      cens_target, tau_quantile, method_id, method_label, n, test_frac, alpha, B
    ) %>%
    summarise(
      n_rep = n(),
      tau_mean = mean(tau, na.rm = TRUE),
      coverage_mean = mean(coverage, na.rm = TRUE),
      coverage_sd = sd(coverage, na.rm = TRUE),
      length_mean = mean(mean_len, na.rm = TRUE),
      length_sd = sd(mean_len, na.rm = TRUE),
      censor_rate_emp_mean = mean(censor_rate_emp, na.rm = TRUE),
      admin_mass_emp_mean = mean(admin_mass_emp, na.rm = TRUE),
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
    arrange(method_label, dgp_id, cens_target, tau_quantile)
}

summarise_overall <- function(rep_df, alpha) {
  target <- 1 - alpha

  rep_df %>%
    group_by(method_id, method_label) %>%
    summarise(
      total_reps = n(),
      coverage_mean = mean(coverage, na.rm = TRUE),
      coverage_sd = sd(coverage, na.rm = TRUE),
      length_mean = mean(mean_len, na.rm = TRUE),
      length_sd = sd(mean_len, na.rm = TRUE),
      abs_cov_error = abs(coverage_mean - target),
      .groups = "drop"
    ) %>%
    arrange(method_label)
}

summarise_paired_differences <- function(rep_df) {
  rep_df %>%
    select(setting_id, rep, dgp_id, dgp_label, cens_target, tau_quantile, method_id, coverage, mean_len) %>%
    pivot_wider(
      names_from = method_id,
      values_from = c(coverage, mean_len)
    ) %>%
    mutate(
      coverage_diff_crossfit_minus_full = coverage_cp_cox_crossfit - coverage_cp_cox_full,
      length_diff_crossfit_minus_full = mean_len_cp_cox_crossfit - mean_len_cp_cox_full
    ) %>%
    group_by(setting_id, dgp_id, dgp_label, cens_target, tau_quantile) %>%
    summarise(
      n_rep = n(),
      coverage_diff_mean = mean(coverage_diff_crossfit_minus_full, na.rm = TRUE),
      coverage_diff_sd = sd(coverage_diff_crossfit_minus_full, na.rm = TRUE),
      length_diff_mean = mean(length_diff_crossfit_minus_full, na.rm = TRUE),
      length_diff_sd = sd(length_diff_crossfit_minus_full, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(dgp_id, cens_target, tau_quantile)
}

run_grid <- function(scenarios, cfg) {
  rep_rows <- list()
  sample_rows <- list()
  ptr <- 1L

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

      out <- run_one_rep_both_methods(
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
        n_folds = cfg$n_folds,
        seed = seed
      )

      rep_rows[[ptr]] <- out$method_rows %>%
        mutate(
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
          train_n = out$train_n,
          test_n = out$test_n
        )

      if (r == 1L) {
        n_take <- min(cfg$sample_n, nrow(out$sample_rows))
        sample_rows[[length(sample_rows) + 1L]] <- out$sample_rows %>%
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
          )
      }

      ptr <- ptr + 1L
    }
  }

  list(
    rep_df = bind_rows(rep_rows),
    sample_df = bind_rows(sample_rows),
    run_minutes = as.numeric(difftime(Sys.time(), start_time, units = "mins"))
  )
}

write_outputs <- function(scenarios, results, cfg, out_dir) {
  raw_dir <- file.path(out_dir, "raw")
  derived_dir <- file.path(out_dir, "derived")
  dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(derived_dir, recursive = TRUE, showWarnings = FALSE)

  by_method <- summarise_by_method(results$rep_df, alpha = cfg$alpha)
  overall <- summarise_overall(results$rep_df, alpha = cfg$alpha)
  paired <- summarise_paired_differences(results$rep_df)

  run_config <- tibble(
    tag = cfg$tag,
    out_dir = out_dir,
    methods = "CP + Cox; CP + Cox (Cross-Fit)",
    n_folds = cfg$n_folds,
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

  write_csv(by_method, file.path(derived_dir, "setting_method_summary.csv"))
  write_csv(overall, file.path(derived_dir, "overall_summary.csv"))
  write_csv(paired, file.path(derived_dir, "paired_setting_differences.csv"))

  writeLines(capture.output(sessionInfo()), con = file.path(raw_dir, "sessionInfo.txt"))
}

main <- function() {
  script_dir <- get_script_dir()
  project_root <- normalizePath(file.path(script_dir, "..", ".."), winslash = "/", mustWork = TRUE)

  source(file.path(project_root, "section-1-simulation", "00_utils.R"))
  source(file.path(project_root, "section-1-simulation", "01_dgp_weibull_ph.R"))
  source(file.path(project_root, "section-1-simulation", "02_conformal_cox_tau.R"))
  source(file.path(script_dir, "00_crossfit_support_fit_methods.R"))

  cfg <- parse_cli(commandArgs(trailingOnly = TRUE))
  scenarios <- build_scenarios(p_base = cfg$p, p_high_dim = cfg$high_dim_p)

  if (is.null(cfg$out_dir) || identical(cfg$out_dir, "")) {
    out_dir <- file.path(project_root, "playground", "crossfit-support-fit-comparison", "results", cfg$tag)
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
