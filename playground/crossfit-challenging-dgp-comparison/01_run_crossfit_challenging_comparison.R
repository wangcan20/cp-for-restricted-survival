# 01_run_crossfit_challenging_comparison.R
#
# Compare the current CP + Cox upper-tau method against the cross-fitted
# variant under harder finite-window DGPs where the Cox working model is
# misspecified.

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

parse_chr_csv <- function(x) {
  vals <- str_split(x, ",", simplify = TRUE)
  vals <- trimws(vals)
  vals[nzchar(vals)]
}

parse_num_csv <- function(x) {
  vals <- str_split(x, ",", simplify = TRUE)
  vals <- trimws(vals)
  vals <- vals[nzchar(vals)]
  as.numeric(vals)
}

parse_cli <- function(args) {
  cfg <- list(
    tag = "matched-v1",
    out_dir = NULL,
    scenario_ids = c("aft_lognormal_linear", "aft_lognormal_nonlinear"),
    n_rep = 30L,
    n = 2000L,
    test_frac = 0.30,
    p = 3L,
    beta_scale = 0.30,
    sigma = 0.70,
    tau_quantiles = c(0.70, 0.80, 0.90, 0.99),
    cens_target = 0.20,
    alpha = 0.10,
    B = 1500L,
    n_folds = 5L,
    seed0 = 20260414L
  )

  for (a in args) {
    if (!startsWith(a, "--") || !str_detect(a, "=")) next
    kv <- str_split_fixed(sub("^--", "", a), "=", 2)
    key <- kv[1]
    val <- kv[2]

    if (key %in% c("tag", "out_dir")) cfg[[key]] <- val
    if (key == "scenario_ids") cfg[[key]] <- parse_chr_csv(val)
    if (key == "tau_quantiles") cfg[[key]] <- parse_num_csv(val)
    if (key %in% c("n_rep", "n", "p", "B", "n_folds", "seed0")) cfg[[key]] <- as.integer(val)
    if (key %in% c("test_frac", "beta_scale", "sigma", "cens_target", "alpha")) {
      cfg[[key]] <- as.numeric(val)
    }
  }

  cfg$tau_quantiles <- sort(unique(as.numeric(cfg$tau_quantiles)))
  cfg
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
                                     scenario_id,
                                     p,
                                     beta_scale,
                                     sigma,
                                     tau_quantile,
                                     cens_target,
                                     alpha,
                                     B,
                                     n_folds,
                                     seed) {
  sim <- simulate_rcensored_challenging_tau(
    n = n,
    scenario_id = scenario_id,
    p = p,
    beta_scale = beta_scale,
    sigma = sigma,
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
    )
  )
}

summarise_by_method <- function(rep_df, alpha) {
  target <- 1 - alpha

  rep_df %>%
    group_by(
      scenario_id, scenario_label, tau_quantile, method_id, method_label,
      n, test_frac, p, beta_scale, sigma, cens_target, alpha, B, n_folds
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
      tcrit = if_else(n_rep > 1, qt(0.975, df = pmax(n_rep - 1, 1)), NA_real_),
      coverage_se = coverage_sd / sqrt(n_rep),
      length_se = length_sd / sqrt(n_rep),
      coverage_ci_low = coverage_mean - tcrit * coverage_se,
      coverage_ci_high = coverage_mean + tcrit * coverage_se,
      length_ci_low = length_mean - tcrit * length_se,
      length_ci_high = length_mean + tcrit * length_se,
      target_coverage = target,
      abs_cov_error = abs(coverage_mean - target)
    ) %>%
    arrange(scenario_id, tau_quantile, method_label)
}

summarise_overall <- function(rep_df, alpha) {
  target <- 1 - alpha

  rep_df %>%
    group_by(scenario_id, scenario_label, method_id, method_label) %>%
    summarise(
      total_reps = n(),
      coverage_mean = mean(coverage, na.rm = TRUE),
      coverage_sd = sd(coverage, na.rm = TRUE),
      length_mean = mean(mean_len, na.rm = TRUE),
      length_sd = sd(mean_len, na.rm = TRUE),
      abs_cov_error = abs(coverage_mean - target),
      .groups = "drop"
    ) %>%
    arrange(scenario_id, method_label)
}

summarise_paired_differences <- function(rep_df) {
  rep_df %>%
    select(scenario_id, scenario_label, tau_quantile, rep, method_id, coverage, mean_len) %>%
    pivot_wider(
      names_from = method_id,
      values_from = c(coverage, mean_len)
    ) %>%
    mutate(
      coverage_diff_crossfit_minus_full = coverage_cp_cox_crossfit - coverage_cp_cox_full,
      length_diff_crossfit_minus_full = mean_len_cp_cox_crossfit - mean_len_cp_cox_full
    ) %>%
    group_by(scenario_id, scenario_label, tau_quantile) %>%
    summarise(
      n_rep = n(),
      coverage_diff_mean = mean(coverage_diff_crossfit_minus_full, na.rm = TRUE),
      length_diff_mean = mean(length_diff_crossfit_minus_full, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(scenario_id, tau_quantile)
}

main <- function() {
  script_dir <- get_script_dir()
  project_root <- normalizePath(file.path(script_dir, "..", ".."), winslash = "/", mustWork = TRUE)

  source(file.path(project_root, "section-1-simulation", "00_utils.R"))
  source(file.path(project_root, "section-1-simulation", "01_dgp_weibull_ph.R"))
  source(file.path(project_root, "section-1-simulation", "02_conformal_cox_tau.R"))
  source(file.path(project_root, "playground", "crossfit-support-fit-comparison", "00_crossfit_support_fit_methods.R"))
  source(file.path(project_root, "playground", "challenging-dgp-comparison", "00_dgp_challenging.R"))

  cfg <- parse_cli(commandArgs(trailingOnly = TRUE))
  catalog <- scenario_catalog() %>% filter(scenario_id %in% cfg$scenario_ids)

  out_dir <- if (is.null(cfg$out_dir) || identical(cfg$out_dir, "")) {
    file.path(project_root, "playground", "crossfit-challenging-dgp-comparison", "results", cfg$tag)
  } else if (is_absolute_path(cfg$out_dir)) {
    cfg$out_dir
  } else {
    file.path(project_root, cfg$out_dir)
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  raw_dir <- file.path(out_dir, "raw")
  derived_dir <- file.path(out_dir, "derived")
  dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(derived_dir, recursive = TRUE, showWarnings = FALSE)

  rep_rows <- list()
  ptr <- 1L

  for (scenario_id in catalog$scenario_id) {
    scenario_label <- catalog$scenario_label[match(scenario_id, catalog$scenario_id)]
    for (tau_q in cfg$tau_quantiles) {
      message(sprintf("Running %s at tau quantile %.2f", scenario_label, tau_q))
      for (r in seq_len(cfg$n_rep)) {
        seed <- cfg$seed0 +
          match(scenario_id, catalog$scenario_id) * 1000000L +
          as.integer(round(1000 * tau_q)) * 10000L + r

        message(sprintf("[%s] scenario=%s tau_q=%.2f rep=%d/%d seed=%d",
                        format(Sys.time(), "%H:%M:%S"),
                        scenario_id, tau_q, r, cfg$n_rep, seed))

        out <- run_one_rep_both_methods(
          n = cfg$n,
          test_frac = cfg$test_frac,
          scenario_id = scenario_id,
          p = cfg$p,
          beta_scale = cfg$beta_scale,
          sigma = cfg$sigma,
          tau_quantile = tau_q,
          cens_target = cfg$cens_target,
          alpha = cfg$alpha,
          B = cfg$B,
          n_folds = cfg$n_folds,
          seed = seed
        )

        rep_rows[[ptr]] <- out$method_rows %>%
          mutate(
            scenario_id = scenario_id,
            scenario_label = scenario_label,
            rep = as.integer(r),
            seed = as.integer(seed),
            tau_quantile = tau_q,
            tau = out$tau,
            n = cfg$n,
            test_frac = cfg$test_frac,
            p = cfg$p,
            beta_scale = cfg$beta_scale,
            sigma = cfg$sigma,
            cens_target = cfg$cens_target,
            alpha = cfg$alpha,
            B = cfg$B,
            n_folds = cfg$n_folds,
            censor_rate_emp = out$censor_rate_emp,
            admin_mass_emp = out$admin_mass_emp
          )
        ptr <- ptr + 1L
      }
    }
  }

  rep_df <- bind_rows(rep_rows)
  summary_df <- summarise_by_method(rep_df, alpha = cfg$alpha)
  overall_df <- summarise_overall(rep_df, alpha = cfg$alpha)
  paired_df <- summarise_paired_differences(rep_df)

  run_config <- tibble(
    tag = cfg$tag,
    out_dir = out_dir,
    scenario_ids = paste(cfg$scenario_ids, collapse = ","),
    n_rep = cfg$n_rep,
    n = cfg$n,
    test_frac = cfg$test_frac,
    p = cfg$p,
    beta_scale = cfg$beta_scale,
    sigma = cfg$sigma,
    tau_quantiles = paste(cfg$tau_quantiles, collapse = ","),
    cens_target = cfg$cens_target,
    alpha = cfg$alpha,
    B = cfg$B,
    n_folds = cfg$n_folds,
    seed0 = cfg$seed0,
    created_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  )

  write_csv(rep_df, file.path(raw_dir, "rep_level_results.csv"))
  write_csv(run_config, file.path(raw_dir, "run_config.csv"))
  write_csv(summary_df, file.path(derived_dir, "setting_method_summary.csv"))
  write_csv(overall_df, file.path(derived_dir, "overall_summary.csv"))
  write_csv(paired_df, file.path(derived_dir, "paired_setting_differences.csv"))
  writeLines(capture.output(sessionInfo()), con = file.path(raw_dir, "sessionInfo.txt"))

  message("Saved outputs to: ", out_dir)
}

if (sys.nframe() == 0L) {
  main()
}
