# 01_run_challenging_model_comparison.R
#
# Preliminary six-method comparison under harder finite-window DGPs where
# the Cox working model is misspecified.

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
  if (length(hit) == 0) return(getwd())
  normalizePath(dirname(sub(key, "", hit[1])), winslash = "/", mustWork = TRUE)
}

is_absolute_path <- function(path) {
  grepl("^(/|[A-Za-z]:[/\\\\])", path)
}

parse_num_csv <- function(x) {
  vals <- str_split(x, ",", simplify = TRUE)
  vals <- trimws(vals)
  vals <- vals[nzchar(vals)]
  as.numeric(vals)
}

parse_chr_csv <- function(x) {
  vals <- str_split(x, ",", simplify = TRUE)
  vals <- trimws(vals)
  vals[nzchar(vals)]
}

pretty_dist_name <- function(dist) {
  switch(
    dist,
    weibull = "Weibull",
    lognormal = "Log-Normal",
    loglogistic = "Log-Logistic",
    exponential = "Exponential",
    str_to_title(dist)
  )
}

parse_cli <- function(args) {
  cfg <- list(
    tag = "prelim-v1",
    out_dir = NULL,
    scenario_ids = c("aft_lognormal_linear", "aft_lognormal_nonlinear"),
    n_rep = 20L,
    n = 2000L,
    test_frac = 0.30,
    p = 3L,
    beta_scale = 0.35,
    sigma = 0.70,
    tau_quantiles = c(0.70, 0.80, 0.90, 0.99),
    cens_target = 0.20,
    alpha = 0.10,
    B = 1500L,
    wrong_dist = "lognormal",
    seed0 = 20260414L
  )

  for (a in args) {
    if (!startsWith(a, "--") || !str_detect(a, "=")) next
    kv <- str_split_fixed(sub("^--", "", a), "=", 2)
    key <- kv[1]
    val <- kv[2]

    if (key %in% c("tag", "out_dir", "wrong_dist")) cfg[[key]] <- val
    if (key == "scenario_ids") cfg[[key]] <- parse_chr_csv(val)
    if (key == "tau_quantiles") cfg[[key]] <- parse_num_csv(val)
    if (key %in% c("n_rep", "n", "p", "B", "seed0")) cfg[[key]] <- as.integer(val)
    if (key %in% c("test_frac", "beta_scale", "sigma", "cens_target", "alpha")) {
      cfg[[key]] <- as.numeric(val)
    }
  }

  cfg$tau_quantiles <- sort(unique(as.numeric(cfg$tau_quantiles)))
  cfg
}

base_method_manifest <- function(wrong_dist) {
  wrong_label <- pretty_dist_name(wrong_dist)

  tibble::tribble(
    ~method_id, ~method_label, ~estimator_family, ~working_model, ~is_conformal,
    "plugin_weibull", "Weibull Plug-in", "plug-in", "Weibull", FALSE,
    "plugin_cox", "Cox Plug-in", "plug-in", "Cox", FALSE,
    paste0("plugin_", wrong_dist), paste0(wrong_label, " Plug-in"), "plug-in", wrong_label, FALSE,
    "cp_weibull", "CP + Weibull", "conformal", "Weibull", TRUE,
    "cp_cox", "CP + Cox", "conformal", "Cox", TRUE,
    paste0("cp_", wrong_dist), paste0("CP + ", wrong_label), "conformal", wrong_label, TRUE
  )
}

fit_all_methods <- function(train_df,
                            test_df,
                            tau,
                            alpha,
                            B,
                            seed,
                            wrong_dist) {
  wrong_plugin_id <- paste0("plugin_", wrong_dist)
  wrong_cp_id <- paste0("cp_", wrong_dist)

  fits <- list(
    plugin_weibull = survreg_plugin_pi_Tstar_tau_upper(
      train_df = train_df, test_df = test_df, tau = tau, alpha = alpha, dist = "weibull"
    ),
    plugin_cox = cox_plugin_pi_Tstar_tau_upper(
      train_df = train_df, test_df = test_df, tau = tau, alpha = alpha
    ),
    cp_weibull = conformal_pi_survreg_Tstar(
      train_df = train_df, test_df = test_df, tau = tau, alpha = alpha,
      B = B, seed = seed + 20000L, dist = "weibull"
    ),
    cp_cox = conformal_pi_cox_Tstar(
      train_df = train_df, test_df = test_df, tau = tau, alpha = alpha,
      B = B, seed = seed + 30000L
    )
  )

  fits[[wrong_plugin_id]] <- survreg_plugin_pi_Tstar_tau_upper(
    train_df = train_df, test_df = test_df, tau = tau, alpha = alpha, dist = wrong_dist
  )
  fits[[wrong_cp_id]] <- conformal_pi_survreg_Tstar(
    train_df = train_df, test_df = test_df, tau = tau, alpha = alpha,
    B = B, seed = seed + 40000L, dist = wrong_dist
  )

  fits
}

collect_method_results <- function(fits, manifest, test_df, scenario_id) {
  method_ids <- names(fits)

  bind_rows(lapply(method_ids, function(id) {
    meta <- manifest %>% filter(method_id == id)
    if (nrow(meta) != 1L) stop("Method manifest mismatch for method_id = ", id)

    fits[[id]]$intervals %>%
      mutate(
        covered = as.integer(T_star >= pi_lower & T_star <= pi_upper),
        length = pmax(pi_upper - pi_lower, 0),
        method_id = meta$method_id,
        method_label = meta$method_label,
        estimator_family = meta$estimator_family,
        working_model = meta$working_model,
        is_conformal = meta$is_conformal,
        spec_status = scenario_spec_status(scenario_id, meta$working_model)
      )
  }))
}

summarise_method_metrics <- function(rep_df, alpha) {
  target <- 1 - alpha

  rep_df %>%
    group_by(
      scenario_id, scenario_label, scenario_description, tau_quantile,
      method_id, method_label, estimator_family, working_model, spec_status,
      is_conformal, n, test_frac, p, beta_scale, sigma, cens_target, alpha, B
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
    arrange(scenario_id, tau_quantile, estimator_family, working_model)
}

summarise_overall_by_method <- function(method_tau_summary) {
  method_tau_summary %>%
    group_by(
      scenario_id, scenario_label, method_id, method_label,
      estimator_family, working_model, spec_status, is_conformal
    ) %>%
    summarise(
      n_tau = n(),
      coverage_mean_avg = mean(coverage_mean, na.rm = TRUE),
      length_mean_avg = mean(length_mean, na.rm = TRUE),
      abs_cov_error_avg = mean(abs_cov_error, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(scenario_id, estimator_family, working_model)
}

write_outputs <- function(rep_df, cfg, out_dir) {
  raw_dir <- file.path(out_dir, "raw")
  derived_dir <- file.path(out_dir, "derived")
  dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(derived_dir, recursive = TRUE, showWarnings = FALSE)

  method_tau_summary <- summarise_method_metrics(rep_df, alpha = cfg$alpha)
  method_overall_summary <- summarise_overall_by_method(method_tau_summary)

  run_config <- tibble(
    tag = cfg$tag,
    out_dir = out_dir,
    interval_form = "[lower, tau] for all six methods",
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
    wrong_dist = cfg$wrong_dist,
    seed0 = cfg$seed0,
    created_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  )

  write_csv(rep_df, file.path(raw_dir, "rep_level_results.csv"))
  write_csv(run_config, file.path(raw_dir, "run_config.csv"))
  write_csv(method_tau_summary, file.path(derived_dir, "method_tau_summary.csv"))
  write_csv(method_overall_summary, file.path(derived_dir, "method_overall_summary.csv"))
  writeLines(capture.output(sessionInfo()), con = file.path(raw_dir, "sessionInfo.txt"))
}

main <- function() {
  script_dir <- get_script_dir()
  project_root <- normalizePath(file.path(script_dir, "..", ".."), winslash = "/", mustWork = TRUE)

  source(file.path(project_root, "section-1-simulation", "00_utils.R"))
  source(file.path(project_root, "section-1-simulation", "01_dgp_weibull_ph.R"))
  source(file.path(project_root, "section-1-simulation", "02_conformal_cox_tau.R"))
  source(file.path(project_root, "section-1-simulation", "09_cox_plugin_Tstar.R"))
  source(file.path(project_root, "playground", "model-comparison", "00_survreg_working_models.R"))
  source(file.path(script_dir, "00_dgp_challenging.R"))

  cfg <- parse_cli(commandArgs(trailingOnly = TRUE))
  catalog <- scenario_catalog() %>% filter(scenario_id %in% cfg$scenario_ids)
  manifest <- base_method_manifest(cfg$wrong_dist)

  out_dir <- if (is.null(cfg$out_dir) || identical(cfg$out_dir, "")) {
    file.path(project_root, "playground", "challenging-dgp-comparison", "results", cfg$tag)
  } else if (is_absolute_path(cfg$out_dir)) {
    cfg$out_dir
  } else {
    file.path(project_root, cfg$out_dir)
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  rep_rows <- list()

  for (scenario_id in catalog$scenario_id) {
    scenario_label <- catalog$scenario_label[match(scenario_id, catalog$scenario_id)]
    scenario_description <- catalog$description[match(scenario_id, catalog$scenario_id)]

    for (tau_q in cfg$tau_quantiles) {
      message(sprintf("Running %s at tau quantile %.2f", scenario_label, tau_q))

      for (r in seq_len(cfg$n_rep)) {
        seed <- cfg$seed0 +
          match(scenario_id, catalog$scenario_id) * 1000000L +
          as.integer(round(1000 * tau_q)) * 10000L + r

        message(sprintf("[%s] scenario=%s tau_q=%.2f rep=%d/%d seed=%d",
                        format(Sys.time(), "%H:%M:%S"),
                        scenario_id, tau_q, r, cfg$n_rep, seed))

        sim <- simulate_rcensored_challenging_tau(
          n = cfg$n,
          scenario_id = scenario_id,
          p = cfg$p,
          beta_scale = cfg$beta_scale,
          sigma = cfg$sigma,
          tau_quantile = tau_q,
          cens_target = cfg$cens_target,
          seed = seed
        )

        dat <- sim$dat
        tau <- sim$tau

        set.seed(seed + 10000L)
        idx <- sample.int(nrow(dat))
        n_test <- floor(cfg$test_frac * nrow(dat))
        test_id <- idx[seq_len(n_test)]
        train_id <- idx[(n_test + 1):nrow(dat)]

        train_df <- dat[train_id, , drop = FALSE]
        test_df <- dat[test_id, , drop = FALSE]

        fits <- fit_all_methods(
          train_df = train_df,
          test_df = test_df,
          tau = tau,
          alpha = cfg$alpha,
          B = cfg$B,
          seed = seed,
          wrong_dist = cfg$wrong_dist
        )

        all_res <- collect_method_results(
          fits = fits,
          manifest = manifest,
          test_df = test_df,
          scenario_id = scenario_id
        )

        censor_rate_emp_train <- mean(train_df$delta == 0)
        admin_mass_emp_train <- mean(train_df$delta == 0 & is_tau(train_df$y, tau))

        metrics <- all_res %>%
          group_by(method_id, method_label, estimator_family, working_model, spec_status, is_conformal) %>%
          summarise(
            coverage = mean(covered),
            mean_len = mean_na(length),
            sd_len = sd_na(length),
            .groups = "drop"
          ) %>%
          mutate(
            scenario_id = scenario_id,
            scenario_label = scenario_label,
            scenario_description = scenario_description,
            rep = r,
            seed = seed,
            tau_quantile = tau_q,
            tau = tau,
            censor_rate_emp = censor_rate_emp_train,
            admin_mass_emp = admin_mass_emp_train,
            n = cfg$n,
            test_frac = cfg$test_frac,
            p = cfg$p,
            beta_scale = cfg$beta_scale,
            sigma = cfg$sigma,
            cens_target = cfg$cens_target,
            alpha = cfg$alpha,
            B = cfg$B
          ) %>%
          select(
            scenario_id, scenario_label, scenario_description, tau_quantile, rep, seed,
            method_id, method_label, estimator_family, working_model, spec_status, is_conformal,
            n, test_frac, p, beta_scale, sigma, cens_target, alpha, B,
            tau, censor_rate_emp, admin_mass_emp, coverage, mean_len, sd_len
          )

        rep_rows[[length(rep_rows) + 1L]] <- metrics
      }
    }
  }

  rep_df <- bind_rows(rep_rows)
  write_outputs(rep_df, cfg, out_dir)
  message("Saved outputs to: ", out_dir)
}

if (sys.nframe() == 0L) {
  main()
}
