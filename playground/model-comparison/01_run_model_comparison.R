# 01_run_model_comparison.R
#
# Expanded working-model comparison under the baseline Weibull PH DGP.
# The comparison includes:
#   - correctly specified Weibull plug-in
#   - correctly specified Cox plug-in
#   - misspecified parametric plug-in
#   - conformal + Weibull working model
#   - conformal + Cox working model
#   - conformal + misspecified parametric working model

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

parse_bool <- function(x) {
  x <- tolower(trimws(x))
  if (x %in% c("true", "t", "1", "yes", "y")) return(TRUE)
  if (x %in% c("false", "f", "0", "no", "n")) return(FALSE)
  stop("Cannot parse logical value from: ", x)
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
    tag = "baseline-full-v5-all-upper",
    out_dir = NULL,
    n_rep = 30L,
    n = 2000L,
    test_frac = 0.30,
    p = 3L,
    beta_scale = 0.30,
    gamma = 1.2,
    kappa = 1.0,
    tau_quantiles = c(0.70, 0.80, 0.90, 0.99),
    cens_target = 0.20,
    alpha = 0.10,
    B = 1500L,
    seed0 = 20260321L,
    sample_n = 140L,
    wrong_dist = "lognormal"
  )

  for (a in args) {
    if (!startsWith(a, "--") || !str_detect(a, "=")) next
    kv <- str_split_fixed(sub("^--", "", a), "=", 2)
    key <- kv[1]
    val <- kv[2]

    if (key %in% c("tag", "out_dir", "wrong_dist")) cfg[[key]] <- val
    if (key == "tau_quantiles") cfg[[key]] <- parse_num_csv(val)
    if (key %in% c("n_rep", "n", "p", "B", "seed0", "sample_n")) cfg[[key]] <- as.integer(val)
    if (key %in% c("test_frac", "beta_scale", "gamma", "kappa", "cens_target", "alpha")) {
      cfg[[key]] <- as.numeric(val)
    }
  }

  cfg$tau_quantiles <- sort(unique(as.numeric(cfg$tau_quantiles)))
  cfg
}

model_manifest <- function(wrong_dist) {
  wrong_label <- pretty_dist_name(wrong_dist)

  tibble::tribble(
    ~method_id, ~method_label, ~estimator_family, ~working_model, ~spec_status, ~is_conformal,
    "plugin_weibull", "Weibull Plug-in", "plug-in", "Weibull", "correct", FALSE,
    "plugin_cox", "Cox Plug-in", "plug-in", "Cox", "correct", FALSE,
    paste0("plugin_", wrong_dist), paste0(wrong_label, " Plug-in"), "plug-in", wrong_label, "misspecified", FALSE,
    "cp_weibull", "CP + Weibull", "conformal", "Weibull", "correct", TRUE,
    "cp_cox", "CP + Cox", "conformal", "Cox", "correct", TRUE,
    paste0("cp_", wrong_dist), paste0("CP + ", wrong_label), "conformal", wrong_label, "misspecified", TRUE
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

collect_method_results <- function(fits, manifest, test_df, beta_true) {
  xnames <- grep("^x\\d+$", names(test_df), value = TRUE)
  eta_true <- as.numeric(as.matrix(test_df[, xnames, drop = FALSE]) %*% beta_true)
  method_ids <- names(fits)

  bind_rows(lapply(method_ids, function(id) {
    meta <- manifest %>% filter(method_id == id)
    if (nrow(meta) != 1L) {
      stop("Method manifest did not contain exactly one row for method_id = ", id)
    }

    fits[[id]]$intervals %>%
      mutate(
        test_row_index = row_number(),
        eta_true = eta_true,
        covered = as.integer(T_star >= pi_lower & T_star <= pi_upper),
        length = pmax(pi_upper - pi_lower, 0),
        method_id = meta$method_id,
        method_label = meta$method_label,
        estimator_family = meta$estimator_family,
        working_model = meta$working_model,
        spec_status = meta$spec_status,
        is_conformal = meta$is_conformal
      )
  }))
}

summarise_method_metrics <- function(rep_df, alpha) {
  target <- 1 - alpha

  rep_df %>%
    group_by(
      tau_quantile, method_id, method_label, estimator_family, working_model,
      spec_status, is_conformal, n, test_frac, p, beta_scale, gamma, kappa,
      cens_target, alpha, B
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
    arrange(tau_quantile, estimator_family, working_model)
}

summarise_overall_by_method <- function(method_tau_summary) {
  method_tau_summary %>%
    group_by(method_id, method_label, estimator_family, working_model, spec_status, is_conformal) %>%
    summarise(
      n_tau = n(),
      coverage_mean_avg = mean(coverage_mean, na.rm = TRUE),
      length_mean_avg = mean(length_mean, na.rm = TRUE),
      abs_cov_error_avg = mean(abs_cov_error, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(estimator_family, working_model)
}

write_outputs <- function(rep_df, sample_df, cfg, out_dir) {
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
    n_rep = cfg$n_rep,
    n = cfg$n,
    test_frac = cfg$test_frac,
    p = cfg$p,
    beta_scale = cfg$beta_scale,
    gamma = cfg$gamma,
    kappa = cfg$kappa,
    tau_quantiles = paste(cfg$tau_quantiles, collapse = ","),
    cens_target = cfg$cens_target,
    alpha = cfg$alpha,
    B = cfg$B,
    wrong_dist = cfg$wrong_dist,
    seed0 = cfg$seed0,
    sample_n = cfg$sample_n,
    created_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  )

  write_csv(rep_df, file.path(raw_dir, "rep_level_results.csv"))
  write_csv(sample_df, file.path(raw_dir, "test_level_samples_rep1_by_tau.csv"))
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
  source(file.path(script_dir, "00_survreg_working_models.R"))

  cfg <- parse_cli(commandArgs(trailingOnly = TRUE))
  manifest <- model_manifest(cfg$wrong_dist)
  beta_true <- rep(cfg$beta_scale, cfg$p)

  out_dir <- if (is.null(cfg$out_dir) || identical(cfg$out_dir, "")) {
    file.path(project_root, "playground", "model-comparison", "results", cfg$tag)
  } else if (is_absolute_path(cfg$out_dir)) {
    cfg$out_dir
  } else {
    file.path(project_root, cfg$out_dir)
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  rep_rows <- list()
  sample_rows <- list()

  for (tau_q in cfg$tau_quantiles) {
    message(sprintf("Running tau quantile %.2f", tau_q))

    for (r in seq_len(cfg$n_rep)) {
      seed <- cfg$seed0 + as.integer(round(1000 * tau_q)) * 10000L + r
      message(sprintf("[%s] tau_q=%.2f rep=%d/%d seed=%d",
                      format(Sys.time(), "%H:%M:%S"), tau_q, r, cfg$n_rep, seed))

      sim <- simulate_rcensored_finite_tau(
        n = cfg$n,
        p = cfg$p,
        beta = beta_true,
        gamma = cfg$gamma,
        kappa = cfg$kappa,
        tau_quantile = tau_q,
        cens_target = cfg$cens_target,
        seed = seed
      )
      dat <- sim$dat
      tau <- sim$tau

      set.seed(seed + 10000L)
      idx <- sample.int(nrow(dat))
      n_test <- floor(cfg$test_frac * nrow(dat))
      test_id <- idx[1:n_test]
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
        beta_true = beta_true
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
          gamma = cfg$gamma,
          kappa = cfg$kappa,
          cens_target = cfg$cens_target,
          alpha = cfg$alpha,
          B = cfg$B
        ) %>%
        select(
          tau_quantile, rep, seed, method_id, method_label, estimator_family, working_model,
          spec_status, is_conformal, n, test_frac, p, beta_scale, gamma, kappa,
          cens_target, alpha, B, tau, censor_rate_emp, admin_mass_emp,
          coverage, mean_len, sd_len
        )

      rep_rows[[length(rep_rows) + 1L]] <- metrics

      if (r == 1L) {
        set.seed(seed + 20000L)
        n_take <- min(cfg$sample_n, nrow(test_df))
        selected_ids <- sort(sample(seq_len(nrow(test_df)), size = n_take, replace = FALSE))

        rank_df <- tibble(
          test_row_index = selected_ids,
          eta_true = all_res %>%
            filter(method_id == "plugin_weibull", test_row_index %in% selected_ids) %>%
            arrange(test_row_index) %>%
            pull(eta_true)
        ) %>%
          arrange(eta_true) %>%
          mutate(obs_rank = row_number())

        sample_rows[[length(sample_rows) + 1L]] <- all_res %>%
          filter(test_row_index %in% selected_ids) %>%
          left_join(rank_df, by = c("test_row_index", "eta_true")) %>%
          mutate(
            tau_quantile = tau_q,
            tau = tau,
            rep = r,
            seed = seed
          ) %>%
          arrange(tau_quantile, method_id, obs_rank)
      }
    }
  }

  rep_df <- bind_rows(rep_rows)
  sample_df <- bind_rows(sample_rows)

  write_outputs(rep_df, sample_df, cfg, out_dir)

  message("Saved outputs to: ", out_dir)
}

if (sys.nframe() == 0L) {
  main()
}
