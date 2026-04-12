# 03_simulation_section1.R

# source("00_utils.R")
# source("01_dgp_weibull_ph.R")
# source("02_conformal_cox_tau.R")

run_one_rep <- function(n = 1000,
                        test_frac = 0.3,
                        p = 3,
                        beta = rep(0.3, 3),
                        gamma = 1.2,
                        kappa = 1.0,
                        tau_quantile = 0.9,
                        cens_target = 0.2,
                        alpha = 0.1,
                        B = 1000,
                        seed = 1,
                        return_fit = FALSE,
                        pi_t0_grid = NULL) {
  sim <- simulate_rcensored_finite_tau(
    n = n, p = p, beta = beta, gamma = gamma, kappa = kappa,
    tau_quantile = tau_quantile, cens_target = cens_target, seed = seed
  )
  dat <- sim$dat
  tau <- sim$tau

  # Split
  set.seed(seed + 10000)
  idx <- sample.int(nrow(dat))
  n_test <- floor(test_frac * nrow(dat))
  test_id <- idx[1:n_test]
  train_id <- idx[(n_test + 1):nrow(dat)]

  train_df <- dat[train_id, ]
  test_df  <- dat[test_id, ]

  # Run conformal
  fit <- conformal_pi_cox_Tstar(
    train_df = train_df,
    test_df = test_df,
    tau = tau,
    alpha = alpha,
    B = B,
    seed = seed + 20000
  )

  xnames <- grep("^x\\d+$", names(fit$intervals), value = TRUE)
  eta_hat <- as.numeric(as.matrix(fit$intervals[, xnames, drop = FALSE]) %*% fit$cox_fit$coefficients)

  res <- fit$intervals %>%
    mutate(
      eta_hat = eta_hat,
      covered = as.integer(T_star >= pi_lower & T_star <= pi_upper),
      length = pmax(pi_upper - pi_lower, 0)
    )

  pi_curve <- NULL
  if (!is.null(pi_t0_grid)) {
    pi_curve <- bind_rows(lapply(seq_len(nrow(res)), function(i) {
      pred <- predict_pi_marginal_from_fit(
        fit_obj = fit,
        z_new = as.numeric(res[i, xnames, drop = TRUE]),
        t0_vec = pi_t0_grid
      )
      pred %>%
        mutate(test_row = i)
    }))
  }

  list(
    tau = tau,
    censor_rate_emp = mean(train_df$delta == 0),
    admin_mass_emp = mean(train_df$delta == 0 & is_tau(train_df$y, tau)),
    coverage = mean(res$covered),
    mean_len = mean_na(res$length),
    sd_len = sd_na(res$length),
    res = res,
    pi_curve = pi_curve,
    fit = if (return_fit) fit else NULL
  )
}

run_simulation <- function(n_rep = 20,
                           n = 1000,
                           test_frac = 0.3,
                           p = 3,
                           beta = rep(0.3, 3),
                           gamma = 1.2,
                           kappa = 1.0,
                           tau_quantile = 0.9,
                           cens_target = 0.2,
                           alpha = 0.1,
                           B = 1000,
                           seed0 = 1) {
  out <- vector("list", n_rep)
  for (r in 1:n_rep) {
    out[[r]] <- run_one_rep(
      n = n, test_frac = test_frac, p = p, beta = beta,
      gamma = gamma, kappa = kappa,
      tau_quantile = tau_quantile, cens_target = cens_target,
      alpha = alpha, B = B,
      seed = seed0 + r
    )
  }
  summ <- tibble(
    rep = 1:n_rep,
    tau = sapply(out, `[[`, "tau"),
    censor_rate_emp = sapply(out, `[[`, "censor_rate_emp"),
    admin_mass_emp = sapply(out, `[[`, "admin_mass_emp"),
    coverage = sapply(out, `[[`, "coverage"),
    mean_len = sapply(out, `[[`, "mean_len"),
    sd_len = sapply(out, `[[`, "sd_len")
  )

  list(
    summary = summ,
    coverage_mean = mean(summ$coverage),
    coverage_sd = sd(summ$coverage),
    length_mean = mean(summ$mean_len),
    length_sd = sd(summ$mean_len),
    details = out
  )
}

# ---- Example run ----
# sim_out <- run_simulation(
#   n_rep = 20,
#   n = 1000,
#   p = 3,
#   beta = rep(0.3, 3),
#   gamma = 1.2,
#   kappa = 1.0,
#   tau_quantile = 0.9,
#   cens_target = 0.2,
#   alpha = 0.1,
#   B = 1000,
#   seed0 = 123
# )
# print(sim_out$summary)
# cat("Mean coverage:", sim_out$coverage_mean, "\n")
# cat("Mean length:", sim_out$length_mean, "\n")
