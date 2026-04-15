# 00_dgp_challenging.R
#
# Harder finite-window survival DGPs for robustness checks.
# These are intentionally more challenging than the thesis baseline Weibull PH
# setup and are designed so that the Cox working model is misspecified.

suppressPackageStartupMessages({
  library(dplyr)
})

scenario_catalog <- function() {
  tibble::tribble(
    ~scenario_id, ~scenario_label, ~description,
    "aft_lognormal_linear", "Log-Normal AFT (Linear)",
    "Linear log-normal accelerated-failure-time DGP; log-normal is correct, Cox and Weibull are misspecified.",
    "aft_lognormal_nonlinear", "Log-Normal AFT (Nonlinear)",
    "Nonlinear log-normal accelerated-failure-time DGP; all three working models are misspecified."
  )
}

scenario_spec_status <- function(scenario_id, working_model) {
  stopifnot(length(scenario_id) == 1L, length(working_model) == 1L)

  if (scenario_id == "aft_lognormal_linear") {
    return(
      dplyr::case_when(
        working_model == "Log-Normal" ~ "correct",
        TRUE ~ "misspecified"
      )
    )
  }

  if (scenario_id == "aft_lognormal_nonlinear") {
    return("misspecified")
  }

  stop("Unknown scenario_id: ", scenario_id)
}

scenario_true_score <- function(Z, scenario_id, beta_scale = 0.35) {
  Z <- as.matrix(Z)
  p <- ncol(Z)
  if (p < 3L) stop("This robustness playground expects at least 3 covariates.")

  if (scenario_id == "aft_lognormal_linear") {
    beta <- rep(beta_scale, p)
    return(as.numeric(Z %*% beta))
  }

  if (scenario_id == "aft_lognormal_nonlinear") {
    score <- beta_scale * (
      0.9 * Z[, 1] +
        0.8 * (Z[, 2]^2 - 1) +
        0.9 * sin(Z[, 3])
    )

    if (p > 3L) {
      score <- score + 0.15 * rowSums(Z[, 4:p, drop = FALSE])
    }

    return(as.numeric(score))
  }

  stop("Unknown scenario_id: ", scenario_id)
}

rchallenging_T <- function(n,
                           Z,
                           scenario_id,
                           beta_scale = 0.35,
                           sigma = 0.70) {
  eps <- rnorm(n)
  eta <- scenario_true_score(Z = Z, scenario_id = scenario_id, beta_scale = beta_scale)

  # Larger eta corresponds to earlier events by shifting log T downward.
  log_T <- -eta + sigma * eps
  T <- exp(log_T)

  list(T = T, eta = eta)
}

simulate_rcensored_challenging_tau <- function(n,
                                               scenario_id,
                                               p = 3,
                                               beta_scale = 0.35,
                                               sigma = 0.70,
                                               tau_quantile = 0.90,
                                               cens_target = 0.20,
                                               seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  catalog <- scenario_catalog()
  if (!scenario_id %in% catalog$scenario_id) {
    stop("Unknown scenario_id: ", scenario_id)
  }

  Z <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(Z) <- paste0("x", seq_len(p))

  draw <- rchallenging_T(
    n = n,
    Z = Z,
    scenario_id = scenario_id,
    beta_scale = beta_scale,
    sigma = sigma
  )
  T <- draw$T
  eta_true <- draw$eta

  tau <- as.numeric(stats::quantile(T, probs = tau_quantile))

  rate_c <- calibrate_censor_rate(T, tau, cens_target = cens_target)
  C_latent <- rexp(n, rate = rate_c)
  C <- pmin(C_latent, tau)

  Y <- pmin(T, C)
  Delta <- as.integer(T <= C)
  T_star <- pmin(T, tau)

  dat <- data.frame(
    y = Y,
    delta = Delta,
    C = C,
    T_true = T,
    T_star = T_star,
    tau = tau,
    eta_true = eta_true
  ) %>%
    bind_cols(as.data.frame(Z))

  list(
    dat = dat,
    tau = tau,
    rate_c = rate_c,
    scenario_id = scenario_id,
    scenario_label = catalog$scenario_label[match(scenario_id, catalog$scenario_id)]
  )
}
