# 09_cox_plugin_Tstar.R

# Cox plug-in predictive interval for T* = min(T, tau).
# This is a semi-parametric comparator with no conformal calibration.

# source("00_utils.R")
# source("02_conformal_cox_tau.R")

cox_tstar_quantile <- function(bh_time, bh_haz, eta, q, tau) {
  q <- min(max(as.numeric(q), 0), 1)
  if (is.na(q)) return(NA_real_)
  if (q <= 0) return(0)
  if (q >= 1) return(tau)

  if (length(bh_time) == 0 || length(bh_haz) == 0) {
    stop("Baseline hazard grid is empty; cannot form Cox plug-in quantiles.")
  }

  H_tau <- bh_haz[length(bh_haz)]
  F_tau <- 1 - exp(-H_tau * exp(eta))

  # The target is the mixed distribution of T* = min(T, tau):
  #   F_{T*}(t | z) = F_T(t | z),  t < tau,
  #   F_{T*}(tau | z) = 1.
  # Hence Q_{T*}(q | z) equals the Cox quantile for T when q < F_T(tau | z),
  # and equals tau once the requested quantile lands in the atom at tau.
  if (q >= F_tau) return(tau)

  target <- -log1p(-q) / exp(eta)

  if (!is.finite(target)) return(tau)

  idx <- which(bh_haz >= target)
  if (length(idx) == 0) return(tau)
  bh_time[min(idx)]
}

cox_plugin_pi_Tstar <- function(train_df,
                                test_df,
                                tau,
                                alpha = 0.1) {
  xnames <- grep("^x\\d+$", names(train_df), value = TRUE)
  if (length(xnames) == 0) stop("No covariate columns named x1, x2, ... found.")

  fit_obj <- fit_cox_working(train_df, xnames)
  cox_fit <- fit_obj$fit
  bh <- fit_obj$bh %>% arrange(time)

  bh_tau <- bh %>% filter(time <= tau)
  if (nrow(bh_tau) == 0) {
    stop("No baseline hazard jumps at or before tau; cannot form Cox plug-in interval.")
  }

  bh_time <- bh_tau$time
  bh_haz <- bh_tau$hazard

  if (tail(bh_time, 1) < tau) {
    bh_time <- c(bh_time, tau)
    bh_haz <- c(bh_haz, tail(bh_haz, 1))
  }

  eta_test <- as.numeric(as.matrix(test_df[, xnames, drop = FALSE]) %*% cox_fit$coefficients)

  q_low <- alpha / 2
  q_high <- 1 - alpha / 2

  lower <- vapply(eta_test, function(eta) {
    cox_tstar_quantile(bh_time = bh_time, bh_haz = bh_haz, eta = eta, q = q_low, tau = tau)
  }, numeric(1))

  upper <- vapply(eta_test, function(eta) {
    cox_tstar_quantile(bh_time = bh_time, bh_haz = bh_haz, eta = eta, q = q_high, tau = tau)
  }, numeric(1))

  out <- test_df %>%
    mutate(
      pi_lower = pmax(lower, 0),
      pi_upper = pmin(upper, tau)
    )

  list(
    intervals = out,
    cox_fit = cox_fit,
    bh_time = bh_time,
    bh_haz = bh_haz,
    tau = tau,
    alpha = alpha
  )
}

cox_plugin_pi_Tstar_tau_upper <- function(train_df,
                                          test_df,
                                          tau,
                                          alpha = 0.1) {
  xnames <- grep("^x\\d+$", names(train_df), value = TRUE)
  if (length(xnames) == 0) stop("No covariate columns named x1, x2, ... found.")

  fit_obj <- fit_cox_working(train_df, xnames)
  cox_fit <- fit_obj$fit
  bh <- fit_obj$bh %>% arrange(time)

  bh_tau <- bh %>% filter(time <= tau)
  if (nrow(bh_tau) == 0) {
    stop("No baseline hazard jumps at or before tau; cannot form Cox plug-in interval.")
  }

  bh_time <- bh_tau$time
  bh_haz <- bh_tau$hazard

  if (tail(bh_time, 1) < tau) {
    bh_time <- c(bh_time, tau)
    bh_haz <- c(bh_haz, tail(bh_haz, 1))
  }

  eta_test <- as.numeric(as.matrix(test_df[, xnames, drop = FALSE]) %*% cox_fit$coefficients)

  lower <- vapply(eta_test, function(eta) {
    cox_tstar_quantile(bh_time = bh_time, bh_haz = bh_haz, eta = eta, q = alpha, tau = tau)
  }, numeric(1))
  upper <- rep(tau, nrow(test_df))

  out <- test_df %>%
    mutate(
      pi_lower = pmax(lower, 0),
      pi_upper = upper
    )

  list(
    intervals = out,
    cox_fit = cox_fit,
    bh_time = bh_time,
    bh_haz = bh_haz,
    tau = tau,
    alpha = alpha
  )
}
