# 00_split_support_fit_methods.R
#
# Exploratory variants of the main CP + Cox upper-tau algorithm that split
# the training sample into separate halves for:
#   (a) support / joint CDF estimation, and
#   (b) Cox working-model fitting.

suppressPackageStartupMessages({
  library(dplyr)
})

split_train_support_fit <- function(train_df, split_frac = 0.5, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  n <- nrow(train_df)
  if (n < 4) stop("train_df is too small to split.")

  n_support <- floor(split_frac * n)
  n_support <- max(2L, min(n_support, n - 2L))

  idx <- sample.int(n)
  support_id <- idx[seq_len(n_support)]
  fit_id <- idx[-seq_len(n_support)]

  list(
    support_df = train_df[support_id, , drop = FALSE],
    fit_df = train_df[fit_id, , drop = FALSE]
  )
}

conformal_pi_cox_Tstar_split_support_fit <- function(train_df,
                                                     test_df,
                                                     tau,
                                                     alpha = 0.1,
                                                     B = 1000,
                                                     seed = NULL,
                                                     split_frac = 0.5,
                                                     split_seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  xnames <- grep("^x\\d+$", names(train_df), value = TRUE)
  if (length(xnames) == 0) stop("No covariate columns named x1, x2, ... found.")

  split_obj <- split_train_support_fit(
    train_df = train_df,
    split_frac = split_frac,
    seed = split_seed
  )
  support_df <- split_obj$support_df
  fit_df <- split_obj$fit_df

  fit_obj <- fit_cox_working(fit_df, xnames)
  cox_fit <- fit_obj$fit
  bh <- fit_obj$bh

  bh_tau <- bh %>% filter(time <= tau)
  if (nrow(bh_tau) == 0) {
    bh_tau <- bh[1, , drop = FALSE]
  }
  bh_time <- bh_tau$time
  bh_haz <- bh_tau$hazard

  sup <- build_weighted_support_Tstar(support_df, tau)
  support_weighted_df <- sup$support_df

  idx_samp <- sample.int(
    nrow(support_weighted_df),
    size = B,
    replace = TRUE,
    prob = support_weighted_df$w
  )

  T_star_b <- support_weighted_df$T_star_support[idx_samp]
  Z_b <- as.matrix(support_weighted_df[idx_samp, xnames, drop = FALSE])

  eta_b <- as.numeric(Z_b %*% cox_fit$coefficients)
  U_b <- pivot_from_bh(bh_time, bh_haz, t = T_star_b, eta = eta_b)

  if (tail(bh_time, 1) < tau) {
    bh_time <- c(bh_time, tau)
    bh_haz <- c(bh_haz, tail(bh_haz, 1))
  }

  u_cut <- compute_upper_tau_u_cut(U_b = U_b, alpha = alpha)

  Z_test <- as.matrix(test_df[, xnames, drop = FALSE])
  eta_test <- as.numeric(Z_test %*% cox_fit$coefficients)

  q_lower <- pmin(pmax(1 - u_cut, 0), 1)
  lower <- vapply(eta_test, function(eta) {
    H_tau <- bh_haz[length(bh_haz)]
    F_tau <- 1 - exp(-H_tau * exp(eta))
    if (q_lower <= 0) return(0)
    if (q_lower >= F_tau) return(tau)

    target <- -log1p(-q_lower) / exp(eta)
    if (!is.finite(target)) return(tau)
    inv_cumhaz_lower(bh_time, bh_haz, target)
  }, numeric(1))
  upper <- rep(tau, nrow(test_df))

  out <- test_df %>%
    mutate(
      pi_lower = pmax(lower, 0),
      pi_upper = upper
    )

  list(
    intervals = out,
    u_cut = u_cut,
    cox_fit = cox_fit,
    bh_time = bh_time,
    bh_haz = bh_haz,
    U_b = U_b,
    T_star_b = T_star_b,
    tau = tau,
    support_n = nrow(support_df),
    fit_n = nrow(fit_df),
    support_fail_n = sum(support_df$delta == 1),
    support_tau_n = sum(support_df$delta == 0 & is_tau(support_df$y, tau))
  )
}
