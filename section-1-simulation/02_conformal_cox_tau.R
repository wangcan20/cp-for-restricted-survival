# 02_conformal_cox_tau.R

# Build IPCW-style weighted support for (T*, Z) under finite tau
# We include:
#  - failures: delta=1 (event observed before tau)
#  - administrative mass at tau: y=tau & delta=0, reweighted by 1/G(tau-)
build_weighted_support_Tstar <- function(train_df, tau) {
  # Reverse KM for censoring survival. With status = (1 - delta):
  # - KM right-continuous curve estimates P(C > t)
  # - IPCW here needs G(t-) = P(C >= t), so evaluate left limits.
  fit_G <- survfit(Surv(y, 1 - delta) ~ 1, data = train_df)

  # Evaluate left limit Ghat(t-) from a survfit curve.
  # This is important at t = tau when C has a point mass at tau.
  eval_G_left <- function(tvec, tol = 1e-10) {
    times <- fit_G$time
    surv <- fit_G$surv

    if (length(times) == 0) return(rep(1, length(tvec)))

    surv_before <- c(1, surv[-length(surv)])
    out <- numeric(length(tvec))

    for (j in seq_along(tvec)) {
      t <- tvec[j]
      if (is.na(t)) {
        out[j] <- NA_real_
        next
      }

      idx <- findInterval(t, times, rightmost.closed = TRUE)
      if (idx <= 0) {
        out[j] <- 1
      } else if (abs(t - times[idx]) <= tol) {
        out[j] <- surv_before[idx]  # left limit at jump time
      } else {
        out[j] <- surv[idx]         # between jumps
      }
    }

    pmax(out, 1e-12)  # avoid division by 0
  }

  # Failure part
  df_fail <- train_df %>% filter(delta == 1)
  if (nrow(df_fail) == 0) stop("No observed failures in training set; cannot build IPCW support.")

  G_fail <- eval_G_left(df_fail$y)
  w_fail <- 1 / G_fail

  # Administrative mass at tau: y==tau and delta==0
  df_tau <- train_df %>% filter(delta == 0, is_tau(y, tau))
  if (nrow(df_tau) > 0) {
    G_tau <- eval_G_left(rep(tau, nrow(df_tau)))
    w_tau <- 1 / G_tau
  } else {
    w_tau <- numeric(0)
  }

  support_df <- bind_rows(
    df_fail %>% mutate(T_star_support = y, w = w_fail),
    df_tau  %>% mutate(T_star_support = tau, w = w_tau)
  )

  # Normalize weights
  support_df$w <- support_df$w / sum(support_df$w)

  list(
    support_df = support_df,
    fit_G = fit_G
  )
}

# Fit Cox and extract cumulative baseline hazard up to tau
fit_cox_working <- function(train_df, xnames) {
  fmla <- as.formula(paste("Surv(y, delta) ~ ", paste(xnames, collapse = "+")))
  fit <- coxph(fmla, data = train_df, x = TRUE)

  # basehaz gives cumulative baseline hazard H0(t)
  bh <- basehaz(fit, centered = FALSE)
  # ensure increasing times
  bh <- bh[order(bh$time), ]
  list(fit = fit, bh = bh)
}

# Compute pivot U = exp(-H0(t) * exp(eta))
pivot_from_bh <- function(bh_time, bh_haz, t, eta) {
  # stepwise cumhaz at t
  # use last hazard with time <= t
  idx <- findInterval(t, bh_time)
  H0 <- numeric(length(idx))
  pos <- which(idx > 0)
  if (length(pos) > 0) {
    H0[pos] <- bh_haz[idx[pos]]
  }
  exp(-H0 * exp(eta))
}

# Estimate fixed-horizon marginal risk score:
#   pi_marg(t0, z0) = P(U >= u0), where
#   U = S(T* | Z) from Monte Carlo draws under Fhat and
#   u0 = S(t0 | z0).
# fit_obj should be the returned object from conformal_pi_cox_Tstar().
predict_pi_marginal_from_fit <- function(fit_obj, z_new, t0_vec) {
  if (is.null(fit_obj$U_b)) {
    stop("fit_obj$U_b is missing. Refit with conformal_pi_cox_Tstar() that returns U_b.")
  }
  if (is.null(fit_obj$cox_fit) || is.null(fit_obj$bh_time) || is.null(fit_obj$bh_haz)) {
    stop("fit_obj does not contain required Cox/baseline hazard objects.")
  }

  beta_hat <- fit_obj$cox_fit$coefficients
  p <- length(beta_hat)

  z_mat <- as.matrix(z_new)
  if (is.vector(z_new)) {
    z_mat <- matrix(as.numeric(z_new), nrow = 1)
  }
  if (ncol(z_mat) != p) {
    stop(sprintf("z_new has %d columns but beta has length %d.", ncol(z_mat), p))
  }

  t0_vec <- as.numeric(t0_vec)
  if (length(t0_vec) == 0) stop("t0_vec must be non-empty.")

  if (nrow(z_mat) == 1 && length(t0_vec) > 1) {
    z_mat <- z_mat[rep(1, length(t0_vec)), , drop = FALSE]
  } else if (length(t0_vec) == 1 && nrow(z_mat) > 1) {
    t0_vec <- rep(t0_vec, nrow(z_mat))
  } else if (nrow(z_mat) != length(t0_vec)) {
    stop("Provide either one z_new with many t0 values, one t0 with many z_new rows, or equal lengths.")
  }

  if (!is.null(fit_obj$tau)) {
    t0_vec <- pmin(pmax(t0_vec, 0), fit_obj$tau)
  }

  eta_hat <- as.numeric(z_mat %*% beta_hat)
  u0 <- pivot_from_bh(
    bh_time = fit_obj$bh_time,
    bh_haz = fit_obj$bh_haz,
    t = t0_vec,
    eta = eta_hat
  )

  pi_marg <- vapply(u0, function(u) mean(fit_obj$U_b >= u), numeric(1))

  data.frame(
    t0 = t0_vec,
    eta_hat = eta_hat,
    u0 = u0,
    pi_marg = pi_marg
  )
}

# Compute one global upper-tail cutoff in pivot space for the upper-tau interval.
compute_upper_tau_u_cut <- function(U_b, alpha) {
  U_b <- as.numeric(U_b)
  if (!length(U_b)) stop("U_b must be non-empty.")
  as.numeric(quantile(U_b, probs = 1 - alpha, na.rm = TRUE, names = FALSE))
}

# Main: conformal PI for T* with upper endpoint fixed at tau
conformal_pi_cox_Tstar <- function(train_df,
                                  test_df,
                                  tau,
                                  alpha = 0.1,
                                  B = 1000,
                                  seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  xnames <- grep("^x\\d+$", names(train_df), value = TRUE)
  if (length(xnames) == 0) stop("No covariate columns named x1, x2, ... found.")

  # 1) Fit Cox
  fit_obj <- fit_cox_working(train_df, xnames)
  cox_fit <- fit_obj$fit
  bh <- fit_obj$bh

  # Restrict baseline hazard grid to <= tau (still fine if last time < tau)
  bh_tau <- bh %>% filter(time <= tau)
  if (nrow(bh_tau) == 0) {
    # if no basehaz time <= tau (rare), keep smallest time
    bh_tau <- bh[1, , drop = FALSE]
  }
  bh_time <- bh_tau$time
  bh_haz  <- bh_tau$hazard

  # 2) Build weighted support for Monte Carlo sampling of (T*, Z)
  sup <- build_weighted_support_Tstar(train_df, tau)
  support_df <- sup$support_df

  # 3) Monte Carlo pivots
  # sample rows from support with prob w
  idx_samp <- sample.int(nrow(support_df), size = B, replace = TRUE, prob = support_df$w)

  T_star_b <- support_df$T_star_support[idx_samp]
  Z_b <- as.matrix(support_df[idx_samp, xnames, drop = FALSE])

  eta_b <- as.numeric(Z_b %*% cox_fit$coefficients)
  U_b <- pivot_from_bh(bh_time, bh_haz, t = T_star_b, eta = eta_b)

  if (tail(bh_time, 1) < tau) {
    bh_time <- c(bh_time, tau)
    bh_haz <- c(bh_haz, tail(bh_haz, 1))
  }

  # 4) One global cutoff in pivot space for [lower(x), tau]
  u_cut <- compute_upper_tau_u_cut(U_b = U_b, alpha = alpha)

  # 5) Predict intervals on test set
  Z_test <- as.matrix(test_df[, xnames, drop = FALSE])
  eta_test <- as.numeric(Z_test %*% cox_fit$coefficients)

  # Solve S(t | x) = u_cut, equivalently F_{T*|x}(t) = 1 - u_cut.
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

  lower <- pmax(lower, 0)

  out <- test_df %>%
    mutate(
      pi_lower = lower,
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
    tau = tau
  )
}
