# 00_survreg_working_models.R
#
# Parametric working-model helpers for T* = min(T, tau).
# These functions are designed to be sourced after:
#   - section-1-simulation/00_utils.R
#   - section-1-simulation/02_conformal_cox_tau.R

suppressPackageStartupMessages({
  library(survival)
  library(dplyr)
})

validate_survreg_dist <- function(dist) {
  allowed <- c("weibull", "lognormal", "loglogistic", "exponential")
  if (!dist %in% allowed) {
    stop("Unsupported survreg distribution: ", dist,
         ". Allowed values: ", paste(allowed, collapse = ", "))
  }
  invisible(dist)
}

fit_survreg_working <- function(train_df, xnames, dist = "weibull") {
  validate_survreg_dist(dist)

  fmla <- as.formula(paste("Surv(y, delta) ~ ", paste(xnames, collapse = "+")))
  fit <- survreg(fmla, data = train_df, dist = dist, x = TRUE, y = TRUE)

  list(
    fit = fit,
    dist = dist,
    xnames = xnames
  )
}

prepare_survreg_newdata <- function(fit_obj, z_new) {
  xnames <- fit_obj$xnames

  if (is.data.frame(z_new)) {
    if (!all(xnames %in% names(z_new))) {
      stop("z_new data frame is missing required covariate columns.")
    }
    nd <- z_new[, xnames, drop = FALSE]
  } else {
    z_mat <- as.matrix(z_new)
    if (is.vector(z_new)) {
      z_mat <- matrix(as.numeric(z_new), nrow = 1)
    }
    if (ncol(z_mat) != length(xnames)) {
      stop(sprintf("z_new has %d columns but expected %d.", ncol(z_mat), length(xnames)))
    }
    nd <- as.data.frame(z_mat)
    names(nd) <- xnames
  }

  nd
}

predict_survreg_lp <- function(fit_obj, z_new) {
  nd <- prepare_survreg_newdata(fit_obj, z_new)
  as.numeric(predict(fit_obj$fit, newdata = nd, type = "lp"))
}

survreg_survival <- function(fit_obj, z_new, t) {
  nd <- prepare_survreg_newdata(fit_obj, z_new)
  lp <- predict_survreg_lp(fit_obj, nd)
  t <- as.numeric(t)

  if (length(t) == 1L) {
    t <- rep(t, nrow(nd))
  }
  if (length(t) != nrow(nd)) {
    stop("Provide either one t value or one t value per row of z_new.")
  }

  out <- numeric(length(t))
  pos <- which(t > 0)

  out[t <= 0] <- 1
  if (length(pos) > 0) {
    out[pos] <- 1 - psurvreg(
      q = t[pos],
      mean = lp[pos],
      scale = fit_obj$fit$scale,
      distribution = fit_obj$dist
    )
  }

  pmax(pmin(out, 1), 0)
}

# Explicit mixed-distribution CDF for T* = min(T, tau).
# For t < tau this equals F_T(t | z), and at t = tau it jumps to 1.
survreg_tstar_cdf <- function(fit_obj, z_new, t, tau) {
  nd <- prepare_survreg_newdata(fit_obj, z_new)
  lp <- predict_survreg_lp(fit_obj, nd)
  t <- as.numeric(t)

  if (length(t) == 1L) {
    t <- rep(t, nrow(nd))
  }
  if (length(t) != nrow(nd)) {
    stop("Provide either one t value or one t value per row of z_new.")
  }

  out <- numeric(length(t))
  out[t <= 0] <- 0

  mid <- which(t > 0 & t < tau)
  if (length(mid) > 0) {
    out[mid] <- psurvreg(
      q = t[mid],
      mean = lp[mid],
      scale = fit_obj$fit$scale,
      distribution = fit_obj$dist
    )
  }

  out[t >= tau] <- 1
  pmax(pmin(out, 1), 0)
}

survreg_tstar_quantile <- function(fit_obj, z_new, q, tau) {
  nd <- prepare_survreg_newdata(fit_obj, z_new)
  lp <- predict_survreg_lp(fit_obj, nd)
  q <- as.numeric(q)

  if (length(q) == 1L) {
    q <- rep(q, nrow(nd))
  }
  if (length(q) != nrow(nd)) {
    stop("Provide either one q value or one q value per row of z_new.")
  }

  q <- pmin(pmax(q, 0), 1)
  out <- numeric(length(q))

  out[q <= 0] <- 0

  # The target distribution is T* = min(T, tau), not T. Its quantile
  # function is piecewise:
  #   Q_{T*}(q | z) = Q_T(q | z), if q < F_T(tau | z),
  #                 = tau,        otherwise.
  # This is algebraically equivalent to min(Q_T(q | z), tau), but the
  # explicit mixed-distribution form makes the target of the plug-in
  # procedure transparent.
  F_tau <- psurvreg(
    q = rep(tau, nrow(nd)),
    mean = lp,
    scale = fit_obj$fit$scale,
    distribution = fit_obj$dist
  )

  mid <- which(q > 0 & q < 1 & q < F_tau)
  if (length(mid) > 0) {
    out[mid] <- qsurvreg(
      p = q[mid],
      mean = lp[mid],
      scale = fit_obj$fit$scale,
      distribution = fit_obj$dist
    )
  }

  out[q >= F_tau] <- tau
  pmax(0, pmin(out, tau))
}

survreg_plugin_pi_Tstar <- function(train_df,
                                    test_df,
                                    tau,
                                    alpha = 0.1,
                                    dist = "weibull") {
  xnames <- grep("^x\\d+$", names(train_df), value = TRUE)
  if (length(xnames) == 0) stop("No covariate columns named x1, x2, ... found.")

  fit_obj <- fit_survreg_working(train_df, xnames = xnames, dist = dist)

  lower <- survreg_tstar_quantile(
    fit_obj = fit_obj,
    z_new = test_df[, xnames, drop = FALSE],
    q = alpha / 2,
    tau = tau
  )
  upper <- survreg_tstar_quantile(
    fit_obj = fit_obj,
    z_new = test_df[, xnames, drop = FALSE],
    q = 1 - alpha / 2,
    tau = tau
  )

  out <- test_df %>%
    mutate(
      pi_lower = lower,
      pi_upper = upper
    )

  list(
    intervals = out,
    fit_obj = fit_obj,
    tau = tau,
    alpha = alpha
  )
}

survreg_plugin_pi_Tstar_tau_upper <- function(train_df,
                                              test_df,
                                              tau,
                                              alpha = 0.1,
                                              dist = "weibull") {
  xnames <- grep("^x\\d+$", names(train_df), value = TRUE)
  if (length(xnames) == 0) stop("No covariate columns named x1, x2, ... found.")

  fit_obj <- fit_survreg_working(train_df, xnames = xnames, dist = dist)

  lower <- survreg_tstar_quantile(
    fit_obj = fit_obj,
    z_new = test_df[, xnames, drop = FALSE],
    q = alpha,
    tau = tau
  )
  upper <- rep(tau, nrow(test_df))

  out <- test_df %>%
    mutate(
      pi_lower = lower,
      pi_upper = upper
    )

  list(
    intervals = out,
    fit_obj = fit_obj,
    tau = tau,
    alpha = alpha
  )
}

conformal_pi_survreg_Tstar <- function(train_df,
                                       test_df,
                                       tau,
                                       alpha = 0.1,
                                       B = 1000,
                                       seed = NULL,
                                       dist = "weibull") {
  if (!is.null(seed)) set.seed(seed)

  xnames <- grep("^x\\d+$", names(train_df), value = TRUE)
  if (length(xnames) == 0) stop("No covariate columns named x1, x2, ... found.")

  fit_obj <- fit_survreg_working(train_df, xnames = xnames, dist = dist)

  sup <- build_weighted_support_Tstar(train_df, tau)
  support_df <- sup$support_df

  idx_samp <- sample.int(nrow(support_df), size = B, replace = TRUE, prob = support_df$w)
  T_star_b <- support_df$T_star_support[idx_samp]
  Z_b <- support_df[idx_samp, xnames, drop = FALSE]

  U_b <- survreg_survival(
    fit_obj = fit_obj,
    z_new = Z_b,
    t = T_star_b
  )

  u_cut <- compute_upper_tau_u_cut(U_b = U_b, alpha = alpha)

  lower <- survreg_tstar_quantile(
    fit_obj = fit_obj,
    z_new = test_df[, xnames, drop = FALSE],
    q = 1 - u_cut,
    tau = tau
  )
  upper <- rep(tau, nrow(test_df))

  out <- test_df %>%
    mutate(
      pi_lower = lower,
      pi_upper = upper
    )

  list(
    intervals = out,
    u_cut = u_cut,
    fit_obj = fit_obj,
    U_b = U_b,
    T_star_b = T_star_b,
    tau = tau,
    alpha = alpha,
    dist = dist
  )
}
