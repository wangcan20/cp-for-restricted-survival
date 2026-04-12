# 01_dgp_weibull_ph.R

# Draw T from Weibull PH model:
# hazard(t|z) = (gamma / kappa) (t/kappa)^(gamma-1) * exp(z^T beta)
# Equivalent sampling:
# T = kappa * { -log(U) / exp(z^T beta) }^(1/gamma)
rweibull_ph <- function(n, Z, beta, gamma = 1.2, kappa = 1.0) {
  linpred <- as.numeric(Z %*% beta)
  U <- runif(n)
  T <- kappa * ((-log(U)) / exp(linpred))^(1 / gamma)
  T
}

# Calibrate exponential censoring rate to achieve an approximate censoring proportion.
# We generate latent censoring C_latent ~ Exp(rate), then define observed
# censoring C = min(C_latent, tau), which induces a point mass at tau.
# We target: P(Delta=0 and Y < tau) ~ cens_target (roughly).
calibrate_censor_rate <- function(T, tau, cens_target = 0.2, max_rate = 10, tol = 1e-3) {
  # Objective: censoring proportion among those with T < tau (rough heuristic)
  obj <- function(rate) {
    C_latent <- rexp(length(T), rate = rate)
    C <- pmin(C_latent, tau)
    Y <- pmin(T, C)
    Delta <- as.integer(T <= C)
    # "non-admin" censoring: censored before tau
    cens_before_tau <- mean(Delta == 0 & Y < tau)
    cens_before_tau - cens_target
  }
  # If already below target even with huge rate, just return huge rate
  if (obj(max_rate) < 0) return(max_rate)
  # If above target even with tiny rate, return tiny rate
  if (obj(1e-6) > 0) return(1e-6)

  uniroot(obj, lower = 1e-6, upper = max_rate, tol = tol)$root
}

# Simulate observed (Y, Delta, Z) with finite horizon tau
simulate_rcensored_finite_tau <- function(n,
                                         p = 3,
                                         beta = rep(0.3, 3),
                                         gamma = 1.2,
                                         kappa = 1.0,
                                         tau_quantile = 0.9,
                                         cens_target = 0.2,
                                         seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  Z <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(Z) <- paste0("x", 1:p)

  T <- rweibull_ph(n, Z, beta = beta, gamma = gamma, kappa = kappa)

  # Choose a fixed tau from this sample's marginal T distribution
  tau <- as.numeric(quantile(T, probs = tau_quantile))

  # Calibrate censoring rate (roughly) and draw censoring times
  rate_c <- calibrate_censor_rate(T, tau, cens_target = cens_target)
  C_latent <- rexp(n, rate = rate_c)
  C <- pmin(C_latent, tau)  # mixed censoring variable with mass at tau

  # Observed:
  # Y = min(T, C)
  # Delta = I(T <= C)
  Y <- pmin(T, C)
  Delta <- as.integer(T <= C)

  # True target for evaluation: T* = min(T, tau)
  T_star <- pmin(T, tau)

  dat <- data.frame(
    y = Y,
    delta = Delta,
    C = C,
    T_true = T,
    T_star = T_star,
    tau = tau
  ) %>%
    bind_cols(as.data.frame(Z))

  list(dat = dat, tau = tau, rate_c = rate_c)
}
