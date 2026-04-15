# 00_crossfit_support_fit_methods.R
#
# Cross-fitted variant of the main CP + Cox upper-tau algorithm.
# Training data are partitioned into K folds. For each fold k:
#   - fold k estimates the weighted support / joint distribution of (T*, Z)
#   - the complement fits the Cox working model
# Pooled fold-specific pivots define one global cutoff, while the final
# inversion is done with a Cox fit on the full training set.

suppressPackageStartupMessages({
  library(dplyr)
})

make_crossfit_folds <- function(train_df, n_folds = 5L, seed = NULL, max_tries = 50L) {
  n <- nrow(train_df)
  if (n_folds < 2L) stop("n_folds must be at least 2.")
  if (n < 2L * n_folds) stop("train_df is too small for the requested number of folds.")

  for (attempt in seq_len(max_tries)) {
    if (!is.null(seed)) set.seed(seed + attempt - 1L)

    idx <- sample.int(n)
    fold_id <- rep(seq_len(n_folds), length.out = n)
    fold_id[idx] <- fold_id

    fail_counts <- vapply(seq_len(n_folds), function(k) {
      sum(train_df$delta[fold_id == k] == 1)
    }, numeric(1))

    if (all(fail_counts > 0)) {
      return(fold_id)
    }
  }

  stop("Failed to create cross-fit folds with at least one observed failure in every support fold.")
}

allocate_mc_draws <- function(B, n_folds) {
  base <- rep(B %/% n_folds, n_folds)
  rem <- B %% n_folds
  if (rem > 0) {
    base[seq_len(rem)] <- base[seq_len(rem)] + 1L
  }
  base
}

conformal_pi_cox_Tstar_crossfit <- function(train_df,
                                            test_df,
                                            tau,
                                            alpha = 0.1,
                                            B = 1000,
                                            seed = NULL,
                                            n_folds = 5L,
                                            fold_seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  xnames <- grep("^x\\d+$", names(train_df), value = TRUE)
  if (length(xnames) == 0) stop("No covariate columns named x1, x2, ... found.")

  fold_id <- make_crossfit_folds(
    train_df = train_df,
    n_folds = n_folds,
    seed = fold_seed
  )
  B_alloc <- allocate_mc_draws(B = B, n_folds = n_folds)

  U_parts <- vector("list", n_folds)
  T_parts <- vector("list", n_folds)

  for (k in seq_len(n_folds)) {
    support_df <- train_df[fold_id == k, , drop = FALSE]
    fit_df <- train_df[fold_id != k, , drop = FALSE]

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

    B_k <- B_alloc[k]
    if (B_k <= 0) next

    idx_samp <- sample.int(
      nrow(support_weighted_df),
      size = B_k,
      replace = TRUE,
      prob = support_weighted_df$w
    )

    T_star_b <- support_weighted_df$T_star_support[idx_samp]
    Z_b <- as.matrix(support_weighted_df[idx_samp, xnames, drop = FALSE])

    eta_b <- as.numeric(Z_b %*% cox_fit$coefficients)
    U_parts[[k]] <- pivot_from_bh(bh_time, bh_haz, t = T_star_b, eta = eta_b)
    T_parts[[k]] <- T_star_b
  }

  U_b <- unlist(U_parts, use.names = FALSE)
  T_star_b <- unlist(T_parts, use.names = FALSE)
  if (!length(U_b)) stop("Cross-fitted pivot sample is empty.")

  u_cut <- compute_upper_tau_u_cut(U_b = U_b, alpha = alpha)

  full_fit_obj <- fit_cox_working(train_df, xnames)
  full_cox_fit <- full_fit_obj$fit
  full_bh <- full_fit_obj$bh

  full_bh_tau <- full_bh %>% filter(time <= tau)
  if (nrow(full_bh_tau) == 0) {
    full_bh_tau <- full_bh[1, , drop = FALSE]
  }

  bh_time <- full_bh_tau$time
  bh_haz <- full_bh_tau$hazard

  if (tail(bh_time, 1) < tau) {
    bh_time <- c(bh_time, tau)
    bh_haz <- c(bh_haz, tail(bh_haz, 1))
  }

  Z_test <- as.matrix(test_df[, xnames, drop = FALSE])
  eta_test <- as.numeric(Z_test %*% full_cox_fit$coefficients)

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

  out <- test_df %>%
    mutate(
      pi_lower = pmax(lower, 0),
      pi_upper = tau
    )

  list(
    intervals = out,
    u_cut = u_cut,
    cox_fit = full_cox_fit,
    bh_time = bh_time,
    bh_haz = bh_haz,
    U_b = U_b,
    T_star_b = T_star_b,
    tau = tau,
    n_folds = n_folds,
    fold_sizes = as.integer(table(fold_id))
  )
}
