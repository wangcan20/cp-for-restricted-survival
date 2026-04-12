# 00_utils.R
suppressPackageStartupMessages({
  library(survival)
  library(dplyr)
})

# Safe equality check for tau (because of floating point)
is_tau <- function(x, tau, tol = 1e-10) {
  abs(x - tau) <= tol
}

# Invert a nondecreasing step function given its knots and values.
# Input:
#   times: increasing time grid (length m)
#   haz  : cumulative hazard evaluated at those times (length m), nondecreasing
# Returns:
#   smallest time t such that H0(t) >= target; if target <= 0, returns 0.
inv_cumhaz_lower <- function(times, haz, target) {
  if (is.na(target) || target <= 0) return(0)
  idx <- which(haz >= target)
  if (length(idx) == 0) return(max(times))
  return(times[min(idx)])
}

# For upper bound we want largest t such that H0(t) <= target.
# If target >= max(haz), return max(times).
inv_cumhaz_upper <- function(times, haz, target) {
  if (is.na(target)) return(NA_real_)
  idx <- which(haz <= target)
  if (length(idx) == 0) return(0)
  return(times[max(idx)])
}

# Summaries
mean_na <- function(x) mean(x, na.rm = TRUE)
sd_na   <- function(x) sd(x, na.rm = TRUE)
