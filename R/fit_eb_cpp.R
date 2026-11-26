
#' Fit the empirical Bayes model
#'
#' Estimates the empirical Bayes parameter and the center of a
#' generalized centered Erdős–Rényi model from a collection of
#' observed edge counts.
#'
#' The parameter is estimated by minimizing the negative log
#' marginal likelihood using golden-section optimization.
#' The center is obtained through a majority-rule criterion.
#'
#' @param counts Integer vector containing the observed edge counts
#'   across the networks or observations.
#' @param n Number of possible edges in each network, corresponding
#'   to the number of trials for each edge count.
#'
#' @return A list containing:
#' \describe{
#'   \item{alpha_hat}{Estimated empirical Bayes parameter,
#'   constrained to the interval (0, 0.5).}
#'   \item{C_hat}{Integer vector representing the estimated
#'   majority-rule center. An entry is equal to 1 if the
#'   corresponding edge count is at least n/2, and 0 otherwise.}
#'   \item{objective}{Value of the negative log marginal likelihood
#'   evaluated at the estimated parameter.}
#' }
#'
#' @details
#' The empirical Bayes estimate is obtained by minimizing the
#' negative log marginal likelihood over (0, 0.5), using
#' golden-section optimization. The optimization tolerance is
#' 1e-11.
#'
#' The majority-rule center is computed as
#' \deqn{\widehat C_j = \mathbb{I}(2\,\text{counts}_j \geq n).}
#'
#' @name fit_eb_cpp
NULL
