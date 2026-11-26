
#' Compute the negative log marginal likelihood under the gCER model
#' using parallel computation
#'
#' Evaluates the negative log marginal likelihood of the generalized
#' centered Erdős–Rényi (gCER) model for a given value of the scale
#' parameter. The computation integrates out the representative
#' graph and evaluates the edge-wise contributions in parallel
#' using RcppParallel.
#'
#' @param alpha Value of the scale parameter, with
#'   \eqn{0 < \alpha < 0.5}.
#' @param n Number of observed networks.
#' @param n_0 Prior effective sample size associated with the
#'   representative graph.
#' @param G_bar Numeric vector containing the edge-wise sample
#'   proportions across the observed networks.
#' @param W Numeric vector containing the prior location parameters
#'   for the representative graph, with entries in \eqn{[0,1]}.
#'   Its length must equal that of \code{G_bar}.
#'
#' @return The negative log marginal likelihood evaluated at
#'   \code{alpha}.
#'
#' @details
#' The marginal likelihood is computed by integrating out the
#' representative graph under the gCER prior. The calculation
#' uses a log-sum-exp representation for numerical stability.
#' Edge-wise contributions are evaluated in parallel using
#' \pkg{RcppParallel}.
#'
#' This function provides a parallel alternative to
#' \code{\link{marg_lik_CER_gCER}}.
#'
#' @references
#' Barile, F., Capozio, S., and Nipoti, B.
#' Efficient Bayesian inference for multiple network data.
#' Accepted for publication in Statistics and Probability Letters.
#'
#' @name marg_lik_CER_gCER_par
NULL
