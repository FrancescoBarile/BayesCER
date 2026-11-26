
#' Simulation results for the generalized CER model
#'
#' A dataset containing the results of a simulation study designed
#' to evaluate and compare empirical Bayes (EB) estimation and
#' Gibbs sampling under the generalized centered Erdős–Rényi
#' (gCER) model.
#'
#' The simulation study considers different numbers of observed
#' networks, network sizes, and values of the model parameter
#' \eqn{\alpha}. For each scenario, 100 independent replications
#' are performed. The data-generating process uses a stochastic
#' block model to generate the center network.
#'
#' @format A data frame in long format, with one row for each
#' combination of simulation scenario, replication, and inference
#' method. It contains the following variables:
#' \describe{
#'   \item{n}{Number of observed networks.}
#'   \item{N}{Number of nodes in each network.}
#'   \item{alpha}{True value of the gCER model parameter.}
#'   \item{rep}{Simulation replication index.}
#'   \item{ecc_true}{True expected global clustering coefficient
#'   under the data-generating model.}
#'   \item{density_true}{True expected network density under the
#'   data-generating model.}
#'   \item{starting_time}{Starting time of the inference procedure.}
#'   \item{ending_time}{Ending time of the inference procedure.}
#'   \item{starting_time_rgen}{Starting time of the random network
#'   generation procedure.}
#'   \item{ending_time_rgen}{Ending time of the random network
#'   generation procedure.}
#'   \item{discrepancy_C_hat}{Normalized Hamming distance between
#'   the estimated and true center networks.}
#'   \item{discrepancy_alpha_hat}{Absolute error in estimating
#'   the parameter \eqn{\alpha}.}
#'   \item{ecc_hat}{Estimated expected global clustering coefficient.}
#'   \item{density}{Estimated expected network density.}
#'   \item{method}{Inference method: \code{"EB"} or \code{"Gibbs"}.}
#' }
#'
#' @details
#' The simulation considers network sizes
#' \eqn{N \in \{10,100,1000,2000\}}, sample sizes
#' \eqn{n \in \{10,100,1000\}}, and parameter values
#' \eqn{\alpha \in \{0.1,0.2,0.3\}}.
#' Each scenario is replicated 100 times.
#'
#' @source Simulation study conducted using the \pkg{BayesCER}
#' package.
#'
#' @examples
#' data(sim_data_result)
#' head(sim_data_result)
"sim_data_result"

