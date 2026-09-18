#' cast a M-dimensional numeric binary vector encoding the lower triangular part of the adjacency matrix associated to a graph into a N \times N matrix.
#' @param low_adj [numeric] M-dimensional numeric binary vector encoding the lower triangular part of the adjacency matrix associated to a graph.
#'
#' @return [matrix] N x N binary matrix encoding the adjacency matrix associated to a graph.
#' @export
from_low_adj_to_full_adj<-function(low_adj){
  M<-length(low_adj)
  N<-(1+sqrt(1+8*M))/2
  m1 <- matrix(NA, N, N)
  m1<-`diag<-`(m1, 0)
  m1[lower.tri(m1, diag=FALSE)] <- low_adj
  m2 <- t(m1)
  m2[lower.tri(m2, diag=FALSE)] <- low_adj
  return(m2)
}

#' cast a N \times N matrix encoding the adjacency matrix associated to a graph.
#' into a M-dimensional numeric binary vector encoding the lower triangular part of the adjacency matrix associated to the graph.
#' @param full_adj [matrix] N \times N matrix encoding the adjacency matrix associated to a graph.
#'
#' @return [numeric] M-dimensional numeric binary vector encoding the lower triangular part of the adjacency matrix associated to a graph.
#' @export

from_full_adj_to_low_adj<-function(full_adj){
  low_adj<-t(full_adj)[lower.tri(t(full_adj))]
  return(low_adj)
}

#' generate a CER-distributed graph.
#' @param C [numeric] M-dimensional numeric binary vector encoding the lower triangular part of the adjacency matrix associated to the representative graph C.
#' @param alpha [numeric] probability an edge differs from that of the representative graph C.
#'
#' @return [numeric] M-dimensional numeric binary vector encoding the lower triangular part of the adjacency matrix associated to generated CER-distributed graph.

G_from_CER <- function(C, alpha) {
  etas<-((C==1)-(C==0))
  p_bin<-(1+(alpha/(1-alpha))^(etas))^(-1)
  G<-rbinom(n=length(C), size=1, prob = p_bin)
  return(G)
}

#' generate a set of n i.i.d. CER-distributed graphs.
#' @param C [numeric] M-dimensional numeric binary vector encoding the lower triangular part of the adjacency matrix associated to the representative graph C.
#' @param alpha [numeric] probability in \eqn{(0,1/2)} that an edge differs from that of the representative graph C.
#' @param n [numeric] number of graphs to generate.
#'
#' @return [matrix] n x M binary matrix with each row encoding the lower triangular part of the adjacency matrix associated to the generated CER-distributed graph.

G_n_from_CER<-function(C, alpha, n) {
  G_group<-t( sapply(1:n, function(i) G_from_CER(C, alpha)) )
  return(G_group)
}

#' generate a SBM-distributed graph.
#' @param N [numeric] number of nodes.
#' @param K [numeric] Number of blocks.
#' @param intra_prob [numeric] Probability of edges within blocks.
#' @param inter_prob [numeric] Probability of edges between blocks.
#'
#' @return [igraph] An igraph object representing a stochastic block model graph.
#' @importFrom igraph sample_sbm
generate_sbm <- function(N, K = 2, intra_prob = 0.8, inter_prob = 0.1) {
  split_blocks <- function(N, K) {
    base <- floor(N / K)
    remainder <- N %% K
    block_sizes <- rep(base, K)
    if (remainder > 0) {
      block_sizes[1:remainder] <- block_sizes[1:remainder] + 1
    }
    return(block_sizes)
  }
  block_sizes <- split_blocks(N, K)
  pref.matrix <- matrix(inter_prob, nrow = K, ncol = K)
  diag(pref.matrix) <- intra_prob
  sample_sbm(N, pref.matrix = pref.matrix, block.sizes = block_sizes,
             directed = FALSE, loops = FALSE)
}

#' generate a SBM-distributed graph.
#' @param N [numeric] number of nodes.
#' @param K [numeric] Number of blocks.
#' @param intra_prob [numeric] Probability of edges within blocks.
#' @param inter_prob [numeric] Probability of edges between blocks.
#'
#' @return [matrix] N \times N binary matrix encoding the adjacency matrix associated to the generated SBM-distributed graph.
#' @export
generate_sbm_fast <- function(N, K = 2, intra_prob = 0.8, inter_prob = 0.1,
                              z = sample(seq_len(K), N, TRUE)) {

  # Create SBM probability matrix
  B <- matrix(inter_prob, nrow = K, ncol = K)
  diag(B) <- intra_prob

  # Node-specific probability matrix
  P <- B[z, z]

  # Generate adjacency matrix
  U <- matrix(runif(N * N), N, N)

  A <- (U < P) * 1

  # Make symmetric
  A[lower.tri(A)] <- t(A)[lower.tri(A)]

  # Remove self-loops
  diag(A) <- 0

  return(
    # list(
    A = A
    # , z = z, B = B
  # )
  )
}


#' Compute the tau transformation of alpha.
#' \deqn{\tau(\alpha,n_0) = \frac{\alpha^{n_0}}
#' {(1-\alpha)^{n_0}+\alpha^{n_0}}}
#' which maps a probability parameter \eqn{\alpha} into the interval
#' \eqn{(0,1/2)}.
#'
#' @param alpha [numeric] A probability parameter in the interval \eqn{(0,1/2)}.
#' @param n_0 [numeric] A positive exponent controlling the concentration
#' of the transformation and interpreted as prior sample size.
#'
#' @return [numeric] The value of \eqn{\tau(\alpha,n_0)}.
#' @export
tau_alpha <- function(alpha, n_0){
  out <- ( (alpha)^n_0 ) / ( (1 - alpha)^n_0 + (alpha)^n_0 )
  return(out)
}

#' Generate random samples from a truncated Beta distribution.
#'
#' This function generates random draws from a Beta distribution with shape
#' parameters \eqn{\alpha} and \eqn{\beta}, truncated to the interval
#' \eqn{(a,b)}.
#' @param n [integer] Number of random samples to generate.
#'
#' @param alpha [numeric] First shape parameter of the TBeta distribution.
#' Must be positive.
#'
#' @param beta [numeric] Second shape parameter of the TBeta distribution.
#' Must be positive.
#'
#' @param a [numeric] Lower truncation bound. Default is \eqn{0}.
#'
#' @param b [numeric] Upper truncation bound. Default is \eqn{1/2}.
#'
#' @return [numeric] A vector of length \code{n} containing independent samples
#' from the truncated Beta distribution on the interval \eqn{(a,b)}.
#' @export
rtbeta <- function(n, alpha, beta, a = 0, b = 1){
  # stopifnot(n > 0 & all(beta > 0) & all(alpha > 0))
  x  <- stats::runif(n)
  Fa <- stats::pbeta(a, alpha, beta)
  Fb <- stats::pbeta(b, alpha, beta)
  y  <- (1 - x) * Fa + x * Fb
  return(stats::qbeta(y, alpha, beta))
}


#' Expected number of triangles incident to a node under a CER model.
#'
#' Computes the expected number of triangles involving node \eqn{i} in a
#' random graph generated from a centered Erdős--Rényi (CER) model based on
#' a representative graph \eqn{C}. The expectation depends on the node degree,
#' the total number of edges in the reference graph, and the number of
#' triangles incident to node \eqn{i} in the reference graph.
#'
#' The expected triangle count is computed as
#' \deqn{
#' E[T_i(G)] =
#' \alpha^3 {N-1 \choose 2}
#' + \alpha^2(1-2\alpha)((N-3)d_i + d_C)
#' + \alpha(1-2\alpha)^2
#' \left[
#' \frac{d_i(d_i-3)}{2}
#' + \sum_{j \neq i} C_{ij}d_j
#' \right]
#' + (1-2\alpha)^3 T_i(C)
#' }
#' where \eqn{d_i} is the degree of node \eqn{i}, \eqn{d_C} is the total
#' number of edges in the reference graph, and \eqn{T_i(C)} is the number of
#' triangles containing node \eqn{i} in \eqn{C}.
#'
#' @param i [numeric] Index of the node for which the expected triangle count is computed.
#' @param alpha [numeric] probability in \eqn{(0,1/2)} that an edge differs from that of the representative graph C.
#' @param C_mat [matrix] A N \times N matrix encoding the adjacency matrix associated to the representative graph of the CER model.
#'
#' @return [numeric] Expected number of triangles incident to node \eqn{i}.
#' @export

expected_count_triangles_node_CER <- function(i, alpha, C_mat){

  N = ncol(C_mat)

  degree_vec = colSums(C_mat)
  # degree_vec = degree(C_igraph, mode = "all", loops = FALSE, normalized = FALSE)
  degree_i   = degree_vec[i]
  degree_C = sum(degree_vec)/2

  T_i_C = diag(C_mat%*%C_mat%*%C_mat)[i]/2
  # T_i_C = count_triangles(C_igraph)[i]

  E_T_i_G =
    alpha^3 * choose(N-1, 2) + alpha^2 * (1-2*alpha)* ( (N-3)*degree_i + degree_C ) +
    alpha * (1-2*alpha)^2 * (0.5* degree_i*(degree_i-3) + sum( C_mat[i,][-i] * degree_vec[-i] ) )+
    # 0.5 * ( degree_i*(degree_i-2) + sum( C_mat[i,][-i] * degree_vec[-i] ) ) +
    (1-2*alpha)^3 * T_i_C

  return(E_T_i_G)
}

#' Expected total number of triangles under a CER model.
#'
#' Computes the expected total number of triangles in a random graph generated
#' from a centered Erdős--Rényi (CER) model with representative graph
#' \eqn{C}. The expectation depends on the Erdős--Rényi contribution,
#' the degree sequence of the representative graph, and the number of triangles
#' in the representative graph.
#'
#' The expected triangle count is computed as
#' \deqn{
#' E[T(G)] =
#' \alpha^3 {N \choose 3}
#' + \alpha^2(1-2\alpha)(N-2)d_C
#' + \frac{1}{2}\alpha(1-2\alpha)^2
#' \left(\sum_i d_i^2 - 2d_C\right)
#' + (1-2\alpha)^3 T(C)
#' }
#' where \eqn{d_i} denotes the degree of node \eqn{i},
#' \eqn{d_C=\sum_i d_i/2} is the number of edges in the reference graph,
#' and \eqn{T(C)} is the total number of triangles in \eqn{C}.
#'
#' @param alpha [numeric] probability in \eqn{(0,1/2)} that an edge differs from that of the representative graph C.
#' @param C_mat [matrix] A N \times N matrix encoding the adjacency matrix associated to the representative graph of the CER model.
#' @return [numeric] Expected total number of triangles in the graph.
#' @export
expected_count_triangles_CER <- function(alpha, C_mat){

  # C_mat = as.matrix(as_adjacency_matrix(C_igraph))
  # T_C = sum(diag(C_mat %*% C_mat %*% C_mat)) / 6

  N = ncol(C_mat)

  # degree_vec = degree(C_igraph, mode = "all", loops = FALSE, normalized = FALSE)
  degree_vec = colSums(C_mat)
  degree_C = sum(degree_vec)/2

  # T_C = sum(count_triangles(C_igraph))/3
  T_C = sum(diag(C_mat%*%C_mat%*%C_mat))/6

  E_T_G = alpha^3 * choose(N,3) + alpha^2 * (1-2*alpha)*(N-2)*degree_C +
    0.5*alpha*(1-2*alpha)^2 * (sum(degree_vec^2) - 2*degree_C) +
    (1-2*alpha)^3 * T_C

  return(E_T_G)
}

#' Expected value of G_i(G_i - 1).
#'
#' Computes the expected value of \eqn{G_i(G_i-1)} for node \code{i}
#' under the centered Erdős–Rényi model.
#'
#' @param i Index of the node for which the expected value of G_i(G_i - 1) is computed.
#' @param alpha probability in \eqn{(0,1/2)} that an edge differs from that of the representative graph C.
#' @param C_mat A N \times N matrix encoding the adjacency matrix associated to the representative graph of the CER model.
#'
#' @return [numeric] Expected value of G_i(G_i - 1).
#'
#' @export

E_G_i_x_G_i_minus_1 <- function(i, alpha, C_mat){
  # C_mat = as.matrix(as_adjacency_matrix(C_igraph))

  N = ncol(C_mat)
  # degree_vec = degree(C_igraph, mode = "all", loops = FALSE, normalized = FALSE)
  degree_vec = colSums(C_mat)
  degree_i   = degree_vec[i]

  out <-
    # (alpha*(N-1))^2 + alpha * (N-1) * (N-2) +
    alpha^2 * (N-1) * (N-2) +
    # alpha^2 *  ( (N-1)^2 - (N-1) ) +
    (1-2*alpha)*degree_i * ( (1-2*alpha)*degree_i + 2*alpha *(N-1) -1 )

  return(out)
}

#' Expected number of non-oriented paths of length two.
#'
#' Computes the expected number of non-oriented paths of length two
#' under the centered Erdős–Rényi model.
#'
#' @param alpha probability in \eqn{(0,1/2)} that an edge differs from that of the representative graph C.
#' @param C_mat A N \times N matrix encoding the adjacency matrix associated to the representative graph of the CER model.
#'
#' @return [numeric] Expected number of non-oriented paths of length two.
#'
#' @export
#'
expected_count_non_oriented_path2_CER <- function(alpha, C_mat){
  # C_mat = as.matrix(as_adjacency_matrix(C_igraph))

  N = ncol(C_mat)

  # degree_vec = degree(C_igraph, mode = "all", loops = FALSE, normalized = FALSE)
  degree_vec = colSums(C_mat)
  degree_C = sum(degree_vec)/2

  out <- 3*alpha^2*choose(N,3) + (1-2*alpha) * ( 2*(N-1)*alpha - 1  ) * degree_C +
    0.5*(1-2*alpha)^2*sum(degree_vec^2)

  return(out)
}

#' Approximate expected transitivity under the CER model.
#'
#' Computes the approximation of the expected transitivity (global clustering
#' coefficient) under the centered Erdős–Rényi model.
#' @param alpha probability in \eqn{(0,1/2)} that an edge differs from that of the representative graph C.
#' @param C_mat Reference graph defining the centered Erdős--Rényi model.
#'
#' @return [numeric] Approximate expected transitivity under the CER model.
#'
#' @export

expected_transitivity_CER_approx <- function(alpha, C_mat){

  N = ncol(C_mat)

  degree_vec = colSums(C_mat)
  degree_C = sum(degree_vec)/2

  T_C = sum(diag(C_mat%*%C_mat%*%C_mat))/6

  E_T_G = alpha^3 * choose(N,3) + alpha^2 * (1-2*alpha)*(N-2)*degree_C +
    0.5*alpha*(1-2*alpha)^2 * (sum(degree_vec^2) - 2*degree_C) +
    (1-2*alpha)^3 * T_C

  num = 3*E_T_G
  den = 3*alpha^2*choose(N,3) + (1-2*alpha) * ( 2*(N-1)*alpha - 1  ) * degree_C +
    0.5*(1-2*alpha)^2*sum(degree_vec^2)

  C_star = num/den

  return(C_star)
}

#' Generate a graph from a Generalized Centered Erdős–Rényi (GCER) model.
#'
#' @param W M-dimensional numeric vector encoding the location parameter W.
#' @param alpha concentration parameter in \eqn{(0,1/2)}.
#'
#' @return [numeric] A binary M-dimensional vector encoding the lower triangular part of the generated graph.
#'
#' @importFrom stats rbinom
G_from_GCER <- function(W, alpha){
  p_ij = ( 1 + (alpha/(1-alpha))^( 2*W-1 ) )^(-1)
  G_new<-rbinom(n=length(W), size=1, prob = p_ij)
  return(G_new)
}

#' @title Generalized CER model inference
#'
#' @param G_bar [numeric vector] M-dimensional vector containing the edge-wise
#' sample means of the observed graphs, where M = N(N - 1)/2. The entries encode
#' the lower triangular part of the edge-by-edge average of the adjacency matrices
#' of the observed graphs of $G^{(1:n)}$,
#' \eqn{\bar{G} = n^{-1}\sum_{i=1}^n G_i}, assuming the observed graphs follow
#' a Centered Erdős–Rényi (CER) distribution.
#' @param n [integer] The sample size.
#' @param prior_list [list] Ordered list of prior hyperparameters:
#' the location parameter \code{W} and \eqn{n_0 > 0} of the GCER prior for
#' \eqn{C}; and the shape parameters \eqn{a,b > 0} of the truncated Beta prior
#' for the scale parameter \eqn{\alpha}.
#' \code{W} is an M-dimensional numeric vector with entries in \eqn{[0,1]}
#' encoding the lower triangular part of the representative graph. The default
#' value is \eqn{1/2} for every edge, reflecting prior ignorance and yielding a
#' semiconjugate model. By default, \eqn{n_0 = a = b = 1}.
#'
#' @param mcmc [logical] If \code{TRUE}, posterior inference is performed via
#' MCMC. If \code{FALSE}, an empirical Bayes inferential approach is used.
#'
#' @param mcmc_options [list] List of MCMC options, including
#' \code{n_iter}, \code{burnin}, and \code{thinning}.
#' @param parallel [logical] If \code{TRUE}, parallel computation is used from
#' the EB method to maximize the marginal likelihood w.r.t. alpha.
#'
#' @return [list] A list containing the inferential output.
#' @importFrom parallel detectCores
#' @export

gCER_inference <- function(G_bar, n,
                           prior_list=list(W=1/2, n_0=1, a=1, b=1),
                           mcmc=FALSE,
                           mcmc_options=list( n_iter=1000, burnin=0, thinning=1  ),
                           parallel=TRUE){

  M = length(G_bar)

  W   = prior_list$W

  if( length(W)==1 ){
    W=rep(W, M)
  }

  n_0 = prior_list$n_0
  a   = prior_list$a
  b   = prior_list$b

  W_prime = (n_0 / (n_0 + n )) * W + (n / (n_0 + n )) * G_bar

  if(mcmc){
    n_iter   = mcmc_options$n_iter
    burnin   = mcmc_options$burnin
    thinning = mcmc_options$thinning

    iterations = burnin+thinning*n_iter

    # Gibbs sampling
    alpha_chain = rep( NA, n_iter )
    expected_clust_coeff_approx = rep( NA, n_iter )
    total_degree = rep( NA, n_iter )
    C_cum = rep(0, M)

    alpha_t = 0.5 # starting value
    count_it = 0
    for( t in 1:iterations ){
      tau_alpha_t = tau_alpha(alpha_t, n_0+n)
      p_ij = ( 1 + (tau_alpha_t/(1-tau_alpha_t))^( 2*W_prime-1 ) )^(-1)

      C_t     = rmvbern(n=1, p_ij)
      C_t_mat = from_low_adj_to_full_adj(C_t)

      # S_par = sum( sapply(1:nrow(G_n), function(i) sum(abs(G_n[i,]-C_t)) ) )
      S_par <- n * sum(abs(G_bar - C_t))

      alpha_t = rtbeta(n = 1 , alpha = a + S_par , beta = b + n*M - S_par, a = 0, b = 1/2)


      if ( t > burnin && ((t - burnin) %% thinning == 0) ){
        count_it = count_it + 1

        alpha_chain[count_it] = alpha_t

        C_cum = C_cum + C_t

        expected_clust_coeff_approx[count_it] = expected_transitivity_CER_approx_cpp(alpha_t, C_t_mat)
        total_degree[count_it] = M * alpha_t + sum(C_t) * ( 1 - 2*alpha_t )

      }
    }

    C_mean = C_cum / n_iter

    out = list(C_mean=C_mean, n_iter = n_iter, alpha_chain = alpha_chain,
               ecc_dist = expected_clust_coeff_approx,
               total_degree_dist = total_degree,
               type="gibbs" )
  }
  else{
    eps=1e-2

    # EB inference
    C_hat = ifelse( W_prime >= 1/2, 1, 0  )

    C_hat_mat = from_low_adj_to_full_adj(C_hat)

    # RcppParallel::setThreadOptions(numThreads = detectCores() - 1 )

    if(parallel==T){
      opt <- optim(par=0.1, fn=marg_lik_CER_gCER_par,
                   lower = 0+eps, upper=0.5-eps,
                   method = "L-BFGS-B",
                   n=n, n_0=n_0, G_bar=G_bar, W=W)
    }
    else{
      opt <- optim(par=0.1, fn=marg_lik_CER_gCER,
                   lower = 0+eps, upper=0.5-eps,
                   method = "L-BFGS-B",
                   n=n, n_0=n_0, G_bar=G_bar, W=W)
    }

    alpha_hat <- opt$par

    ecc = expected_transitivity_CER_approx_cpp(alpha_hat, C_hat_mat)

    total_degree = M * alpha_hat + sum(C_hat) * ( 1 - 2*alpha_hat )

    out = list(C_hat=C_hat, alpha_hat=alpha_hat, ecc=ecc, total_degree=total_degree, type="EB" )
  }

  return(out)
}

