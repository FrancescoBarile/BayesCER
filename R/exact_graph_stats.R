
#' Compute exact graph statistics
#'
#' Computes exact graph statistics from the upper-triangular
#' vector representation of an undirected, unweighted network.
#' The function returns the number of edges, the sum of squared
#' node degrees, and the number of triangles.
#'
#' @param C Integer vector containing the upper-triangular entries
#'   of the network adjacency matrix, stored column-wise.
#'   Entries equal to 1 indicate the presence of an edge; entries
#'   equal to 0 indicate its absence.
#'
#' @return A numeric vector of length three containing:
#' \describe{
#'   \item{edges}{Number of edges in the network.}
#'   \item{d2}{Sum of squared node degrees.}
#'   \item{triangles}{Number of triangles in the network.}
#' }
#'
#' @details
#' The number of nodes is inferred from the length of `C`, which
#' must equal N(N-1)/2 for some integer N. The function assumes
#' that `C` represents the upper triangular part of a symmetric
#' adjacency matrix, excluding the diagonal.
#'
#' @name exact_graph_stats
NULL
