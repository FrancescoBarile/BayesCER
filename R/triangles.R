
#' Count the number of triangles in a graph
#'
#' Computes the number of triangles in an undirected, unweighted
#' graph represented by a bit-packed adjacency matrix. The
#' computation uses bitwise operations to identify common
#' neighbours of adjacent vertices.
#'
#' @param a Vector of unsigned 64-bit integers representing the
#'   bit-packed adjacency matrix. Each row occupies `nw`
#'   consecutive entries.
#' @param N Number of vertices in the graph.
#' @param nw Number of 64-bit words used to represent each row
#'   of the adjacency matrix, typically `ceiling(N / 64)`.
#'
#' @return The number of triangles in the graph.
#'
#' @details
#' Each triangle is counted once for each of its three edges
#' through the common-neighbour calculation. The result is
#' therefore divided by three to obtain the total number of
#' distinct triangles.
#'
#' @name triangles
NULL
