#' FastPotts: Loopy Belief Propagation for the Potts Model
#'
#' Provides efficient inference algorithms for the Potts model, a classical
#' graphical model used for spatial labeling problems. Implements loopy belief
#' propagation (LBP) for computing marginal probabilities on undirected graphs.
#'
#' The primary use case is cell type segmentation in spatial transcriptomics data,
#' where individual transcripts are assigned to cells based on their spatial locations
#' and gene expression signatures.
#'
#' @section Main Functions:
#' - [run_crf()] - Performs cell type assignment via Potts model CRF
#' - [hilbert_index()] - Computes Hilbert curve indices for spatial ordering
#' - [potts_lbp_parallel_cpp()] - Parallel loopy belief propagation inference
#' - [build_potts_lbp_graph()] - Converts sparse matrices to LBP graph format
#'
#' @section C++ Backend:
#' The package uses Rcpp to call optimized C++ implementations of:
#' - Loopy belief propagation message passing
#' - Numerically stable log-domain computations
#' - Hilbert curve index calculations
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom dplyr arrange
#' @importFrom dplyr filter
#' @importFrom dplyr mutate
#' @importFrom dplyr row_number
#' @importFrom magrittr %>%
#' @importFrom Rcpp evalCpp
#' @importFrom Rcpp sourceCpp
#' @importFrom rlang .data
#' @importFrom rlang sym
#' @useDynLib FastPotts, .registration = TRUE
## usethis namespace: end
NULL
  
