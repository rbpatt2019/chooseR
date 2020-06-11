#' Find clusters in Seurat object
#'
#' This function combines the necessary and sequential calls
#'  to Seurat::FindNeighbors and Seurat::FindClusters
#'
#' @param obj Seurat object
#' @param reduction String naming the reduction for \code{Seurat::FindNeighbors}
#' @param npcs Integer of number of dimensions for \code{Seurat::FindNeighbors}
#' @param assay String naming the assay for \code{Seurat::FindNeighbors}
#'  This assay is also used to define features in the default settings
#' @param features Character vector containing features for \code{Seurat::FindNeighbors}
#' @param resolution Numeric of resolution to use for \code{Seurat::FindClusters}
#' @param verbose Boolean controlling verbosity
#'  If \code{TRUE}, print updates to terminal
#'
#' @return A Seurat object with the new clusters stored in Idents and the
#'  old clusters stored in \code{seurat_clusters}
#'
#' @section Warning:
#' Each call to this function overwrites \code{seurat_clusters}
#'
#' @seealso \url{https://www.rdocumentation.org/packages/Seurat/versions/3.1.0/topics/FindClusters}
#'  and\url{https://www.rdocumentation.org/packages/Seurat/versions/3.1.0/topics/FindNeighbors}
#'  for further documentation on the component functions.
#'
#' @examples
#' data <- find_clusters(obj)
#' data <- find_clusters(obj, npcs = 30, verbose = TRUE)
#' \dontrun{
#' data <- find_clusters()
#' }
#' @export
find_clusters <- function(
                          obj,
                          reduction = "pca",
                          npcs = 100,
                          assay = "integrated",
                          features = NULL,
                          resolution = 0.8,
                          verbose = FALSE) {
  obj <- Seurat::FindNeighbors(
    obj,
    reduction = reduction,
    dims = 1:npcs,
    assay = integtrated,
    features = features,
    verbose = verbose,
    graph.name = paste(reduction, assay, sep = ".")
  )
  obj <- Seurat::FindClusters(
    obj,
    resolution = resolution,
    graph.name = paste(reduction, assay, sep = "."),
    verbose = verbose
  )
  return(obj)
}

#' Generate n sub-samples
#'
#' Wrapper function around rplicate and split to generate
#'  \code{n} subsamples of size \code{size} from \code{input}
#'
#' @param n Integer number of replicates
#' @param input Vector to be subset
#' @param size Numeric fraction of input to take for subsample
#' @param replace Boolean whether to sample with or without replacement
#' @param simplify Whether to simplify output to array or matrix, if possible
#'
#' @return A list of subsamples of \code{input}.
#'  If \code{simplify}, then in a simplfied array or matrix
#'
#' @seealso \code{\link[base]{replicate}} and \code{\link[base]{sample}}
#'  for documentation on the component functions.
#
#' @examples
#' n_samples(3, c(1, 3, 5, 2))
#' n_samples(5, c('a', 'b', 'c'), size = 0.5, replace = TRUE)
#'
#' \dontrun{
#' n_samples(c(1, 2, 3))
#' }
#' @export
n_samples <- function(
                     n,
                     input,
                     size = 0.8,
                     replace = FALSE,
                     simplify = FALSE) {
  splits <- replicate(
    n,
    sample(
      input,
      as.integer(length(input) * size),
      replace = replace
    ),
    simplify = simplify
  )
}
