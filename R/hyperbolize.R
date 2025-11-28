#' Hyperbolic cluster mapping
#'
#' This function takes a data matrix, maps both cluster centers and observations into the
#' 2D hyperbolic disk based on given cluster labels.
#'
#' @param D Numeric matrix, with observations in rows and variables in columns.
#' @param numk Integer. Number of clusters.
#' @param preprocess Either a character string specifying the clustering method
#'   (\code{"kmeansWrap"}, \code{"kmeans"}, \code{"kmedians"},
#'   \code{"hddc"}, \code{"pgmm"}), or directly a clustering result / label vector:
#' @param scale Character. Scaling choice for the data:
#'   \itemize{
#'     \item \code{"std"}: standardize using \code{scale(D)}.
#'     \item \code{"norm"}: center and rescale to a fixed average norm (paper sec 3.2).
#'   }
#' @param use_input_order Logical. If \code{FALSE} (default), cluster centers are
#'   reordered using \code{order_clusters_tsp}. If \code{TRUE}, the input cluster order is kept.
#'
#' @return A list with components:
#' \describe{
#'   \item{points}{\code{n × 2} matrix. Final 2D hyperbolic embedding of all samples.}
#'   \item{centers}{\code{numk × 2} matrix. 2D positions of cluster centroids.}
#'   \item{cluster}{Integer vector of length \code{n}. Cluster label of each sample.}
#'   \item{analysis}{Model fit object returned by the chosen clustering method
#'    , or \code{NULL} if cluster labels were supplied directly.}
#' }
#'
#' @examples
#' \dontrun{
#'   set.seed(100000)
#'   D <- scale(iris[, 1:4])
#'   out <- hyperbolize(D, numk = 3,scale="std", preprocess = "kmeansWrap")
#'   plot_disk(out,3)
#' }

#' @export
hyperbolize <- function(D,
                        numk,
                        preprocess,
                        scale=c("std", "norm"),
                        use_input_order=F
) {
  if (scale == "norm") {
    D <- PrepareData(D)
  } else if (scale == "std") {
    D <- scale(D)
  }
  if (is.character(preprocess)){
    if (preprocess == "kmeansWrap") {
      cl <- kmeansWrap(D, numk); analysis <- cl
      C <- cl$centers; memv <- cl$cluster
    } else if ((preprocess == "kmeans") || (preprocess == "kmedians")) {
      cl <- flexclust::kcca(D, numk, family=flexclust::kccaFamily(preprocess))
      analysis <- cl
      C <- flexclust::parameters(cl)
      memv <- flexclust::clusters(cl)
    } else if (preprocess == "hddc") {
      prms <- HDclassif::hddc(D, K=numk); analysis <- prms
      C <- prms$mu; C <- C[, ]; memv <- prms$class
    } else if (preprocess == "pgmm") {
      E <- scale(D)
      cl <- pgmm::pgmmEM(E, rG=numk); analysis <- cl
      memv <- cl$map
      C <- CentroidsFromMembership(D, memv); C <- C[, ]
    }
  }else {
    ## ----------when cluster labels are provided directly
    analysis <- NULL
    cl <- preprocess

    if (is.vector(cl) || is.factor(cl)) {
      memv <- as.integer(cl)
      C <- CentroidsFromMembership(D, memv)
    } else if (is.list(cl)) {
      memv <- as.integer(cl$cluster)
      C <- CentroidsFromMembership(D, memv)
    }
  }
  if (!use_input_order){
    ord_res <- order_clusters_tsp(C, memv)
    memv <- ord_res$memv_new
    C<-ord_res$centers
  }
  hypC <- CentersToDisk(C)
  W <- ComputeWeightList(D, C)
  st <- star(C)
  hx <- matrix(c(0, 0), nrow=1)
  sz <- nrow(D)
  for (i in 1:sz) {
    wt <- as.vector(W[i, ])
    pt <- matrix(Map(wt, st), nrow=1)
    hx <- rbind(hx, pt)
  }


  return(list(points=hx, centers=hypC, cluster=memv, analysis=analysis))
}
