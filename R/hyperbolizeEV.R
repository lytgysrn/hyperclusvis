#' Hyperbolic cluster mapping with multiple anchors
#'
#' This function takes a data matrix, maps both cluster centers and observations into the
#' 2D hyperbolic disk based on given cluster labels (the extended version in the paper).
#'
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
#' @param local_model Character. Clustering model used to define local
#'   sub-centers within each global cluster (e.g. \code{"kmeansWrap"}).
#' @param kg_per_class Integer. Number of local sub-centers per cluster.
#' @param lambda Numeric. Ridge parameter used in the local projection step.
#' @param spread_scale Numeric. Controls how far local sub-centers are spread
#'   from the global centroid in the tangent plane (corresponds to \eqn{\delta}
#'   in the paper).
#' @param use_input_order Logical. If \code{FALSE} (default), cluster centers are
#'   reordered using \code{order_clusters_tsp}. If \code{TRUE}, the input cluster order is kept.
#'
#'
#'
#' @return A list with components:
#' \describe{
#'   \item{points}{\code{n × 2} matrix. Final 2D hyperbolic embedding of all samples.}
#'   \item{centers}{\code{numk × 2} matrix. 2D positions of cluster centroids.}
#'   \item{cluster}{Integer vector of length \code{n}. Cluster label of each sample.}
#'   \item{analysis}{Model fit object returned by the chosen clustering method
#'    , or \code{NULL} if cluster labels were supplied directly.}
#'   \item{local_info}{List of length \code{numk}. For each global cluster \eqn{g}
#'     it contains local mapping objects: indices of observations, centroid in
#'     original space and tangent space, fitted local sub-centers, their 2D
#'     coordinates, and both tangent-plane and disk embeddings of samples.}
#'   \item{params}{List of the main tuning parameters used in the call.}
#' }
#'
#' @examples
#' \dontrun{
#' ## Iris demo
#' set.seed(42)
#' out <- hyperbolizeEV(
#'   iris[, 1:4],
#'   numk       = 3,
#'   preprocess = "kmeansWrap",
#'   scale      = "std"
#' )
#' plot_disk(out, 3)
#'
#' ## Yale faces demo (from PPCI)
#' if (!requireNamespace("PPCI", quietly = TRUE)) {
#'   message("Package 'PPCI' is not installed")
#' } else {
#'   set.seed(42)
#'   library(PPCI)
#'   out <- hyperbolizeEV(
#'     PPCI::yale$x,
#'     numk         = 10,
#'     preprocess   = PPCI::yale$c,
#'     scale        = "norm",
#'     local_model  = "kmeansWrap",
#'     kg_per_class = 10,
#'     spread_scale = 0.6
#'   )
#'   plot_disk(out, 10)
#' }
#' }
#'
#' @export
hyperbolizeEV <- function(D,
                          numk,
                          preprocess,
                          scale=c("std","norm"),
                          local_model="kmeansWrap",
                          kg_per_class=5,
                          lambda=1e-6,
                          spread_scale=1.0,
                          use_input_order=F) {
  #input
  #D: data set
  #numK: num clusters
  #scale: scaling method
  #local_model: clustering models for defining local sub-centers
  #kg_per_class: num local sub-centers
  #lambda: ridge parameter
  #spread_scale：delta in the paper
  #use_input_order: if just use the input labels' order

  scale <- match.arg(scale)
  n_obs=nrow(D)

  if (scale == "norm") {
    Dp <- PrepareData(D)
  } else {
    Dp <- scale(D)
  }
  if (is.character(preprocess)){
    if (preprocess == "kmeansWrap") {
      cl <- kmeansWrap(Dp, numk); analysis <- cl
      C <- cl$centers; memv <- cl$cluster
    } else if ((preprocess == "kmeans") || (preprocess == "kmedians")) {
      cl <- flexclust::kcca(Dp, numk, family=flexclust::kccaFamily(preprocess))
      analysis <- cl
      C <- flexclust::parameters(cl)
      memv <- flexclust::clusters(cl)
    } else if (preprocess == "hddc") {
      prms <- HDclassif::hddc(Dp, K=numk); analysis <- prms
      C <- prms$mu; C <- C[, ]; memv <- prms$class
    } else if (preprocess == "pgmm") {
      E <- scale(Dp)
      cl <- pgmm::pgmmEM(E, rG=numk); analysis <- cl
      memv <- cl$map
      C <- CentroidsFromMembership(Dp, memv); C <- C[, ]
    }
  }else {
    ## ----------when cluster labels are provided directly
    analysis <- NULL
    cl <- preprocess

    if (is.vector(cl) || is.factor(cl)) {
      memv <- as.integer(cl)
      C <- CentroidsFromMembership(Dp, memv)
    } else if (is.list(cl)) {
      memv <- as.integer(cl$cluster)
      C <- CentroidsFromMembership(Dp, memv)
    }
  }

  if (!use_input_order){
    ord_res <- order_clusters_tsp(C, memv)
    memv <- ord_res$memv_new
    C<-ord_res$centers
  }

  ## AS OV
  C2d_star <- star(C)              # centroids in T_{c0}
  C2d_disk <- CentersToDisk(C)     # centroids in disk

  #for each cluster, define local sub-centers and map to 2d tangent plane
  G <- sort(unique(memv))
  local_info <- vector("list", length(G))
  names(local_info) <- paste0("G", G)

  for (idx in seq_along(G)) {
    g    <- G[idx]
    indg <- which(memv == g)
    Xg   <- Dp[indg, , drop = FALSE]

    cg      <- C[idx, , drop = FALSE]  # centroid in original space
    cg2d    <- C2d_star[idx, ]         # its representation in T_{c0}
    #local sub-centers
    kg <- min(kg_per_class, max(1, nrow(unique(Xg))))
    sub <- fit_local_subcenters(Xg, kg = kg, local_model = local_model)
    Mg  <- sub$centers

    # treat Mg as data，use C to obtain weights
    Wg <- ComputeWeightList(Mg, C)
    # local sub_centers map to the tangent plane ṽ_{kg}
    kg_now <- nrow(Mg)
    Mg_global2d <- matrix(0, nrow = kg_now, ncol = 2)
    for (jj in 1:kg_now) {
      wt_j <- as.vector(Wg[jj, ])
      Mg_global2d[jj, ] <- Map_tan(wt_j, C2d_star)
    }

    #    u_{kg}(s) = cg2d + s ( ṽ_{kg} - cg2d )
    diff_global <- sweep(Mg_global2d, 2, as.numeric(cg2d), FUN = "-")
    M2d         <- sweep(spread_scale * diff_global, 2,
                         as.numeric(cg2d), FUN = "+")  # k_g x 2

    # 4) relative coordinates（cg2d as origin）for latter project_by_coeffs
    Dg2_t <- t( sweep(M2d, 2, as.numeric(cg2d), FUN = "-") )  # 2 x k_g

    #mapping in 2d tangent plane
    proj <- project_by_coeffs(
      Mg    = Mg,
      cg    = cg,
      Xg    = Xg,
      Dg2   = Dg2_t,
      cg2d  = cg2d,
      lambda = lambda
    )

    local_info[[idx]] <- list(
      class      = g,
      ind        = indg,
      cg         = cg,
      cg2d       = cg2d,
      subcenters = Mg,
      submodel   = sub,
      centers2d  = M2d,
      Dg2_t      = Dg2_t,
      samples2d  = proj
    )
  }

  ## to Poincaré disk

  map_to_disk <- function(M) {
    map_to_disk_point <- function(pt2) {
      D <- sqrt(sum(pt2^2))
      if (D < 1e-12) return(c(0,0))
      r <- convert(D)
      (r / D) * pt2
    }
    if (is.null(M) || length(M)==0) return(M)
    t(apply(M, 1, map_to_disk_point))
  }


  for (i in seq_along(local_info)) {
    li <- local_info[[i]]
    local_info[[i]]$centersDisk <- map_to_disk(li$centers2d)
    local_info[[i]]$samplesDisk <- map_to_disk(li$samples2d)
    local_info[[i]]$cgDisk      <- map_to_disk(matrix(li$cg2d, nrow = 1))
  }

  Xdisk    <- matrix(0, nrow = n_obs, ncol = 2)
  C2d_disk <- matrix(0, nrow = numk,  ncol = 2)

  for (i in seq_along(local_info)) {
    li  <- local_info[[i]]
    ind <- li$ind
    g   <- li$class
    Xdisk[ind, ]  <- li$samplesDisk
    C2d_disk[g, ] <- li$cgDisk[1, ]
  }

  list(points=Xdisk,
       centers=C2d_disk,
       cluster=memv,
       analysis=analysis,
       local_info=local_info,
       params=list(numk=numk, kg_per_class=kg_per_class,
                   preprocess=preprocess, scale=scale,
                   local_model=local_model,
                   lambda=lambda, spread_scale=spread_scale))
}
