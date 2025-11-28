
#' @keywords internal
######Define local sub-centers###############
fit_local_subcenters <- function(Xg, kg=5, local_model=c("kmeansWrap","kmedian")) {
  #Xg: n_g x d  points in this Cluster.
  #kg: num of local sub-centers.
  #local_model: clustering models for defining local sub-centers
  if (local_model == "kmeansWrap") {
    km <- kmeansWrap(Xg, kg)
    list(centers=km$centers, model=km, type="kmeans")
  } else if (local_model=="kmedian"){
    cl <- flexclust::kcca(Xg, kg, family=flexclust::kccaFamily('kmedian'))
    C <- flexclust::parameters(cl)
    C<-C[,]
    list(centers=C,model=cl,type='kmedian')
  }
}


#projection to the tangent plane for each cluster.
project_by_coeffs <- function(Mg, cg, Xg, Dg2, cg2d, lambda = 1e-6) {
  #   Mg   : k_g x d local sub-centers
  #   cg   : 1 x d centroid for this group
  #   Xg   : n_g x d  points in this Cluster.
  #   Dg2  : 2 x k_g  sub-centers coordinates in the 2d tangent plane
  #   cg2d : 1 x 2  centroid in 2d tangent plane
  #   lambda : ridge parameter
  #
  # return：
  #   Yg2_abs : n_g x 2 points coordinates in the 2d tangent plane
  k <- nrow(Mg); d <- ncol(Mg); n <- nrow(Xg)

  # directional vector from sub-centers to cg
  Dg <- t( t(Mg) - as.vector(cg) )   # k x d
  Dg <- t(Dg)                        # d x k

  # points' directional vector towards cg
  Yg <- t( t(Xg) - as.vector(cg) )   # n x d
  Yg <- t(Yg)                        # d x n

  #Lagrange multiplier with constraint sum a=1
  G   <- crossprod(Dg)                   # k x k
  L   <- G + diag(lambda, k)            # k x k
  ones <- matrix(1, nrow = k, ncol = 1) # k x 1

  #  (k+1) x (k+1) linear system
  Kmat <- rbind(
    cbind(L,    ones),
    cbind(t(ones), 0)
  )

  # right sight (k+1) x n ：with first k row  Dg^T Y, last row all 1
  RHS <- rbind(
    crossprod(Dg, Yg),                  # k x n
    matrix(1, nrow = 1, ncol = n)       # 1 x n
  )
  #solve this to obtain coefficients
  sol <- solve(Kmat, RHS)
  A   <- sol[1:k, , drop = FALSE]



  # -------- then use those coefficients in the tangent plane：y^(2) = Dg2 %*% a --------
  Yg2_rel <- t(Dg2 %*% A)            # n x 2

  # translate to the global position
  Yg2_abs <- sweep(Yg2_rel, 2, as.numeric(cg2d), `+`)  # n x 2

  return(Yg2_abs)
}

#map to the tangent plane
Map_tan <- function(wts, star) {
  drop(wts %*% star)
}
