
#' @keywords internal

################ Hyperbolic distance

vlength <- function(v) {
  ret <- sqrt(crossprod(v)[1, 1])
  return(ret)
}

vlengthsq <- function(v) {
  ret <- crossprod(v)[1, 1]
  return(ret)
}

delta <- function(u, v) {
  numer <- vlengthsq(u - v)#square norm of u-v
  denom <- (1 - vlengthsq(u)) * (1 - vlengthsq(v))
  ret <- 2 * numer / denom
  return(ret)
}


#' Hyperbolic distance in the Poincaré disk
#'
#' Computes the hyperbolic distance between two points inside the
#' 2-dimensional Poincaré disk model.
#'
#' @param u Numeric vector of length 2. Must satisfy \code{||u|| < 1}.
#' @param v Numeric vector of length 2. Must satisfy \code{||v|| < 1}.
#'
#' @return A numeric value giving the hyperbolic distance between \code{u} and \code{v}.
#'
#' @examples
#' hypdist(c(0.1, 0), c(0.2, 0.3))
#'
#' @export
hypdist <- function(u, v) {
  ret <- acosh(1 + delta(u, v))
  return(ret)
}

#to polar form
convert <- function(d) {
  ret <- (cosh(d) - 1) / (cosh(d) + 1)
  return(sqrt(ret))
}

################ Orthogonal projection

on <- function(B) {
  # Computes an orthonormal basis.
  #
  # Args:
  #   B: Matrix with linearly independent rows.
  #
  # Returns:
  #   Matrix with dimensions matching those of B.
  #   The rows are orthonormal, and the span of the rows
  #   is the span of the rows of B.
  Q <- qr.Q(qr(t(B))) #t(B)=QR
  return(t(Q))
}

orthproj <- function(v, B) {
  # Computes an orthogonal projection.
  #
  # Args:
  #   v: Numeric vector.
  #   B: Matrix with linearly independent rows.
  #      Require length(v) = ncol(B).
  # Returns:
  #   Numeric vector equal to the orthogonal projection of v
  #   onto the subspace spanned by the rows of B.
  A <- on(B)
  coeffs <- A %*% v # the expression in the new coordiate system
  n <- nrow(B)
  pr <- rep(0, ncol(B))
  for (i in 1:n) {
    pr <- pr + coeffs[i] * as.vector(A[i, ])
  }
  return(pr)
}

#probably be a simpler version
orthproj2 <- function(v, B) {
  A <- on(B)
  as.vector(t(A) %*% (A %*% v))  # = (A^T A) v
}

################ Convex combinations

CentroidsToBasis <- function(C) {
  n = nrow(C)
  B = matrix(C[1, ], nrow=n, ncol=ncol(C), byrow=TRUE)
  A <- C - B
  A <- A[2:n, ]
  if (!is.matrix(A)) {
    A <- matrix(A, nrow = 1)
  }else {
    ret <- A
  }
  return(ret)
}

ConvexWeights <- function(p, B) {
  A <- B %*% t(B)
  b <- B %*% p
  x <- as.vector(qr.solve(A, b))
  val <- 1 - sum(x)
  return(c(val, x))
}

#3.5
ComputeWeightList <- function(data, C) {
  bp <- as.vector(C[1, ])
  B <- CentroidsToBasis(C)#move all points to the coordinate system with C1 as origin
  N <- nrow(data) # size of data
  TT <- matrix(rep(bp, N), nrow=N, byrow=TRUE)
  data <- data - TT
  pjM <- rbind(rep(0, ncol(data)))
  WL <- rep(0, nrow(C))
  for (i in 1:N) {
    dat <- as.vector(data[i, ])
    pj <- orthproj(dat, B)
    wts <- ConvexWeights(pj, B)
    WL <- rbind(WL, wts)
    pjM <- rbind(pjM, pj)
  }
  pjM <- pjM[2:(N + 1), ]
  return(WL[2:(N + 1), ])
}

################ Construct radial simplex in the tangent space

## 2D tangent space version

Center <- function(C) {
  return(as.vector(colMeans(C)))
}

angle <- function(u, v) {
  cosVal <- sum(u * v) / sqrt(sum(u * u) * sum(v * v))
  return(acos(cosVal))
}

CentralAngles <- function(C) {

  c <- Center(C)
  n=nrow(C)
  V   <- sweep(C, 2, c, FUN = "-")
  ang <- numeric(n)
  for (i in 1:n) {
    j <- if (i == n) 1L else (i + 1L)
    v <- V[i, ]
    w <- V[j, ]
    ang[i] <- angle(v, w)
  }
  return(ang)
}

## This is the intrinsically 2D part,
## in that in 2D there's a circle's worth of angles.
# theta in 2.4
AngleSeq <- function(E) {
  angs <- CentralAngles(E)
  total <- sum(angs)
  for (i in 2:length(angs)) {
    angs[i] <- angs[i] + angs[i - 1]
  }
  return(angs * 2 * pi / total)
}

################ Euclidean distances

spray <- function(E) {
  n <- nrow(E)
  c <- Center(E)
  dists <- rep(0, n)
  for (i in 1:n) {
    dists[i] <- vlength(c - as.vector(E[i, ]))
  }
  return(dists)
}

star <- function(E) {
  dists <- spray(E)
  angs <- c(0, AngleSeq(E)) # last angle is only used implicitly (in total of AngleSeq procedure)
  n <- nrow(E)
  vts <- matrix(rep(0, 2 * n), nrow=n)
  for (i in 1:n) {
    vts[i, 1] <- dists[i] * cos(angs[i])
    vts[i, 2] <- dists[i] * sin(angs[i])
  }
  return(vts)
}

################ Map

Map <- function(wts, star) {
  tanv <- c(0, 0)
  for (i in 1:nrow(star)) {
    tanv <- tanv + wts[i] * as.vector(star[i, ])
  }
  D = sqrt(tanv[1] ^ 2 + tanv[2] ^ 2)
  r = convert(D) #sqrt(cosh(R)-1/cos(R)+1)
  return((r/D) * tanv) #scale to rcos rsin
}

################ Automated toy example building


CentersToDisk <- function(C) {
  st <- star(C) #image in tangent space
  n <- nrow(C)
  HC <- c(0, 0) # default initialization
  wt <- rep(0, n) # holds current weights
  for (i in 1:n) {
    wt[i] <- 1
    hc <- Map(wt, st)#disk coordinate
    HC <- rbind(HC, hc)
    wt[i] <- 0
  }
  return(HC[2:(n + 1), ])
}



################


################

## Center and scale the data

CenterData <- function(D) {
  c <- colSums(D) / (nrow(D))
  TT <- matrix(rep(c, nrow(D)), byrow=TRUE, nrow=nrow(D))
  return(D - TT)
}


AvgNorm <- function(D) {
  sum <- 0
  N <- nrow(D)
  d <- 0
  for (i in 1:N) {
    d <- sqrt(D[i, ] %*% D[i, ])
    sum <- sum + d
  }
  return(sum / N)
}

# perform only on CENTERED data
ScaleData <- function(D) {
  kF <- 0.5  # scaling constants
  l <- kF * AvgNorm(D)
  l <- as.numeric(l)
  A <- D / l
  return(A)
}

#3.2
PrepareData <- function(D) {
  D <- scale(D, scale=FALSE)
  D <- ScaleData(D)
  return(D)
}

################

## preclustering

SumSqs <- function(D, C, memv) {
  # Sums squares of distances from data to centers.
  #
  # Args:
  #   D   : Data matrix. Each row is a data point.
  #   C   : Matrix. Each row is a center.
  #   memv: Membership vector.
  #
  # Returns:
  #   Sum of squares of distances from each data point to its
  #   nearest center.
  ret <- 0
  inc <- 0
  N <- nrow(D)
  for (i in 1:N) {
    dat <- as.vector(D[i, ])
    cen <- as.vector(C[memv[i], ])
    inc <- sum((cen - dat)^2)
    ret <- ret + inc
  }
  return(ret)
}

kmeansWrap <- function(D, numk, r=5) {
  cl <- stats::kmeans(D, numk)
  for (i in 2:r) {
    tmp <- stats::kmeans(D, numk)
    if (tmp$tot.withinss < cl$tot.withinss) {
      cl <- tmp
    }
  }
  return(cl)
}

################

CentroidsFromMembership <- function(D, memv) {
  cl <- sort(unique(memv))
  k <- length(cl)
  C <- rbind(rep(0, ncol(D)))
  for (i in 1:k) {
    centroid <- rep(0, ncol(D))
    ind <- which(memv == i)
    for (j in 1:length(ind)) {
      centroid <- centroid + as.vector(D[ind[j], ])
    }
    centroid <- centroid / length(ind)
    C <- rbind(C, centroid)
  }
  return(as.matrix(C[-1, ]))
}

#sort centroids based on greedy tsp
order_clusters_tsp <- function(C, memv) {
  #C: centroids
  #memv: clustering membership
  labs <- sort(unique(memv))
  K <- length(labs)

  dist_mat <- as.matrix(stats::dist(C))
  # define c1
  start <- which.min(rowSums(dist_mat))

  diag(dist_mat) <- Inf


  visited <- rep(FALSE, K)
  order_idx <- integer(K)
  cur <- start
  for (i in 1:K) {
    order_idx[i] <- cur
    visited[cur] <- TRUE
    drow <- dist_mat[cur, ]
    drow[visited] <- Inf
    cur <- if (i == K) start else which.min(drow)
  }

  # order map to memv
  labs_order <- labs[order_idx]
  new_ids <- seq_len(K)
  names(new_ids) <- labs_order
  memv_new <- unname(new_ids[as.character(memv)])

  C_sorted <- C[order_idx, , drop = FALSE]
  list(
    memv_new   = memv_new,
    centers    = C_sorted
  )
}
