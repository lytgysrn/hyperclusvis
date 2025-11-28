# Hyperbolic cluster mapping with multiple anchors

This function takes a data matrix, maps both cluster centers and
observations into the 2D hyperbolic disk based on given cluster labels
(the extended version in the paper).

## Usage

``` r
hyperbolizeEV(
  D,
  numk,
  preprocess,
  scale = c("std", "norm"),
  local_model = "kmeansWrap",
  kg_per_class = 5,
  lambda = 1e-06,
  spread_scale = 1,
  use_input_order = F
)
```

## Arguments

- D:

  Numeric matrix, with observations in rows and variables in columns.

- numk:

  Integer. Number of clusters.

- preprocess:

  Either a character string specifying the clustering method
  (`"kmeansWrap"`, `"kmeans"`, `"kmedians"`, `"hddc"`, `"pgmm"`), or
  directly a clustering result / label vector:

- scale:

  Character. Scaling choice for the data:

  - `"std"`: standardize using `scale(D)`.

  - `"norm"`: center and rescale to a fixed average norm (paper sec
    3.2).

- local_model:

  Character. Clustering model used to define local sub-centers within
  each global cluster (e.g. `"kmeansWrap"`).

- kg_per_class:

  Integer. Number of local sub-centers per cluster.

- lambda:

  Numeric. Ridge parameter used in the local projection step.

- spread_scale:

  Numeric. Controls how far local sub-centers are spread from the global
  centroid in the tangent plane (corresponds to \\\delta\\ in the
  paper).

- use_input_order:

  Logical. If `FALSE` (default), cluster centers are reordered using
  `order_clusters_tsp`. If `TRUE`, the input cluster order is kept.

## Value

A list with components:

- points:

  `n × 2` matrix. Final 2D hyperbolic embedding of all samples.

- centers:

  `numk × 2` matrix. 2D positions of cluster centroids.

- cluster:

  Integer vector of length `n`. Cluster label of each sample.

- analysis:

  Model fit object returned by the chosen clustering method , or `NULL`
  if cluster labels were supplied directly.

- local_info:

  List of length `numk`. For each global cluster \\g\\ it contains local
  mapping objects: indices of observations, centroid in original space
  and tangent space, fitted local sub-centers, their 2D coordinates, and
  both tangent-plane and disk embeddings of samples.

- params:

  List of the main tuning parameters used in the call.

## Examples

``` r
if (FALSE) { # \dontrun{
## Iris demo
set.seed(42)
out <- hyperbolizeEV(
  iris[, 1:4],
  numk       = 3,
  preprocess = "kmeansWrap",
  scale      = "std"
)
plot_disk(out, 3)

## Yale faces demo (from PPCI)
if (!requireNamespace("PPCI", quietly = TRUE)) {
  message("Package 'PPCI' is not installed")
} else {
  set.seed(42)
  library(PPCI)
  out <- hyperbolizeEV(
    PPCI::yale$x,
    numk         = 10,
    preprocess   = PPCI::yale$c,
    scale        = "norm",
    local_model  = "kmeansWrap",
    kg_per_class = 10,
    spread_scale = 0.6
  )
  plot_disk(out, 10)
}
} # }
```
