# Hyperbolic cluster mapping

This function takes a data matrix, maps both cluster centers and
observations into the 2D hyperbolic disk based on given cluster labels.

## Usage

``` r
hyperbolize(D, numk, preprocess, scale = c("std", "norm"), use_input_order = F)
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

## Examples

``` r
if (FALSE) { # \dontrun{
  set.seed(100000)
  D <- scale(iris[, 1:4])
  out <- hyperbolize(D, numk = 3,scale="std", preprocess = "kmeansWrap")
  plot_disk(out,3)
} # }
```
