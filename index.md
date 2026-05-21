# hyperclusvis

Hyperbolic Cluster Mapping for Visualization of High-Dimensional
Clustering Results

## Overview

`hyperclusvis` implements the hyperbolic cluster mapping method
described in:

> Wang, S., Cappadocia, C. & McNicholas, P.D. “Hyperbolic Cluster
> Mapping: A Visualization Method for Assessing Clustering Quality.”

The method maps high-dimensional clustered data into the 2D Poincaré
disk, producing interpretable visualizations where well-separated
clusters appear as distinct regions and poorly separated clusters
overlap visibly. Two versions are provided:

- **OV (Original Version)** : uses one centroid per cluster as the
  anchor for the mapping.
- **EV (Extended Version)** : uses multiple local sub-centers per
  cluster, capturing within-cluster structure more faithfully.

## Installation

``` r

install.packages("devtools")
devtools::install_github("lytgysrn/hyperclusvis")
```

## Quick start

### OV: single-centroid mapping

``` r

library(hyperclusvis)

set.seed(42)
out <- hyperbolize(
  iris[, 1:4],
  numk       = 3,
  preprocess = "kmeansWrap",
  scale      = "std"
)
plot_disk(out, 3)
```

### EV: multi-anchor mapping

``` r

set.seed(42)
out <- hyperbolizeEV(
  iris[, 1:4],
  numk         = 3,
  preprocess   = "kmeansWrap",
  scale        = "std",
  kg_per_class = 5,
  spread_scale = 1.0,
  lambda       = 0
)
plot_disk(out, 3)
```

### Using your own cluster labels

Both
[`hyperbolize()`](https://lytgysrn.github.io/hyperclusvis/reference/hyperbolize.md)
and
[`hyperbolizeEV()`](https://lytgysrn.github.io/hyperclusvis/reference/hyperbolizeEV.md)
accept a pre-computed label vector (or factor) through the `preprocess`
argument instead of a method name:

``` r

labels <- my_clustering_function(data)
out <- hyperbolizeEV(data, numk = 4, preprocess = labels, scale = "std")
plot_disk(out, 4)
```

## Main functions

| Function | Description |
|----|----|
| [`hyperbolize()`](https://lytgysrn.github.io/hyperclusvis/reference/hyperbolize.md) | OV mapping — one centroid per cluster |
| [`hyperbolizeEV()`](https://lytgysrn.github.io/hyperclusvis/reference/hyperbolizeEV.md) | EV mapping — multiple local sub-centers per cluster |
| [`plot_disk()`](https://lytgysrn.github.io/hyperclusvis/reference/plot_disk.md) | Visualize the embedding on the Poincaré disk |
| [`hypdist()`](https://lytgysrn.github.io/hyperclusvis/reference/hypdist.md) | Compute the hyperbolic distance between two points in the disk |
| [`PrepareData()`](https://lytgysrn.github.io/hyperclusvis/reference/PrepareData.md) | Center and rescale data to a fixed average norm |

## Parameters

### Shared parameters

| Parameter | Description |
|----|----|
| `D` | Numeric data matrix (observations in rows, variables in columns) |
| `numk` | Number of clusters |
| `preprocess` | Clustering method (`"kmeansWrap"`, `"kmeans"`, `"kmedians"`, `"hddc"`, `"pgmm"`) or a label vector |
| `scale` | `"std"` for column standardization; `"norm"` for center-and-rescale to fixed average norm |
| `use_input_order` | If `FALSE` (default), centroids are reordered by a greedy TSP heuristic for better layout |

### EV-specific parameters

| Parameter | Paper notation | Default | Description |
|----|----|----|----|
| `kg_per_class` | *G* | 5 | Number of local sub-centers per cluster |
| `spread_scale` | δ | 1.0 | Contraction factor; values below 1 pull sub-centers toward the global centroid |
| `lambda` | λ | 1e-6 | Ridge parameter in the local projection step; larger values produce more compact layouts |

## Dependencies

`stats`, `flexclust`, `HDclassif`, `pgmm`, `ggplot2`, `RColorBrewer`,
`grDevices`, `rlang`

## Examples

Representative examples used in the manuscript are provided in
example_codes/run_hyper_vis.R.

## License

MIT
