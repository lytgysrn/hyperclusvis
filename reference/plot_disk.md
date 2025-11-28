# Plot hyperbolic embedding on the Poincaré disk

Visualize the 2D hyperbolic embedding returned by
[`hyperbolizeEV()`](https://lytgysrn.github.io/hyperclusvis/reference/hyperbolizeEV.md)
on the Poincaré disk, with points colored by cluster and global
centroids highlighted.

## Usage

``` r
plot_disk(out, numk, label = F)
```

## Arguments

- out:

  A list, typically the output of
  [`hyperbolizeEV()`](https://lytgysrn.github.io/hyperclusvis/reference/hyperbolizeEV.md)
  or
  [`hyperbolize()`](https://lytgysrn.github.io/hyperclusvis/reference/hyperbolize.md).

- numk:

  Integer. Number of clusters .

- label:

  Logical. If `TRUE`, show a color legend; otherwise hide it.

## Value

A `ggplot` object representing the Poincaré disk visualization.
