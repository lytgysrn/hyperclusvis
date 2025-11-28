# Hyperbolic distance in the Poincaré disk

Computes the hyperbolic distance between two points inside the
2-dimensional Poincaré disk model.

## Usage

``` r
hypdist(u, v)
```

## Arguments

- u:

  Numeric vector of length 2. Must satisfy `||u|| < 1`.

- v:

  Numeric vector of length 2. Must satisfy `||v|| < 1`.

## Value

A numeric value giving the hyperbolic distance between `u` and `v`.

## Examples

``` r
hypdist(c(0.1, 0), c(0.2, 0.3))
#> [1] 0.6689373
```
