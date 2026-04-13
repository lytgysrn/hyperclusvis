#' Plot hyperbolic embedding on the Poincaré disk
#'
#' Visualize the 2D hyperbolic embedding returned by \code{hyperbolizeEV()}
#' on the Poincaré disk, with points colored by cluster and global centroids
#' highlighted.
#'
#' @param out A list, typically the output of \code{hyperbolizeEV()} or
#'   \code{hyperbolize()}.
#' @param numk Integer. Number of clusters.
#' @param label Logical. If \code{TRUE}, show a color legend; otherwise hide it.
#'
#' @return A \code{ggplot} object representing the Poincaré disk visualization.
#' @importFrom rlang .data
#' @export
plot_disk <- function(out, numk, label = FALSE) {
  if (numk > 10) {
    get_colors <- function(numk) {
      base_cols <- RColorBrewer::brewer.pal(12, "Set3")
      grDevices::colorRampPalette(base_cols)(numk)
    }
    cluster_colors <- get_colors(numk)
  } else {
    cluster_colors <- c(
      "#E41A1C", # Red
      "#377EB8", # Blue
      "#4DAF4A", # Green
      "#984EA3", # Purple
      "#FF7F00", # Orange
      "#A65628", # Brown
      "#F781BF", # Pink
      "#999999", # Grey
      "#FFD92F", # Yellow
      "#A6CEE3"  # Light Blue
    )[seq_len(numk)]
  }

  if (nrow(out$points) > length(out$cluster)) {
    out$points <- out$points[-1, , drop = FALSE]
  }

  plot_df <- data.frame(
    x = out$points[, 1],
    y = out$points[, 2],
    cluster = factor(out$cluster)
  )

  centers_df <- data.frame(
    x = out$centers[, 1],
    y = out$centers[, 2],
    cluster = factor(seq_len(numk))
  )

  theta <- seq(0, 2 * pi, length.out = 720)
  circle_df <- data.frame(
    x = cos(theta),
    y = sin(theta)
  )

  p <- ggplot2::ggplot() +
    ggplot2::geom_path(
      data = circle_df,
      ggplot2::aes(x = .data$x, y = .data$y),
      linewidth = 0.8,
      color = "grey20"
    ) +
    ggplot2::annotate(
      "polygon",
      x = circle_df$x,
      y = circle_df$y,
      alpha = 0.05,
      fill = "grey50"
    ) +
    ggplot2::geom_point(
      data = plot_df,
      ggplot2::aes(x = .data$x, y = .data$y, color = .data$cluster),
      size = 2,
      alpha = 0.9
    ) +
    ggplot2::geom_point(
      data = centers_df,
      ggplot2::aes(x = .data$x, y = .data$y),
      shape = 8,
      size = 3.2,
      stroke = 1.1,
      color = "black"
    ) +
    ggplot2::scale_color_manual(
      values = cluster_colors,
      guide = if (isTRUE(label)) "legend" else "none"
    ) +
    ggplot2::coord_equal(
      xlim = c(-1.02, 1.02),
      ylim = c(-1.02, 1.02),
      expand = FALSE
    ) +
    ggplot2::labs(
      title = NULL,
      x = NULL,
      y = NULL
    ) +
    ggplot2::theme_void()

  return(p)
}
