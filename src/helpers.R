
#'
#' This function plots a tile map
#'
plot_map <- function(map, title = NULL) {
  ggplot(
    data = map,
    mapping = aes(x = x, y = y, fill = z)
  ) +
    geom_tile() +
    guides(
      fill = guide_colorbar(title.vjust = 0.8)
    ) +
    scale_x_continuous(
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      expand = c(0, 0)
    ) +
    coord_fixed() +
    labs(
      title = title
    ) +
    theme(
      legend.position = 'bottom',
      axis.line = element_blank(),
      panel.grid = element_blank()
    )
}