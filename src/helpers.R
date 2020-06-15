
#'
#' This function plots a tile map
#'
plot_map <- function(map, title = NULL, var = "z") {
  ggplot() +
    geom_tile(
      data = map,
      mapping = aes_string(x = "x", y = "y", fill = var)
    ) +
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

#'
#' This function plots the empirical variaogram.
#'
plot_variogram <- function(empirical_variogram) {
  ggplot() +
    geom_point(
      data = empirical_variogram,
      mapping = aes(x = x, y = vario_mean)
    ) +
    labs(
      x = "Distance",
      y = latex2exp::TeX("$\\frac{1}{2N_h}\\sum_{i = 1}^{N_h}(Z(s_{i+h}) - Z(s_i))^2$"),
      title = "Empirical variogram"
    )
}

#'
#' This function plots the distribution of the semi-variance as boxplot for the different groups of
#' distances.
#'
plot_variogram_box <- function(point_variogram, empirical_variogram) {
  ggplot() +
    geom_boxplot(
      data = point_variogram,
      mapping = aes(x = as.character(x), y = semi_variance),
      varwidth = TRUE
    ) +
    geom_jitter(
      width = 0.2
    ) +
    geom_label(
      data = empirical_variogram,
      mapping = aes(x = as.character(x), y = -.8, label = count)
    ) +
    labs(
      x = "Distance",
      y = latex2exp::TeX("$\\frac{(Z(s_i) - Z(s_j))^2}{2}$"),
      title = "Distribution of the semi-variance"
    )
}


#'
#' This function extracts the center of class obtained from cut function.
#'
define_center <- function(x) {
  x %>%
    as.character() %>%
    stri_replace_all_regex("(\\(|\\[|\\])", "") %>%
    stri_split_fixed(",") %>%
    sapply(
      FUN = function(x) sum(as.numeric(x)) / 2
    )
}
