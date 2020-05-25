
# Spatial interpolation with Inverse Distance Weighting. The method is applied on a simple simulated
# case study and on a more challenging one.

# Set-up --------------------------------------------------------------------------------------

library("data.table") # Fast dataset manipulation
import::from("fields", "rdist") # rdist function
library("ggplot2") # Data visualisations using the Grammar of Graphics
library("gstat") # Geostatistical modelling
library("ipdw") # Inverse path distance weighting
library("sf") # Spatial data manipulation
source("src/helpers.R") # Load custom functions

# Set default ggplot theme
theme_set(
  theme_light(
    base_size = 20
  ) +
    theme(
      text = element_text(family = "Gibson", colour = "gray10"),
      panel.border = element_blank(),
      axis.line = element_line(colour = "gray50", size = .5),
      axis.ticks = element_blank(),
      strip.background = element_rect(colour = "gray50", fill = "transparent", size = .7),
      strip.text.x = element_text(colour = "gray10"),
      strip.text.y = element_text(colour = "gray10"),
      legend.key.size = unit(1.5, "cm")
    )
)

# Set default scales
scale_colour_continuous <- function(...) ggplot2::scale_colour_viridis_c(..., option = "viridis")
scale_colour_discrete <- function(...) ggplot2::scale_colour_viridis_d(..., option = "viridis")
scale_fill_continuous <- function(...) ggplot2::scale_fill_viridis_c(..., option = "viridis")
scale_fill_discrete <- function(...) ggplot2::scale_fill_viridis_d(..., option = "viridis")


# Data import ---------------------------------------------------------------------------------

# Import toy data for the simple case and the more advanced ones.

# Load the map
simulated_grid <- readRDS("data/grf_spheric_01x01.rds")

# Sample 100 points
set.seed(123)
sample_points <- simulated_grid[sample(1:10201, 100), ]

# Visualise the points and the map
plot_map(simulated_grid, "The 100 sample points on the true map") +
  geom_point(
    data = sample_points,
    mapping = aes(x = x, y = y)
  )


# Variogram step by step ----------------------------------------------------------------------

# Build the variogram. First define the empirical variogram by varying the distances classes then
# adjust a paremeter variogram.

### Distance matrix

# Distance matrix between all points
dist_matrix <- dist(sample_points[, .(x, y)])

# Coerce matrix to long table so it is easier to plot later on (with the variogram)
dist_matrix <- dist_matrix %>%
  as.matrix() %>%
  as.data.table() %>%
  .[, point_a := 1:nrow(sample_points)] %>%
  melt(
    id.vars = "point_a",
    variable.name = "point_b",
    value.name = "distance"
  ) %>%
  .[, point_b := as.numeric(point_b)]

### Semi-variance

# Compute the variogram (variance of z(si) - z(sj) since mean is supposed to be the same by
# stationarity)
vario_matrix <- outer(
  X = sample_points$z,
  Y = sample_points$z,
  FUN = function(x, y) (x - y)^2 / 2
)

# Coerce to long table
vario_matrix <- vario_matrix %>%
  as.data.table() %>%
  setNames(as.character(1:nrow(sample_points))) %>%
  .[, point_a := 1:nrow(sample_points)] %>%
  melt(
    id.vars = "point_a",
    variable.name = "point_b",
    value.name = "semi_variance"
  ) %>%
  .[, point_b := as.numeric(point_b)]

# Build one table from distance and variogram
point_variogram <- merge(
  x = dist_matrix,
  y = vario_matrix,
  by = c("point_a", "point_b"),
  all = TRUE,
  sort = FALSE
)

# Visualisation of the variogram scatter plot
ggplot(
  data = point_variogram,
  mapping = aes(x = distance, y = semi_variance)
) +
  geom_point(
    alpha = .3
  ) +
  labs(
    x = "Distance",
    y = latex2exp::TeX("$\\frac{(Z(s_i) - Z(s_j))^2}{2}$"),
    title = "Semi-variance by pair of points"
  )

### Empirical variogram

# Create classes of distance on which average semi variance will be computed
point_variogram[, group1 := cut(distance, 15)]

# Compute the avegare semi-variance by distance group
empirical_variogram <- point_variogram[, .(count = .N, vario_mean = mean(semi_variance)), by = group1]

# Visualisation of the empirical variogram
plot_variogram(empirical_variogram, "group1")

# Visualisation of the distribution of the semi-variogram by group of distances
plot_variogram_box(point_variogram, empirical_variogram, "group1")

# Adjust the empirical variogram by defining better classes and specifying a max.dist
# Create evenly populated classes of distance on which average semi variance will be computed
point_variogram[, group2 := cut_number(distance, 15)]

# Compute the avegare semi-variance by distance group
empirical_variogram <- point_variogram[, .(count = .N, vario_mean = mean(semi_variance)), by = group2]

# Visualisation of the empirical variogram
plot_variogram(empirical_variogram, "group2")

# Visualisation of the distribution of the semi-variogram by group of distances
plot_variogram_box(point_variogram, empirical_variogram, "group2")

# Create evenly populated classes of distance on which average semi variance will be computed
point_variogram[distance < .8, group3 := cut_number(distance, 10)]

# Compute the avegare semi-variance by distance group
empirical_variogram <- point_variogram[, .(count = .N, vario_mean = mean(semi_variance)), by = group3]

# Visualisation of the empirical variogram
plot_variogram(empirical_variogram[!is.na(group3)], "group3")

# Visualisation of the distribution of the semi-variogram by group of distances
plot_variogram_box(point_variogram[!is.na(group3)], empirical_variogram[!is.na(group3)], "group3")

# Smaller groups at the beginning
point_variogram[distance < .8, group4 := cut(distance, breaks = quantile(distance, c(seq(0, .15, by = .05), seq(.2, 1, length.out = 9))))]

# Compute the avegare semi-variance by distance group
empirical_variogram <- point_variogram[, .(count = .N, vario_mean = mean(semi_variance)), by = group4]

# Visualisation of the empirical variogram
plot_variogram(empirical_variogram[!is.na(group4)], "group4")

# Visualisation of the distribution of the semi-variogram by group of distances
plot_variogram_box(point_variogram[!is.na(group4)], empirical_variogram[!is.na(group4)], "group4")

### Parametric variogram

#

# Kriging step by step ------------------------------------------------------------------------


# Coordinates of point to predict
point_pred <- c(.5, .5)

# Add row names to easily spot points
rownames(sample_points) <- paste0("s", 1:nrow(sample_points))

# Distance matrix of the points in the sample and the point to predict
matrix_dist <- as.matrix(sample_points[, 1:2]) %>% # transform the coordinates into a matrix
  rbind("s_pred" = point_pred) %>% # Add the point to predict
  wordspace::dist.matrix(method = "euclidean") # compute the distance between all points

# Extract the distance between point to predict and sample point
predict_dist <- matrix_dist[1:100, 101]

# Remove distance with point to predict from matrix_dist
echantillon_dist <- matrix_dist[1:100, 1:100]

# Function to compute the variogram
vario_func <- function(h, vario_fit) {
  return(vario_fit$nugget + vario_fit$cov.pars[1] * (1 - exp(-h / vario_fit$cov.pars[2])))
}

# Variogram matrix
variogram_matrix <- apply(
  X = echantillon_dist,
  MARGIN = c(1, 2),
  FUN = function(x) vario_func(x, vario_fit_exp)
)

# Add the extra row and column of 1 (0 on diag)
variogram_matrix <- rbind(variogram_matrix, 1)
variogram_matrix <- cbind(variogram_matrix, 1)
variogram_matrix[nrow(variogram_matrix), ncol(variogram_matrix)] <- 0

# Inverse variogram
variogram_matrix_inv <- solve(variogram_matrix)

# Variogram vector
variogram_vec <- vario_func(predict_dist, vario_fit_exp)

# Add extra 1
variogram_vec <- c(variogram_vec, 1)

# Compute the vector of weights
weights <- variogram_matrix_inv %*% variogram_vec

# Compute the value at the point pred and compare with the one obtained with the Kriginig model
sum(weights[1:100] * sample_points$z)
krige$predict[point_pred[1] + 1 + point_pred[2] * 101]

# Compute the variance
weights[101] + sum(weights[1:100] * variogram_vec[1:100])
krige$krige.var[point_pred[1] + 1 + point_pred[2] * 101]

