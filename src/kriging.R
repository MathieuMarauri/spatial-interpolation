
# Spatial interpolation with Inverse Distance Weighting. The method is applied on a simple simulated
# case study and on a more challenging one.

# Set-up --------------------------------------------------------------------------------------

library("data.table") # Fast dataset manipulation
import::from("fields", "rdist") # rdist function
library("ggplot2") # Data visualisations using the Grammar of Graphics
library("gstat") # Geostatistical modelling
library("ipdw") # Inverse path distance weighting
library("sf") # Spatial data manipulation
library("stringi") # String manipulation
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
      legend.key.size = unit(1, "cm")
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
point_variogram[, x := define_center(group1)]

# Compute the avegare semi-variance by distance group
empirical_variogram <- point_variogram[, .(count = .N, vario_mean = mean(semi_variance)), by = x]

# Visualisation of the empirical variogram
plot_variogram(empirical_variogram)

# Visualisation of the distribution of the semi-variogram by group of distances
plot_variogram_box(point_variogram, empirical_variogram)

# Adjust the empirical variogram by defining better classes and specifying a max.dist
# Create evenly populated classes of distance on which average semi variance will be computed
point_variogram[, group2 := cut_number(distance, 15)]
point_variogram[, x := define_center(group2)]

# Compute the avegare semi-variance by distance group
empirical_variogram <- point_variogram[, .(count = .N, vario_mean = mean(semi_variance)), by = x]

# Visualisation of the empirical variogram
plot_variogram(empirical_variogram)

# Visualisation of the distribution of the semi-variogram by group of distances
plot_variogram_box(point_variogram, empirical_variogram)

# Create evenly populated classes of distance on which average semi variance will be computed
point_variogram[distance < .8, group3 := cut_number(distance, 10)]
point_variogram[, x := define_center(group3)]

# Compute the avegare semi-variance by distance group
empirical_variogram <- point_variogram[, .(count = .N, vario_mean = mean(semi_variance)), by = x]

# Visualisation of the empirical variogram
plot_variogram(empirical_variogram[!is.na(x)])

# Visualisation of the distribution of the semi-variogram by group of distances
plot_variogram_box(point_variogram[!is.na(x)], empirical_variogram[!is.na(x)])

# Smaller groups at the beginning
point_variogram[distance < .8, group4 := cut(distance, breaks = quantile(distance, c(seq(0, .15, by = .05), seq(.2, 1, length.out = 9))))]
point_variogram[, x := define_center(group4)]

# Compute the avegare semi-variance by distance group
empirical_variogram <- point_variogram[, .(count = .N, vario_mean = mean(semi_variance)), by = x]

# Visualisation of the empirical variogram
plot_variogram(empirical_variogram[!is.na(x)])

# Visualisation of the distribution of the semi-variogram by group of distances
plot_variogram_box(point_variogram[!is.na(x)], empirical_variogram[!is.na(x)])

setorderv(empirical_variogram, "x")

### Parametric variogram

# Define the function to optimise
loss_function <- function(params, x, vario_emp, cov_model) {
  sigma <- params[1]
  phi <- params[2]
  nugget <- params[3]
  if (cov_model == "spherical") {
    vario_mod <- ifelse(
      test = x > phi,
      yes = nugget + sigma,
      no = nugget + sigma * (3 / 2 * (x / phi) - 1 / 2 * (x / phi)^3)
    )
  } else if (cov_model == "gaussian") {
    vario_mod <- ifelse(
      test = x == 0,
      yes = 0,
      no = nugget + sigma * (1 - exp(-x^2 / phi^2))
    )
  } else if (cov_model == "exponential") {
    vario_mod <- ifelse(
      test = x == 0,
      yes = 0,
      no = nugget + sigma * (1 - exp(-x / phi))
    )
  }
  loss <- sum((vario_emp - vario_mod)^2)
  return(loss)
}

# Numeric optimisation to find the best parameters
gausian_param <- optim(
  par = c(2.2, .4, 0),
  fn = loss_function,
  x =  empirical_variogram$x[!is.na(empirical_variogram$x)],
  vario_emp = empirical_variogram$vario_mean[!is.na(empirical_variogram$x)],
  cov_model = "gaussian"
)
spherical_param <- optim(
  par = c(2.2, .4, 0),
  fn = loss_function,
  x =  empirical_variogram$x[!is.na(empirical_variogram$x)],
  vario_emp = empirical_variogram$vario_mean[!is.na(empirical_variogram$x)],
  cov_model = "spherical"
)
exponential_param <- optim(
  par = c(2.2, .4, 0),
  fn = loss_function,
  x =  empirical_variogram$x[!is.na(empirical_variogram$x)],
  vario_emp = empirical_variogram$vario_mean[!is.na(empirical_variogram$x)],
  cov_model = "exponential"
)

# Visualise the parametric variogram and the empirical one
vario_fit_data <- data.frame(
  x = rep(seq(0, .8, length.out = 1000), 3),
  type = rep(c("Exponential", "Gaussian", "Spherical"), each = 1000),
  y = c(
    exponential_param$par[3] + exponential_param$par[1] * (1 - exp(-seq(0, .8, length.out = 1000) / exponential_param$par[2])),
    gausian_param$par[3] + gausian_param$par[1] * (1 - exp(-seq(0, .8, length.out = 1000)^2 / gausian_param$par[2]^2)),
    ifelse(
      test = seq(0, .8, length.out = 1000) > spherical_param$par[2],
      yes = spherical_param$par[3] + spherical_param$par[1],
      no = spherical_param$par[3] + spherical_param$par[1] * (3 / 2 * (seq(0, .8, length.out = 1000) / spherical_param$par[2]) - 1 / 2 * (seq(0, .8, length.out = 1000) / spherical_param$par[2])^3)
    )
  )
)
plot_variogram(empirical_variogram[!is.na(x)]) +
  geom_line(
    data = vario_fit_data,
    mapping = aes(x = x, y = y, color = type)
  ) +
  guides(
    colour = guide_legend(override.aes = list(size = 2))
  ) +
  labs(
    y = "Semi-variogram",
    color = "Parametric model"
  ) +
  theme(
    legend.position = c(.8, .2)
  )


# Kriging step by step ------------------------------------------------------------------------

# Pediction is made step by step for a single point first, then for the grid.

### For one point

# Coordinates of point to predict
point_pred <- c(.5, .5)

# Add row names to easily spot points
rownames(sample_points) <- paste0("s", 1:nrow(sample_points))

# Distance matrix of the points in the sample and the point to predict
matrix_dist <- as.matrix(sample_points[, c("x", "y")]) %>% # transform the coordinates into a matrix
  rbind("s_pred" = point_pred) %>% # Add the point to predict
  wordspace::dist.matrix(method = "euclidean") # compute the distance between all points

# Extract the distance between point to predict and sample points
predict_dist <- matrix_dist[1:100, 101]

# Remove distance with point to predict from matrix_dist
sample_dist <- matrix_dist[1:100, 1:100]

# Function to compute the variogram
vario_func <- function(x) {
  return(ifelse(
    test = x > spherical_param$par[2],
    yes = spherical_param$par[3] + spherical_param$par[1],
    no = spherical_param$par[3] + spherical_param$par[1] * (3 / 2 * (x / spherical_param$par[2]) - 1 / 2 * (x / spherical_param$par[2])^3)
  ))
}

# Variogram matrix
variogram_matrix <- apply(
  X = sample_dist,
  MARGIN = c(1, 2),
  FUN = function(x) vario_func(x)
)

# Add the extra row and column of 1 (0 on diag)
variogram_matrix <- rbind(variogram_matrix, 1)
variogram_matrix <- cbind(variogram_matrix, 1)
variogram_matrix[nrow(variogram_matrix), ncol(variogram_matrix)] <- 0

# Inverse variogram
variogram_matrix_inv <- solve(variogram_matrix)

# Variogram vector
variogram_vec <- vario_func(predict_dist)

# Add extra 1
variogram_vec <- c(variogram_vec, 1)

# Compute the vector of weights
weights <- variogram_matrix_inv %*% variogram_vec

# Compute the value at the point pred and compare with the one obtained with the Kriginig model
sum(weights[1:100] * sample_points$z)

### For the entire grid

# Make the prediction
pb <- progress_bar$new(total = nrow(simulated_grid), format = "[:bar] :current/:total (:percent) :elapsedfull")
grid_pred <- sapply(
  X = seq_len(nrow(simulated_grid)),
  FUN = function(i) {
    point_pred <- as.numeric(simulated_grid[i, c("x", "y")])
    matrix_dist <- as.matrix(sample_points[, c("x", "y")]) %>% # transform the coordinates into a matrix
      rbind("s_pred" = point_pred) %>% # Add the point to predict
      wordspace::dist.matrix(method = "euclidean") # compute the distance between all points
    predict_dist <- matrix_dist[1:100, 101]
    sample_dist <- matrix_dist[1:100, 1:100]
    vario_func <- function(x) {
      return(ifelse(
        test = x > spherical_param$par[2],
        yes = spherical_param$par[3] + spherical_param$par[1],
        no = spherical_param$par[3] + spherical_param$par[1] * (3 / 2 * (x / spherical_param$par[2]) - 1 / 2 * (x / spherical_param$par[2])^3)
      ))
    }
    variogram_matrix <- apply(
      X = sample_dist,
      MARGIN = c(1, 2),
      FUN = function(x) vario_func(x)
    )
    variogram_matrix <- rbind(variogram_matrix, 1)
    variogram_matrix <- cbind(variogram_matrix, 1)
    variogram_matrix[nrow(variogram_matrix), ncol(variogram_matrix)] <- 0
    variogram_matrix_inv <- solve(variogram_matrix)
    variogram_vec <- vario_func(predict_dist)
    variogram_vec <- c(variogram_vec, 1)
    weights <- variogram_matrix_inv %*% variogram_vec
    pb$tick()
    return(sum(weights[1:100] * sample_points$z))
  }
)

# Visualise the result
plot_map(cbind(simulated_grid[, c("x", "y")], "z" = grid_pred), "")




