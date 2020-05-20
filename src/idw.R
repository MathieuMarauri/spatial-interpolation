
# Spatial interpolation with Inverse Distance Weighting. The method is applied on a simple simulated
# case study and on a more challenging one.

# Set-up --------------------------------------------------------------------------------------

library("data.table") # Fast dataset manipulation
library("ggplot2") # Data visualisations using the Grammar of Graphics
library("gstat") # Geostatistical modelling
import::from("fields", "rdist") # rdist function
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


# Simple case ---------------------------------------------------------------------------------

# Perform IDW on a simple case by computing everything with base functions.

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

# First compute the distance matrix between all points of the grid and all sample points
distance <- rdist(simulated_grid[, c("x", "y")], sample_points[, c("x", "y")])

# Define the power order
p <- 5

# Take the inverse of every element
weights <- apply(
  X = distance,
  MARGIN = 1:2,
  FUN = function(x) 1 / x^p
)

# Divide (FUN = "/") the weights by the row sum (STATS = rowSums(weights)) to every element of the
# weights matrix (x = weights) row wise (MARGIN = 1)
weights <- sweep(
  x = weights,
  MARGIN = 1,
  STATS = rowSums(weights),
  FUN = "/"
)

# Multiply the weights matrix with the z value of the sample points
pred <- weights %*% sample_points$z

# Add value of sample points
pred[is.na(pred)] <- sample_points$z

# Plot the result alongside the true data
data.table(simulated_grid[, c("x", "y")], "pred" = pred[, 1], "true" = simulated_grid$z) %>%
  merge(
    y = sample_points,
    by = c("x", "y"),
    all.x = TRUE
  ) %>%
  .[is.na(pred), pred := z] %>%
  .[, z := NULL] %>%
  melt(
    id.vars = c("x", "y"),
    variable.name = "type",
    value.name = "z"
  ) %>%
plot_map("Prediction with IDW and real data") +
  facet_grid(cols = vars(type)) +
  theme(
    strip.text = element_text(size = 12),
    panel.spacing = unit(2.5, "lines")
  )


# Choice of P ---------------------------------------------------------------------------------

# Build leave-one-out estimator of error to select p
rmse_errors <- sapply(
  X = 1:20,
  FUN = function(p) {
    pred_sample <- numeric(length = 100)
    for (i in length(pred_sample)) {
      distance <- rdist(sample_points[i, c("x", "y")], sample_points[-i, c("x", "y")])
      weights <- apply(
        X = distance,
        MARGIN = 1:2,
        FUN = function(x) 1 / x^p
      )
      weights <- sweep(
        x = weights,
        MARGIN = 1,
        STATS = rowSums(weights),
        FUN = "/"
      )
      pred_sample[i] <- weights %*% sample_points[-i, ]$z %>% as.numeric()
    }
    return(sqrt(sum((pred_sample - sample_points$z)^2) / nrow(sample_points)))
  }
)
which.min(rmse_errors)

# Visualise the evolution of the error
ggplot(
  data = data.frame(x = 1:20, y = rmse_errors),
  mapping = aes(x = x, y = y)
) +
  geom_line() +
  labs(
    x = "Power order",
    y = "Leave-one-out RMSE",
    title = "RMSE for different orders p"
  )


# Using specific functions --------------------------------------------------------------------

# Directly compute the idw using function idw from gstat package.
pred_gstat <- idw(
  formula = z ~ 1,
  locations = ~ x + y,
  data = sample_points,
  newdata = simulated_grid[, c("x", "y")],
  idp = 5
)
names(pred_gstat) <- c("x", "y", "z", "var")

data.table(simulated_grid[, c("x", "y")], "idw_base" = pred[, 1], "idw_gstat" = pred_gstat$z) %>%
  merge(
    y = sample_points,
    by = c("x", "y"),
    all.x = TRUE
  ) %>%
  .[is.na(idw_base), idw_base := z] %>%
  .[, z := NULL] %>%
  melt(
    id.vars = c("x", "y"),
    variable.name = "type",
    value.name = "z"
  ) %>%
  plot_map("Prediction with IDW and gstat IDW") +
  facet_grid(cols = vars(type)) +
  theme(
    strip.text = element_text(size = 12),
    panel.spacing = unit(2.5, "lines")
  )


# Interpolation with spatial constraints ------------------------------------------------------


