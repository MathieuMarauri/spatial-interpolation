
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

# Import the points and barriers that will be used.
barriers <- readRDS("data/barriers.rds")
points <- readRDS("data/points.rds")

# Interpolation with IDW
pred_idw <- idw(
  formula = z ~ 1,
  locations = ~ x + y,
  data = points,
  newdata = simulated_grid[, c("x", "y")],
  idp = 2
)
names(pred_idw) <- c("x", "y", "z", "var")

# Visualise the results
plot_map(pred_idw, "Interpolation with IDW") +
  geom_sf(
    data = barriers
  )

# First compute the distance matrix between all points of the grid and all sample points
distance <- rdist(simulated_grid[, c("x", "y")], points[, c("x", "y")])

# Define the power order
p <- 5

# Take the inverse of every element
weights <- apply(
  X = distance,
  MARGIN = 1:2,
  FUN = function(x) 1 / x^p
)

# Remove connection that are not in "line of sight". For each pair of point build a line and check
# wether it crosses a barrier
to_remove <- matrix(FALSE, nrow = nrow(weights), ncol = ncol(weights))
for (i in seq_len(nrow(simulated_grid))) {
  for (j in seq_len(nrow(points))) {
    to_remove[i, j] <- rbind(as.numeric(points[j, c("x", "y")]), as.numeric(simulated_grid[i, c("x", "y")])) %>%
      st_linestring() %>%
      st_intersects(barriers) %>%
      .[[1]] %>%
      length() > 0
  }
}
weights[to_remove] <- 0


all_pairs <- merge.data.table(
  x = as.data.table(points[, c("x", "y")])[, id := 1],
  y = as.data.table(simulated_grid[, c("x", "y")])[, id := 1],
  by = c("id"),
  suffixes = c("_point", "_grid"),
  allow.cartesian = TRUE
)

lines <- lapply(
  X = 1:nrow(all_pairs),
  FUN = function(i) {
    st_linestring(rbind(
      as.numeric(all_pairs[i, .(x_point, y_point)]),
      as.numeric(all_pairs[i, .(x_grid, y_grid)])
    ))
  }
)

a <- st_sfc(lines)
b <- st_intersects(a, barriers)
c <- lengths(b)




# Divide (FUN = "/") the weights by the row sum (STATS = rowSums(weights)) to every element of the
# weights matrix (x = weights) row wise (MARGIN = 1)
weights <- sweep(
  x = weights,
  MARGIN = 1,
  STATS = rowSums(weights),
  FUN = "/"
)

# Multiply the weights matrix with the z value of the sample points
pred <- weights %*% points$z

# Plot the result alongside the IDW prediction
data.table(simulated_grid[, c("x", "y")], "pred" = pred[, 1], "idw" = pred_idw$z) %>%
  merge(
    y = points,
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
  plot_map("Prediction with IDW and constrained IDW") +
  geom_sf(data = barriers) +
  facet_grid(cols = vars(type)) +
  theme(
    strip.text = element_text(size = 12),
    panel.spacing = unit(2.5, "lines")
  )
