
# Generate sample data that will be used as the base case study for the different methods.

# Set-up --------------------------------------------------------------------------------------

library("data.table") # Fast dataset manipulation
library("geoR") # Analysis of geo-statistical data
options(geoR.messages = FALSE)
library("ggplot2") # Data visualisations using the Grammar of Graphics
library("magrittr") # Pipe operators
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
      # axis.ticks = element_blank(),
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


# Simple data ---------------------------------------------------------------------------------

# Data is generated using gaussian random fields

# Map on the [0,1]x[0,1] grid
base_map <- grf(
  n = 10201,
  grid = expand.grid(x = seq(0, 1, length.out = 101), y = seq(0, 1, length.out = 101)),
  cov.model = "spherical",
  cov.pars = c(2, .4)
)
base_map <- data.table(base_map$coords, "z" = base_map$data)

# Visualise the generated map
plot_map(base_map, "Simulated data")

# Export simulated data
saveRDS(base_map, "data/grf_spheric_01x01.rds")


# Map on the [0,1]x[0,1] grid
base_map <- grf(
  n = 10201,
  grid = expand.grid(x = seq(0, 100, length.out = 101), y = seq(0, 100, length.out = 101)),
  cov.model = "gaussian",
  cov.pars = c(2, 25)
)
base_map <- data.table(base_map$coords, "z" = base_map$data)

# Visualise the generated map
plot_map(base_map, "Simulated data")

# Export simulated data
saveRDS(base_map, "data/grf_gauss_100x100.rds")


# Barrier data --------------------------------------------------------------------------------

# Generate sample data points on a region with barriers.

# Generate some barriers
barrier_1 <- rbind(c(0, .3), c(.3, .3), c(.3, .4), c(0, .4), c(0, .3))
barrier_2 <- rbind(c(.5, .4), c(.55, .4), c(.75, .75), c(.7, .75), c(.5, .4))
barrier_3 <- rbind(c(.7, .8), c(.7, 1), c(.3, 1), c(.3, .8), c(.7, .8))
barriers <- st_multipolygon(list(list(barrier_1), list(barrier_2), list(barrier_3))) %>%
  st_sfc() %>%
  st_as_sf()

# Generate some points
points <- data.frame(
  x = c(.1, .3, .5, .8, .55, .2, .8, .9, .05, .9, .1, .2, .05, .62),
  y = c(0, .1, .2, .1, .6, .6, .3, .7, .65, .9, .9, .42, .43, .7),
  z = c(-2, -1, 0, -1, 2, 2, 0, -1, 2, -2, 1, 2, 2, 1)
)

# Export the 2 objects
saveRDS(barriers, "data/barriers.rds")
saveRDS(points, "data/points.rds")




