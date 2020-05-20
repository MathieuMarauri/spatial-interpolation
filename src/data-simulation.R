
# Generate sample data that will be used as the base case study for the different methods.

# Set-up --------------------------------------------------------------------------------------

library("data.table") # Fast dataset manipulation
library("geoR") # Analysis of geostatistical data
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


# Generation ----------------------------------------------------------------------------------

# Data is generated using gaussian random fileds

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
