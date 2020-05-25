
# Spatial interpolation with Triangulated Irregular Network. The method is applied on a simple simulated
# case study and on a more challenging one.

# Set-up --------------------------------------------------------------------------------------

library("data.table") # Fast dataset manipulation
import::from("fields", "rdist") # rdist function
library("ggplot2") # Data visualisations using the Grammar of Graphics
library("spatstat") # Spatial data analysis
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

# Perform a Delaunay triangulation to define the triangulated surface needed for interpolation
test <- delaunay(as.ppp(sample_points[, c("x", "y")], W = owin(xrange = c(0, 1), yrange = c(0, 1))))
plot(test)
