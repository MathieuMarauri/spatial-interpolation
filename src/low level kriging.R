
# Objective: perform Kriging on simulated data with lox level functions and then kriging on raster
# data with appropriate functions to evaluate performance.

# Set-up --------------------------------------------------------------------------------------

library("data.table") # dataset manipulation
# library("dplyr") # dataset manipulation
library("geoR") # geostatistical analysis
library("ggplot2") # data visualisation
library("tidyr") # data tidying

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
      legend.key.size = unit(2, "cm")
    )
)

# Set default scales
scale_colour_continuous <- function(...) ggplot2::scale_colour_viridis_c(..., option = "viridis")
scale_colour_discrete <- function(...) ggplot2::scale_colour_viridis_d(..., option = "viridis")
scale_fill_continuous <- function(...) ggplot2::scale_fill_viridis_c(..., option = "viridis")
scale_fill_discrete <- function(...) ggplot2::scale_fill_viridis_d(..., option = "viridis")


# Simulated data ------------------------------------------------------------------------------

# Sample points from a grid 100 x 100.

# Import des données
grid <- fread("data/simulation.txt")

# Visualisation of the data simulated
ggplot(
  data = grid,
  mapping = aes(x = x, y = y, fill = z)
) +
  geom_tile() +
  guides(
    fill = guide_colorbar(title.vjust = 0.8)
  ) +
  labs(
    title = "Initial data"
  ) +
  theme_void() +
  theme(
    legend.position = 'bottom'
  )

# Sample 100 points
set.seed(1234)
sample_points <- grid[sample(1:nrow(grid), 100)]

# Visualisation of the data sampled
ggplot(
  data = simulation,
  mapping = aes(x = x, y = y, fill = z)
) +
  geom_tile() +
  geom_point(
    data = sample_points,
    mapping = aes(x = x, y = y)
  ) +
  guides(
    fill = guide_colorbar(title.vjust = 0.8)
  ) +
  labs(
    title = "Sampled points on the initial map"
  ) +
  theme_void() +
  theme(
    legend.position = 'bottom'
  )





# 1.3 Variogramme -----------------------------------------------------------------------------

# Le but est de définir le variogramme empirique en faisant varier différents réglages et d'ajuster
# un modèle de variogramme paramétrique.

# Arguments par défaut (max.dist est la distance maximale entre 2 points, pas de pépite) et sortie
# sous forme de nuage de points (n*(n-1)/2 points qui correspondent à chaque paires de points)
variogram_c <- variog(geodata, option = "cloud")
ggplot(
  data = data.frame(x = variogram_c$u, y = variogram_c$v),
  mapping = aes(x = x, y = y)) +
  geom_point() +
  labs(
    x = "Distance entre les points",
    y = "Semi variance",
    title = "Nuage de points variographique"
  )

# Arguments par défaut (max.dist est la distance maximale entre 2 points, pas de pépite, 13
# classes) et sortie sous forme de moyenne par classe
variogram <- variog(geodata)
ggplot(
  data = data.frame(x = variogram$u, y = variogram$v),
  mapping = aes(x = x, y = y)) +
  geom_point() +
  labs(
    x = "Distance entre les points",
    y = "Semi variance",
    title = "Variogramme empirique"
  )

# La distance maximale est fixée à 70 et le nombre de classe à 8 (breaks donne les limites)
variogram <- variog(
  geodata = geodata, # les points sur lesquels calculer le variogramme
  max.dist = 70, # la distance maximale jusqu'à laquelle des paires de points sont considérées
  breaks = seq(0, 80, 10) # les limites des classes de distance
)
ggplot(
  data = data.frame(x = variogram$u, y = variogram$v),
  mapping = aes(x = x, y = y)) +
  geom_point() +
  labs(
    x = "Distance entre les points",
    y = "Semi variance",
    title = "Variogramme empirique"
  )

# Les classes sont modifiées, plus resserrées pour les faibles distance puis plus larges
variogram <- variog(
  geodata = geodata,
  max.dist = 70,
  breaks = c(seq(0, 30, 5), seq(30, 70, 10))
)
ggplot(
  data = data.frame(x = variogram$u, y = variogram$v),
  mapping = aes(x = x, y = y)) +
  geom_point() +
  labs(
    x = "Distance entre les points",
    y = "Semi variance",
    title = "Variogramme empirique"
  )

# Visualisation de la distribution de la semi variance par classe
variogram_c <- variog(
  geodata = geodata,
  max.dist = 70,
  breaks = seq(0, 80, 10),
  bin.cloud = TRUE # garder chaque valeur indivudelle par classe
)
ggplot(
  data = data.frame(x = as.factor(rep(variogram_c$u, variogram_c$n)), y = unlist(variogram_c$bin.cloud)),
  mapping = aes(x = x, y = y)) +
  geom_boxplot() +
  labs(
    x = "Centre des classes",
    y = "Semi variance",
    title = "Box-plot sur le variogramme empirique"
  )

# Fit d'un modèle exponentiel sur le variogram
vario_fit_exp <- variofit(
  vario = variogram, # le variogramme empirique
  cov.model = "exponential", # la nature du modèle qu'on souhaite ajuster
  ini.cov.pars = c(1.5, 30) # les paramtères initiaux nécessaires pour l'optimisation (l'ajustement) sigma et a dans le cours
)
vario_fit_exp
vario_fit_data <- data.frame(
  x = seq(0, 70, length.out = 1000),
  y = vario_fit_exp$nugget + vario_fit_exp$cov.pars[1] * (1 - exp(-seq(0, 70, length.out = 1000) / vario_fit_exp$cov.pars[2]))
)
ggplot(
  data = data.frame(x = variogram$u, y = variogram$v),
  mapping = aes(x = x, y = y)
) +
  geom_point() +
  geom_line(
    data = vario_fit_data,
    mapping = aes(x = x, y = y)
  ) +
  labs(
    x = "Distance entre les points",
    y = "Semi variance",
    title = "Variogramme empirique et exponentiel"
  )

# Fit d'un modèle exponentiel sur le variogram avec pépite
vario_fit_exp_nug <- variofit(
  vario = variogram, # le variogramme empirique
  cov.model = "exponential", # la nature du modèle qu'on souhaite ajuster
  ini.cov.pars = c(1.5, 30), # les paramtères initiaux nécessaires pour l'optimisation (l'ajustement) sigma et a dans le cours
  nugget = 0.1, # valeur initiale de la pépipte
  fix.nugget = TRUE # la pépite n'est pas évaluée, la valeur intiale est conservée
)
vario_fit_exp_nug
vario_fit_data <- data.frame(
  x = seq(0, 70, length.out = 1000),
  y = vario_fit_exp$nugget + vario_fit_exp$cov.pars[1] * (1 - exp(-seq(0, 70, length.out = 1000) / vario_fit_exp$cov.pars[2]))
)
ggplot(
  data = data.frame(x = variogram$u, y = variogram$v),
  mapping = aes(x = x, y = y)
) +
  geom_point() +
  geom_line(
    data = vario_fit_data,
    mapping = aes(x = x, y = y)
  ) +
  labs(
    x = "Distance entre les points",
    y = "Semi variance",
    title = "Variogramme empirique et exponentiel"
  )

# Fit d'un modèle gaussien sur le variogram
vario_fit_norm <- variofit(
  vario = variogram, # le variogramme empirique
  cov.model = "gaussian", # la nature du modèle qu'on souhaite ajuster
  ini.cov.pars = c(1.5, 30) # les paramtères initiaux nécessaires pour l'optimisation (l'ajustement) sigma et a dans le cours
)
vario_fit_norm
vario_fit_data <- data.frame(
  x = seq(0, 70, length.out = 1000),
  y = vario_fit_exp$nugget + vario_fit_exp$cov.pars[1] * (1 - exp(-seq(0, 70, length.out = 1000)^2 / vario_fit_exp$cov.pars[2]^2))
)
ggplot(
  data = data.frame(x = variogram$u, y = variogram$v),
  mapping = aes(x = x, y = y)
) +
  geom_point() +
  geom_line(
    data = vario_fit_data,
    mapping = aes(x = x, y = y)
  ) +
  labs(
    x = "Distance entre les points",
    y = "Semi variance",
    title = "Variogramme empirique et gaussien"
  )

# Visualisation des 2 modèles
vario_fit_data <- data.frame(
  x = rep(seq(0, 70, length.out = 1000), 2),
  type = rep(c("Exponentiel", "Gaussien"), each = 1000),
  y = c(vario_fit_exp$nugget + vario_fit_exp$cov.pars[1] * (1 - exp(-seq(0, 70, length.out = 1000) / vario_fit_exp$cov.pars[2])),
        vario_fit_exp$nugget + vario_fit_exp$cov.pars[1] * (1 - exp(-seq(0, 70, length.out = 1000)^2 / vario_fit_exp$cov.pars[2]^2)))
)
ggplot(
  data = data.frame(x = variogram$u, y = variogram$v),
  mapping = aes(x = x, y = y)
) +
  geom_point() +
  geom_line(
    data = vario_fit_data,
    mapping = aes(x = x, y = y, col = type)
  ) +
  scale_color_manual(
    name = "Modèle",
    values = c("Exponentiel" = "#3D6CE8", "Gaussien" = "#EA619D")
  ) +
  labs(
    x = "Distance entre les points",
    y = "Semi variance",
    title = "Variogramme empirique et paramétrique"
  ) +
  theme(
    legend.position = c(.8, .25),
    legend.key.size = unit(1.5, 'lines')
  )


# 1.4 Krigeage --------------------------------------------------------------------------------

# Définition des points pour lesquels la prédiction sera faite
grille <- expand.grid(x = 0:100, y = 0:100)
ggplot(
  data = grille,
  mapping = aes(x = x, y = y)
) +
  geom_point(
    size = .1
  )

# Construction du modèle de Kriging avec le variogramme gaussien
krige_control <- krige.control(
  type.krige = "ok", # type de modèle, ordinary, simple, universal
  obj.model = vario_fit_norm # variogramme utilisé
)
krige <- krige.conv(
  geodata = geodata, # observations
  locations = grille, # points pour lesquels faire une prédiction
  krige = krige_control # paramètres de contrôle
)

# Erreur de prédiction
sum(abs(krige$predict - simulation$z))

# Visualisation des prédictions
cbind(grille, Estimation = krige$predict, Réalité = simulation$z) %>%
  pivot_longer(cols = Estimation:Réalité, names_to = "type", values_to = "Z") %>%
  ggplot(
    mapping = aes(x = x, y = y, fill = Z)
  ) +
  geom_tile() +
  geom_tile(
    data = echantillon,
    mapping = aes(x = x, y = y),
    fill = "black"
  ) +
  guides(
    fill = guide_colorbar(title.vjust = 0.8)
  ) +
  facet_grid(cols = vars(type)) +
  labs(
    title = paste0("Comparaison modèle Krigeage et réalité : erreur absolue = ", sum(abs(krige$predict - simulation$z)))
  ) +
  scale_x_continuous(
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    expand = c(0, 0)
  ) +
  theme_void() +
  theme(
    legend.position = 'bottom',
    strip.text = element_text(size = 12)
  )

# Construction du modèle de Kriging avec le variogramme exponentiel
krige_control <- krige.control(
  type.krige = "ok", # type de modèle, ordinary, simple, universal
  obj.model = vario_fit_exp # variogramme utilisé
)
krige <- krige.conv(
  geodata = geodata, # observations
  locations = grille, # points pour lesquels faire une prédiction
  krige = krige_control # paramètres de contrôle
)

# Erreur de prédiction
sum(abs(krige$predict - simulation$z))

# Visualisation des prédictions
cbind(grille, Estimation = krige$predict, Réalité = simulation$z) %>%
  pivot_longer(cols = Estimation:Réalité, names_to = "type", values_to = "Z") %>%
  ggplot(
    mapping = aes(x = x, y = y, fill = Z)
  ) +
  geom_tile() +
  geom_tile(
    data = echantillon,
    mapping = aes(x = x, y = y),
    fill = "black"
  ) +
  guides(
    fill = guide_colorbar(title.vjust = 0.8)
  ) +
  facet_grid(cols = vars(type)) +
  labs(
    title = "Comparaison modèle Krigeage et réalité"
  ) +
  theme_void() +
  theme(
    legend.position = 'bottom'
  )

# Visualisation des écart-types
ggplot(
  data = cbind(grille, val = sqrt(krige$krige.var)),
  mapping = aes(x = x, y = y, fill = val)
) +
  geom_tile() +
  geom_tile(
    data = echantillon,
    mapping = aes(x = x, y = y),
    fill = "grey20"
  ) +
  labs(
    title = "Ecart-types avec Krigeage sur variogramme gaussien"
  ) +
  theme(
    panel.grid = element_blank()
  )

# Visualisation de la carte initiale
ggplot(
  data = cbind(grille, val = simulation$z),
  mapping = aes(x = x, y = y, fill = val)
) +
  geom_tile() +
  geom_tile(
    data = echantillon,
    mapping = aes(x = x, y = y),
    fill = "grey20"
  ) +
  labs(
    title = "Données initiales"
  ) +
  theme(
    panel.grid = element_blank()
  )


# 1.5 Simulations conditionnelles -------------------------------------------------------------

# Définition de la grille pour laquelle les simulations seront faites
grille <- expand.grid(x = seq(0, 100, 2), y = seq(0, 100, 2))

# Modèle de Kriging avec simulations pour obtenir les probas
out_control <- output.control(
  n.predictive = 100, # nombre de simulations
  simul = TRUE,
  threshold = 2 # seuil avec lequel calculer les porbabilités de non dépassement
)
krige_prob <- krige.conv(
  geodata = geodata, # observations
  locations = grille, # points pour lesquels faire une prédiction
  krige = krige_control, # paramètres de contrôle
  output = out_control # paramètres de contrôle de la sortie, nécessaire pour avoir toutes les simulations
)

# Visualisation de quelques simulations
cbind(grille, krige_prob$simulations[, sample(1:100, 4)]) %>%
  pivot_longer(cols = -x:-y, names_to = "type", values_to = "Z") %>%
  ggplot(
    mapping = aes(x = x, y = y, fill = Z)
  ) +
  geom_tile() +
  geom_tile(
    data = echantillon,
    mapping = aes(x = x, y = y),
    fill = "black"
  ) +
  guides(
    fill = guide_colorbar(title.vjust = 0.8)
  ) +
  facet_wrap(facets = vars(type), nrow = 2) +
  labs(
    title = "Exemples de simulations"
  ) +
  theme_void() +
  theme(
    legend.position = 'bottom'
  )

# Visualisation de la carte de probabilité que la valeur soit strictement supérieure au seuil de 2
ggplot(
  data = cbind(grille, val = 1 - krige_prob$probabilities.simulations),
  mapping = aes(x = x, y = y, fill = val)
) +
  geom_tile() +
  geom_tile(
    data = echantillon,
    mapping = aes(x = x, y = y),
    fill = "grey20"
  ) +
  labs(
    title = "Probabilité de dépassement d'un seuil",
    fill = "P(Z>2)"
  ) +
  theme(
    panel.grid = element_blank()
  )

# Nettoyer l'environnement
rm(list = setdiff(ls(), ls(pattern = "^scale_")))


# Pollution de l'air --------------------------------------------------------------------------

### Import des données

# Importer les données
donnees <- read.table("data/stationsKm.txt", header = TRUE)
grille <- read.table("data/grilleKm.txt", header = TRUE)
geodata <- as.geodata(donnees)

# Aperçu des données
ggplot(
  data = grille,
  mapping = aes(x = x, y = y, fill = z)
) +
  geom_tile() +
  labs(
    title = "Modèle initial"
  ) +
  theme(
    panel.grid = element_blank()
  )
ggplot(
  data = donnees,
  mapping = aes(x = z)
) +
  geom_histogram(bins = 15)


### Variogramme empirique et paramétrique

# Nuage variogrpahique
variogram_c <- variog(geodata, option = "cloud")
ggplot(
  data = data.frame(x = variogram_c$u, y = variogram_c$v),
  mapping = aes(x = x, y = y)) +
  geom_point() +
  labs(
    x = "Distance entre les points",
    y = "Semi variance",
    title = "Nuage de points variographique"
  )

# Création du variogramme
max_dist <- 100
variogram <- variog(
  geodata = geodata, # les points sur lesquels calculer le variogramme
  max.dist = max_dist, # la distance maximale jusqu'à laquelle des paires de points sont considérées
  breaks = seq(10, max_dist, by = 20), # les limites des classes de distance
  pairs.min = 2, # nombre minimal de paires par classe
)

# Visualiser le variogramme
ggplot(
  data = data.frame(x = variogram$u, y = variogram$v),
  mapping = aes(x = x, y = y)) +
  geom_point() +
  labs(
    x = "Distance entre les points",
    y = "Semi variance",
    title = "Variogramme empirique"
  )

# Ajuster le variogramme
vario_fit <- variofit(
  vario = variogram,
  cov.model = "gaussian",
  nugget = 10,
  fix.nugget = TRUE
)
vario_fit_data <- data.frame(
  x = seq(0, max_dist, length.out = 1000),
  y = vario_fit$nugget + vario_fit$cov.pars[1] * (1 - exp(-seq(0, max_dist, length.out = 1000)^2 / vario_fit$cov.pars[2]^2))
)
ggplot(
  data = data.frame(x = variogram$u, y = variogram$v),
  mapping = aes(x = x, y = y)
) +
  geom_point() +
  geom_line(
    data = vario_fit_data,
    mapping = aes(x = x, y = y)
  ) +
  labs(
    x = "Distance entre les points",
    y = "Semi variance",
    title = "Variogramme empirique et gaussien",
    subtitle = paste0("portee = ", round(vario_fit$cov.pars[2], 2),
                      ", palier = ", round(vario_fit$cov.pars[1], 2))
  )


### Krigeage

# Modèle de krigeage à partir du variogramme précédent pour une prédiction en chaque point de la grille
krige <- krige.conv(
  geodata = geodata,
  locations = grille[, 1:2],
  krige = krige.control(type.krige = "ok", obj.model = vario_fit)
)
cbind(grille, val = krige$predict) %>%
  ggplot(
    mapping = aes(x = x, y = y)
  ) +
  geom_tile(
    mapping = aes(fill = val)
  ) +
  geom_contour(
    mapping = aes(z = val),
    color = "black",
    binwidth = 50
  ) +
  metR::geom_text_contour(
    mapping = aes(z = val),
    binwidth = 50,
    stroke = .2
  ) +
  geom_point(
    data = donnees,
    mapping = aes(x = x, y = y),
    color = "black"
  ) +
  guides(
    fill = guide_colorbar(title.vjust = 0.8)
  ) +
  labs(
    title = "Modèle de Krigeage"
  ) +
  theme_void() +
  theme(
    legend.position = 'bottom'
  )

# Visualisation des écart-types
cbind(grille, val = sqrt(krige$krige.var)) %>%
  ggplot(
    mapping = aes(x = x, y = y)
  ) +
  geom_tile(
    mapping = aes(fill = val)
  ) +
  geom_contour(
    mapping = aes(z = val),
    color = "black",
  ) +
  metR::geom_text_contour(
    mapping = aes(z = val),
    stroke = .2
  ) +
  geom_point(
    data = donnees,
    mapping = aes(x = x, y = y),
    color = "black"
  ) +
  guides(
    fill = guide_colorbar(title.vjust = 0.8)
  ) +
  labs(
    title = "Ecart-types du modèle"
  ) +
  theme_void() +
  theme(
    legend.position = 'bottom'
  )

### Comparaison modèle / Krigeage

# Visualisation des 2 cartes
cbind(grille, Krigeage = krige$predict) %>%
  setNames(c("x", "y", "Modèle", "Krigeage")) %>%
  pivot_longer(
    -x:-y,
    names_to = "type",
    values_to = "val"
  ) %>%
  ggplot(
    mapping = aes(x = x, y = y)
  ) +
  geom_tile(
    mapping = aes(fill = val)
  ) +
  geom_contour(
    mapping = aes(z = val),
    color = "black",
    binwidth = 50
  ) +
  metR::geom_text_contour(
    mapping = aes(z = val),
    binwidth = 50,
    stroke = .2
  ) +
  facet_wrap(vars(type)) +
  guides(
    fill = guide_colorbar(title.vjust = 0.8)
  ) +
  labs(
    title = "Comparaison modèle et Krigeage"
  ) +
  scale_x_continuous(
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    expand = c(0, 0)
  ) +
  theme_void() +
  theme(
    legend.position = 'bottom',
    strip.text = element_text(size = 12)
  )

# Comparaison des valeurs mesurées et modélisées aux stations
ggplot(
  data = donnees,
  mapping = aes(x = x, y = y, color = mod - z)
) +
  geom_point(
    size = 5
  ) +
  labs(
    title = "Différence entre le modèle et les mesures",
    color = "Modèle - observations"
  ) +
  guides(
    color = guide_colorbar(title.vjust = 0.8)
  ) +
  # theme_void() +
  theme(
    legend.position = 'bottom'
  )

### Correction du modèle

# Construction des observations comme la différence entre le modéle est la réalité
zerr <- z - donnees[, 4]
geodata <- donnees %>%
  transmute(x, y, err = z - mod) %>%
  as.geodata()

# Création du variogramme
max_dist <- 100
variogram <- variog(
  geodata = geodata, # les points sur lesquels calculer le variogramme
  max.dist = max_dist, # la distance maximale jusqu'à laquelle des paires de points sont considérées
  breaks = seq(10, max_dist, by = 20), # les limites des classes de distance
  pairs.min = 2, # nombre minimal de paires par classe
  messages.screen = FALSE
)

# Visualiser le variogramme
ggplot(
  data = data.frame(x = variogram$u, y = variogram$v),
  mapping = aes(x = x, y = y)) +
  geom_point() +
  labs(
    x = "Distance entre les points",
    y = "Semi variance",
    title = "Variogramme empirique"
  )

# Ajustement du variogramme
vario_fit <- variofit(
  vario = variogram,
  cov.model = "exponential",
  ini.cov.pars = c(500, 50),
  fix.nugget = TRUE
)
vario_fit_data <- data.frame(
  x = seq(0, max_dist, length.out = 1000),
  y = vario_fit$nugget + vario_fit$cov.pars[1] * (1 - exp(-seq(0, max_dist, length.out = 1000) / vario_fit$cov.pars[2]))
)
ggplot(
  data = data.frame(x = variogram$u, y = variogram$v),
  mapping = aes(x = x, y = y)
) +
  geom_point() +
  geom_line(
    data = vario_fit_data,
    mapping = aes(x = x, y = y)
  ) +
  labs(
    x = "Distance entre les points",
    y = "Semi variance",
    title = "Variogramme empirique et exponentiel",
    subtitle = paste0("portee = ", round(vario_fit$cov.pars[2], 2),
                      ", palier = ", round(vario_fit$cov.pars[1], 2))
  )

# Construction du modèle de Krigeage sur les erreurs de modélisation
krige <- krige.conv(
  geodata = geodata,
  locations = grille[, 1:2],
  krige = krige.control(type.krige = "ok", obj.model = vario_fit)
)

# Visualisation de la carte obtenue après correction du modèle
cbind(grille, val = grille$z + krige$predict) %>%
  ggplot(
    mapping = aes(x = x, y = y)
  ) +
  geom_tile(
    mapping = aes(fill = val)
  ) +
  geom_contour(
    mapping = aes(z = val),
    color = "black",
  ) +
  metR::geom_text_contour(
    mapping = aes(z = val),
    stroke = .2
  ) +
  geom_point(
    data = donnees,
    mapping = aes(x = x, y = y),
    color = "black"
  ) +
  guides(
    fill = guide_colorbar(title.vjust = 0.8)
  ) +
  labs(
    title = "Modéle corrigé"
  ) +
  theme_void() +
  theme(
    legend.position = 'bottom'
  )

# Visualisation des écart-types du modèle
cbind(grille, val = krige$krige.var) %>%
  ggplot(
    mapping = aes(x = x, y = y)
  ) +
  geom_tile(
    mapping = aes(fill = val)
  ) +
  geom_contour(
    mapping = aes(z = val),
    color = "black",
  ) +
  metR::geom_text_contour(
    mapping = aes(z = val),
    stroke = .2
  ) +
  geom_point(
    data = donnees,
    mapping = aes(x = x, y = y),
    color = "black"
  ) +
  guides(
    fill = guide_colorbar(title.vjust = 0.8)
  ) +
  labs(
    title = "Ecart-types"
  ) +
  theme_void() +
  theme(
    legend.position = 'bottom'
  )
