rm(list = ls())

setwd("C:/Users/GabrielNascimento/Downloads/HGP-INLA-novo")
# getwd()

source("blockNNGPrgeneric.R")
source("Irregblock.R")
source("utils.R")

library(terra)
library(sf)
library(rnaturalearth)
library(dplyr)
library(ggplot2)
library(patchwork) # Para combinar múltiplos gráficos
library(tidyr)
library(INLA)
library(stats)
library(knitr)
require(Matrix)

alpha = 0.8
n.partition <- 4
n.blocks <- 2^n.partition
num.nb <- 8

# ======== Aereal data ========
# Get world PM2.5 concentration

crs_proj = "EPSG:3310"
# world_data <- c(
#   rast("./CA_data/aereal_PM25_CA_2010/aereal_PM25_CA_2010.tif"),
#   rast("./CA_data/aereal_PM25_CA_2011/aereal_PM25_CA_2011.tif"),
#   rast("./CA_data/aereal_PM25_CA_2012/aereal_PM25_CA_2012.tif")
# )
# plot(world_data)

# Get california shapefile
# shapefile <- ne_states(
#   country = "United States of America",
#   returnclass = "sf"
# )

# shapefile <- shapefile[shapefile$name == "California", ]

# ======== Point data ========
fileList <- list.files(
  path = "./CA_data/point_PM25_CA_2012", #
  pattern = "\\.csv$",
  full.names = TRUE
)

# Get California point Data
pointData <- lapply(fileList, function(file) {
  f <- st_read(
    file,
    options = c(
      "X_POSSIBLE_NAMES=Site.Longitude",
      "Y_POSSIBLE_NAMES=Site.Latitude"
    ),
    quiet = TRUE
  )
  return(f)
})

world_raster = terra::rast(
  "./CA_data/global-annual-avg-pm2-5/PM25_201001_201212.tif"
) |>
  project(crs_proj)

CA_shp = st_read(
  "./CA_data/ca_state/",
  quiet = TRUE
)
# shapefile <- st_transform(CA_shp, crs(world_raster)) # sync CRS
CA_shp = st_transform(CA_shp, crs = crs_proj)
plot(CA_shp)
# === Transformar em grid
quad_tam = sqrt(101.95 * 10)
# CA_raster = mask(world_raster, CA_shp)
aereal_raster_proj <- mask(crop(world_raster, CA_shp), CA_shp)
aereal_data <- as.polygons(aereal_raster_proj) %>%
  st_as_sf() %>%
  rename(mean = PM25_201001_201212)

aereal_data_new = st_transform(aereal_data, "EPSG:3310")

aereal_data_grid <- st_make_grid(aereal_data_new, cellsize = 10195) %>%
  st_sf() %>%
  st_intersection(aereal_data_new) %>%
  st_cast("MULTIPOLYGON") %>%
  mutate(id = row_number())


plot(aereal_raster_proj)
# ====

plot(aereal_data_grid)

aerealData_sf = aereal_data_grid
#  ! use pipe operator

#  Get mean values
#  Cada estação possuí várias medidas ao longo do ano, o que faremos é tirar a média de todas as medidas dentro de um ano, e então tirar a média de todos os anos, assim, cada estação terá apenas uma medida, que é a média total

#  "Daily.Mean.PM2.5.Concentration", "Site.Latitude", "Site.Longitude"
annual_means <- lapply(pointData, function(df) {
  df$Daily.Mean.PM2.5.Concentration <- as.numeric(as.character(
    df$Daily.Mean.PM2.5.Concentration
  ))

  # Calcular médias por estação
  result <- df %>%
    group_by(Site.Latitude, Site.Longitude, Site.ID) %>%
    summarise(
      Annual_Mean = mean(Daily.Mean.PM2.5.Concentration, na.rm = TRUE),
      Count = sum(!is.na(Daily.Mean.PM2.5.Concentration)),
      .groups = "drop"
    )

  # Converter para objeto sf com geometria
  if (inherits(df, "sf")) {
    result <- st_as_sf(
      result,
      coords = c("Site.Longitude", "Site.Latitude"),
      crs = st_crs(df)
    )
  } else {
    # Se não for objeto sf, criar geometria manualmente
    result <- st_as_sf(
      result,
      coords = c("Site.Longitude", "Site.Latitude"),
      crs = 4326
    )
  }

  return(result)
})


all_stations <- do.call(rbind, annual_means) |> st_transform(crs_proj)

#  Get total mean
pointData_df <- all_stations %>%
  group_by(Site.ID) %>%
  summarise(
    mean = mean(Annual_Mean, na.rm = TRUE),
    .groups = "drop"
  )


mixedData <- bind_rows(
  aerealData = select(aerealData_sf, mean, geometry),
  pointData = select(pointData_df, mean, geometry),
  .id = "type"
)

# ggplot() +
#   geom_sf(data = mixedData, aes(color = mean, fill = mean)) +
#   labs(title = "") +
#   theme_minimal()

# plot(annual_means[[1]]$Count)

# cat(
#   "Qtd. de bservações: \n Ano 1:",
#   nrow(annual_means[[1]]),
#   "\nAno 2:",
#   nrow(annual_means[[2]]),
#   "\nAno 3:",
#   nrow(annual_means[[3]]),
#   "\n"
# )
# cat(
#   "Há",
#   length(unique(pointData_df$geometry)),
#   "estações de coleta e",
#   nrow(aerealData_sf),
#   "níveis de dados de área ->",
#   nrow(mixedData)
# )

# ======== Inferência ========
sf = st_transform(st_sf(mixedData$geometry), 4326)
sf <- st_make_valid(sf)
n = nrow(sf)

# n <- nrow(sf)
HGPdata <- get_HGPdata(
  sf,
  n.blocks,
  num.nb,
  alpha = 0.8,
  priors = list(a = 0, b = 796.88)
)


y <- mixedData$mean[HGPdata$order]

# View(mixedData)

# 1. Organizar os dados
data1 <- data.frame(y = y)
data1$idx <- 1:nrow(data1)

# 2. Modelo de efeito aleatório espacial
f.spatial <- y ~ 1 + f(idx, model = HGPdata$model)

# 3. Rodar INLA

exec.time = system.time({
  resf <- inla(
    f.spatial,
    data = as.data.frame(data1),
    family = "gaussian",
    control.compute = list(dic = TRUE, waic = TRUE)
  )
})


summary(resf)
# Extraction of hyperparameters
tau.est <- inla.emarginal(function(x) 1 / x, resf$marginals.hyperpar[[1]])
sigmasq.est <- inla.emarginal(function(x) exp(-x), resf$marginals.hyperpar[[2]])
phi.est <- inla.emarginal(
  function(x) 1 - 1 / (1 + exp(x)),
  resf$marginals.hyperpar[[3]]
)
w.est <- resf$summary.random$idx$mean
beta.est <- resf$summary.fixed[, "mean"]

summary.theta <- c(beta.est, tau.est, sigmasq.est, phi.est)

print(c("B0", "tau.est", "sigmasq.est", "phi.est"))
print(c(beta.est, tau.est, sigmasq.est, phi.est))

exec_time = exec.time
show_results(resf) # COLOCAR LOOIC
summarize_inla_hgp_CA(resf, alpha)


unorder_idx <- order(HGPdata$order)
results_df <- mixedData
results_df$fitted_mean <- resf$summary.fitted.values$mean[unorder_idx]
results_df$fitted_sd <- resf$summary.fitted.values$sd[unorder_idx]
results_df$spatial_effect <- resf$summary.random$idx$mean[unorder_idx]

# ======== ANÁLISE DOS BLOCOS ========

# --- 1. Adicionar o ID do Bloco ao dataframe de resultados ---
# Garante que cada observação (ponto ou polígono) saiba a qual bloco pertence.
unorder_idx <- order(HGPdata$order)
results_df$block_id <- HGPdata$blocks[unorder_idx]

# --- 2. Obter o Tamanho e Outras Informações dos Blocos ---
# Para calcular a área geográfica, primeiro transformamos para um sistema de coordenadas
# projetado apropriado para a Califórnia (EPSG:3310 - California Albers).
california_albers_crs <- 3310

block_summary_table <- results_df %>%
  # Transforma para o CRS apropriado para cálculo de área
  st_transform(crs = california_albers_crs) %>%
  # Agrupa por ID do bloco
  group_by(block_id) %>%
  # Calcula as estatísticas resumidas
  summarise(
    n_observacoes = n(),
    area_km2 = as.numeric(st_area(st_union(geometry))) / 1e6, # Une geometrias e calcula área em km²
    efeito_espacial_medio = mean(spatial_effect, na.rm = TRUE)
  ) %>%
  # Remove a coluna de geometria para uma tabela mais limpa
  st_drop_geometry() %>%
  arrange(block_id)

cat("--- Tabela Resumo por Bloco (Tamanho e Efeito Espacial) ---\n")
print(kable(
  block_summary_table,
  digits = 2,
  caption = "Tamanho (nº de observações e área) e efeito espacial médio por bloco."
))


# --- 3. Mapa 1: Plotar cada Bloco com uma Cor Diferente ---
# Este mapa serve para visualizar a partição espacial que o modelo criou.

mapa_particao_blocos <- ggplot(results_df) +
  geom_sf(aes(fill = as.factor(block_id)), color = "white", linewidth = 0.1) +
  # Usamos uma paleta com cores bem distintas para dados categóricos
  scale_fill_viridis_d(option = "turbo") +
  labs(
    title = "Visualização da Partição em 16 Blocos Espaciais",
    subtitle = "Cada cor representa um bloco diferente",
    fill = "ID do Bloco"
  ) +
  theme_void()

print(mapa_particao_blocos)


# --- 4. Mapa 2: Plotar o Efeito Espacial Médio por Bloco ---
# Este é o mapa choropleth que mostra as tendências espaciais.

# Primeiro, preparamos os dados (unir geometrias e calcular a média)
block_effect_sf <- results_df %>%
  group_by(block_id) %>%
  mutate(
    # O summarise em um objeto sf une as geometrias
    efeito_espacial_medio = mean(spatial_effect, na.rm = TRUE)
  )

mapa_efeito_blocos <- ggplot(data = block_effect_sf) +
  # Camada 1: Plota os polígonos dos blocos (usa o dado principal com 16 linhas)
  geom_sf(aes(fill = efeito_espacial_medio), color = "black", linewidth = 0.3) +

  # Camada 2: Adiciona os rótulos
  geom_sf_text(
    # CORREÇÃO: Especificamos o mesmo dado aqui para evitar ambiguidade
    data = block_effect_sf,
    aes(label = block_id),
    color = "white",
    fontface = "bold",
    size = 3.5
  ) +

  # O resto do código permanece o mesmo
  scale_fill_gradient2(
    low = "#00539C", # Azul para valores negativos
    mid = "white", # Branco para valores próximos a zero
    high = "#EEA47F", # Laranja/Vermelho para valores positivos
    midpoint = 0
  ) +
  labs(
    title = "Efeito Espacial Médio por Bloco",
    subtitle = "Média do componente espacial latente (w)",
    fill = "Efeito Médio"
  ) +
  theme_void()

# Exibir o gráfico corrigido
print(mapa_efeito_blocos)
# ======== VIZUALIZAÇÃO ========
library(ggplot2)


# 3. Adicionar as colunas do resultado do INLA, aplicando a reordenação

# 4. Calcular os resíduos (Observado - Ajustado)
results_df$residuals <- results_df$mean - results_df$fitted_mean

plot_obs_vs_fit <- ggplot(results_df, aes(x = mean, y = fitted_mean)) +
  geom_point(alpha = 0.5, color = "dodgerblue") +
  geom_abline(
    intercept = 0,
    slope = 1,
    color = "red",
    linetype = "dashed",
    size = 1
  ) +
  labs(
    title = "Valores Observados vs. Ajustados pelo Modelo",
    subtitle = "Ajuste do modelo HGP para dados de PM2.5",
    x = "PM2.5 Observado (y)",
    y = "PM2.5 Ajustado (Fitted)"
  ) +
  theme_minimal() +
  coord_fixed()

print(plot_obs_vs_fit)

library(ggplot2)

plot_residuals_hist <- ggplot(results_df, aes(x = residuals)) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 40,
    fill = rgb(207, 224, 252, maxColorValue = 255),
    color = rgb(66, 135, 245, maxColorValue = 255)
  ) +
  geom_density(
    aes(color = "Densidade dos Resíduos", linetype = "Densidade dos Resíduos"),
    linewidth = 1,
    key_glyph = "path"
  ) +
  stat_function(
    aes(color = "Curva Normal Teórica", linetype = "Curva Normal Teórica"),
    fun = dnorm,
    args = list(
      mean = mean(results_df$residuals),
      sd = sd(results_df$residuals)
    ),
    linewidth = 1
  ) +
  scale_color_manual(
    name = NULL,
    values = c(
      "Densidade dos Resíduos" = "blue",
      "Curva Normal Teórica" = rgb(230, 143, 117, maxColorValue = 255)
    )
  ) +
  scale_linetype_manual(
    name = NULL, # Remove o título da legenda
    values = c(
      "Densidade dos Resíduos" = "solid",
      "Curva Normal Teórica" = "dashed"
    )
  ) +
  labs(
    title = "Distribuição dos Resíduos do Modelo",
    x = "Resíduo (Observado - Ajustado)",
    y = "Densidade"
  ) +
  theme_minimal() +
  theme(legend.position = "top")

print(plot_residuals_hist)


# Resíduos vs. Valores Ajustados
plot_residuals_fit <- ggplot(results_df, aes(x = fitted_mean, y = residuals)) +
  geom_point(alpha = 0.5, color = "dodgerblue") +
  geom_hline(yintercept = 0, color = "gray", linetype = "dashed", size = 1) +
  labs(
    title = "Resíduos vs. Valores Ajustados",
    x = "PM2.5 Ajustado",
    y = "Resíduo"
  ) +
  theme_minimal()

polygons_sf <- results_df %>% filter(type == "aerealData")
points_sf <- results_df %>% filter(type == "pointData")

range_observed <- range(results_df$mean, na.rm = TRUE)
range_fitted <- range(results_df$fitted_mean, na.rm = TRUE)
range_effect <- range(results_df$spatial_effect, na.rm = TRUE)
range_residuals <- range(results_df$residuals, na.rm = TRUE)

# --- Mapa dos valores observados ---
map_observed <- ggplot() +
  geom_sf(data = polygons_sf, aes(fill = mean), color = NA) +
  geom_sf(data = points_sf, aes(color = mean), size = 2.5) +
  scale_fill_viridis_c(option = "plasma", limits = range_observed) +
  scale_color_viridis_c(option = "plasma", limits = range_observed) +
  labs(title = "PM2.5 Observado", fill = "PM2.5", color = "PM2.5") +
  theme_void()


# --- Mapa dos valores ajustados pelo modelo ---
map_fitted <- ggplot() +
  geom_sf(data = polygons_sf, aes(fill = fitted_mean), color = NA) +
  geom_sf(data = points_sf, aes(color = fitted_mean), size = 2.5) +
  scale_fill_viridis_c(option = "plasma", limits = range_fitted) +
  scale_color_viridis_c(option = "plasma", limits = range_fitted) +
  labs(title = "PM2.5 Ajustado", fill = "PM2.5", color = "PM2.5") +
  theme_void()

# --- Mapa do Efeito Aleatório Espacial ---
map_spatial_effect <- ggplot() +
  geom_sf(data = polygons_sf, aes(fill = spatial_effect), color = NA) +
  geom_sf(data = points_sf, aes(color = spatial_effect), size = 2.5) +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    limits = range_effect
  ) +
  scale_color_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    limits = range_effect
  ) +
  labs(
    title = "Efeito Espacial Latente (w)",
    fill = "Efeito",
    color = "Efeito"
  ) +
  theme_void()

# --- Mapa dos resíduos ---
map_residuals <- ggplot() +
  geom_sf(data = polygons_sf, aes(fill = residuals), color = NA) +
  geom_sf(data = points_sf, aes(color = residuals), size = 2.5) +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    limits = range_residuals
  ) +
  scale_color_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    limits = range_residuals
  ) +
  labs(
    title = "Resíduos (Observado - Ajustado)",
    fill = "Resíduo",
    color = "Resíduo"
  ) +
  theme_void()


print(plot_residuals_hist)
print(plot_residuals_fit)

print(map_observed | map_fitted)

print(map_spatial_effect)
print(map_residuals)


plot_obs_vs_fit <- ggplot(results_df, aes(x = fitted_mean, y = mean)) +
  geom_point(alpha = 0.5, color = "dodgerblue") +
  geom_hline(yintercept = 0, color = "gray", linetype = "dashed", size = 1) +
  labs(
    title = "Observado vs. Ajustado",
    x = "Ajustado",
    y = "Observado"
  ) +
  theme_minimal()

# 4. Calcular os resíduos (Observado - Ajustado)
results_aereal = results_df |> filter(type == "aerealData")
results_pt = results_df |> filter(type == "aerealData")

# # --- Mapa dos valores ajustados pelo modelo ---
# map_fitted_pt <- ggplot() +
#   geom_sf(data = points_sf, aes(color = fitted_mean), size = 2.5) +
#   scale_fill_viridis_c(option = "plasma", limits = range_fitted) +
#   scale_color_viridis_c(option = "plasma", limits = range_fitted) +
#   labs(
#     title = "Valores ajustados",
#     fill = "PM2.5",
#     color = "PM2.5"
#   ) +
#   theme_void()
# map_fitted_poly <- ggplot() +
#   geom_sf(data = polygons_sf, aes(fill = fitted_mean), color = NA) +
#   scale_fill_viridis_c(option = "plasma", limits = range_fitted) +
#   scale_color_viridis_c(option = "plasma", limits = range_fitted) +
#   labs(
#     title = "Valores ajustados",
#     fill = "PM2.5",
#     color = "PM2.5"
#   ) +
#   theme_void()

# map_residuals_pt <- ggplot() +
#   geom_sf(data = points_sf, aes(color = residuals), size = 2.5) +
#   scale_fill_gradient2(
#     low = "blue",
#     mid = "white",
#     high = "red",
#     midpoint = 0,
#     limits = range_residuals
#   ) +
#   scale_color_gradient2(
#     low = "blue",
#     mid = "white",
#     high = "red",
#     midpoint = 0,
#     limits = range_residuals
#   ) +
#   labs(
#     title = "Resíduos (Observado - Ajustado)",
#     fill = "Resíduo",
#     color = "Resíduo"
#   ) +
#   theme_void()
# map_residuals_poly <- ggplot() +
#   geom_sf(data = polygons_sf, aes(fill = residuals), color = NA) +
#   scale_fill_gradient2(
#     low = "blue",
#     mid = "white",
#     high = "red",
#     midpoint = 0,
#     limits = range_residuals
#   ) +
#   scale_color_gradient2(
#     low = "blue",
#     mid = "white",
#     high = "red",
#     midpoint = 0,
#     limits = range_residuals
#   ) +
#   labs(
#     title = "Resíduos (Observado - Ajustado)",
#     fill = "Resíduo",
#     color = "Resíduo"
#   ) +
#   theme_void()

# print(map_fitted_poly | map_residuals_poly)
# print(map_fitted_pt | map_residuals_pt)
