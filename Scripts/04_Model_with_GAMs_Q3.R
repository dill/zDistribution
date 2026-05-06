#### 3D Distribution
#### Models testing Q3: Depth distribution across xy
#### Fall 2025
#### EKJ&AVC

library(mgcv)
library(tidyverse)
library(PNWColors)
library(marmap)
library(terra)

load("./ProcessedData/detect_data.RData")
load("./ProcessedData/detect_data_clean.RData")
detect_data <- detect_data %>% mutate(BestTaxon = as.factor(BestTaxon))
mmEcoEvo <- read.csv("./Data/MM_metadata.csv")

### Get bathymetry for pred grid -----------------------------------------------

bathy <- getNOAA.bathy(lon1 = min(detect_data_clean$lon), 
                       lon2 = max(detect_data_clean$lon), 
                       lat1 = min(detect_data_clean$lat),  
                       lat2 = max(detect_data_clean$lat),
                       resolution = 1)

bathy_raster <- marmap::as.raster(bathy)
bathy_r <- rast(bathy_raster)

### Get density estimates from Becker et al. 2020 tech memo --------------------

density_files <- list.files(path = "./Data/CCE_model_run", pattern = "predgrid_yearly", full.names = TRUE)

density_data <- density_files %>%
  set_names(basename(.)) %>%
  map_dfr(read_csv, .id = "file") %>% 
  filter(year == 2014) %>% 
  separate(file, into = c(NA,NA,"species",NA), sep = "_")

#subset detect data to target the three main species
detect_data_sub <- detect_data_clean %>% 
  filter(BestTaxon %in% c("Berardius bairdii", 
                          "Lagenorhynchus obliquidens",
                          "Megaptera novaeangliae"))

#join data using a spatial join
library(sf)

detections_sf <- detect_data_sub %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326)

density_sf <- density_data %>%
  st_as_sf(coords = c("mlon", "mlat"), crs = 4326)

detect_data_dens <- detections_sf %>%
  group_by(BestTaxon) %>%
  group_modify(~ {density_sub <- density_sf %>% filter(species == unique(.y$BestTaxon))
    st_join(.x, density_sub, join = st_nearest_feature)}) %>%
  ungroup() %>% 
  select(-geometry)


# Q3: Does depth distribution of detections vary across xy spatial distribution or species?
# H3.0: Depth distribution of detections does not vary across xy spatial distribution or species.
  #H3.0a: xy variability + z variability
  #H3.0b: xyz covariability
  #H3.0c: xyz covariability by species
  #H3.1: z and density covariability by species

### H3.0a: Depth distribution plus xy space distribution -----------------------

m3.0a <- bam(Detected ~ s(depth) + s(utm.lat, utm.lon),
             family = "binomial",
             data = detect_data_dens,
             method = "fREML",
             discrete = FALSE,
             nthreads = 4)
summary(m3.0a)

# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq  p-value    
# s(depth)            2.763  3.269  23.80 4.74e-05 ***
#   s(utm.lat,utm.lon) 14.600 18.878  80.77  < 2e-16 ***
#   ---
#   
# R-sq.(adj) =  0.0229   Deviance explained = 7.88%

gam.check(m3.0a)
# k'   edf k-index p-value
# s(depth)            9.00  2.76     0.9    0.25
# s(utm.lat,utm.lon) 29.00 14.60     0.9    0.38

AIC(m3.0a)
#1629

### H3.0b: Depth distribution smoothed by xy space distribution ----------------

m3.0b <- bam(Detected ~ te(depth, utm.lat, utm.lon),
             family = "binomial",
             data = detect_data_dens,
             method = "fREML",
             discrete = TRUE,
             nthreads = 4)
summary(m3.0b)
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value    
# te(depth,utm.lat,utm.lon) 28.72  36.26  142.1  <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.0261   Deviance explained = 11.2%
AIC(m3.0b)
#1596

### m3.0b predictions ----------------------------------------------------------

m3.0b_pred_grid <- expand_grid(depth = seq(0,500, by = 100),
                               utm.lat = seq(min(detect_data$utm.lat, na.rm = TRUE),
                                         max(detect_data$utm.lat, na.rm = TRUE),
                                         by = 100),
                               utm.lon = seq(min(detect_data$utm.lon, na.rm = TRUE),
                                         max(detect_data$utm.lon, na.rm = TRUE),
                                         by = 100))


m3.0bpreds <- predict.gam(m3.0b, m3.0b_pred_grid, se.fit = TRUE)

m3.0b_sePreds <- data.frame(m3.0b_pred_grid,
                            mu   = binomial()$linkinv(m3.0bpreds$fit),
                            low  = binomial()$linkinv(m3.0bpreds$fit - 1.96 * m3.0bpreds$se.fit),
                            high = binomial()$linkinv(m3.0bpreds$fit + 1.96 * m3.0bpreds$se.fit))

### H3.0c: Depth smoothed over xy with shape and intercept variable by species -

m3.0c <-
  bam(Detected ~ 
        # main effects of space, depth, taxon
        # ti(utm.lon, utm.lat,
        #    d=2,
        #    k=20,
        #    bs="tp")+
        #  ti(depth,
        #     k=5,
        #     bs="ts")+
        ti(BestTaxon,
           k=3,
           bs="re")+
        # interaction between *everything*
        ti(utm.lon, utm.lat, depth, BestTaxon,
           d=c(2,1,1),
           k=c(20, 5, 3),
           bs=c("tp","ts", "re"))+
        # space-taxon effect
        ti(utm.lon, utm.lat, BestTaxon,
           d=c(2,1),
           k=c(20,3),
           bs=c("tp","re"))+
        # depth-taxon effect
        ti(depth, BestTaxon,
           k=c(5,3),
           bs=c("ts","re")),
      family = "binomial",
      method = "fREML",
      data = detect_data_dens,
      discrete = TRUE)

summary(m3.0c)
# Approximate significance of smooth terms:
#   edf  Ref.df  Chi.sq  p-value    
# ti(utm.lon,utm.lat)                 6.931e+00   8.479  12.553 0.152585    
# ti(depth)                           1.809e-05   4.000   0.000 0.090641 .  
# ti(BestTaxon)                       1.624e+00   2.000   7.415 0.002093 ** 
#   ti(BestTaxon,depth,utm.lon,utm.lat) 3.734e+01 228.000 107.762  < 2e-16 ***
#   ti(BestTaxon,utm.lon,utm.lat)       8.151e+00  55.000  15.608 0.000707 ***
#   ti(BestTaxon,depth)                 7.332e+00  12.000  35.902  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.0927   Deviance explained =   26%

### WITHOUT DEPTH OR LAT/LON

# Approximate significance of smooth terms:
#   edf Ref.df  Chi.sq  p-value    
# ti(BestTaxon)                        1.677      2   9.396 0.000499 ***
#   ti(BestTaxon,depth,utm.lon,utm.lat) 40.895    228 108.583  < 2e-16 ***
#   ti(BestTaxon,utm.lon,utm.lat)       17.077     57  41.507  < 2e-16 ***
#   ti(BestTaxon,depth)                  6.948     12  36.440  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.0934   Deviance explained = 26.1%
# fREML = 7856.2  Scale est. = 1         n = 7740

AIC(m3.0c)
# 1413 with all terms
# 1428 with non-significant terms (depth and lat/lon) removed

gam.check(m3.0c)
# k'    edf k-index p-value
# ti(BestTaxon)                         3.00   1.68      NA      NA
# ti(BestTaxon,depth,utm.lon,utm.lat) 228.00  40.89      NA      NA
# ti(BestTaxon,utm.lon,utm.lat)        57.00  17.08      NA      NA
# ti(BestTaxon,depth)                  12.00   6.95      NA      NA

#mean squared Pearson residual dispersion parameter
sum(residuals(m3.0c, type = "pearson")^2) / df.residual(m3.0c)

detect_data_dens$resid_response <- residuals(m3.0c,
                                             type = "response")

detect_data_dens$resid_pearson <- residuals(m3.0c, type = "pearson")

ggplot(detect_data_dens,
       aes(utm.lon, utm.lat,
           color = resid_pearson), size = 2) +
  geom_point(alpha = 0.7) +
  facet_wrap(~BestTaxon) +
  scale_color_gradient2(midpoint = 0) +
  coord_equal()

### m3.0c predictions ----------------------------------------------------------

m3.0c_pred_grid <- expand_grid(depth = seq(from = 0, to = 500, by = 10),
                               utm.lat = seq(min(detect_data_dens$utm.lat, na.rm = TRUE),
                                         max(detect_data_dens$utm.lat, na.rm = TRUE),
                                         by = 100),
                               utm.lon = seq(min(detect_data_dens$utm.lon, na.rm = TRUE),
                                         max(detect_data_dens$utm.lon, na.rm = TRUE),
                                         by = 100),
                               BestTaxon = as.factor(c("Lagenorhynchus obliquidens",
                                                        "Megaptera novaeangliae",
                                                        "Berardius bairdii")))

# mask depths that aren't real
m3.0c_pred_grid$seafloor <- terra::extract(bathy_r,
                                                 m3.0c_pred_grid[, c("utm.lon", "utm.lat")])[,2]

m3.0c_pred_grid_trimmed <- m3.0c_pred_grid[is.na(m3.0c_pred_grid$seafloor) | 
                                                         m3.0c_pred_grid$depth <= pmax(abs(m3.0c_pred_grid$seafloor), 0),]


# response predictions
m3.0cpreds <- predict.bam(m3.0c, m3.0c_pred_grid_trimmed,
                          se.fit = TRUE)

m3.0c_sePreds <- data.frame(m3.0c_pred_grid,
                            mu   = binomial()$linkinv(m3.0cpreds$fit),
                            low  = binomial()$linkinv(m3.0cpreds$fit - 1.96 * m3.0cpreds$se.fit),
                            high = binomial()$linkinv(m3.0cpreds$fit + 1.96 * m3.0cpreds$se.fit),
                            low50  = binomial()$linkinv(m3.0cpreds$fit - 0.674 * m3.0cpreds$se.fit),
                            high50 = binomial()$linkinv(m3.0cpreds$fit + 0.674 * m3.0cpreds$se.fit))


#Marginal effect of depth
ref_lon <- mean(detect_data_dens$utm.lon, na.rm = TRUE)
ref_lat <- mean(detect_data_dens$utm.lat, na.rm = TRUE)

depth_grid <- expand.grid(
  depth = seq(from = 0, to = 500, by = 10),
  utm.lon = ref_lon,
  utm.lat = ref_lat,
  BestTaxon = as.factor(c("Lagenorhynchus obliquidens",
                          "Megaptera novaeangliae",
                          "Berardius bairdii")))

pred <- predict.bam(m3.0c, newdata = depth_grid, se.fit = TRUE)

depth_grid$fit <- binomial()$linkinv(pred$fit)

depth_grid$low <- binomial()$linkinv(pred$fit - 1.96 * pred$se.fit)

depth_grid$high <- binomial()$linkinv(pred$fit + 1.96 * pred$se.fit)

ggplot(depth_grid, aes(depth, fit, color = BestTaxon, fill = BestTaxon)) +
  geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1.2) +
  labs(y = "Predicted detection probability",
    x = "Depth") +
  facet_wrap(~BestTaxon, scales = "free")

### Save -----------------------------------------------------------------------

save(m3.0a, m3.0b, m3.0c, m3.0b_sePreds, m3.0cpreds, m3.0c_sePreds,
     file = "./ProcessedData/m3.0models_preds_0.05degree.Rdata")

save(m3.0a, m3.0b, m3.0c, m3.0b_sePreds, m3.0cpreds, m3.0c_sePreds,
     file = "./ProcessedData/m3.0models_preds_0.05degree_depthmask.Rdata")

### H3.1: Detection smoothed over depth + density by species ------------------

m3.1 <-
  bam(Detected ~ 
        # main effects of density, depth, taxon
        # ti(depth,
        #    k=5,
        #    bs="ts")+
        ti(BestTaxon,
           k=3,
           bs="re")+
        ti(D,
           k=5,
           bs="ts")+
        # interaction between *everything*
        ti(D, depth, BestTaxon,
           d=c(1,1,1),
           k=c(5, 5, 3),
           bs=c("ts","ts", "re")) +
        # depth-taxon effect
        ti(depth, BestTaxon,
           k=c(5,3),
           bs=c("ts","re")) +
        # # density-taxon effect
        ti(D, BestTaxon,
           k=c(5,3),
           bs=c("ts","re")),
      family = "binomial",
      method = "fREML",
      data = detect_data_dens,
      discrete = TRUE)

summary(m3.1)
# Approximate significance of smooth terms:
#   edf Ref.df  Chi.sq  p-value    
# ti(depth)             3.706e-01      4   0.000 0.999861    
# ti(BestTaxon)         1.669e+00      2 273.060 0.000111 ***
#   ti(D)                 3.177e-05      4   2.093  < 2e-16 ***
#   ti(D,BestTaxon,depth) 1.152e+01     46  69.786  < 2e-16 ***
#   ti(BestTaxon,depth)   5.830e+00      8  29.189  < 2e-16 ***
#   ti(BestTaxon,D)       3.061e+00      8  30.901  < 2e-16 ***

# R-sq.(adj) =  -0.415   Deviance explained =  -29%
# fREML = 8267.8  Scale est. = 1         n = 7740

### WITHOUT DEPTH 

# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq  p-value    
# ti(BestTaxon)         1.4163      2 41.659 0.024028 *  
#   ti(D)                 0.1461      4  0.314 0.081469 .  
# ti(D,BestTaxon,depth) 4.5703     46 14.399 0.000269 ***
#   ti(BestTaxon,depth)   7.5964     12 40.239  < 2e-16 ***
#   ti(BestTaxon,D)       3.3138      8 52.233  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.0523   Deviance explained = 15.2%
# fREML = 7874.4  Scale est. = 1         n = 7740

AIC(m3.1)
# 1507 with all terms
# 1507 with non-significant term (depth) removed

#mean squared Pearson residual dispersion parameter
sum(residuals(m3.1, type = "pearson")^2) / df.residual(m3.1)

detect_data_dens$resid_pearson <- residuals(m3.1, type = "pearson")

ggplot(detect_data_dens,
       aes(utm.lon, utm.lat,
           color = resid_pearson), size = 2) +
  geom_point(alpha = 0.7) +
  facet_wrap(~BestTaxon) +
  scale_color_gradient2(midpoint = 0) +
  coord_equal()

### m3.1 predictions ----------------------------------------------------------

m3.1_pred_grid <- expand_grid(depth = seq(from = 0, to = 500, by = 10),
                               D = seq(from = min(detect_data_dens$D),
                                       to = max(detect_data_dens$D), 
                                       by = 0.05),
                               BestTaxon = as.factor(c("Lagenorhynchus obliquidens",
                                                       "Megaptera novaeangliae",
                                                       "Berardius bairdii")))
# response predictions
m3.1preds <- predict.bam(m3.1, m3.1_pred_grid,
                          se.fit = TRUE)

m3.1_sePreds <- data.frame(m3.1_pred_grid,
                            mu   = binomial()$linkinv(m3.1preds$fit),
                            low  = binomial()$linkinv(m3.1preds$fit - 1.96 * m3.1preds$se.fit),
                            high = binomial()$linkinv(m3.1preds$fit + 1.96 * m3.1preds$se.fit),
                            low50  = binomial()$linkinv(m3.1preds$fit - 0.674 * m3.1preds$se.fit),
                            high50 = binomial()$linkinv(m3.1preds$fit + 0.674 * m3.1preds$se.fit))

###Some plots for m3.1 ---------------------------------------------------------
#Detection rate by depth and density estimate (3D interactive)
library(plotly)

plot_ly(
  m3.1_sePreds,
  x = ~depth,
  y = ~D,
  z = ~mu,
  color = ~BestTaxon,
  type = "scatter3d",
  mode = "markers"
)

#Detection rate by depth and density estimate (2D)
g <- m3.1_sePreds %>%
  arrange(BestTaxon, depth, D)

g_list <- split(g, g$BestTaxon)

plots <- lapply(names(g_list), function(sp) {
  
  df <- g_list[[sp]]
  
  ggplot(df, aes(depth, D, fill = mu)) +
    geom_raster() +
    scale_fill_viridis_c(limits = range(df$mu)) +
    labs(title = sp) +
    theme(legend.position = "bottom")
})

patchwork::wrap_plots(plots)

#Marginal effect of density
D_grid <- expand.grid(
  depth = median(detect_data_dens$depth, na.rm = TRUE),
  D = seq(min(detect_data_dens$D, na.rm = TRUE),
          max(detect_data_dens$D, na.rm = TRUE),
          length.out = 100),
  BestTaxon = as.factor(c("Lagenorhynchus obliquidens",
                                      "Megaptera novaeangliae",
                                      "Berardius bairdii")))

pred_D <- predict.bam(m3.1, newdata = D_grid, se.fit = TRUE)

D_grid$fit <- binomial()$linkinv(pred_D$fit)
D_grid$low <- binomial()$linkinv(pred_D$fit - 1.96 * pred_D$se.fit)
D_grid$high <- binomial()$linkinv(pred_D$fit + 1.96 * pred_D$se.fit)

ggplot(D_grid, aes(D, fit)) +
  geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.2) +
  geom_line(linewidth = 1) +
  facet_wrap(~BestTaxon) +
  labs(y = "Predicted detection probability")

#Marginal effect of depth
depth_grid <- expand.grid(
  depth = seq(0, 500, by = 5),
  D = median(detect_data_dens$D, na.rm = TRUE),
  BestTaxon = as.factor(c("Lagenorhynchus obliquidens",
                          "Megaptera novaeangliae",
                          "Berardius bairdii")))

pred_depth <- predict.bam(m3.1, newdata = depth_grid, se.fit = TRUE)

depth_grid$fit <- binomial()$linkinv(pred_depth$fit)
depth_grid$low <- binomial()$linkinv(pred_depth$fit - 1.96 * pred_depth$se.fit)
depth_grid$high <- binomial()$linkinv(pred_depth$fit + 1.96 * pred_depth$se.fit)

ggplot(depth_grid, aes(depth, fit)) +
  geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.2) +
  geom_line(linewidth = 1) +
  facet_wrap(~BestTaxon) +
  labs(y = "Predicted detection probability")
