#### 3D Distribution
#### Models testing Q3: Depth distribution across xy
#### Fall 2025
#### EKJ&AVC

library(mgcv)
library(tidyverse)
library(PNWColors)
library(marmap)
library(terra)
library(rstan)

load("../ProcessedData/detect_data.Rdata")
load("../ProcessedData/detect_data_clean.RData")
detect_data <- detect_data %>% mutate(BestTaxon = as.factor(BestTaxon))
mmEcoEvo <- read.csv("../Data/MM_metadata.csv")


load("../ProcessedData/m3.0models_preds_0.05degree.Rdata")

load("../ProcessedData/m3.0models_preds_0.05degree_depthmask.Rdata")


### Get bathymetry for pred grid -----------------------------------------------

bathy <- getNOAA.bathy(lon1 = min(detect_data_clean$lon), 
                       lon2 = max(detect_data_clean$lon), 
                       lat1 = min(detect_data_clean$lat),  
                       lat2 = max(detect_data_clean$lat),
                       resolution = 1)

bathy_raster <- marmap::as.raster(bathy)
bathy_r <- rast(bathy_raster)

w <- matrix(1, 5, 5)


### m3.0c_clean predictions ----------------------------------------------------

m3.0c_clean_pred_grid <- expand_grid(depth = seq(from = 0, to = 500, by = 10),
                               lat = seq(min(detect_data_clean$lat, na.rm = TRUE),
                                         max(detect_data_clean$lat, na.rm = TRUE),
                                         by = 0.05),
                               lon = seq(min(detect_data_clean$lon, na.rm = TRUE),
                                         max(detect_data_clean$lon, na.rm = TRUE),
                                         by = 0.05),
                               BestTaxon = as.factor(c("Lagenorhynchus obliquidens",
                                                       "Megaptera novaeangliae",
                                                       "Berardius bairdii")))

# mask depths that aren't real
m3.0c_clean_pred_grid$seafloor <- terra::extract(bathy_r,
                              m3.0c_clean_pred_grid[, c("lon", "lat")])[,2]

m3.0c_clean_pred_grid_trimmed <- m3.0c_clean_pred_grid[is.na(m3.0c_clean_pred_grid$seafloor) | 
                                                         m3.0c_clean_pred_grid$depth <= pmax(abs(m3.0c_clean_pred_grid$seafloor), 0),]

#m3.0c_clean_pred_grid[is.na(m3.0c_clean_pred_grid$seafloor),]

# response predictions
ind <- 1:100
m3.0c_clean_preds <- predict.bam(m3.0c_clean, m3.0c_clean_pred_grid_trimmed[ind,],
                          se.fit = TRUE)
mt <- predict(m3.0c_clean, m3.0c_clean_pred_grid_trimmed[ind,], discrete=FALSE, type="response")

m3.0c_clean_sePreds <- data.frame(m3.0c_clean_pred_grid_trimmed[ind,],
                            mu   = binomial()$linkinv(m3.0c_clean_preds$fit),
                            low  = binomial()$linkinv(m3.0c_clean_preds$fit - 1.96 * m3.0c_clean_preds$se.fit),
                            high = binomial()$linkinv(m3.0c_clean_preds$fit + 1.96 * m3.0c_clean_preds$se.fit),
                            low50  = binomial()$linkinv(m3.0c_clean_preds$fit - 0.674 * m3.0c_clean_preds$se.fit),
                            high50 = binomial()$linkinv(m3.0c_clean_preds$fit + 0.674 * m3.0c_clean_preds$se.fit))



# model is BIG so need to do this in a memory efficient way
# gam.mh will break as will gratia::fitted_samples
source("../Scripts/XX_binomial_sampler.R")
aa <- binomial_sampler(m3.0c_clean, m3.0c_clean_pred_grid_trimmed[ind,], iter=500, warmup=0)


source("../Scripts/XX_binomial_sampler_lap.R")
aa2 <- binomial_sampler(m3.0c_clean, m3.0c_clean_pred_grid_trimmed[ind,])


library(posterior)
yr <- aa2$summary(variables = "y_rep",
       posterior::default_summary_measures()[1:4],
       quantiles = ~ quantile2(., probs = c(0.025, 0.975)),
       posterior::default_convergence_measures()
       )

plot(yr$mean, mt)

plot(yr$mean, m3.0c_clean_sePreds$mu)
points(yr$q5, m3.0c_clean_sePreds$low)
points(yr$q95, m3.0c_clean_sePreds$l)


plot(yr$sd, m3.0c_clean_preds$se.fit)


