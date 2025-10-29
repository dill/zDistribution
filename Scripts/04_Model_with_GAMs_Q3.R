#### 3D Distribution
#### Models testing Q3: Depth distribution across xy
#### Fall 2025
#### EKJ&AVC

library(mgcv)
library(tidyverse)
library(PNWColors)

load("./ProcessedData/detect_data.RData")
detect_data <- detect_data %>% mutate(BestTaxon = as.factor(BestTaxon))
mmEcoEvo <- read.csv("./Data/MM_metadata.csv")

# Q3: Does depth distribution of detections vary across xy spatial distribution?
# H3.0: Depth distribution of detections does not vary across xy spatial distribution.
  #H3.0a: xy variability + z variability
  #H3.0b: xyz covariability
  #H3.0c: xyz covariability by species

### H3.0a: Depth distribution plus xy space distribution -----------------------

m3.0a <- bam(Detected ~ s(depth) + s(utm.lat, utm.lon),
             family = "binomial",
             data = detect_data,
             method = "fREML",
             discrete = TRUE,
             nthreads = 4)
summary(m3.0a)
gam.check(m3.0a)
#both depth and xy are significant
AIC(m3.0a)
#4511

### H3.0b: Depth distribution smoothed by xy space distribution ----------------

m3.0b <- bam(Detected ~ te(depth, lat, lon),
             family = "binomial",
             data = detect_data,
             method = "fREML",
             discrete = TRUE,
             nthreads = 4)
summary(m3.0b)
#te(depth, utm.lat, utm.lon) is significant
AIC(m3.0b)
#4477

### m3.0b predictions ----------------------------------------------------------

m3.0b_pred_grid <- expand_grid(depth = seq(0,500, by = 100),
                               lat = seq(min(detect_data$lat, na.rm = TRUE),
                                         max(detect_data$lat, na.rm = TRUE),
                                         by = 0.1),
                               lon = seq(min(detect_data$lon, na.rm = TRUE),
                                         max(detect_data$lon, na.rm = TRUE),
                                         by = 0.1))


m3.0bpreds <- predict.gam(m3.0b, m3.0b_pred_grid, se.fit = TRUE)

m3.0b_sePreds <- data.frame(m3.0b_pred_grid,
                            mu   = exp(m3.0bpreds$fit),
                            low  = exp(m3.0bpreds$fit - 1.96 * m3.0bpreds$se.fit),
                            high = exp(m3.0bpreds$fit + 1.96 * m3.0bpreds$se.fit))

### H3.0c: Depth smoothed over xy with shape and intercept variable by species -

m3.0c <-
  bam(Detected ~ 
        # main effects of space, depth, taxon
        ti(lon, lat,
           d=2,
           k=20,
           bs="tp")+
        ti(depth,
           k=5,
           bs="ts")+
        ti(BestTaxon,
           k=16,
           bs="re")+
        # interaction between *everything*
        ti(lon, lat, depth, BestTaxon,
           d=c(2,1,1),
           k=c(20, 5, 16),
           bs=c("tp","ts", "re"))+
        # space-taxon effect
        ti(lon, lat, BestTaxon,
           d=c(2,1),
           k=c(10,16),
           bs=c("tp","re"))+
        # depth-taxon effect
        ti(depth, BestTaxon,
           k=c(10,16),
           bs=c("ts","re")),
      family = "binomial",
      method = "fREML",
      data = detect_data,
      discrete = TRUE)

summary(m3.0c)
# Approximate significance of smooth terms:
#                               edf       Ref.df  Chi.sq   p-value    
#   ti(lon,lat)                 1.048e+01   13.27  16.75   0.242    
#   ti(depth)                   2.557e-05    4.00   0.00   0.736    
#   ti(BestTaxon)               1.255e+01   15.00 122.49  <2e-16 ***
#   ti(BestTaxon,depth,lon,lat) 4.904e+01 1216.00 102.55  <2e-16 ***
#   ti(lon,lat,BestTaxon)       4.169e+01  142.00  95.53  <2e-16 ***
#   ti(depth,BestTaxon)         2.591e+01  144.00  90.24  <2e-16 ***

# R-sq.(adj) =  0.082   Deviance explained = 23.7%
# fREML =  25017  Scale est. = 1         n = 25088

AIC(m3.0c)
#3810 with all terms
#4947 without non-significant terms

### m3.0c predictions ----------------------------------------------------------

m3.0c_pred_grid <- expand_grid(depth = seq(from = 0, to = 500, by = 10),
                               lat = seq(min(detect_data$lat, na.rm = TRUE),
                                         max(detect_data$lat, na.rm = TRUE),
                                         by = 0.05),
                               lon = seq(min(detect_data$lon, na.rm = TRUE),
                                         max(detect_data$lon, na.rm = TRUE),
                                         by = 0.05),
                               BestTaxon = as.factor(c("Lagenorhynchus obliquidens",
                                                        "Megaptera novaeangliae",
                                                        "Berardius bairdii")))
# response predictions
m3.0cpreds <- predict.bam(m3.0c, m3.0c_pred_grid, 
                          type = "response",
                          se.fit = TRUE)

m3.0c_sePreds_response <- data.frame(m3.0c_pred_grid,
                            mu   = exp(m3.0cpreds$fit),
                            low  = exp(m3.0cpreds$fit - 1.96 * m3.0cpreds$se.fit),
                            high = exp(m3.0cpreds$fit + 1.96 * m3.0cpreds$se.fit),
                            low50  = exp(m3.0cpreds$fit - 0.674 * m3.0cpreds$se.fit),
                            high50 = exp(m3.0cpreds$fit + 0.674 * m3.0cpreds$se.fit))

# link predictions
m3.0cpredsL <- predict.bam(m3.0c, m3.0c_pred_grid, 
                          se.fit = TRUE)

m3.0c_sePreds_link <- data.frame(m3.0c_pred_grid,
                                     mu   = exp(m3.0cpredsL$fit),
                                     low  = exp(m3.0cpredsL$fit - 1.96 * m3.0cpredsL$se.fit),
                                     high = exp(m3.0cpredsL$fit + 1.96 * m3.0cpredsL$se.fit),
                                     low50  = exp(m3.0cpredsL$fit - 0.674 * m3.0cpredsL$se.fit),
                                     high50 = exp(m3.0cpredsL$fit + 0.674 * m3.0cpredsL$se.fit))

### Save -----------------------------------------------------------------------

save(m3.0a, m3.0b, m3.0c, m3.0b_sePreds, m3.0c_sePreds_response, 
     m3.0c_sePreds_link, m3.0cpredsL, m3.0cpreds,
     file = "./ProcessedData/m3.0models_preds_0.05degree.Rdata")

### DLM explanation of ti vs. te -----------------------------------------------

# te can be decomposed into multiple ti terms
# "te produces a full tensor product smooth, while ti produces a tensor product 
# interaction, appropriate when the main effects (and any lower interactions) 
# are also present."
# Using UTM proj coordinates
# d groups lat and lon into one smooth, so 2D smooth on lat/lon
# this is cheaper to set up, and you get a 2D spline
# 1D smooth on depth, and 1D on BestTaxon
# changed depth to k = 5 because there are really only 5 common values
# thin plate for lat/lon and for depth, BestTaxon is a random effect
# first term is the full interaction
# also term for just interaction between lat/lon and BestTaxon
# also term for just interaction between depth and BestTaxon
# this allows you to see which bits are important
# individual ti terms could be estimated as zero, indicating
# that they are not important (hypothesis testing!)

# tensor between spline and random effect is the same as a factor-smooth
# model (Pedersen et al., 2019), where each level of the factor gives a
# related set of basis functions (a spline), but the corresponding
# coefficients are shrunk so we have more parameter efficiency (contrast
# to by= where the smooths are independent)

# I think for ti() terms you always need to include the "main effects"
# that is, if you have ti(x,y,z) you also must include ti(x)+ti(y)+ti(z)

### NOT RUN:

### H3.1a: Depth and oceanographic variables -----------------------------------

# m3.1a <- bam(Detected ~ te(depth, SST, bottom_depth,
#                           k = c(10,5,5)),
#                 family = "binomial",
#                 data = "REML",
#                 discrete = TRUE,
#                 nthreads = 40)
# 
# summary(m3.1a)
# #
# AIC(m3.1a)
# #
# 
# save(m3.1a, file = "./ProcessedData/m3.1a.RData")

### H3.1b: Depth and oceanographic variables by species ------------------------ 

# m3.1b <- bam(Detected ~ te(depth, SST, bottom_depth,     
#                            k=c(10,5,5),  
#                            by = as.factor(BestTaxon)), 
#              family = "binomial",   
#              data = detect_data,  
#              discrete = TRUE,
#              nthreads = 40)   
# 
# save(m3.1b, file = "m3.1b.RData")   
# 
# summary(m3.1b)
# #
# AIC(m3.1b)
# #
# 
# ### M3.1b predictions ----------------------------------------------------------
# 
# m3.1b_pred_grid <- expand_grid(depth = seq(0,500, by = 100),
#                                SST = seq(min(detect_data$SST, na.rm = TRUE),
#                                          max(detect_data$SST, na.rm = TRUE),
#                                          by = 0.1),
#                                lon = seq(min(detect_data$lon, na.rm = TRUE),
#                                          max(detect_data$lon, na.rm = TRUE),
#                                          by = 0.1),
#                                BestTaxon = as.factor(unique(detect_data$BestTaxon)))
# 
# 
# m3.1bpreds <- predict.bam(m3.1b, m3.1b_pred_grid,
#                           se.fit = TRUE,
#                           discrete = TRUE,
#                           n.threads = 40)
# 
# m3.1b_sePreds <- data.frame(m3.1b_pred_grid,
#                             mu   = exp(m3.1bpreds$fit),
#                             low  = exp(m3.1bpreds$fit - 1.96 * m3.1bpreds$se.fit),
#                             high = exp(m3.1bpreds$fit + 1.96 * m3.1bpreds$se.fit))
# 
# save(m3.1bpreds,m3.1b_sePreds, file = "m3.1b_preds.Rdata")
