#### 3D Distribution
#### Models testing Q1: detection probability across depth
#### January2025
#### EKJ&AVC

library(mgcv)
library(tidyverse)
library(PNWColors)

load("./ProcessedData/detect_data.RData")
mmEcoEvo <- read.csv("./Data/MM_metadata.csv")
detect_data$BestTaxon <- as.factor(detect_data$BestTaxon)

# Q1: Does the probability of detecting cetaceans in eDNA samples vary with sample depth?
# H0: Probability of detection does not vary with depth.
# H1: Probability of detection varies across depth agnostic to species or functional group.
# H2: Probability of detection varies across depth according to species or functional group.

### H1: POD by depth alone -----------------------------------------------------
# basic model with no species-specific terms
m1.0 <- gam(Detected ~ s(depth, k = 5), family = "binomial", data = detect_data, method="REML") 
summary(m1.0)
# depth p-value = 2.6e-06
AIC(m1.0)
# AIC 6124

# compare gam to jags version of the same model
# m1.0_sePreds$mu_jags <- exp(predict(jam,newdata=pd, scale = "response"))
# ggplot(m1.0_sePreds) +
#   geom_line(aes(x=depth, y = mu))+
#   geom_line(aes(x=depth, y = mu_jags), color = "blue")+
#   ylab("P(Detection)")+
#   xlab("Depth")+
#   theme_bw()

### H2: POD by depth across species --------------------------------------------
# this model will have a different intercept for each species, but spline will be same shape
m1.1 <- gam(Detected ~ s(depth) + BestTaxon, family = "binomial", data = detect_data, method="REML")
summary(m1.1)

AIC(m1.1)
# AIC = 5604

# separate the depth effect from the taxon effect from the depth-taxon effect. 
m1.2 <- gam(Detected ~ 
              ti(depth, k=5, bs="ts")+
              ti(BestTaxon, k=16, bs="re")+
              ti(depth, BestTaxon, k=c(5, 16), bs=c("ts","re")),
            family = "binomial", data = detect_data,
            method = "REML")

summary(m1.2)

AIC(m1.2)
# AIC 5310

### H2a: POD by depth across taxonomic family ----------------------------------
detect_data$Family <- as.factor(detect_data$Family)
m1.2a <- gam(Detected ~ 
               ti(depth, k=5, bs="ts")+
               ti(Family, k=6, bs="re")+
               ti(depth, Family, k=c(5, 6), bs=c("ts","re")),
             family = "binomial", data = detect_data, method = "REML")
summary(m1.2a)

AIC(m1.2a)
#AIC 5712



### H2b:POD by depth across prey category --------------------------------------
detect_data$Prey.family <- as.factor(detect_data$Prey.family)
m1.2b <-  gam(Detected ~ 
                     ti(depth, k=5, bs="ts")+
                     ti(Prey.family, k=3, bs="re")+
                     ti(depth, Prey.family, k=c(5, 3), bs=c("ts","re")),
                   family = "binomial", data = detect_data,  method = "REML")
summary(m1.2b)
#significant for all three types
AIC(m1.2b)
#AIC 5823

### H2c: POD by time-at-depth --------------------------------------------------
# 
# m1.2c <- gam(Detected ~ s(time_10m), 
#              family = "binomial", data = detect_species_divetime, method="REML")
# summary(m1.2c)
# #p<2e-16
# AIC(m1.2c)
# 
# 
# m1.2c_predictions <- data.frame(time_10m = seq(min(detect_species_divetime$time_10m),
#                                                  max(detect_species_divetime$time_10m), 
#                                                  by = 5))
# 
# m1.2cpreds <- predict(m1.2c, m1.2c_predictions, type = "response", se.fit = TRUE)
# 
# m1.2c_sePreds <- data.frame(m1.2c_predictions,
#                            mu   = exp(m1.2cpreds$fit),
#                            low  = exp(m1.2cpreds$fit - 1.96 * m1.2cpreds$se.fit),
#                            high = exp(m1.2cpreds$fit + 1.96 * m1.2cpreds$se.fit))
# 
# ggplot(m1.2c_sePreds, aes(x = time_10m, y = mu)) +
#   geom_smooth(aes(ymin = low, ymax = high, y = mu), stat = "identity", 
#               alpha = 0.2, color = "#74677e", fill = "#74677e") +
#   ylab("POD") +
#   xlab("Time at depth") +
#   theme_minimal()

### H2d: POD by time-at-depth across Family ------------------------------------
# 
# m1.2d <- gam(Detected ~ s(time_10m, by = as.factor(Family)), 
#              family = "binomial", data = detect_species_divetime, method="REML")
# summary(m1.2d)
# #significant for all but bowhead and grey whale
# AIC(m1.2d)
# #4299 - within ~2 of depth-by-species
# 
# m1.2d_predictions <- expand_grid(time_10m = seq(min(detect_species_divetime$time_10m),
#                                                 max(detect_species_divetime$time_10m), 
#                                                 by = 5),
#                                  Family = as.factor(unique(detect_species_divetime$Family)))
# 
# m1.2dpreds <- predict(m1.2d, m1.2d_predictions, type = "response", se.fit = TRUE)
# 
# m1.2d_sePreds <- data.frame(m1.2d_predictions,
#                             mu   = exp(m1.2dpreds$fit),
#                             low  = exp(m1.2dpreds$fit - 1.96 * m1.2dpreds$se.fit),
#                             high = exp(m1.2dpreds$fit + 1.96 * m1.2dpreds$se.fit))
# 
# ggplot(m1.2d_sePreds, aes(x = time_10m, y = mu, color = Family, fill = Family)) +
#   geom_smooth(aes(ymin = low, ymax = high, y = mu), stat = "identity", 
#               alpha = 0.2) +
#   
#   ylab("POD") +
#   xlab("Time at depth") +
#   scale_fill_manual(values = c(pnw_palette("Bay",8, type = "continuous"))) +
#   scale_color_manual(values = c(pnw_palette("Bay",8, type = "continuous"))) +
#   facet_wrap(~Family, scales = "free_y") +
#   theme_minimal()
# 
# ### H2e: POD by time-at-depth across species -----------------------------------
# ## This one takes a really long time to run.
# 
# m1.2e <- gam(Detected ~ s(time_10m, by = as.factor(BestTaxon)), 
#              family = "binomial", data = detect_species_divetime, method="REML")
# summary(m1.2e)
# #10m bin summary:
# #sig: Bphy, Lobl, Lbor, Mnov, Ppho, Pdal, Scoe
# #<0.09: Zcav, Ggri, Erob, Bbai, Bmus, Bacu, Bmys
# #no sig: Mste, Oorc
# #13.2% deviance explained
# 
# AIC(m1.2e)
# #10m bins = 4075
# #equal bins = 4082
# 
# m1.2e_predictions <- expand_grid(time_10m = seq(min(detect_species_divetime$time_10m),
#                                                   max(detect_species_divetime$time_10m), 
#                                                   by = 5),
#                                  BestTaxon = as.factor(unique(detect_species_divetime$BestTaxon)))
# 
# m1.2epreds <- predict(m1.2e, m1.2e_predictions, type = "response", se.fit = TRUE)
# 
# m1.2e_sePreds <- data.frame(m1.2e_predictions,
#                             mu   = exp(m1.2epreds$fit),
#                             low  = exp(m1.2epreds$fit - 1.96 * m1.2epreds$se.fit),
#                             high = exp(m1.2epreds$fit + 1.96 * m1.2epreds$se.fit)) %>% 
#   left_join(mmEcoEvo, by = c("BestTaxon" = "Species"))
# 
# ggplot(m1.2e_sePreds, aes(x = time_10m, color = Broad_taxa, fill = Broad_taxa)) +
#   geom_smooth(aes(ymin = low, ymax = high, y = mu), stat = "identity", 
#               alpha = 0.2) +
#   ylab("POD") +
#   xlab("Prop time at depth") +
#   facet_wrap(~abbrev, scales = "free_y") +
#   scale_fill_manual(values = c(pnw_palette("Cascades",2, type = "continuous"),
#                                pnw_palette("Sunset",2, type = "continuous"))) +
#   scale_color_manual(values = c(pnw_palette("Cascades",2, type = "continuous"),
#                                 pnw_palette("Sunset",2, type = "continuous"))) +
#   geom_rug(data = detect_species_divetime, aes(x=time_10m), color = "grey")+
#   geom_rug(data = filter(detect_species_divetime, Detected == 1), aes(x=time_10m))+
#   theme_minimal() +
#   theme(legend.position = "bottom")
#   
# 
# save(m1.2e, file = "./ProcessedData/m1.2e.RData")

### H2f: POD by depth + time-at-depth across species ---------------------------
## This one takes a really long time to run.
# m1.2f <- gam(Detected ~ s(time_per_m, by = as.factor(BestTaxon)) +
#                s(depth, by = as.factor(BestTaxon)), 
#              family = "binomial", data = detect_species_divetime)
# summary(m1.2f)
# 
# AIC(m1.2f)
# # 3477.885
# 
# m1.2f_predictions <- expand_grid(depth = 0:500, time_per_m = min(detect_species_divetime$time_per_m):max(detect_species_divetime$time_per_m),
#                                  BestTaxon = as.factor(unique(detect_species_divetime$BestTaxon)))
# m1.2f_predictions$pred <- predict.gam(m1.2f, m1.2f_predictions, type = "response")
# 
# ggplot(m1.2f_predictions) +
#   geom_line(aes(x=time_per_m, y = pred, group = BestTaxon, color = BestTaxon)) +
#   xlab("Time at depth")+
#   ylab("POD")+
#   theme_bw()
# 
# m1.2fpreds <- predict(m1.2f, m1.2f_predictions, se.fit = TRUE)
# 
# m1.2f_sePreds <- data.frame(m1.2f_predictions,
#                             mu   = exp(m1.2fpreds$fit),
#                             low  = exp(m1.2fpreds$fit - 1.96 * m1.2fpreds$se.fit),
#                             high = exp(m1.2fpreds$fit + 1.96 * m1.2fpreds$se.fit)) %>% 
#   left_join(mmEcoEvo, by = c("BestTaxon" = "Species"))
# 
# ggplot(m1.2f_sePreds, aes(x = depth, y = mu, color = BestTaxon, fill = BestTaxon)) +
#   geom_smooth(aes(ymin = low, ymax = high, y = mu), stat = "identity", 
#               alpha = 0.2) +
#   
#   ylab("POD") +
#   xlab("Depth") +
#   scale_fill_manual(values = c(pnw_palette("Bay",24, type = "continuous"))) +
#   scale_color_manual(values = c(pnw_palette("Bay",24, type = "continuous"))) +
#   facet_wrap(Broad_taxa~BestTaxon, scales = "free") +
#   theme_minimal() +
#   theme(legend.position = "none")
# 
# ### H2g: POD by time-at-depth across depths ------------------------------------
# 
# bins <- c(0, 100, 300, 500)
# detect_species_depthbin <- detect_species_divetime %>% 
#   mutate(depth_bin = cut(depth, bins, include.lowest = TRUE))
# 
# m1.2g <- gam(Detected ~ s(time_per_m, by = as.factor(depth_bin)), 
#              family = "binomial", data = detect_species_depthbin)
# summary(m1.2g)
# 
# AIC(m1.2g)
# # 4273
# 
# m1.2g2 <- gam(Detected ~ s(time_per_m, by = depth), 
#              family = "binomial", data = detect_species_divetime)
# summary(m1.2g2)
# 
# AIC(m1.2g2)
# # 3991.547
# m1.2g_predictions <- expand_grid(time_per_m = min(detect_species_depthbin$time_per_m, na.rm = TRUE):max(detect_species_depthbin$time_per_m, na.rm = TRUE),
#                                  depth_bin = as.factor(unique(detect_species_depthbin$depth_bin)))
# m1.2g_predictions$pred <- predict.gam(m1.2g, m1.2g_predictions, type = "response")
# 
# 
# m1.2gpreds <- predict(m1.2g, m1.2g_predictions, se.fit = TRUE)
# 
# m1.2g_sePreds <- data.frame(m1.2g_predictions,
#                             mu   = exp(m1.2gpreds$fit),
#                             low  = exp(m1.2gpreds$fit - 1.96 * m1.2gpreds$se.fit),
#                             high = exp(m1.2gpreds$fit + 1.96 * m1.2gpreds$se.fit))
# 
# #very large uncertainty bands
# ggplot(m1.2g_sePreds, aes(x = time_per_m, y = mu, color = depth_bin, fill = depth_bin)) +
#   geom_smooth(aes(ymin = low, ymax = high, y = mu), stat = "identity", 
#               alpha = 0.2) +
#   
#   ylab("POD") +
#   xlab("Time at depth") +
#   scale_fill_manual(values = c(pnw_palette("Bay",5, type = "continuous"))) +
#   scale_color_manual(values = c(pnw_palette("Bay",5, type = "continuous"))) +
#   facet_wrap(~depth_bin, scales = "free_y") +
#   coord_cartesian(ylim = c(0, 1)) +
#   theme_minimal()

### Aggregate model AIC --------------------------------------------------------

modelAIC <- AIC(m1.0, m1.1, m1.2, m1.2a, m1.2b) %>% 
  rownames_to_column(var = "model")

######## save all models -------------------------------------------------------

save(m1.0, m1.1, m1.2, m1.2a, m1.2b, modelAIC,
     file = "./ProcessedData/H1models_binomial.RData")
# 1.2e
# 1.1
# 1.2