library(MCMCvis)
library(boot)
library(tidyverse)
library(mcmcplots)
library(ggplot2)
library(ggdist)
library(dplyr)
library(tidyr)
library(viridis)
library(ggokabeito)
library(patchwork)
library(nimble)

source("./Scripts/attach.nimble_v2.R")

# Import the data
load("./ProcessedData/detect_data_allcet.RData")
mm.data <- detect_data_allcet
mm.data$EKJPrimer <- mm.data$primer
for (i in 1:nrow(mm.data)){
mm.data$EKJPrimer[i] <- ifelse(mm.data$primer[i] == "DL" & mm.data$Thaw[i] > 1, "DL2", mm.data$primer[i])
}

# collapse data so it's the probability of detecting *any* cetacean

mm.data.collapse <- mm.data %>%
  group_by()

# because I've removed some of the unique biosample reference numbers with the above 
# methods removal, I need to renumber the unique biosamples so that they are consecutive
# and start with 1. I do this using the factor trick. 
mm.data$unique_biorep_numeric <- as.numeric(as.factor(mm.data$NWFSCsampleID))

# create a site variable that can hold multiple depths
mm.data$site.numeric <- as.numeric(factor(paste0(mm.data$utm.lat, mm.data$utm.lon)))

# create a primer variable
mm.data$primer.numeric <- as.numeric(factor(mm.data$EKJPrimer))

# pull unique info for each bio sample (can add to these for enviro covariates)
biosamp_dat <- mm.data %>%
  group_by(unique_biorep_numeric) %>%
  slice(1) %>%
  ungroup()

# make data, index vectors and constants for nimble work
biosamp_station_index <- biosamp_dat$site.numeric
#biosamp_method_index <- biosamp_dat$Collection_method_numeric
Y_biosamp_index <- mm.data$unique_biorep_numeric
biosamp_Volume_filt_mL <- as.numeric(biosamp_dat$volume) - mean(as.numeric(biosamp_dat$volume)) #centered water volumes
biosamp_depth_Depth_m <- as.numeric(biosamp_dat$depth) - mean(as.numeric(biosamp_dat$depth)) #centered sample depths
Y_primer_index <- mm.data$primer.numeric

N <- dim(mm.data)[1]
n_sites <- length(unique(mm.data$site.numeric))
n_biosamples <- length(unique(mm.data$unique_biorep_numeric))
#n_methods <- length(unique(mm.data$Collection_method_numeric ))
n_primers <- length(unique(mm.data$primer.numeric))

Y <- mm.data$DetectAny

#############
# SITE x DEPTH MODEL: Site-depth specific occupancy states
#############

# Create site-depth combinations for occupancy states
biosamp_dat$site_depth_combo <- paste(biosamp_dat$site.numeric, 
                                      round(biosamp_dat$depth, 1), 
                                      sep = "_")
unique_site_depth <- unique(biosamp_dat$site_depth_combo)
n_site_depth_states <- length(unique_site_depth)

# Create mapping from site-depth combo to numeric index
site_depth_lookup <- data.frame(
  combo = unique_site_depth,
  index = 1:n_site_depth_states
)

# Extract site and depth for each site-depth state
site_depth_lookup$site <- as.numeric(sapply(strsplit(site_depth_lookup$combo, 
                                                     "_"), `[`, 1))
site_depth_lookup$depth <- as.numeric(sapply(strsplit(site_depth_lookup$combo, 
                                                      "_"), `[`, 2))

# Create index vector for biosamples to site-depth states
biosamp_dat$site_depth_index <- match(biosamp_dat$site_depth_combo, 
                                      site_depth_lookup$combo)
biosamp_site_depth_index <- biosamp_dat$site_depth_index

# Center depths for site-depth states
site_depth_depths <- site_depth_lookup$depth - mean(site_depth_lookup$depth)

# Define the site x depth NIMBLE model (UNCHANGED)
edna_code_site_depth_occupancy <- nimbleCode({
  # Shared intercept for all sites
  intercept ~ dnorm(0, 1)
  
  # Site-depth state occupancy probabilities
  for (i in 1:n_site_depth_states) {
    logit(prob_site_depth_occurrence[i]) <- intercept + 
      b_depth_occ * site_depth_depths[i]
    site_depth_occurrence[i] ~ dbern(prob_site_depth_occurrence[i])
  }
  
  # Capture probability parameters
  cap_prob_hat ~ dnorm(0, 1.7)
  cap_prob_SD ~ dexp(1)
  
  # Site-level capture probability variation
  for (i in 1:n_sites) {
    cap_prob_logit_site[i] ~ dnorm(cap_prob_hat, cap_prob_SD)
  }
  
  # Biosample-level capture probabilities
  for (i in 1:n_biosamples) {
    logit(prob_capture[i]) <- cap_prob_logit_site[biosamp_station_index[i]] +
      b_vol * biosamp_Volume_filt_mL[i] 
    
    bio_capture[i] ~ dbern(site_depth_occurrence[biosamp_site_depth_index[i]] *
                             prob_capture[i])
  }
  
  # Likelihood of detecting in each lab sample
  for (i in 1:N) {
    Y[i] ~ dbern(bio_capture[Y_biosamp_index[i]] * prob_detection[Y_primer_index[i]]) 
  }
  
  # Priors
  b_depth_occ ~ dnorm(0, 1)
  b_vol ~ dnorm(0, 1)
  prob_detection[1] ~ dbeta(1, 1)
  prob_detection[2] ~ dbeta(1, 1)
  prob_detection[3] ~ dbeta(1, 1)
  prob_detection[4] ~ dbeta(1, 1)

})

# Define constants for the site x depth model
# UPDATED: Removed site_depth_sites as it was not used in the nimble model
constants_site_depth <- list(
  n_sites = n_sites,
  n_biosamples = n_biosamples,
  n_site_depth_states = n_site_depth_states,
  N = N,
  biosamp_station_index = biosamp_station_index,
  biosamp_site_depth_index = biosamp_site_depth_index,
  Y_biosamp_index = Y_biosamp_index,
  Y_primer_index = Y_primer_index,
  site_depth_depths = site_depth_depths
)

data_site_depth <- list(
  Y = Y,
  biosamp_Volume_filt_mL = biosamp_Volume_filt_mL
)

# UPDATED: Using a function for initial values for the second model as well.
inits_fn_sitedepth <- function(){
  # Determine if a biosample had any positive detections
  max_obs_per_biosample <- tapply(Y, Y_biosamp_index, max)
  bio_capture_init <- rep(0, n_biosamples)
  bio_capture_init[as.numeric(names(max_obs_per_biosample))] <- max_obs_per_biosample
  
  # Determine if a site-depth combo had any positive biosamples
  site_depth_occurrence_init_raw <- tapply(bio_capture_init, 
                                           biosamp_site_depth_index, 
                                           max)
  site_depth_occurrence_init <- rep(0, n_site_depth_states)
  site_depth_occurrence_init[as.numeric(names(site_depth_occurrence_init_raw))] <- site_depth_occurrence_init_raw
  
  list(
    # Latent states (from before)
    bio_capture = bio_capture_init,
    site_depth_occurrence = site_depth_occurrence_init,
    
    # Model parameters (previously incomplete)
    intercept = 0,
    b_depth_occ = 0,
    b_vol = 0,
    cap_prob_hat = 0,
    cap_prob_SD = 1,
    cap_prob_logit_site = rep(0, n_sites),
#    b_meth = c(NA, 0, 0), # NA for the first fixed element, 0 for the others
    prob_detection = rep(0.5, n_primers)
  )
}


# Run the site x depth NIMBLE model
edna_site_depth_occupancy.run <- nimbleMCMC(code = edna_code_site_depth_occupancy, 
                                            constants = constants_site_depth, 
                                            data = data_site_depth, 
                                            inits = inits_fn_sitedepth, # Using the function here
                                            niter = 200000, 
                                            nburnin = 10000, 
                                            thin = 100, 
                                            nchains = 3,
                                            summary = TRUE,
                                            samplesAsCodaMCMC = TRUE,
                                            WAIC = TRUE)
# Diagnostics
MCMCsummary(edna_site_depth_occupancy.run$samples)
mcmcplot(edna_site_depth_occupancy.run$samples, 
         parms = c("b_depth_occ", "b_vol", "intercept"))

# --- Site x Depth Model 4-Panel Plot ---
attach.nimble(edna_site_depth_occupancy.run$samples)
n.post.new <- length(b_depth_occ)

# Depth effect on occupancy
depth_pred_seq <- seq(from = min(site_depth_lookup$depth), 
                      to = max(site_depth_lookup$depth), 
                      length.out = 100)
depth_pred_centered <- depth_pred_seq - mean(site_depth_lookup$depth)

occupancy_pred <- matrix(NA, n.post.new, length(depth_pred_seq))
for (i in 1:length(depth_pred_seq)) {
  for (j in 1:n.post.new) {
    occupancy_pred[j, i] <- inv.logit(intercept[j] + 
                                        b_depth_occ[j] * depth_pred_centered[i])
  }
}

plot_data_occupancy_depth <- data.frame(
  depth = rep(depth_pred_seq, each = n.post.new),
  occupancy = as.vector(occupancy_pred)
)

# Volume effect on capture
vol_pred_seq <- seq(from = min(as.numeric(biosamp_dat$volume)), 
                    to = max(as.numeric(biosamp_dat$volume)), 
                    length.out = 100)
vol_pred_centered <- vol_pred_seq - mean(as.numeric(biosamp_dat$volume))

capture_pred_vol <- matrix(NA, n.post.new, length(vol_pred_seq))
for (i in 1:length(vol_pred_seq)) {
  for (j in 1:n.post.new) {
    capture_pred_vol[j, i] <- inv.logit(cap_prob_hat[j] + 
                                          b_vol[j] * vol_pred_centered[i])
  }
}

plot_data_capture_vol <- data.frame(
  volume = rep(vol_pred_seq, each = n.post.new),
  capture = as.vector(capture_pred_vol)
)


# Detection comparison
df_detection_new <- data.frame(
  Probability = c(prob_detection[, 1], prob_detection[, 2],
                  prob_detection[, 3], prob_detection[, 4]),
  Method = rep(levels(factor(mm.data$EKJPrimer)), each = nrow(prob_detection))
)

# Create 4-panel plot for site x depth model
p1_new <- ggplot(plot_data_occupancy_depth, 
                 aes(x = depth, y = occupancy)) +
  stat_lineribbon(alpha = 0.25, fill = "#4CAF50", color = "#2E7D32", 
                  .width = c(0.25, 0.5, 0.75)) +
  labs(x = "Depth (m)", y = "Probability of Occupancy") +
  theme_minimal()

ggsave(plot = p1_new, file = "./Figures/PatinEtAl_OccupancyByDepth_wThaw.png",
       width = 5, height = 3, units = "in")

p2_new <- ggplot(plot_data_capture_vol, 
                 aes(x = volume, y = capture)) +
  stat_lineribbon(alpha = 0.25, fill = "#EE7AE9", color = "#DA70D6", 
                  .width = c(0.25, 0.5, 0.75)) +
  labs(x = "Volume Filtered (mL)", y = "Probability of Capture") +
  theme_minimal()

ggsave(plot = p2_new, file = "./Figures/PatinEtAl_CaptureByVolume_wThaw.png",
       width = 5, height = 3, units = "in")

p4_new <- ggplot(df_detection_new, 
                 aes(x = Probability, fill = Method, color = Method)) +
  geom_histogram(aes(y = after_stat(density)), 
                 alpha = 0.3, position = "identity", bins = 30) +
  geom_density(size = 1.2) +
  scale_fill_viridis(discrete = TRUE, alpha = 0.3) +
  scale_color_viridis(discrete = TRUE) +
  labs(x = "Detection Probability", y = "Density") +
  theme_bw() + theme(panel.grid = element_blank(), legend.position = "bottom")

ggsave(plot = p4_new, file = "./Figures/PatinEtAl_DetectionByPrimer_wThaw.png",
       width = 5, height = 3, units = "in")


