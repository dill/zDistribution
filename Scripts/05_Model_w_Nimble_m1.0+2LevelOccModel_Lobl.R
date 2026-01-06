# simple nimble model with just a spline on depth (for Lobl)
# using the jagam object 

library(MCMCvis)
library(boot)
library(tidyverse)
library(ggplot2)
library(ggdist)
library(dplyr)
library(tidyr)
library(viridis)
library(patchwork)
library(nimble)
library(mgcv)
library(PNWColors)

# Import the jagam object we created previously
load("ProcessedData/jagam_m1.0.RData")

# Import the data
load("./ProcessedData/detect_data_lobl.RData")
mm.data <- detect_data_lobl

m1.0 <- gam(Detected ~ s(depth, k = 5, bs = "bs"),  
            diagonalize = TRUE, 
            family = "binomial", data = detect_data_lobl, method="REML", 
            select = TRUE) # to create two lambdas

# because I've removed some of the unique biosample reference numbers with the above 
# methods removal, I need to renumber the unique biosamples so that they are consecutive
# and start with 1. I do this using the factor trick. 
mm.data$unique_biorep_numeric <- as.numeric(as.factor(mm.data$NWFSCsampleID))

# create a site variable that can hold multiple depths
mm.data$site.numeric <- as.numeric(factor(paste0(mm.data$utm.lat, mm.data$utm.lon)))
# there are 182 stations

# create a primer variable
mm.data$primer.numeric <- as.numeric(factor(mm.data$primer))

# pull unique info for each bio sample (can add to these for enviro covariates)
biosamp_dat <- mm.data %>%
  group_by(unique_biorep_numeric) %>%
  slice(1) %>%
  ungroup()
# there are 572 station-depth combos

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

Y <- mm.data$Detected

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

m1_newX <- predict(m1.0, newdata = site_depth_lookup, type = "lpmatrix")

# Create index vector for biosamples to site-depth states
biosamp_dat$site_depth_index <- match(biosamp_dat$site_depth_combo, 
                                      site_depth_lookup$combo)
biosamp_site_depth_index <- biosamp_dat$site_depth_index

### MODEL

m1.0_nimble <- nimbleCode({
  # Priors
  ## Parametric effect priors CHECK tau=1/11^2 is appropriate!
  for (i in 1:1) { b_depth[i] ~ dnorm(0,0.0087) }
  ## prior for s(depth)... 
  for (i in 2:4) { b_depth[i] ~ dnorm(0, lambda[1]) }
  for (i in 5) { b_depth[i] ~ dnorm(0, lambda[2]) }
  ## smoothing parameter priors CHECK...
  for (i in 1:2) {
    lambda[i] ~ dgamma(.05,.005)
    rho[i] <- log(lambda[i])
  }
  
  # primer detection probability priors
  prob_detection[1] ~ dbeta(1, 1)
  prob_detection[2] ~ dbeta(1, 1)
  prob_detection[3] ~ dbeta(1, 1)
  
  # Linear predictor, effect of depth
  # eta has dimensions of # site-depth combos
  eta[1:n_site_depth_states] <- X[1:n_site_depth_states, 1:5] %*% b_depth[1:5] 
    
  # Depth-level occurrence
  for (i in 1:n_site_depth_states) { 
    prob_site_depth_occurrence[i] <- ilogit(eta[i]) # probability of occurrence at each depth
    site_depth_occurrence[i] ~ dbern(prob_site_depth_occurrence[i]) # OCCUPANCY (UNOBSERVED)
  } # end for n_site_depth_states
  
  # Biosample-level capture 
  for (i in 1:n_biosamples) {
    bio_capture[i] <- site_depth_occurrence[biosamp_site_depth_index[i]]
  } # end for n_biosamples
  
  # Replicate-level likelihood of detecting in each of N lab samples
  for (i in 1:N) {
    Y[i] ~ dbern(bio_capture[Y_biosamp_index[i]] * prob_detection[Y_primer_index[i]]) 
  } # end for N
   
}) # end model definition

# Data
data <- list(X = m1_newX,
             Y = Y)

# Constants
constants <- list(  n_sites = n_sites,
                    n_biosamples = n_biosamples,
                    n_site_depth_states = n_site_depth_states,
                    N = q1Model_m1.0$jags.data$n,
                    biosamp_station_index = biosamp_station_index,
                    biosamp_site_depth_index = biosamp_site_depth_index,
                    Y_biosamp_index = Y_biosamp_index,
                    Y_primer_index = Y_primer_index)


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
#    bio_capture = bio_capture_init,
    site_depth_occurrence = site_depth_occurrence_init,
    
    # Model parameters 
#    b_vol = 0,
#    cap_prob_hat = 0,
#    cap_prob_SD = 1,
#    cap_prob_logit_site = rep(0, n_sites),
  #  prob_capture = rep(0.5, n_biosamples),
    prob_detection = rep(0.5, n_primers),
    lambda = m1.0$sp, # 2 vec
    #rho = log(q1Model_m1.0$jags.ini$lambda),
    b_depth = m1.0$coefficients
  )
}

# Run NIMBLE model
nimbleOut_m1.0_2LevelOcc <- nimbleMCMC(code = m1.0_nimble, 
                             data = data, 
                             inits = inits_fn_sitedepth(),
                             constants = constants,
                             niter = 250000, 
                             nburnin = 225000, 
                             thin = 10, 
                             nchains = 4,
                             summary=TRUE,
                             samplesAsCodaMCMC = TRUE,
                             WAIC = TRUE)

save(nimbleOut_m1.0_2LevelOcc, file = "./Results/nimbleOut_m1.0_2LevelOcc_Lobl.RData")

# Gelman-Rubin diagnostic
MCMCsummary(nimbleOut_m1.0_2LevelOcc$samples)

# Visualize MCMC chains
#mcmcplot(nimbleOut_m1.0_2LevelOcc$samples)
#mcmc.output_coda <- as.mcmc.list(lapply(nimbleOut_m1.0_2LevelOcc, as.mcmc))

n.post <- 10000
post.samples <- rbind.data.frame(nimbleOut_m1.0_2LevelOcc$samples$chain1,
                                 nimbleOut_m1.0_2LevelOcc$samples$chain2,
                                 nimbleOut_m1.0_2LevelOcc$samples$chain3,
                                 nimbleOut_m1.0_2LevelOcc$samples$chain4)

post.samples$chain <- c(rep(1, 2500), rep(2, 2500), rep(3, 2500), rep(4, 2500))
post.samples$sample <- rep(1:2500, times = 4)

ggplot(post.samples) +
  geom_line(aes(x = sample, y = log(`lambda[2]`), color = as.factor(chain))) +
  scale_color_discrete()

# check spline values against mgcv fit

post.samples_spline <- post.samples %>%
  select(`b_depth[1]`, `b_depth[2]`, `b_depth[3]`, `b_depth[4]`, `b_depth[5]`,
         `lambda[1]`, `lambda[2]`) %>%
  pivot_longer(cols = 1:7, names_to = "Parameter", values_to = "Est") %>%
  group_by(Parameter) %>%
  summarize(Med = median(Est))

View(data.frame(cbind(post.samples_spline$Parameter, 
                 "MGCV" = c(m1.0$coefficients, m1.0$sp), 
                 "Nimble" = post.samples_spline$Med)))

# smooth of depth

# make a new design matrix just for the depths we want
depth_vals <- data.frame(depth=seq(0, 500, by=10))
m1_newX_pred <- predict(m1.0, newdata = depth_vals, type = "lpmatrix")

eta.post <- m1_newX_pred %*% t(as.matrix(post.samples[,1:5]))
mu.post <- ilogit(eta.post)

# the mean so so smoooth
mu.post.med <- apply(mu.post, 1, mean)

# look at the avg only
goblin <- depth_vals
goblin$mp <- mu.post.med
plot(goblin$depth, goblin$mp, type="l", ylim=c(0,1))


mu.post.long <- cbind.data.frame(Depth = depth_vals$depth,
                                 mu.post) %>%
  pivot_longer(cols = 2:(ncol(mu.post)+1), 
               names_to = "Sample", values_to = "PDetect")

mu.post.median <- mu.post.long %>%
  group_by(Depth) %>%
  # mean is smoother
  summarize(Median = median(PDetect),
            LCI = quantile(PDetect, 0.025),
            UCI = quantile(PDetect, 0.975))

p <- ggplot() +
  geom_ribbon(data = mu.post.median, 
            aes(x=Depth, ymin= LCI, ymax = UCI), fill = "lightgrey") +
  geom_line(data = mu.post.median, aes(x=Depth, y = Median))+
  ylab("Posterior Median P(Occ/Capture)") +
  ylim(c(0,1))+
  theme_bw()

ggsave(plot = p, file = "./Figures/m1.0_nimble2LevOccModel_Lobl.png", 
       width = 4, height = 4, units = "in")

# plot detectability results

post.samples_detectability <- post.samples %>%
  select(`prob_detection[1]`, `prob_detection[2]`, 
         `prob_detection[3]`) %>%
  pivot_longer(cols = 1:3, names_to = "Parameter", values_to = "Est") %>%
  mutate(Primer = case_when(Parameter == "prob_detection[1]" ~ "DL1",
                            Parameter == "prob_detection[2]" ~ "MFU",
                            Parameter == "prob_detection[3]" ~ "MV1"))

save(post.samples_detectability, file = "./Results/post.samples_detectability_m1.0+2LevelOcc_Lobl.RData")

Q2_Detectability <- ggplot(post.samples_detectability) +
  geom_density(aes(x=Est, fill = Primer, color = Primer), alpha = 0.75) +
  scale_fill_manual(values = c(pnw_palette("Cascades",5, type = "discrete")[c(2, 3, 5)],
                                pnw_palette("Sunset",1, type = "discrete"))) +
  scale_color_manual(values = c(pnw_palette("Cascades",5, type = "discrete")[c(2, 3, 5)],
                               pnw_palette("Sunset",1, type = "discrete"))) +
  theme_bw() +
  xlab("P(Detection)") +
  ylab("Posterior Density")

save(Q2_Detectability, file = "./Figures/Q2_DetectabilityByPrimer_Lobl.RData")

ggsave(plot = last_plot(), file = "./Figures/m1.0+2LevelOcc_PDetection_Lobl.png", 
       width = 6, height = 4, units = "in")  

