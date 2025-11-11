# simple nimble model with just a spline on depth
# using the jagam object 

#library(MCMCvis)
#library(boot)
#library(tidyverse)
#library(mcmcplots)
#library(ggplot2)
#library(ggdist)
library(dplyr)
#library(tidyr)
#library(viridis)
#library(patchwork)
library(nimble)
library(mgcv)

setwd("..")
load("./ProcessedData/detect_data.Rdata")
detect_data$BestTaxon <- as.factor(detect_data$BestTaxon)
load("./ProcessedData/jagam_m1.2.RData")

m1.2 <- gam(Detected ~ 
              ti(depth, k=5, bs="ts")+
              ti(BestTaxon, k=16, bs="re")+
              ti(depth, BestTaxon, k=c(5, 16), bs=c("ts","re")),
            family = "binomial", data = detect_data,
            method = "REML")


detect_data <- as.data.frame(ungroup(detect_data))

m1.2jg <- jagam(Detected ~ 
              ti(depth, k=5, bs="ts")+
              ti(BestTaxon, k=16, bs="re")+
              ti(depth, BestTaxon, k=c(5, 16), bs=c("ts","re")),
            family = "binomial", data = detect_data,
            diagonalize=TRUE, file="test.jags")

mm.data <- detect_data
mm.data$EKJPrimer <- mm.data$primer
for (i in 1:nrow(mm.data)){
  mm.data$EKJPrimer[i] <- ifelse(mm.data$primer[i] == "DL" & mm.data$Thaw[i] > 1, "DL2", mm.data$primer[i])
}

# need these to be consecutive
mm.data$sample_numeric <- as.numeric(as.factor(mm.data$NWFSCsampleID))

# create a primer variable
mm.data$primer_numeric <- as.numeric(factor(mm.data$EKJPrimer))

# create a numeric species variable
mm.data$species_numeric <- as.numeric(factor(mm.data$BestTaxon))

# create a variable to hold the species-sample info
mm.data$sp_sample <- paste(mm.data$species_numeric, 
                           mm.data$sample_numeric,
                           sep = "_")
mm.data$sp_sample_numeric <- as.numeric(factor(mm.data$sp_sample))
# 9152 sp-sample combos


# make data, index vectors and constants for nimble work

Y_sample_index <- mm.data$sample_numeric
Y_sp_sample_index <- mm.data$sp_sample_numeric
Y_primer_index <- mm.data$primer_numeric

N <- dim(mm.data)[1]
n_species <- length(unique(mm.data$BestTaxon))
n_samples <- length(unique(mm.data$sample_numeric))
n_primers <- length(unique(mm.data$primer_numeric))

Y <- mm.data$Detected

#############
# SITE x DEPTH MODEL: Site-depth specific occupancy states
#############

sample_info <- mm.data %>%
                group_by(sample_numeric) %>%
                slice(1) %>%
                ungroup() %>%
                select(sample_numeric, depth)

sp_sample_info <- expand.grid("sample_numeric" = sample_info$sample_numeric,
                              "species_numeric" = unique(mm.data$species_numeric))

# get corresponding depth info
sp_sample_info$depth <- sp_sample_info$depth[sp_sample_info$sample_numeric]

# Create sample-species combos with species-depth-site combinations for occupancy states
sp_sample_info$sp_sample_combo <- paste(sp_sample_info$species_numeric, 
                                    sp_sample_info$sample_numeric,
                                    sep = "_")
unique_samples <- sample_info$sample_numeric
unique_sp_samples <- sp_sample_info$sp_sample_combo
n_sp_sample_states <- nrow(sp_sample_info)

# Create mapping from site-depth combo to numeric index
sp_sample_lookup <- data.frame(
  combo = unique_sp_samples,
  index = 1:n_sp_sample_states
)

# Extract site and depth for each site-depth state
sp_sample_lookup$species_numeric <- as.numeric(sapply(strsplit(sp_sample_lookup$combo, 
                                                      "_"), `[`, 1))
sp_sample_lookup$sample_numeric <- as.numeric(sapply(strsplit(sp_sample_lookup$combo, 
                                                    "_"), `[`, 2))

sp_sample_lookup$depth <- sample_info$depth[sp_sample_lookup$sample_numeric]
sp_sample_lookup$BestTaxon <- levels(mm.data$BestTaxon)[sp_sample_lookup$species]

m1.2_newX <- predict(m1.2, newdata = sp_sample_lookup, type = "lpmatrix")


# define the model
m1.2_nimble <- nimbleCode({
  
  # Linear predictor, effect of depth
  # eta has dimensions of # site-depth-species combos
  eta[1:n_sp_sample_states] <- m1.2_newX[1:n_sp_sample_states, 1:81] %*% b[1:81] 
  
  # species-depth probability of occurrence
  prob_sp_depth_occurrence[1:n_sp_sample_states] <- ilogit(eta[1:n_sp_sample_states]) # probability of occurrence at each depth

  # species-depth-site occupancy
  # NOTE need to change inits function to look for any detections BY SPECIES
  for (i in 1:n_sp_sample_states){
    sp_sample_occurrence[i] ~ dbern(prob_sp_depth_occurrence[i]) # OCCUPANCY (UNOBSERVED)
  }

  # Replicate-level likelihood of detecting in each of N lab samples
  for (i in 1:N) {
    Y[i] ~ dbern(sp_sample_occurrence[Y_sp_sample_index[i]] * prob_detection[Y_primer_index[i]]) 
  } # end for N

  ## Parametric effect priors CHECK tau=1/11^2 is appropriate!
  for (i in 1:1) { b[i] ~ dnorm(0,0.0087) }
  ## prior for ti(depth)... 
  for (i in c(2:5)) { b[i] ~ dnorm(0, lambda[1]) }
  ## prior for ti(BestTaxon)... 
  for (i in c(6:21)) { b[i] ~ dnorm(0, lambda[2]) }
  ## prior for ti(depth,BestTaxon)... 
  K3[1:60,1:60] <- S3[1:60,1:60] * lambda[3]  + S3[1:60,61:120] * lambda[4]
  b[22:81] ~ dmnorm(zero[22:81], prec=K3[1:60, 1:60])
  ## smoothing parameter priors CHECK...
  for (i in 1:4) {
    lambda[i] ~ dgamma(.05,.005)
    rho[i] <- log(lambda[i])
  }
  
  # primer detection probability priors
  prob_detection[1] ~ dbeta(1, 1)
  prob_detection[2] ~ dbeta(1, 1)
  prob_detection[3] ~ dbeta(1, 1)
  prob_detection[4] ~ dbeta(1, 1)
  
    
})

# Data
data <- list(Y = q1Model_m1.2$jags.data$y,
             m1.2_newX = m1.2_newX)

# Constants
constants <- list(N = q1Model_m1.2$jags.data$n, # number of data points
                  S3 = q1Model_m1.2$jags.data$S3,
                  zero = q1Model_m1.2$jags.data$zero,
                n_sp_sample_states = n_sp_sample_states,
                Y_primer_index = Y_primer_index,
                Y_sp_sample_index = Y_sp_sample_index)


# 
# m1.2_predictions <- expand_grid(depth = 0:500, BestTaxon = as.factor(unique(detect_data$BestTaxon)))

# response type prediction
#m1.2preds <- predict(m1.2, m1.2_predictions, type = "response", se.fit = TRUE)

#m1.2_predictions$pred <- m1.2preds$fit

#ggplot(m1.2_predictions) +
#  geom_line(aes(x=depth, y = pred))+
#  facet_wrap(~BestTaxon) +
#  theme_bw()

# UPDATED: Using a function for initial values for the second model as well.
inits_fn_sitedepth <- function(){
  # Determine if a biosample had any positive detections by species
  
#  sp_sample_occurrence_raw <- mm.data %>% 
#    ungroup() %>%
#    select(species_numeric, sample_numeric, Detected) %>%
#    group_by(species_numeric, sample_numeric) %>%
#    summarize(DTR = sum(Detected)) 
  
#  sp_sample_occurrence <- ifelse(sp_sample_occurrence_raw$DTR > 0, 1, 0) # 9152
  

  list(
    # Latent states 
    sp_sample_occurrence = rep(1, n_sp_sample_states),
    # Model parameters 
    prob_detection = rep(0.5, n_primers),
    lambda = q1Model_m1.2$jags.ini$lambda,
    b = q1Model_m1.2$jags.ini$b
  )
}

# Run NIMBLE model
nimbleOut_m1.2 <- nimbleMCMC(code = m1.2_nimble, 
                             data = data, 
                             inits = inits_fn_sitedepth,
                             constants = constants,
                             niter = 100000, 
                             nburnin = 75000, 
                             thin = 10, 
                             nchains = 4,
                             summary=TRUE,
                             samplesAsCodaMCMC = TRUE,
                             WAIC = TRUE)

save(nimbleOut_m1.2, file = "./Results/nimbleOut_m1.2+2LevOcc.RData")
# load("./Results/nimbleOut_m1.2+2LevOcc.RData")
#  
# samps <- do.call(rbind, nimbleOut_m1.2$samples)[,1:81]
#  
# bt <- unique(detect_data$BestTaxon)
# 
# alldat <- c()
# 
# # need to install mcmcr
# ss <- list(b=mcmcr::as.mcarray(nimbleOut_m1.2$samples[[1]][,1:81]),
#             rho=array(nimbleOut_m1.2$samples[[1]][,82:85], dim=c(250,4,1)))
#  fakey_bacon <- sim2jam(ss, q1Model_m1.2$pregam)
#  
# for(tax in bt){
#   pg <- data.frame(depth=seq(0, 500, by=10),
#                    BestTaxon=tax)
#   Xp <- predict(fakey_bacon, pg, type="lpmatrix")
# 
#   lp <- Xp %*% t(samps)
# 
#   pg$p <- ilogit(apply(lp, 1, mean))
#   pg$lci <- ilogit(apply(lp, 1, quantile, p=0.025))
#   pg$uci <- ilogit(apply(lp, 1, quantile, p=0.975))
# 
#   alldat <- rbind(alldat, pg)
# }
# 
# ggplot(alldat) +
#   geom_ribbon(aes(x=depth, ymin= lci, ymax = uci), fill = "lightgrey") +
#   geom_line(aes(x=depth, y = p)) +
#   facet_wrap(~BestTaxon) +
#   theme_bw()
# 
# 
# 
# library(bayesplot)
# bayesplot_theme_set(new=theme_minimal())
# 
# mcmc_trace(nimbleOut_m1.2$samples, pars=paste0("prob_detection[",1:3,"]"), transformations=log)
# 
# mcmc_trace(nimbleOut_m1.2$samples, pars=paste0("lambda[",1:4,"]"), transformations=log)
# 
# mcmc_pairs(nimbleOut_m1.2$samples, pars=paste0("lambda[",1:4,"]"), transformations=log)
# 
