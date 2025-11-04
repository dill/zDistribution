# simple nimble model with just a spline on depth
# using the jagam object 

library(MCMCvis)
library(boot)
library(tidyverse)
library(mcmcplots)
library(ggplot2)
library(ggdist)
library(dplyr)
library(tidyr)
library(viridis)
library(patchwork)
library(nimble)
library(mgcv)

load("../ProcessedData/detect_data.Rdata")
detect_data$BestTaxon <- as.factor(detect_data$BestTaxon)
load("../ProcessedData/jagam_m1.2.RData")

# check dimensions
X = q1Model_m1.2$jags.data$X
nrow(X)
ncol(X)

# define the model
m1.2_nimble <- nimbleCode({
  eta[1:44496] <- X[1:44496, 1:81] %*% b[1:81] ## linear predictor (b is beta)
  for (i in 1:n) { mu[i] <-  ilogit(eta[i]) } ## expected response
  for (i in 1:n) { y[i] ~ dbin(mu[i],w[i]) } ## response 
  ## Parametric effect priors CHECK tau=1/11^2 is appropriate!
  for (i in 1:1) { b[i] ~ dnorm(0,0.0087) }
  ## prior for ti(depth)... 
  for (i in c(2:5)) { b[i] ~ dnorm(0, lambda[1]) }
  ## prior for ti(BestTaxon)... 
  for (i in c(6:21)) { b[i] ~ dnorm(0, lambda[2]) }
  ## prior for ti(depth,BestTaxon)... 
  K3[1:60,1:60] <- S3[1:60,1:60] * lambda[3]  + S3[1:60,61:120] * lambda[4]
  b[22:81] ~ dmnorm(zero[22:81],K3[1:60, 1:60]) 
  ## smoothing parameter priors CHECK...
  for (i in 1:4) {
    lambda[i] ~ dgamma(.05,.005)
    rho[i] <- log(lambda[i])
  }
  
  
})

# Data
data <- list(y = q1Model_m1.2$jags.data$y,
             X = q1Model_m1.2$jags.data$X)

# Constants
constants <- list(n = q1Model_m1.2$jags.data$n, # number of data points
                  w = q1Model_m1.2$jags.data$w,
                #  S1 = q1Model_m1.2$jags.data$S1,
                  S3 = q1Model_m1.2$jags.data$S3,
                  zero = q1Model_m1.2$jags.data$zero)

# could use fitted values from mgcv as starting values for the betas
m1.2 <- gam(Detected ~ 
              ti(depth, k=5, bs="ts")+
              ti(BestTaxon, k=16, bs="re")+
              ti(depth, BestTaxon, k=c(5, 16), bs=c("ts","re")),
            family = "binomial", data = detect_data,
            method = "REML")
# summary(m1.2)
# 
# m1.2_predictions <- expand_grid(depth = 0:500, BestTaxon = as.factor(unique(detect_data$BestTaxon)))

# response type prediction
#m1.2preds <- predict(m1.2, m1.2_predictions, type = "response", se.fit = TRUE)

#m1.2_predictions$pred <- m1.2preds$fit

#ggplot(m1.2_predictions) +
#  geom_line(aes(x=depth, y = pred))+
#  facet_wrap(~BestTaxon) +
#  theme_bw()

# Initial values
inits <- list(lambda = q1Model_m1.2$jags.ini$lambda, # 2 vec
             # rho = log(q1Model_m1.0$jags.ini$lambda),
              b = q1Model_m1.2$jags.ini$b)

# Run NIMBLE model
nimbleOut_m1.2 <- nimbleMCMC(code = m1.2_nimble, 
                             data = data, 
                             inits = inits,
                             constants = constants,
                             niter = 100000, 
                             nburnin = 75000, 
                             thin = 100, 
                             nchains = 4,
                             summary=TRUE,
                             samplesAsCodaMCMC = TRUE,
                             WAIC = TRUE)

save(nimbleOut_m1.2, file = "./ProcessedData/nimbleOut_m1.2.RData")
load("./ProcessedData/nimbleOut_m1.2.RData")


samps <- do.call(rbind, nimbleOut_m1.2$samples)[,1:81]

bt <- unique(detect_data$BestTaxon)

alldat <- c()

# need to install mcmcr
ss <- list(b=mcmcr::as.mcarray(nimbleOut_m1.2$samples[[1]][,1:81]),
           rho=array(nimbleOut_m1.2$samples[[1]][,82:85], dim=c(250,4,1)))
fakey_bacon <- sim2jam(ss, q1Model_m1.2$pregam)

for(tax in bt){
  pg <- data.frame(depth=seq(0, 500, by=10),
                   BestTaxon=tax)
  Xp <- predict(fakey_bacon, pg, type="lpmatrix")

  lp <- Xp %*% t(samps)

  pg$p <- ilogit(apply(lp, 1, mean))
  pg$lci <- ilogit(apply(lp, 1, quantile, p=0.025))
  pg$uci <- ilogit(apply(lp, 1, quantile, p=0.975))

  alldat <- rbind(alldat, pg)
}

ggplot(alldat) +
  geom_ribbon(aes(x=depth, ymin= lci, ymax = uci), fill = "lightgrey") +
  geom_line(aes(x=depth, y = p)) +
  facet_wrap(~BestTaxon, scale = "free_y") +
  theme_bw()

ggsave(filename = "./Figures/m1.2_splinesonly.png", plot = last_plot(), 
       width = 10, height = 8, units = "in")


## Gelman-Rubin diagnostic
#MCMCsummary(nimbleOut_m1.2$samples)
#
## Visualize MCMC chains
#mcmcplot(nimbleOut_m1.2$samples)
#
#n.post <- 250
#post.samples <- rbind.data.frame(nimbleOut_m1.0$samples$chain1,
#                                 nimbleOut_m1.0$samples$chain2,
#                                 nimbleOut_m1.0$samples$chain3,
#                                 nimbleOut_m1.0$samples$chain4)
#
## reconstruct the model expectation
#
#mu.post <- matrix(rep(0, q1Model_m1.0$jags.data$n*nrow(post.samples)), 
#                  nrow = q1Model_m1.0$jags.data$n)
#
## create a new lp matrix
#
##predict.gam(object = m1.0, type = "lpmatrix")
#
## make a new design matrix just for the depths we want
## NEED TO ADD BestTaxon here
#depth_vals <- data.frame(depth=seq(0, 500, by=10))
#m1_newX_pred <- predict(m1.0, newdata = depth_vals, type = "lpmatrix")
#
## NEED TO MODIFY
#eta.post <- m1_newX_pred %*% t(as.matrix(post.samples[,1:5]))
#mu.post <- ilogit(eta.post)
#
## look at the avg only
#depth_vals$mp <- apply(mu.post, 1, mean)
#plot(depth_vals$depth, depth_vals$mp, type="l", ylim=c(0,1))
#
#mu.post.long <- cbind.data.frame(Depth = depth_vals$depth,
#                                 mu.post) %>%
#  pivot_longer(cols = 2:(ncol(mu.post)+1), 
#               names_to = "Sample", values_to = "PDetect")
#
#mu.post.med <- mu.post.long %>%
#  group_by(Depth) %>%
#  # mean is smoother
#  summarize(Med = mean(PDetect),
#            LCI = quantile(PDetect, 0.025),
#            UCI = quantile(PDetect, 0.975))
#
#p <- ggplot() +
#  geom_ribbon(data = mu.post.med, 
#            aes(x=Depth, ymin= LCI, ymax = UCI), fill = "lightgrey") +
#  geom_line(data = mu.post.med, aes(x=Depth, y = Med)) +
#  theme_bw()
#
#ggsave(plot = p, file = "./Figures/m1.0_nimble.png", width = 4, height = 4, units = "in")
#
## this next bit doesn't work -- not sure what m1.0_sePreds is/was?
#
#names(m1.0_sePreds)[1] <- "Depth"
#m1.0_compare <- left_join(m1.0_sePreds, mu.post.med, by = "Depth")
#
#ggplot(m1.0_compare) +
#  geom_line(aes(x=Depth, y = mu))+
#  geom_line(aes(x=Depth, y = mu_jags), color = "blue")+
#  geom_point(aes(x=Depth, y = Med), color = "green")+
#  ylab("P(Detection)")+
#  xlab("Depth")+
#  theme_bw()
