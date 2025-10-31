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

load("./ProcessedData/detect_data_allcet.RData")
load("ProcessedData/jagam_m1.0.RData")

# check dimensions
X = q1Model_m1.0$jags.data$X
nrow(X)
ncol(X)

# define the model
m1.0_nimble <- nimbleCode({
    eta[1:2781] <- X[1:2781, 1:5] %*% b_depth[1:5] ## linear predictor (b is beta)
    for (i in 1:n) { mu[i] <-  ilogit(eta[i]) } ## expected response
    for (i in 1:n) { y[i] ~ dbin(mu[i],w[i]) } ## response 
    ## Parametric effect priors CHECK tau=1/7.3^2 is appropriate!
    for (i in 1:1) { b_depth[i] ~ dnorm(0,0.019) }
    ## prior for s(depth)... 
    for (i in c(2:4)) { b_depth[i] ~ dnorm(0, lambda[1]) }
    for (i in c(5)) { b_depth[i] ~ dnorm(0, lambda[2]) }
    ## smoothing parameter priors CHECK...
    for (i in 1:2) {
      lambda[i] ~ dgamma(.05,.005)
      rho[i] <- log(lambda[i])
    }
})

# Data
data <- list(y = q1Model_m1.0$jags.data$y,
             X = q1Model_m1.0$jags.data$X)

# Constants
constants <- list(n = q1Model_m1.0$jags.data$n, # number of data points
                  w = q1Model_m1.0$jags.data$w)

# use fitted values from mgcv as starting values for the betas

m1.0 <- gam(DetectAny ~ s(depth, k = 5, bs = "bs"),  
            diagonalize = TRUE, 
            family = "binomial", data = detect_data_allcet, method="REML", 
            select = TRUE) # to create two lambdas
summary(m1.0)
nd <- data.frame(depth = seq(0, 500, by = 10))
newy <- predict.gam(m1.0, newdata = nd, type = "response")
plot(nd$depth, newy, type = "l", ylim = c(0,1))

# Initial values
inits <- list(lambda = q1Model_m1.0$jags.ini$lambda, # 2 vec
             # rho = log(q1Model_m1.0$jags.ini$lambda),
              b_depth = q1Model_m1.0$jags.ini$b)

# Run NIMBLE model
nimbleOut_m1.0 <- nimbleMCMC(code = m1.0_nimble, 
                             data = data, 
                             inits = inits,
                             constants = constants,
                             niter = 50000, 
                             nburnin = 10000, 
                             thin = 100, 
                             nchains = 4,
                             summary=TRUE,
                             samplesAsCodaMCMC = TRUE,
                             WAIC = TRUE)

save(nimbleOut_m1.0, file = "./ProcessedData/nimbleOut_m1.0.RData")

# Gelman-Rubin diagnostic
MCMCsummary(nimbleOut_m1.0$samples)

# Visualize MCMC chains
mcmcplot(nimbleOut_m1.0$samples)

n.post <- 2000
post.samples <- rbind.data.frame(nimbleOut_m1.0$samples$chain1,
                                 nimbleOut_m1.0$samples$chain2,
                                 nimbleOut_m1.0$samples$chain3,
                                 nimbleOut_m1.0$samples$chain4)

# reconstruct the model expectation

mu.post <- matrix(rep(0, q1Model_m1.0$jags.data$n*nrow(post.samples)), 
                  nrow = q1Model_m1.0$jags.data$n)

# create a new lp matrix

#predict.gam(object = m1.0, type = "lpmatrix")

# make a new design matrix just for the depths we want
depth_vals <- data.frame(depth=seq(0, 500, by=10))
m1_newX_pred <- predict(m1.0, newdata = depth_vals, type = "lpmatrix")

eta.post <- m1_newX_pred %*% t(as.matrix(post.samples[,1:5]))
mu.post <- ilogit(eta.post)

# the mean so so smoooth
mu.post.med <- 

# look at the avg only
depth_vals$mp <- apply(mu.post, 1, mean)
plot(depth_vals$depth, depth_vals$mp, type="l", ylim=c(0,1))

mu.post.long <- cbind.data.frame(Depth = depth_vals$depth,
                                 mu.post) %>%
  pivot_longer(cols = 2:(ncol(mu.post)+1), 
               names_to = "Sample", values_to = "PDetect")

mu.post.med <- mu.post.long %>%
  group_by(Depth) %>%
  # mean is smoother
  summarize(Med = mean(PDetect),
            LCI = quantile(PDetect, 0.025),
            UCI = quantile(PDetect, 0.975))

p <- ggplot() +
  geom_ribbon(data = mu.post.med, 
            aes(x=Depth, ymin= LCI, ymax = UCI), fill = "lightgrey") +
  geom_line(data = mu.post.med, aes(x=Depth, y = Med)) +
  theme_bw()

ggsave(plot = p, file = "./Figures/m1.0_nimble.png", width = 4, height = 4, units = "in")

# this next bit doesn't work -- not sure what m1.0_sePreds is/was?

names(m1.0_sePreds)[1] <- "Depth"
m1.0_compare <- left_join(m1.0_sePreds, mu.post.med, by = "Depth")

ggplot(m1.0_compare) +
  geom_line(aes(x=Depth, y = mu))+
  geom_line(aes(x=Depth, y = mu_jags), color = "blue")+
  geom_point(aes(x=Depth, y = Med), color = "green")+
  ylab("P(Detection)")+
  xlab("Depth")+
  theme_bw()
