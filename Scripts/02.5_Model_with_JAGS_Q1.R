#### 3D Distribution
#### Export best Q1 model to BUGS
#### February2025
#### EKJ&AVC

library(mgcv)
library(tidyverse)
library(rjags)

load("./ProcessedData/detect_data_allcet.RData")
load("./ProcessedData/H1models.Rdata")

### Aggregate model AIC --------------------------------------------------------

#modelAIC <- AIC(m1.0, m1.1, m1.2, m1.2a, m1.2b, m1.2c, m1.2d, m1.2e, m1.2g)

### Build jags model -----------------------------------------------------------

# simplest model to start with
# NOTE b splines are probably better than thin plate splines
# use eigen effect, more spread out
# bs = bs (b-splines better in Bayesian, bc more local, sampling would be more efficient)
q1Model_m1.0 <- jagam(DetectAny ~ s(depth, k = 5, bs = "bs"), 
                      diagonalize = TRUE, 
                      family = "binomial", data = detect_data_allcet,
                      file = "./ProcessedData/m1.0.jag")

save(q1Model_m1.0, file = "ProcessedData/jagam_m1.0.RData")

# depth by species
# q1Model_m1.2 <- jagam(Detected ~ s(depth, by = as.factor(BestTaxon)), 
#                            family = "binomial", data = detect_species_meta,
#                            file = "./ProcessedData/m1.2.jag")

### Run jags model -------------------------------------------------------------

# m1.0
jm <-jags.model("./ProcessedData/m1.0.jag",
                data=q1Model_m1.0$jags.data,
                inits=q1Model_m1.0$jags.ini,
                n.chains=4) # changed to 4 so we can assess convergence

list.samplers(jm)
# simplest model has b and lambda
sam <- jags.samples(jm,c("b","lambda"),n.iter=10000,thin=10) # takes aaages
jam <- sim2jam(sam, q1Model_m1.0$pregam)
plot(jam)  # this is just the spline on depth

pd <- data.frame(depth = 0:500)
fv <- predict(jam,newdata=pd, scale = "response")
plot(pd$depth, exp(fv), type = "l", ylim = c(0,1))
# 
# # m1.2 (not run yet as it seriously takes ages)
# jm <-jags.model("./ProcessedData/m1.2.jag",
#                 data=q1Model_m1.2$jags.data,
#                 inits=q1Model_m1.2$jags.ini,
#                 n.chains=4, n.adapt = 2500) # changed to 4 so we can assess convergence
# 
# list.samplers(jm)
# # simplest model has b and lambda
# sam <- jags.samples(jm,c("b","lambda"), n.iter=25000, thin=10) # takes aaages
# jam <- sim2jam(sam, q1Model_m1.0$pregam)
# plot(jam)  # this is just the spline on depth
