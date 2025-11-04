
# check whether p(detection) across primers is correlated with p(occ/capture) at depth
# using m1.1 2-level occ model
# EKJ Nov 2025

# create model to use for predictions
load("./ProcessedData/detect_data_allcet.RData")
mm.data <- detect_data_allcet

m1.0 <- gam(DetectAny ~ s(depth, k = 5, bs = "bs"),  
            diagonalize = TRUE, 
            family = "binomial", data = detect_data_allcet, method="REML", 
            select = TRUE) # to create two lambdas

load("./Results/nimbleOut_m1.0_2LevelOcc.RData")

n.post <- 10000
post.samples <- rbind.data.frame(nimbleOut_m1.0_2LevelOcc$samples$chain1,
                                 nimbleOut_m1.0_2LevelOcc$samples$chain2,
                                 nimbleOut_m1.0_2LevelOcc$samples$chain3,
                                 nimbleOut_m1.0_2LevelOcc$samples$chain4)

post.samples$chain <- c(rep(1, 2500), rep(2, 2500), rep(3, 2500), rep(4, 2500))
post.samples$sample <- rep(1:2500, times = 4)

# for each sample of each chain, calculate the p(occ/capture) at set depths

# make a new design matrix just for the depths we want
depth_vals <- data.frame(depth=seq(0, 500, by=10))
m1_newX_pred <- predict(m1.0, newdata = depth_vals, type = "lpmatrix")

eta.post <- m1_newX_pred %*% t(as.matrix(post.samples[i,1:5]))
mu.post <- ilogit(eta.post)

df.out <- data.frame("post.sample.index" = 1,
                   "p(detection)[1]" = post.samples$`prob_detection[1]`[1],
                   "p(detection)[2]" = post.samples$`prob_detection[1]`[1],
                   "p(detection)[3]" = post.samples$`prob_detection[1]`[1],
                   "depth" = depth_vals,
                   'p(occ/cap)' = as.numeric(mu.post))

for (i in 2:10000){
  
  eta.post <- m1_newX_pred %*% t(as.matrix(post.samples[i,1:5]))
  mu.post <- ilogit(eta.post)
  
  df <- data.frame("post.sample.index" = i,
                   "p(detection)[1]" = post.samples$`prob_detection[1]`[i],
                   "p(detection)[2]" = post.samples$`prob_detection[1]`[i],
                   "p(detection)[3]" = post.samples$`prob_detection[1]`[i],
                   "depth" = depth_vals,
                   "p(occ/cap)" = as.numeric(mu.post))
  
  df.out <- rbind.data.frame(df.out, df)
}

ggplot(filter(df.out, depth == 250)) +
  geom_point(aes(x= p.detection..1., y = p.occ.cap.))
