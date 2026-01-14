
# script to plot # of tech reps and p(detection in at least one tech rep) by primer

library(dplyr)
library(ggplot2)
library(PNWColors)

load("./Results/post_pdetect_m1.0+2LevelOcc_wThaw.RData")

mean.detectability <- post_pdetect %>%
  filter(Thaw == 1) %>%
  group_by(Primer) %>%
  summarize(Mean = mean(PDetect))

set.seed(20251112)

rdf <- data.frame()

for (p in 1:3){
  
  # simulate tech reps
  
  One <- rbinom(n = 10000, size = 1, mean.detectability$Mean[p])
  Two <- rbinom(n = 10000, size = 1, mean.detectability$Mean[p])
  Three <- rbinom(n = 10000, size = 1, mean.detectability$Mean[p])
  Four <- rbinom(n = 10000, size = 1, mean.detectability$Mean[p])
  Five <- rbinom(n = 10000, size = 1, mean.detectability$Mean[p])
  Six <- rbinom(n = 10000, size = 1, mean.detectability$Mean[p])
  
  # store info
  
  df <- data.frame("Primer" = unique(mean.detectability$Primer[p]), 
                   "One" = One, 
                   "Two" = Two, 
                   "Three" = Three, 
                   "Four" = Four,
                   "Five" = Five, 
                   "Six" = Six)
  
  rdf <- rbind.data.frame(rdf, df)
  
  } # end for p

rdf <- rdf %>%
  mutate(SOne = One) %>%
  mutate(STwo = case_when(One + Two >= 1 ~ 1,
                          One + Two == 0 ~ 0)) %>%
  mutate(SThree = case_when(One + Two + Three >= 1 ~ 1,
                            One + Two + Three == 0 ~ 0)) %>%
  mutate(SFour = case_when (One + Two + Three + Four >= 1 ~ 1,
                            One + Two + Three + Four == 0 ~ 0)) %>%
  mutate(SFive = case_when(One + Two + Three + Four + Five >= 1 ~ 1,
                           One + Two + Three + Four + Five == 0 ~ 0)) %>%
  mutate(SSix = case_when(One + Two + Three + Four + Five + Six >= 1 ~ 1,
                           One + Two + Three + Four + Five + Six == 0 ~ 0)) %>%
  group_by(Primer) %>%
  summarize(POne = sum(SOne)/length(SOne),
            PTwo = sum(STwo)/length(STwo),
            PThree = sum(SThree)/length(SThree),
            PFour = sum(SFour)/length(SFour),
            PFive = sum(SFive)/length(SFive),
            PSix = sum(SSix)/length(SSix)) %>%
  pivot_longer(cols = 2:7, names_to = "NReps", values_to = "P(Detection)") %>%
  mutate(TechReps = case_when(NReps == "POne" ~ 1,
                              NReps == "PTwo" ~ 2,
                              NReps == "PThree" ~ 3,
                              NReps == "PFour" ~ 4,
                              NReps == "PFive" ~ 5, 
                              NReps == "PSix" ~ 6))


Q2_TechRepsXPrimer <- ggplot(rdf) +
  geom_line(aes(x=TechReps, y = `P(Detection)`, color = Primer)) +
  scale_color_manual(values = c(pnw_palette("Cascades",5, type = "discrete")[c(2, 3, 5)],
                                pnw_palette("Sunset",1, type = "discrete"))) +
  theme_bw() +
  xlab("Number of Tech Reps") +
  ylab("Probability of Detection")

ggsave(plot = Q2_TechRepsXPrimer, file = "./Figures/Q2_TechRepsXPrimer.png",
       width = 6, height = 4, units = "in")

save(Q2_TechRepsXPrimer, file = "./Figures/Q2_TechRepsXPrimer.RData")
  
