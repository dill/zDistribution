
# script to plot # of tech reps and p(detection in at least one tech rep) by primer

library(dplyr)
library(ggplot2)
library(PNWColors)

load("./Results/post_pdetect_m1.0+2LevelOcc_wThaw.RData")

# calculation approach

post.detectability <- post_pdetect %>%
  filter(Thaw == 1) %>%
  group_by(Primer) %>%
  summarize(Mean = mean(PDetect),
            LCI = quantile(PDetect, 0.025),
            UCI = quantile(PDetect, 0.975)) %>%
  pivot_longer(cols = 2:4, names_to = "Type", values_to = "P") %>%
  expand_grid("TechReps" = 1:8) %>%
  mutate(POD = 1 - (1-P)^TechReps) %>%
  select(!P) %>%
  pivot_wider(names_from = Type, values_from = POD)

Q2_TechRepsXPrimer <- ggplot(post.detectability) +
  geom_ribbon(aes(x = TechReps, ymin =LCI, ymax = UCI, fill = Primer), alpha = 0.25) +
  geom_line(aes(x=TechReps, y = Mean, color = Primer)) +
  scale_color_manual(values = c(pnw_palette("Cascades",5, type = "discrete")[c(2, 3, 5)],
                                pnw_palette("Sunset",1, type = "discrete"))) +
  scale_fill_manual(values = c(pnw_palette("Cascades",5, type = "discrete")[c(2, 3, 5)],
                                pnw_palette("Sunset",1, type = "discrete"))) +
  theme_bw() +
  geom_hline(yintercept = 0.90, lty = "dashed", color = "grey")+
  xlab("Number of Tech Reps") +
  ylab("Probability of Detection")

ggsave(plot = Q2_TechRepsXPrimer, file = "./Figures/Q2_TechRepsXPrimer.png",
       width = 6, height = 4, units = "in")

save(Q2_TechRepsXPrimer, file = "./Figures/Q2_TechRepsXPrimer.RData")



