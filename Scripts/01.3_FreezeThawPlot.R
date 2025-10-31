# explore effect of freeze-thaw on probability of detection

load("./ProcessedData/detect_data_allcet.Rdata")

summarize_thaw <- detect_data_allcet %>%
  filter(!is.na(Thaw)) %>%
  group_by(primer, Thaw) %>%
  summarize(POD = sum(DetectAny)/length(DetectAny), n = length(DetectAny))

ggplot(summarize_thaw) +
  geom_col(aes(x= Thaw, y = POD, fill = n)) +
  facet_wrap(~primer, nrow = 3) +
  theme_bw()+
  scale_fill_gradient2(low = "lightgrey",
                        midpoint = 327,
                        mid = "darkgrey",
                        high = "black",
                        space="Lab")

ggsave(plot = last_plot(), file = "./Figures/detectThaw.png", 
       width = 6, height = 6, units = "in")
