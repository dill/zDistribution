# explore effect of freeze-thaw on probability of detection

library(tidyverse)

load("./ProcessedData/detect_data_allcet.Rdata")

summarize_thaw <- detect_data_allcet %>%
  filter(!is.na(Thaw)) %>%
  group_by(primer, Thaw, Plate) %>%
  summarize(POD = sum(DetectAny)/length(DetectAny), n = length(DetectAny))

ggplot(summarize_thaw) +
  geom_col(aes(x= Thaw, y = POD, fill = n)) +
  geom_text(aes(x=Thaw, y = 0.2, label = paste("n = ", n)))+
  facet_grid(Plate~primer, scales = "free") +
  theme_bw()+
  scale_fill_gradient2(low = "lightgrey",
                        midpoint = 327,
                        mid = "darkgrey",
                        high = "black",
                        space="Lab")

ggsave(plot = last_plot(), file = "./Figures/detectThaw.png", 
       width = 8, height = 6, units = "in")

detect_data_allcet %>% filter(primer == "MFU" & Thaw == 5) %>% pull(Sample_name)

### Repeat with "clean" dataset ------------------------------------------------

load("./ProcessedData/detect_data_clean.RData")

detect_data_allcet_clean <- detect_data_clean %>%
  select(SampleUID, Sample_name, run, primer, NWFSCsampleID,
         dilution, techRep, seqRep, Detected, Plate, Thaw,
         diluti0n, PopID, sample, station, Niskin, depth, volume,
         date, year, month, day, transect, lat, lon, water.depth,
         bathy.bottom.depth, bottom.depth.consensus, utm.lon, utm.lat) %>%
  group_by(SampleUID, Sample_name, run, primer, NWFSCsampleID,
           dilution, techRep, seqRep, Plate, Thaw,
           diluti0n, PopID, sample, station, Niskin, depth, volume,
           date, year, month, day, transect, lat, lon, water.depth,
           bathy.bottom.depth, bottom.depth.consensus, utm.lon, utm.lat) %>%
  summarize(DetectAny = ifelse(sum(Detected) == 0, 0, 1))


summarize_thaw_clean <- detect_data_allcet_clean %>%
  filter(!is.na(Thaw)) %>%
  group_by(primer, Thaw, Plate) %>%
  summarize(POD = sum(DetectAny)/length(DetectAny), n = length(DetectAny))

ggplot(summarize_thaw_clean) +
  geom_col(aes(x= Thaw, y = POD, fill = n)) +
  geom_text(aes(x=Thaw, y = 0.2, label = paste("n = ", n)))+
  facet_grid(~primer, nrow = 3, scales = "free_y") +
  theme_bw()+
  scale_fill_gradient2(low = "lightgrey",
                       midpoint = 327,
                       mid = "darkgrey",
                       high = "black",
                       space="Lab")

ggsave(plot = last_plot(), file = "./Figures/detectThaw_clean.png", 
       width = 8, height = 6, units = "in")

