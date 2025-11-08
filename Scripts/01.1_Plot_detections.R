### 3D distribution
### Density distribution by depth
### January 2025
### AVC

#### Set up environment --------------------------------------------------------

library(tidyverse)
library(PNWColors)
library(ggridges)

load("./ProcessedData/detect_data.RData")
metadata <- read.csv("./Data/Hake_2019_metadata.csv")
mmEcoEvo <- read.csv("./Data/MM_metadata.csv")

#### Collapse by station/species -----------------------------------------------

detect_by_station <- detect_data %>% 
  group_by(station, depth, BestTaxon) %>% 
  mutate(totReads = sum(nReads)) %>% 
  mutate(detect = case_when(totReads > 0 ~ 1,
                            TRUE ~ 0)) %>% 
  slice_head() %>% 
  ungroup() %>% 
  filter(detect == 1) %>% 
  mutate(depth = case_when(depth > 350~500,
                           depth > 160 & depth < 351~300,
                           depth > 60 & depth < 161~150,
                           depth > 10 & depth < 61~50,
                           depth < 11~0,
                           TRUE~depth))

#### Bubbleplot ----------------------------------------------------------------
detectDepth_bubble <- ggplot(detect_by_station, aes(y = common_name, x = depth, 
                              fill = Broad_taxa, color = Broad_taxa)) +
  geom_count(alpha = 0.7) +
  scale_size_area(max_size = 11) +
  facet_wrap(~Broad_taxa, scales = "free_x", ncol = 1) +
  theme_minimal() + 
  coord_flip(clip = "off") +
  scale_x_reverse() +
  xlab("Sample Depth (m)")+
  ylab("")+
  theme(axis.text.x = element_text(angle = 15, vjust = 1, hjust=0.6)) +
  scale_fill_manual(values = c(pnw_palette("Cascades",5, type = "continuous")[4:5],
                               pnw_palette("Sunset",1, type = "continuous"))) +
  scale_color_manual(values = c(pnw_palette("Cascades",5, type = "continuous")[4:5],
                                pnw_palette("Sunset",1, type = "continuous"))) +
  theme(legend.position = "none",
        strip.text = element_text(size = 12),
        axis.title = element_text(size = 12))

detectDepth_bubble
#### Ridgeplot -----------------------------------------------------------------

lowdetect_subset <- detect_by_station %>% 
  group_by(BestTaxon) %>% 
  mutate(nDetect = n()) %>% 
  filter(nDetect < 3)

detectDepth_ridge <- ggplot(detect_by_station, aes(y = BestTaxon, x = depth, 
                              fill = BestTaxon, color = BestTaxon)) +
  geom_density_ridges(scale = 1, 
                      bandwidth = 100, 
                      jittered_points = TRUE,
                      point_alpha = 1,
                      point_shape = 21,
                      alpha = 0.6) +
  geom_point(data = lowdetect_subset, aes()) +
  theme_minimal() + 
  coord_flip() +
  scale_x_reverse() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_fill_manual(values = c(pnw_palette("Cascades",12, type = "continuous"),
                               pnw_palette("Sunset",12, type = "continuous")[1:12])) +
  scale_color_manual(values = c(pnw_palette("Cascades",12, type = "continuous"),
                                pnw_palette("Sunset",12, type = "continuous")[1:12])) +
  theme(legend.position = "none")

#### Identity ridgeplot --------------------------------------------------------

sumDetect_by_station <- detect_by_station %>% 
  group_by(BestTaxon, depth) %>% 
  mutate(total = sum(detect))

detectDepth_iridge <- ggplot(sumDetect_by_station, aes(y = BestTaxon, x = depth, 
                                                   fill = Broad_taxa, color = Broad_taxa)) +
  geom_ridgeline(data = sumDetect_by_station, aes(height = total/10)) +
  geom_point(data = lowdetect_subset, aes()) +
  facet_wrap(~Broad_taxa, scales = "free") +
  theme_minimal() + 
  coord_flip() +
  scale_x_reverse() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_fill_manual(values = c(pnw_palette("Cascades",3, type = "continuous"),
                               pnw_palette("Sunset",3, type = "continuous")[1:3])) +
  scale_color_manual(values = c(pnw_palette("Cascades",3, type = "continuous"),
                                pnw_palette("Sunset",3, type = "continuous")[1:3])) +
  theme(legend.position = "none")

#### Save figures -------------------------------------------------------------

save(detectDepth_bubble, detectDepth_ridge, file = "./Figures/detectDepth_plots.Rdata")

pdf(file = "./Figures/detectDepth_bubble.pdf")
detectDepth_bubble
dev.off()
