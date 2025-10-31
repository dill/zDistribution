#### Potential Shark Prey
#### AVC March 2025

library(tidyverse)
library(PNWColors)
library(ggOceanMaps)

load("ProcessedData/detect_data.RData")
metadata <- read.csv("./Data/Hake_2019_metadata.csv")
mmEcoEvo <- read.csv("./Data/MM_metadata.csv")
all_detect_data <- read.csv("Data/M3_compiled_taxon_table_wide.csv")
tax_data <- read.csv("Data/MFU_database.csv") %>% 
  bind_rows(read.csv("Data/MV1_database.csv")) %>% 
  distinct(BestTaxon, .keep_all = TRUE)

#### Look @ bottom depth of deep detections ------------------------------------

detect_data_deep <- detect_data %>% 
  filter(depth > 250) %>% 
  filter(Detected == 1) %>% 
  mutate(distToBottom = bottom.depth.consensus - depth) %>% 
  group_by(station, depth, BestTaxon) %>% 
  slice_head() %>% 
  filter(Broad_taxa != "Beaked whale")

deep_summary <- detect_data_deep %>% 
  group_by(Broad_taxa) %>% 
  summarise(nObs = n())

ggplot(detect_data_deep, aes(x= distToBottom, fill = Broad_taxa)) +
  geom_histogram(binwidth = 50) +
  facet_wrap(~Broad_taxa, scales = "free_x") +
  coord_cartesian(xlim = c(0,2000)) +
  theme_minimal() +
  scale_fill_manual(values = c(pnw_palette("Cascades",3, type = "continuous"),
                               pnw_palette("Sunset",3, type = "continuous")[1:3]))

# Plot the location of likely whalefall (dist to bottom < 100m)
potential_whalefall <- detect_data_deep %>% 
  filter(distToBottom < 100 & distToBottom > 0)

whalefall_summary <- potential_whalefall %>% 
  group_by(Broad_taxa) %>% 
  summarize(nObs = n())

basemap(limits = c(min(potential_whalefall$lon)-0.9,
                   max(potential_whalefall$lon)+0.25,
                   min(potential_whalefall$lat)-0.25,
                   max(potential_whalefall$lat)+0.4),
        bathy.style = "rcb", crs = 4236,
        rotate = FALSE) +
  ggspatial::geom_spatial_point(data = potential_whalefall, 
                                aes(x = lon, y = lat, color = BestTaxon),
                                size = 4, alpha = 0.8, fill = "transparent") +
  guides(fill = "none") +
  scale_color_manual(values = c(pnw_palette("Sunset",6, type = "continuous"))) +
  theme(legend.key = element_blank()) +
  theme_minimal() +
  geom_label(data = potential_whalefall, aes(x = lon-0.5, y = lat+0.25, label = round(distToBottom, digits = 2)))

# Look for killer whale prey

detect_data_deep_detect_stations <- detect_data %>% 
  filter(station %in% detect_data_deep$station) %>% 
  group_by(station) %>% 
  filter(any(BestTaxon == "Orcinus orca")) %>% 
  distinct(station, depth, BestTaxon)

deepDetect_noFall <- detect_data_deep %>% 
  filter(distToBottom > 100) 

detectKW <- detect_data %>% 
  ungroup() %>% 
  filter(Detected == 1) %>% 
  filter(BestTaxon == "Orcinus orca") %>% 
  distinct(BestTaxon, station) %>% 
  left_join(deepDetect_noFall, by = "station", multiple = "all") %>% 
  filter(!is.na(BestTaxon.y)) %>% 
  filter(BestTaxon.y != "Orcinus orca")

basemap(limits = c(min(deepDetect_noFall$lon)-0.25,
                   max(deepDetect_noFall$lon)+0.25,
                   min(deepDetect_noFall$lat)-0.25,
                   max(deepDetect_noFall$lat)+0.25),
        bathy.style = "rcb", crs = 4236,
        rotate = FALSE) +
  ggspatial::geom_spatial_point(data = deepDetect_noFall, 
                                aes(x = lon, y = lat, color = BestTaxon),
                                size = 4, alpha = 0.8, fill = "transparent") +
  ggspatial::geom_spatial_point(data = detectKW, 
                                aes(x = lon, y = lat),
                                size = 4, alpha = 0.8, color = "red", fill = "transparent") +
  guides(fill = "none") +
  scale_color_manual(values = c(pnw_palette("Sunset",12, type = "continuous"))) +
  theme(legend.key = element_blank()) +
  theme_minimal() 

#### Get shark data-------------------------------------------------------------

all_detect_data_long <- all_detect_data %>% 
  pivot_longer(cols = -BestTaxon, names_to = "Sample_name", values_to = "nReads") %>% 
  filter(nReads > 1) %>% 
  left_join(tax_data, by = c("BestTaxon" = "BestTaxon"))

shark_detect <- all_detect_data_long %>% 
  filter(Class == "Chondrichthyes") %>% 
  separate(Sample_name, into = c("Plate", "primer", "pop", "sample_no", "dilution", "techRep"), sep = "\\.") %>% 
  separate(techRep, into = c("techRep", NA), sep = "_") %>% 
  unite(col = "Sample_name", primer:techRep, sep = "-") %>% 
  filter(BestTaxon %in% c("Carcharodon carcharias", "Galeocerdo cuvier", 
                          "Carcharhinus obscurus", "Carcharhinus leucas",
                          "Isurus oxyrinchus", "Prionace glauca")) %>% 
  separate(Sample_name, into = c("primer", "pop", "sample_no","dilution", "techRep"), sep = "-") %>% 
  unite("SampleID", pop:sample_no, sep = "-") %>% 
  distinct(SampleID, BestTaxon)

#### sharks with deep detections -----------------------------------------------

shark_prey <- deepDetect_noFall %>% 
  separate(Sample_name, into = c("primer", "pop", "sample_no","dilution", "techRep"), sep = "-") %>% 
  unite("SampleID", pop:sample_no, sep = "-") %>% 
  left_join(shark_detect, by = "SampleID") %>% 
  filter(!(is.na(BestTaxon.y))) %>% 
  filter(Broad_taxa != "Beaked whale") %>% 
  left_join(metadata, by = c("SampleID" = "sampleID"))

basemap(limits = c(min(deepDetect_noFall$lon)-0.25,
                   max(deepDetect_noFall$lon)+0.25,
                   min(deepDetect_noFall$lat)-0.25,
                   max(deepDetect_noFall$lat)+0.25),
        bathy.style = "rcb", crs = 4236,
        rotate = FALSE) +
  ggspatial::geom_spatial_point(data = deepDetect_noFall, 
                                aes(x = lon, y = lat, color = BestTaxon),
                                size = 4, alpha = 0.8, fill = "transparent") +
  ggspatial::geom_spatial_point(data = shark_prey, 
                                aes(x = lon.x, y = lat.x),
                                size = 4, alpha = 0.8, color = "red", fill = "transparent") +
  guides(fill = "none") +
  scale_color_manual(values = c(pnw_palette("Sunset",12, type = "continuous"))) +
  theme(legend.key = element_blank()) +
  theme_minimal() 

#### Filter out likely whalefalls from full species x station df ---------------

# sum(detect_species_meta$Detected)
# 
# detect_species_meta <- detect_species_meta %>% 
#   mutate(distToBottom = bottom.depth.consensus - depth) %>% 
#   mutate(Detected = ifelse(distToBottom < 100 & depth > 400, 0, Detected))
# 
# sum(detect_species_meta$Detected)

#### Resave data with whalefall filtered out -----------------------------------

save(potential_whalefall,
     file = "./ProcessedData/detect_species_meta.RData")
