#### 3D Distribution
#### wrangling data for presence/absence modelling
#### Summer 2025
#### EKJ&AVC

library(tidyverse)

## Get data --------------------------------------------------------------------

metadata <- read.csv("./Data/Hake_2019_metadata.csv")
timeAtDepth <- read.csv("./Data/MM_dive_time_expand.csv")
mmEcoEvo <- read.csv("./Data/MM_metadata.csv")
freezethaw <- read.csv("./Data/HAKE2019_miseq_runs_thaw_v2.csv") # MRS: note I updated this
sample_locs <- read.csv("./Data/HAKE2019_sample_locations_v2.csv") # MRS: note I updated this

detect_data_raw <- read.csv("./Data/M3_compiled_taxon_table_wide.csv") %>% 
  pivot_longer(-c(BestTaxon, Class), names_to = "SampleUID", values_to = "nReads") %>% 
  group_by(SampleUID) %>% 
  mutate(totalReads = sum(nReads)) %>% 
  separate(SampleUID, into = c("Sample_name", NA), remove = FALSE, sep = "_") %>% 
  separate(Sample_name, into = c("run", "primer", "pop", "sample", "dilution", "techRep", "seqRep"), remove = FALSE, sep = "\\.") %>% 
  unite(pop:sample, col = "NWFSCsampleID", sep = "-") %>% 
  mutate(techRep = as.numeric(techRep)) %>% 
  mutate(run = gsub("a","",run)) %>% 
  mutate(run = gsub("b","",run)) %>% 
  mutate(run = gsub("c","",run)) %>% 
  mutate(Detected = ifelse(nReads>0, 1, 0)) %>% 
  filter(!(primer %in% c("MFU", "MV1") & totalReads == 0)) %>% 
  filter(Class == "Mammalia") %>% 
  filter(!BestTaxon %in% c("Moschus", "Equus caballus"))

## Filter out DLL1, C16 primer, plate 309, and DL/DLL1 from plate 314 ----------

detect_data_filt <- detect_data_raw %>% 
  filter(run != "MURI309") %>% 
  filter(!(primer %in% c("DLL1N", "C16", "DLL1"))) %>% 
  filter(!(primer == "DL" & run == "MURI314")) %>% 
  ungroup()

## Add freeze/thaw info --------------------------------------------------------

# # take sample info and change names to match other data streams
# MRS: the re-aliquoted plate (2A) is causing problems, so doing separately and then joining

# remove plate 2A
sampleinfo2 <- sample_locs %>%
  filter(Plate != "2.0A")

detect_data_plate2 <- detect_data_filt %>%
  filter(!run %in% c("MURI342", "MURI343", "MURI344")) %>% # plate 2.0A run on MURI342-MURI344
  left_join(sampleinfo2, by = "NWFSCsampleID") 

# now include plate 2A
sampleinfo2A <- sample_locs %>%
  filter(Plate != "2")

detect_data_plate2A <- detect_data_filt %>%
  filter(run %in% c("MURI342", "MURI343", "MURI344")) %>%
  left_join(sampleinfo2A, by = "NWFSCsampleID") 

# join
detect_data_plate <- rbind(detect_data_plate2, detect_data_plate2A) %>%
  mutate(Plate = as.character(Plate)) %>%
  mutate(Plate = case_when(Plate == "2" ~ "2.0", TRUE ~ Plate))

# join raw detection data with sample info

freezethaw_mod <- freezethaw %>%
  mutate(Plate = as.character(Plate)) %>%
  mutate(Plate = case_when(Plate == "2" ~ "2.0", TRUE ~ Plate)) %>%
  mutate("run" = paste0("MURI", RunNo)) %>%
  select(primer, Plate, Thaw, run)

detect_data_thaw <- detect_data_plate %>%
  left_join(freezethaw_mod, by = c("run", "Plate", "primer"))

rm(freezethaw_mod, detect_data_plate, detect_data_plate2, detect_data_plate2A,
   sampleinfo2, sampleinfo2A, freezethaw)

# MRS: I stopped here! And did multiple checks along the way and I think we are in good shape! 
  
## Reduce sequencing reps ------------------------------------------------------

detect_data_1seq <- detect_data_thaw %>% 
  group_by(run, primer, NWFSCsampleID, dilution, techRep, seqRep) %>% 
  mutate(totReads = sum(nReads)) %>% 
  ungroup() %>% 
  group_by(primer, NWFSCsampleID, dilution, techRep) %>% 
  filter(totReads == max(totReads)) %>% 
  ungroup() %>% 
  mutate(seqRep = replace_na(seqRep, "sr1")) %>% 
  filter(!(totReads == 0 & seqRep %in% c("sr2", "sr3"))) %>% 
  select(-totReads)
  

## Reduce dilutions ------------------------------------------------------------

detect_data_1dil <- detect_data_1seq %>% 
  group_by(run, primer, NWFSCsampleID, dilution, techRep) %>% 
  mutate(totReads = sum(nReads)) %>% 
  ungroup() %>% 
  group_by(primer, NWFSCsampleID, techRep) %>% 
  filter(totReads == max(totReads)) %>% 
  mutate(dilution = substr(dilution, 2, nchar(dilution))) %>% 
  mutate(diluti0n = as.numeric(dilution)) %>% 
  arrange(dilution, .by_group = TRUE) %>% 
  slice_head(n = 22) %>% 
  select(-totReads) 
          
## Reduce tech reps ------------------------------------------------------------
## EKJ Note I think we don't want to do this anymore, since
## each tech rep will have a different number of freeze-thaw cycles
## associated with it.
# detect_data_1rep <- detect_data_1dil %>% 
#   group_by(primer, NWFSCsampleID, BestTaxon) %>% 
#   mutate(nReps = n()) %>% 
#   mutate(Detected = max(Detected)) %>% 
#   slice_head()
  
## Count number of samples -----------------------------------------------------

nSamps_primer <- detect_data_1dil %>%
  group_by(primer, NWFSCsampleID) %>%
  n_groups()

length(unique(detect_data_1dil$NWFSCsampleID))

## Add metadata ----------------------------------------------------------------

detect_data_meta <- detect_data_1dil %>% 
  left_join(metadata, by = c("NWFSCsampleID" = "sampleID")) %>%
  left_join(mmEcoEvo, by = c("BestTaxon" = "Species"))

## Check species are all marine mammals
unique(detect_data_meta$BestTaxon)

## Remove Delphinidae family and all pinniped species, add common names
detect_data <- detect_data_meta %>% 
  filter(!(BestTaxon %in% c('Delphinidae', "Callorhinus ursinus",
                            "Eumetopias jubatus", "Phoca vitulina",
                            "Zalophus californianus", "Mirounga angustirostris")))  

  
## Check species are all cetaceans
unique(detect_data$BestTaxon)
unique(detect_data$common_name)

## Remove delphinid and baleen detections <100m from bottom (likely whalefall) -

detect_data_nowf <- detect_data %>%
  mutate(dist_to_bottom = bottom.depth.consensus - depth) %>%
  mutate(Detected = case_when(dist_to_bottom < 100 &
                  Detected == 1 &
                  bottom.depth.consensus > 200 &
                  Broad_taxa %in% c("Baleen whale", "Dolphin/Porpoise")~0,
                  TRUE~Detected))

## count number of detections by species ---------------------------------------

detect_per_species <- detect_data %>% 
  group_by(BestTaxon) %>% 
  summarize(nDetect = sum(Detected))

detect_per_family <- detect_data %>% 
  group_by(Family) %>% 
  summarize(nDetect = sum(Detected))

detect_per_primer_species <- detect_data %>% 
  group_by(BestTaxon, primer) %>% 
  summarize(nDetect = sum(Detected)) %>% 
  pivot_wider(names_from = "primer", values_from = "nDetect")

## add time at depth per species -----------------------------------------------

# max dive depth by species, per the Navy reports
maxDepth_species <- read.csv("./Data/MM_dive_time_expand.csv") %>% 
  group_by(Species) %>% 
  summarize(maxDepth = max(depth))

detect_species_divetime <- detect_data %>% 
  mutate(common_name = case_when(common_name == "killer whale"~"mammal eating killer whale",
                                 TRUE~common_name)) %>% 
  mutate(depth = case_when(depth == 0~1,
                           TRUE~depth)) %>% 
  left_join(timeAtDepth, by = c("common_name" = "Species", "depth" = "depth")) %>% 
  left_join(maxDepth_species, by = c("common_name" = "Species")) %>% 
  mutate(time_per_m = case_when(depth > maxDepth~0,
                                TRUE~time_per_m)) %>% 
  mutate(time_10m = case_when(depth > maxDepth~0,
                              TRUE~time_10m)) %>% 
  mutate(time_equalBin = case_when(depth > maxDepth~0,
                                   TRUE~time_equalBin))

save(detect_data, detect_species_divetime,
     detect_per_species, detect_per_family, 
     detect_per_primer_species, mmEcoEvo,
     maxDepth_species, file = "./ProcessedData/detect_data.Rdata")

## collapse data to detection/non-detection across all cetaceans
## (this is for the occupancy models where occupancy is 0/1
## for all cetaceans)

detect_data_allcet <- detect_data %>%
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

save(detect_data_allcet, file = "./ProcessedData/detect_data_allcet.RData")
