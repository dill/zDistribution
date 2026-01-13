### Calcuate POD for each sample in the study

library(tidyverse)

load("Results/post_pdetect_m1.0+2LevelOcc_wThaw.RData")
load("ProcessedData/detect_data_clean.RData")

mean_post_pdetect <- post_pdetect %>% 
  group_by(Primer, Thaw) %>% 
  summarize(meanPOD = mean(PDetect)) %>% 
  mutate(Primer = case_when(Primer == "DL1"~"DL",
                            TRUE~Primer)) %>% 
  ungroup()

POD_per_samp <- detect_data_clean %>% 
  left_join(mean_post_pdetect, by = c("primer" = "Primer", "Thaw" = "Thaw")) %>% 
  distinct(Sample_name, .keep_all = TRUE) %>% 
  group_by(NWFSCsampleID) %>% 
  summarise(POD = sum(meanPOD), nReps = n())

POD_per_primer <- detect_data_clean %>% 
  left_join(mean_post_pdetect, by = c("primer" = "Primer", "Thaw" = "Thaw")) %>% 
  distinct(Sample_name, .keep_all = TRUE) %>% 
  group_by(primer) %>% 
  summarise(POD = mean(meanPOD), nSamps = n())

save(POD_per_primer, POD_per_samp, file = "./ProcessedData/POD_est_per_samp.Rdata")
