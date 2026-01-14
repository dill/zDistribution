### Calcuate POD for each sample in the study

library(tidyverse)
library(betareg)

load("Results/post_pdetect_m1.0+2LevelOcc_wThaw.RData")
load("ProcessedData/detect_data_clean.RData")

mean_post_pdetect <- post_pdetect %>% 
  group_by(Primer, Thaw) %>% 
  summarize(meanPOD = mean(PDetect)) %>% 
  mutate(Primer = case_when(Primer == "DL1"~"DL",
                            TRUE~Primer)) %>% 
  ungroup()

# EKJ checking for trends in detection probability 
POD_per_rep <- detect_data_clean %>%
  left_join(mean_post_pdetect, by = c("primer" = "Primer", "Thaw" = "Thaw")) %>%
  select(Sample_name, primer, Thaw, utm.lon, utm.lat, depth, meanPOD) %>%
  distinct()

hist(POD_per_rep$meanPOD)

# check depth
plot(POD_per_rep$depth, POD_per_rep$meanPOD)
m_depth <- betareg(meanPOD ~ depth, data = POD_per_rep)
summary(m_depth)
df <- data.frame(depth = 0:500)
df$pred <- predict(m_depth, newdata = df)
plot(df$depth, df$pred)

# check lat
plot(POD_per_rep$utm.lat, POD_per_rep$meanPOD)
m_lat <- betareg(meanPOD ~ utm.lat, data = POD_per_rep)
summary(m_lat)
df <- data.frame(utm.lat = 800:2000)
df$pred <- predict(m_lat, newdata = df)
plot(df$utm.lat, df$pred)

# check lon
plot(POD_per_rep$utm.lon, POD_per_rep$meanPOD)
m_lon <- betareg(meanPOD ~ utm.lon, data = POD_per_rep)
summary(m_lon)
df <- data.frame(utm.lon = 700:1000)
df$pred <- predict(m_lon, newdata = df)
plot(df$utm.lon, df$pred)






# check lat

# check lon

POD_per_samp <- detect_data_clean %>% 
  left_join(mean_post_pdetect, by = c("primer" = "Primer", "Thaw" = "Thaw")) %>% 
  distinct(Sample_name, .keep_all = TRUE) %>% 
  group_by(NWFSCsampleID) %>% 
  mutate(POnotD = 1-meanPOD) %>% 
    summarise(POD = 1-prod(POnotD), nReps = n())

POD_per_primer <- detect_data_clean %>% 
  left_join(mean_post_pdetect, by = c("primer" = "Primer", "Thaw" = "Thaw")) %>% 
  distinct(Sample_name, .keep_all = TRUE) %>% 
  group_by(primer) %>% 
  summarise(POD = mean(meanPOD), nSamps = n())

save(POD_per_primer, POD_per_samp, file = "./ProcessedData/POD_est_per_samp.Rdata")
