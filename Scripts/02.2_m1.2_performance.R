### Check m1.2 model performance with AUC
### November 2025

library(mgcv)
library(pROC)
library(tidysdm)

load("./ProcessedData/detect_data.RData")

#### group data ----------------------------------------------------------------

random_samples <- detect_data %>% 
  ungroup() %>% 
  select(NWFSCsampleID) %>% 
  distinct(NWFSCsampleID) %>% 
  mutate(groupNo = sample(1:10, size = nrow(.), replace = TRUE))
  

detect_data_AUC <- detect_data %>% 
  filter(BestTaxon %in% c(detect_per_species %>% 
           filter(nDetect > 10) %>% 
           pull(BestTaxon))) %>% 
  mutate(BestTaxon = as.factor(BestTaxon)) %>% 
  left_join(random_samples, by = "NWFSCsampleID")

#### AUC and TSS test ----------------------------------------------------------

AUC_list <- list()
TSS_list <- list()

for (i in 1:10){
  
  detect_data_train <- detect_data_AUC %>% 
    filter(!(groupNo == i))
  
  detect_data_test <- detect_data_AUC %>% 
    filter(groupNo == i)
  
  m1.2train <- gam(Detected ~ 
                ti(depth, k=5, bs="ts")+
                ti(BestTaxon, k=16, bs="re")+
                ti(depth, BestTaxon, k=c(5, 16), bs=c("ts","re")),
              family = "binomial", data = detect_data_train,
              method = "REML")
  
  m1.2preds <- predict(m1.2train, detect_data_test, type = "response")
  
  roc_object <- roc(detect_data_test$Detected, m1.2preds)
  
  AUC_list[[i]] <- auc(roc_object)
  
  threshold <- seq(0.005,max(m1.2preds),by=0.001)
  
  TSS_df <- list()
  
  for (j in 1:length(threshold)){
  TSS_df[[j]] <- detect_data_test %>% 
    as.data.frame() %>% 
    bind_cols(as.data.frame(m1.2preds)) %>% 
    mutate(Detected = as.factor(Detected)) %>% 
    mutate(m1.2preds = case_when(m1.2preds > threshold[j]~1,
                                 TRUE~0)) %>% 
    mutate(m1.2preds = as.factor(m1.2preds)) %>% 
    tss(., Detected, m1.2preds)
  
  }
  
  TSS_list[[i]] <- as.data.frame(unlist(TSS_df)) %>% 
    rename("TSS" = 1) %>% 
    filter(!(TSS %in% c("tss", "binary"))) %>% 
    mutate("threshold" = threshold) %>% 
    mutate("test" = i)
}

AUC_df <-do.call(rbind.data.frame, AUC_list) %>% 
  rename("AUC" = 1)

TSS_df <- as.data.frame(do.call(rbind.data.frame, TSS_list)) %>% 
  mutate(TSS = as.numeric(TSS)) %>% 
  group_by(threshold) %>% 
  summarize(mean = mean(TSS), sd = sd(TSS), nSamp = n()) %>% 
  mutate(se = sd/(nSamp^(1/2))) %>% 
  mutate(low95 = mean - 1.96 * se,
         high95 = mean + 1.96 * se)

TSS_plot <- ggplot(TSS_df, aes(x = threshold)) +
  geom_smooth(aes(ymin = low95, ymax = high95, y = mean), 
              stat = "identity") +
  theme_minimal() +
  ylab("TSS")

#### Repeat per species --------------------------------------------------------

species <- detect_per_species %>% 
  filter(nDetect > 20) %>% 
  pull(BestTaxon)

AUC_by_species <- data.frame()
TSS_by_species <- data.frame()
TSS_raw_species <- data.frame()

for (a in 1:length(species)) {
 
  AUC_list <- list()
  TSS_list <- list()
  
  for (i in 1:10){
    
    detect_data_train <- detect_data_AUC %>% 
      filter(BestTaxon == species[a]) %>% 
      filter(!(groupNo == i))
    
    detect_data_test <- detect_data_AUC %>% 
      filter(BestTaxon == species[a]) %>% 
      filter(groupNo == i)
    
    if(sum(detect_data_test$Detected) >0) {
    
    m1.2train <- gam(Detected ~ 
                       ti(depth, k=5, bs="ts"),
                       #ti(BestTaxon, k=16, bs="re")+
                       #ti(depth, BestTaxon, k=c(5, 16), bs=c("ts","re")),
                     family = "binomial", data = detect_data_train,
                     method = "REML")
    
    m1.2preds <- predict(m1.2train, detect_data_test, type = "response")
    
    roc_object <- roc(detect_data_test$Detected, m1.2preds)
    
    AUC_list[[i]] <- auc(roc_object)
    
    threshold <- seq(min(m1.2preds),max(m1.2preds),by=0.001)
    
    TSS_temp <- list()
    
    for (j in 1:length(threshold)){
      TSS_temp[[j]] <- detect_data_test %>% 
        as.data.frame() %>% 
        bind_cols(as.data.frame(m1.2preds)) %>% 
        mutate(Detected = as.factor(Detected)) %>% 
        mutate(m1.2preds = case_when(m1.2preds > threshold[j]~1,
                                     TRUE~0)) %>% 
        mutate(m1.2preds = as.factor(m1.2preds)) %>% 
        tss(., Detected, m1.2preds)
      
    }
    
    TSS_list[[i]] <- as.data.frame(unlist(TSS_temp)) %>% 
      rename("TSS" = 1) %>% 
      filter(!(TSS %in% c("tss", "binary"))) %>% 
      mutate("threshold" = threshold) %>% 
      mutate("test" = i)
    }
  
  AUC_dftemp <-do.call(rbind.data.frame, AUC_list) %>% 
    rename("AUC" = 1) %>% 
    mutate(species = species[a])
  
  }
  
  AUC_by_species <- rbind(AUC_by_species, AUC_dftemp)
  
  breaks <- seq(0,0.1, by = 0.005)
  names(breaks) <- seq(1:length(breaks))
  
  TSS_dftemp <- as.data.frame(do.call(rbind.data.frame, TSS_list)) %>% 
    mutate(TSS = as.numeric(TSS)) %>% 
    mutate(bins = cut(threshold, breaks = breaks,
                      labels = FALSE, include.lowest = TRUE,
                      ordered_result = TRUE)) %>% 
    group_by(bins) %>% 
    summarize(mean = mean(TSS), sd = sd(TSS), nSamp = n()) %>% 
    mutate(se = sd/(nSamp^(1/2))) %>% 
    mutate(low95 = mean - 1.96 * se,
           high95 = mean + 1.96 * se) %>% 
    mutate(species = species[a]) %>% 
    mutate(bin_label = breaks[names(breaks) == bins])
  
  TSS_by_species <- rbind(TSS_by_species, TSS_dftemp)
  TSS_raw_species <- rbind(TSS_raw_species, as.data.frame(do.call(rbind.data.frame, TSS_list)) %>%
                                            mutate(species = species[a]))
}


TSS_species_plot <- ggplot(TSS_by_species, 
                           aes(x = bin_label, color = species, fill = species)) +
 geom_smooth(aes(ymin = low95, ymax = high95, y = mean), 
              stat = "identity", alpha = 0.2) +
  theme_minimal() +
  ylab("TSS") +
  ylab("POD threshold") +
  scale_color_viridis_d() +
  scale_fill_viridis_d()

meanAUC_by_species <- AUC_by_species %>% 
  group_by(species) %>% 
  summarize(meanAUC = mean(AUC))

TSS_optim_species <- TSS_raw_species %>% 
  mutate(TSS = as.numeric(TSS)) %>% 
  group_by(species,test) %>% 
  arrange(desc(TSS)) %>% 
  slice_head() %>% ungroup() %>% 
  group_by(species) %>% 
  summarize(meanTSS = mean(TSS), meanThreshold = mean(threshold))

species_metrics <- meanAUC_by_species %>% 
  left_join(TSS_optim_species, by = "species")

### save test metrics ----------------------------------------------------------

save(AUC_df, TSS_df, 
     AUC_by_species,
     TSS_by_species, 
     meanAUC_by_species, 
     TSS_species_plot,
     TSS_plot,
     TSS_optim_species,
     species_metrics, file = "./ProcessedData/m1.2performance.Rdata")

### Repeat with "clean" dataset ------------------------------------------------

random_samples_clean <- detect_data_clean %>% 
  ungroup() %>% 
  select(NWFSCsampleID) %>% 
  distinct(NWFSCsampleID) %>% 
  mutate(groupNo = sample(1:10, size = nrow(.), replace = TRUE))


detect_data_AUC_clean <- detect_data_clean %>% 
  filter(BestTaxon %in% c(detect_per_species %>% 
                            filter(nDetect > 10) %>% 
                            pull(BestTaxon))) %>% 
  mutate(BestTaxon = as.factor(BestTaxon)) %>% 
  left_join(random_samples_clean, by = "NWFSCsampleID")

#### AUC and TSS test ----------------------------------------------------------

AUC_list_clean <- list()
TSS_list_clean <- list()

for (i in 1:10){
  
  detect_data_train <- detect_data_AUC_clean %>% 
    filter(!(groupNo == i))
  
  detect_data_test <- detect_data_AUC_clean %>% 
    filter(groupNo == i)
  
  m1.2train <- gam(Detected ~ 
                     ti(depth, k=5, bs="ts")+
                     ti(BestTaxon, k=16, bs="re")+
                     ti(depth, BestTaxon, k=c(5, 16), bs=c("ts","re")),
                   family = "binomial", data = detect_data_train,
                   method = "REML")
  
  m1.2preds <- predict(m1.2train, detect_data_test, type = "response")
  
  roc_object <- roc(detect_data_test$Detected, m1.2preds)
  
  AUC_list_clean[[i]] <- auc(roc_object)
  
  threshold <- seq(0.005,max(m1.2preds),by=0.001)
  
  TSS_df <- list()
  
  for (j in 1:length(threshold)){
    TSS_df[[j]] <- detect_data_test %>% 
      as.data.frame() %>% 
      bind_cols(as.data.frame(m1.2preds)) %>% 
      mutate(Detected = as.factor(Detected)) %>% 
      mutate(m1.2preds = case_when(m1.2preds > threshold[j]~1,
                                   TRUE~0)) %>% 
      mutate(m1.2preds = as.factor(m1.2preds)) %>% 
      tss(., Detected, m1.2preds)
    
  }
  
  TSS_list_clean[[i]] <- as.data.frame(unlist(TSS_df)) %>% 
    rename("TSS" = 1) %>% 
    filter(!(TSS %in% c("tss", "binary"))) %>% 
    mutate("threshold" = threshold) %>% 
    mutate("test" = i)
}

AUC_df_clean <- do.call(rbind.data.frame, AUC_list_clean) %>% 
  rename("AUC" = 1)

TSS_df_clean <- as.data.frame(do.call(rbind.data.frame, TSS_list_clean)) %>% 
  mutate(TSS = as.numeric(TSS)) %>% 
  group_by(threshold) %>% 
  summarize(mean = mean(TSS), sd = sd(TSS), nSamp = n()) %>% 
  mutate(se = sd/(nSamp^(1/2))) %>% 
  mutate(low95 = mean - 1.96 * se,
         high95 = mean + 1.96 * se)

TSS_plot_clean <- ggplot(TSS_df_clean, aes(x = threshold)) +
  geom_smooth(aes(ymin = low95, ymax = high95, y = mean), 
              stat = "identity") +
  theme_minimal() +
  ylab("TSS")

### save "clean" test metrics --------------------------------------------------

save(AUC_df_clean, TSS_df_clean,
     TSS_plot_clean, file = "./ProcessedData/m1.2clean_performance.Rdata")
