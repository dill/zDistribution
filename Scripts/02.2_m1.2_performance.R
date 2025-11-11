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

ggplot(TSS_df, aes(x = threshold)) +
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
  
  TSS_dftemp <- as.data.frame(do.call(rbind.data.frame, TSS_list)) %>% 
    mutate(TSS = as.numeric(TSS)) %>% 
    mutate(bins = cut(threshold, breaks = seq(0,0.1, by = 0.005),
                      labels = FALSE)) %>% 
    group_by(bins) %>% 
    summarize(mean = mean(TSS), sd = sd(TSS), nSamp = n()) %>% 
    mutate(se = sd/(nSamp^(1/2))) %>% 
    mutate(low95 = mean - 1.96 * se,
           high95 = mean + 1.96 * se) %>% 
    mutate(species = species[a])
  
  TSS_by_species <- rbind(TSS_by_species, TSS_dftemp)
  
}

ggplot(TSS_by_species, aes(x = bins, color = species, fill = species)) +
  geom_smooth(aes(ymin = low95, ymax = high95, y = mean), 
              stat = "identity") +
  theme_minimal() +
  ylab("TSS")

meanAUC_by_species <- AUC_by_species %>% 
  group_by(species) %>% 
  summarize(meanAUC = mean(AUC))

### save test metrics ----------------------------------------------------------

save(AUC_df, TSS_df, 
     AUC_by_species,
     TSS_by_species, 
     meanAUC_by_species, file = "./ProcessedData/m1.2performance.Rdata")
