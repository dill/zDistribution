#### z Distribution
#### Figures for Q1 models
#### October 2025
#### AVC

library(tidyverse)
library(PNWColors)

load("./ProcessedData/H1models.RData")
load("./ProcessedData/detect_data.RData")
mmEcoEvo <- read.csv("./Data/MM_metadata.csv")
timeAtDepth <- read.csv("./Data/MM_dive_time_expand.csv")

zval <- 1.96
# 95% = 1.96
# 90% = 1.645
# 75% = 1.15

family_taxon <- mmEcoEvo %>% 
  distinct(Family, .keep_all = TRUE) %>% 
  select(Family, Broad_taxa)

### m1.0 figure ----------------------------------------------------------------

m1.0_predictions <- expand_grid(depth = 0:500)

m1.0preds <- predict(m1.0, m1.0_predictions, se.fit = TRUE)

m1.0_sePreds <- data.frame(m1.0_predictions,
                                mu   = binomial()$linkinv((m1.0preds$fit)),
                                low95  = binomial()$linkinv((m1.0preds$fit - zval * m1.0preds$se.fit)),
                                high95 = binomial()$linkinv(m1.0preds$fit + zval * m1.0preds$se.fit))

m1.0POD <- ggplot(m1.0_sePreds, aes(x = depth, 
                                    color = c(pnw_palette("Cascades",1, type = "continuous")), 
                                    fill = c(pnw_palette("Cascades",1, type = "continuous")))) +
  geom_smooth(aes(ymin = low95, ymax = high95, y = mu), stat = "identity") +
  ylab("POD") +
  geom_rug(data = detect_data, aes(x=depth), color = "grey30")+
  geom_rug(data = filter(detect_data, Detected == 1), aes(x=depth),
           sides = "t")+
  theme_minimal() +
  theme(legend.position = "none")


png(file = "./Figures/m1.0POD.png")
m1.0POD
dev.off()

### m1.1 figure ----------------------------------------------------------------

m1.1_predictions <- expand_grid(depth = 0:500, BestTaxon = as.factor(unique(detect_data$BestTaxon)))

m1.1preds <- predict(m1.1, m1.1_predictions, se.fit = TRUE)
m1.1_sePreds <- data.frame(m1.1_predictions,
                           mu   = binomial()$linkinv(m1.1preds$fit),
                           low  = binomial()$linkinv(m1.1preds$fit - zval * m1.1preds$se.fit),
                           high = binomial()$linkinv(m1.1preds$fit + zval * m1.1preds$se.fit)) %>% 
  left_join(mmEcoEvo, by = c("BestTaxon" = "Species"))

m1.1POD <- ggplot(m1.1_sePreds, aes(x = depth, color = Broad_taxa, fill = Broad_taxa)) +
  geom_smooth(aes(ymin = low, ymax = high, y = mu), stat = "identity") +
  ylab("POD") +
  scale_fill_manual(values = c(pnw_palette("Cascades",2, type = "continuous"),
                               pnw_palette("Sunset",2, type = "continuous"))) +
  scale_color_manual(values = c(pnw_palette("Cascades",2, type = "continuous"),
                                pnw_palette("Sunset",2, type = "continuous"))) +
  facet_wrap(~abbrev, scales = "free_y") +
  geom_rug(data = filter(detect_data, BestTaxon %in% (detect_per_species %>% 
                                                      filter(nDetect >= 10) %>% 
                                                        pull(BestTaxon))), 
           aes(x=depth), color = "grey30")+
  geom_rug(data = (detect_data %>% filter(Detected == 1) %>% 
                   filter(BestTaxon %in% (detect_per_species %>% 
                                          filter(nDetect >= 10) %>% 
                                          pull(BestTaxon))) %>%
                  mutate(depth_jittered = depth + runif(n(), -5, 30))),
           aes(x=depth),
           sides = "t")+
  theme_minimal() +
  theme(legend.position = "bottom")

png(file = "./Figures/m1.1POD.png")
m1.1POD
dev.off()

### m1.2 figure ----------------------------------------------------------------

m1.2_predictions <- expand_grid(depth = 0:500, BestTaxon = as.factor(unique(detect_data$BestTaxon)))


# link type prediction
m1.2predsL <- predict(m1.2, m1.2_predictions, se.fit = TRUE)
m1.2_sePreds_link <- data.frame(m1.2_predictions,
                           mu   = binomial()$linkinv((m1.2predsL$fit)),
                           low95  = binomial()$linkinv((m1.2predsL$fit - zval * m1.2predsL$se.fit)),
                           high95 = binomial()$linkinv(m1.2predsL$fit + zval * m1.2predsL$se.fit),
                           low75  = binomial()$linkinv((m1.2predsL$fit - 1.15 * m1.2predsL$se.fit)),
                           high75 = binomial()$linkinv(m1.2predsL$fit + 1.15 * m1.2predsL$se.fit),
                           low50  = binomial()$linkinv((m1.2predsL$fit - 0.674 * m1.2predsL$se.fit)),
                           high50 = binomial()$linkinv(m1.2predsL$fit + 0.674 * m1.2predsL$se.fit)) %>% 
  left_join(mmEcoEvo, by = c("BestTaxon" = "Species")) %>% 
  mutate(common_name = case_when(common_name == "killer whale"~"mammal eating killer whale",
                                 TRUE~common_name)) %>% 
  mutate(depth = case_when(depth == 0~1,
                           TRUE~depth)) %>%
  left_join(timeAtDepth, by = c("common_name" = "Species", "depth" = "depth")) %>% 
  left_join(maxDepth_species, by = c("common_name" = "Species")) %>% 
  mutate(time_10m = case_when(depth > maxDepth~0,
                              TRUE~time_10m)) %>% 
  filter(BestTaxon %in% (detect_per_species %>% 
                           filter(nDetect >= 10) %>% 
                           pull(BestTaxon)))

# transform time_10m to plot with POD by depth
m1.2_scaled <- m1.2_sePreds_link %>% 
  group_by(abbrev) %>%
  mutate(time_scaled = time_10m/545) %>% 
  mutate(
    a = diff(range(mu, na.rm = TRUE)) / diff(range(time_10m, na.rm = TRUE)),
    b = min(mu, na.rm = TRUE) - a * min(time_10m, na.rm = TRUE),
    time_scaled = a * time_10m + b)

# trim CI to ranges with detections per species
# det_ranges <- detect_data %>%
#   group_by(BestTaxon) %>%
#   filter(Detected == 1) %>%
#   summarise(min_depth = min(depth, na.rm = TRUE),
#             max_depth = max(depth, na.rm = TRUE))
# 
# m1.2_scaled_trimmed <- m1.2_scaled %>%
#   left_join(det_ranges, by = "BestTaxon") %>%
#   mutate(low  = ifelse(depth >= min_depth & depth <= max_depth, low, NA),
#          high = ifelse(depth >= min_depth & depth <= max_depth, high, NA))

# Plot
m1.2POD <- ggplot(m1.2_scaled, aes(x = depth, color = Broad_taxa, fill = Broad_taxa)) +
  geom_line(aes(y = mu)) +
  geom_smooth(aes(ymin = low95, ymax = high95, y = mu), stat = "identity") +
  geom_ribbon(aes(ymin = low75, ymax = high75), stat = "identity", alpha = 0.3, color = NA) +
  geom_ribbon(aes(ymin = low50, ymax = high50), stat = "identity", alpha = 0.3, color = NA) +
  geom_line(aes(y = time_scaled), linetype = 2) +
  scale_fill_manual(values = c(pnw_palette("Cascades",5, type = "continuous")[4:5],
                               pnw_palette("Sunset",1, type = "continuous"))) +
  scale_color_manual(values = c(pnw_palette("Cascades",5, type = "continuous")[4:5],
                                pnw_palette("Sunset",1, type = "continuous"))) +
  facet_wrap(~abbrev, scales = "free_y") +
  geom_rug(data = filter(detect_data, BestTaxon %in% (detect_per_species %>% 
                                                        filter(nDetect >= 10) %>% 
                                                        pull(BestTaxon))),
                         aes(x=depth), color = "grey30")+
  geom_rug(data = (detect_data %>% 
                     filter(Detected == 1) %>% 
                     filter(BestTaxon %in% (detect_per_species %>% 
                                              filter(nDetect >= 10) %>% 
                                              pull(BestTaxon))) %>%
                     mutate(depth_jittered = depth + runif(n(), -5, 30))), 
           aes(x=depth_jittered),
           sides = "t") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank()) +
  ylab("POD") +
  xlab("Depth")

png(file = "./Figures/m1.2POD.png")
m1.2POD
dev.off()

### m1.2a figure ---------------------------------------------------------------

m1.2a_predictions <- expand_grid(depth = 0:500, Family = as.factor(unique(detect_data$Family)))

m1.2apreds <- predict(m1.2a, m1.2a_predictions, se.fit = TRUE)

m1.2a_sePreds <- data.frame(m1.2a_predictions,
                            mu   = binomial()$linkinv(m1.2apreds$fit),
                            low  = binomial()$linkinv(m1.2apreds$fit - zval * m1.2apreds$se.fit),
                            high = binomial()$linkinv(m1.2apreds$fit + zval * m1.2apreds$se.fit)) %>% 
  left_join(family_taxon, by = c("Family" = "Family")) %>% 
  filter(Family %in% (detect_per_family %>% 
                           filter(nDetect >= 10) %>% 
                           pull(Family)))

m1.2aPOD <- ggplot(m1.2a_sePreds, aes(x = depth, color = Broad_taxa, fill = Broad_taxa)) +
  geom_smooth(aes(ymin = low, ymax = high, y = mu), stat = "identity") +
  scale_fill_manual(values = c(pnw_palette("Cascades",5, type = "continuous")[4:5],
                               pnw_palette("Sunset",1, type = "continuous"))) +
  scale_color_manual(values = c(pnw_palette("Cascades",5, type = "continuous")[4:5],
                                pnw_palette("Sunset",1, type = "continuous"))) +
  facet_wrap(~Family, scales = "free_y") +
  geom_rug(data = filter(detect_data, BestTaxon %in% (detect_per_species %>% 
                                                        filter(nDetect >= 10) %>% 
                                                        pull(BestTaxon))),
           aes(x=depth), color = "grey30") +
  geom_rug(data = (detect_data %>% 
                     filter(Detected == 1) %>%
                     filter(Family %in% (detect_per_family %>% 
                                        filter(nDetect >= 10) %>% 
                                        pull(Family))) %>%
                     mutate(depth_jittered = depth + runif(n(), -5, 30))),
           aes(x=depth_jittered),
           sides = "t") +
  theme_minimal() +
  ylab("POD") +
  xlab("Depth") +
  theme(legend.position = "bottom", legend.title = element_blank())

png(file = "./Figures/m1.2aPOD.png")
m1.2aPOD
dev.off()


### m1.2b figure ---------------------------------------------------------------

m1.2b_predictions <- expand_grid(depth = 0:500, Prey.family = as.factor(unique(detect_data$Prey.family)))

m1.2bpreds <- predict(m1.2b, m1.2b_predictions, se.fit = TRUE)

m1.2b_sePreds <- data.frame(m1.2b_predictions,
                            mu   = binomial()$linkinv(m1.2bpreds$fit),
                            low  = binomial()$linkinv(m1.2bpreds$fit - zval * m1.2bpreds$se.fit),
                            high = binomial()$linkinv(m1.2bpreds$fit + zval * m1.2bpreds$se.fit))

m1.2bPOD <- ggplot(m1.2b_sePreds, aes(x = depth, color = Prey.family, fill = Prey.family)) +
  geom_smooth(aes(ymin = low, ymax = high, y = mu), stat = "identity") +
  scale_fill_manual(values = c(pnw_palette("Mushroom",5, type = "continuous")[c(1,3,4)])) +
  scale_color_manual(values = c(pnw_palette("Mushroom",5, type = "continuous")[c(1,3,4)])) +
  #guides(color = "none", fill = "none") +
  facet_wrap(~Prey.family, scales = "free_y", ncol = 1,
             labeller = as_labeller(c(fish = "Bony Fishes", invert = "Invertebrates", squid = "Cephalopods"))) +
  geom_rug(data = detect_data, aes(x=depth), color = "grey30")+
  geom_rug(data = filter(detect_data, Detected == 1) %>%
             mutate(depth_jittered = depth + runif(n(), -5, 30)), aes(x=depth), sides = "t")+
  theme_minimal() +
  theme_minimal() +
  xlab("Depth")+
  ylab("POD") +
  theme(legend.position = "none", legend.title = element_blank())

png(file = "./Figures/m1.2bPOD.png")
m1.2bPOD
dev.off()


save(m1.0POD, m1.1POD, m1.2POD, m1.2aPOD, m1.2bPOD, 
     file = paste0("./Figures/H1PODplots.Rdata"))
