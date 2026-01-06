# explore effect of freeze-thaw on probability of detection

library(tidyverse)

### Get detections fo all cetaceans per sample ---------------------------------

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
  summarize(DetectAny = ifelse(sum(Detected) == 0, 0, 1)) %>% 
  #filter(!(primer == "MFU" & Thaw > 4)) %>% 
  mutate(primer = as.factor(primer)) %>% 
  ungroup()


## Get POD per primer/thaw -----------------------------------------------------
summarize_thaw_clean <- detect_data_allcet_clean %>%
  filter(!is.na(Thaw)) %>%
  group_by(primer, Thaw) %>%
  summarize(POD = sum(DetectAny)/length(DetectAny), n = length(DetectAny))

ggplot(summarize_thaw_clean) +
  geom_col(aes(x= Thaw, y = POD, fill = n)) +
  geom_text(aes(x=Thaw, y = 0.2, label = paste("n = ", n)))+
  facet_wrap(~primer, nrow = 3, scales = "free_y") +
  theme_bw()+
  scale_fill_gradient2(low = "lightgrey",
                       midpoint = 327,
                       mid = "darkgrey",
                       high = "black",
                       space="Lab")

ggsave(plot = last_plot(), file = "./Figures/detectThaw_clean.png", 
       width = 8, height = 6, units = "in")

#### Model effect of freeze/thaw on detection ----------------------------------

m <- glm(DetectAny ~ Thaw * primer,
           family = binomial, data = detect_data_allcet_clean)

summary(m)

# plot expectation

new.df <- expand.grid("primer" = unique(detect_data_allcet_clean$primer),
                      "Thaw" = 1:6)
new.df$pred <- predict(m, newdata = new.df, se.fit = TRUE)$fit
new.df$se <- predict(m, newdata = new.df, se.fit = TRUE)$se
new.df$lci <- inv.logit(new.df$pred + )
### TODO finish calculating confidence intervals, plot

ggplot(new.df) +
  geom_line(aes(x=Thaw, y = pred, color = primer)) +
  theme_bw()

# plot expectation on top of observed 

# Function: get model coefficients per primer
get_primer_coefs <- function(model) {
  coefs <- coef(model)
  
  intercepts <- c(MFU = coefs["(Intercept)"],
                  DL = coefs["(Intercept)"] + coefs["primerDL"],
                  MV = coefs["(Intercept)"] + coefs["primerMV1"])
  
  slopes <- c(MFU = coefs["Thaw"],
              DL = coefs["Thaw"] + coefs["Thaw:primerDL"],
              MV1 = coefs["Thaw"] + coefs["Thaw:primerMV1"])
  
  data.frame(Primer = names(intercepts),
             Intercept = intercepts,
             Slope = slopes)
}

# Bootstrap
set.seed(42)

n_boot <- 10000  

boot_results <- replicate(n_boot, {
  
  boot_data <- detect_data_allcet_clean[sample(nrow(detect_data_allcet_clean), replace = TRUE), ]
  
  boot_model <- glm(DetectAny ~ Thaw * primer,
                    family = binomial,
                    data = boot_data)
  
  get_primer_coefs(boot_model)
  
}, simplify = FALSE)

boot_df <- bind_rows(boot_results, .id = "Bootstrap") %>% 
  crossing(Thaw = 1:6) %>% 
  mutate(POD = plogis(Intercept + Slope * Thaw)) 

boot_df$Bootstrap <- as.integer(boot_df$Bootstrap)


ggplot(boot_df, aes(x = Thaw, y = POD, group = Bootstrap)) +
  geom_line(alpha = 0.1, linewidth = 0.3) +
  facet_wrap("Primer", scales = "free_y", nrow = 3) +
  theme_minimal()

#Get CI from bootstrap

ci_df_thaw <- boot_df %>%
  mutate(Primer = case_when(Primer == "DL.(Intercept)"~"DL",
                            Primer == "MFU.(Intercept)"~"MFU",
                            TRUE~"MV1")) %>% 
  group_by(Primer) %>%
  summarise(Intercept_lower = quantile(Intercept, 0.025, na.rm = TRUE),
            Intercept_upper = quantile(Intercept, 0.975, na.rm = TRUE),
            Slope_lower = quantile(Slope, 0.025, na.rm = TRUE),
            Slope_upper = quantile(Slope, 0.975, na.rm = TRUE),
            Intercept_mean = mean(Intercept, na.rm = TRUE),
            Slope_mean = mean(Slope, na.rm = TRUE))


## Plot

newdat_thaw <- expand.grid(Thaw = seq(1,7,1),
                      primer = unique(detect_data_allcet_clean$primer)) %>% 
  left_join(ci_df_thaw %>% select(Primer, Intercept_mean, Slope_mean,
                             Intercept_lower, Intercept_upper,
                             Slope_lower, Slope_upper),
            by = c("primer" = "Primer")) %>%
  mutate(fit = plogis(Intercept_mean + Slope_mean * Thaw),
         lower = plogis(Intercept_lower + Slope_lower * Thaw),
         upper = plogis(Intercept_upper + Slope_upper * Thaw))

PODthaw_plot <- ggplot(newdat_thaw, aes(x = Thaw, y = fit, color = primer, fill = primer)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.5, color = NA) +
  geom_line(linewidth = 1.2) +
  labs(x = "Thaw cycles", y = "Detection probability") +
  theme_minimal() +
  facet_wrap(~primer, scales = "free_y", nrow = 3) +
  ylim(0,1) +
  scale_color_manual(values = c(pnw_palette("Cascades",5, type = "discrete")[c(3)],
                                pnw_palette("Sunset",7, type = "discrete")[c(3,5)])) +
  scale_fill_manual(values = c(pnw_palette("Cascades",5, type = "discrete")[c(3)],
                                pnw_palette("Sunset",7, type = "discrete")[c(3,5)]))

save(ci_df_thaw, PODthaw_plot, file = "ProcessedData/PODxThawxPrimer.Rdata")

## Estimate POD average across grid -------------------------------------------

POD_samples_exp <- summarize_thaw_clean %>% 
  rename("PODobs" = "POD") %>% 
  left_join(newdat_thaw %>% 
              dplyr::select(fit, lower, upper, Thaw, primer), 
            by = c("Thaw", "primer")) %>%
  group_by(primer) %>% 
  mutate(fit_n = n*fit, lower_n = n*lower, upper_n = n*upper) %>% 
  summarize(meanPOD = sum(fit_n)/516, maxPOD = sum(upper_n)/516, minPOD = sum(lower_n)/516) %>% 
  ungroup() %>% 
  summarise(sum(meanPOD), sum(maxPOD), sum(minPOD))


