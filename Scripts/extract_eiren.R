###

library(dplyr)
library(tidyr)
library(stringr)

ff <- list.files(path = "./Data/ddf", pattern=".Rdata", full.names = TRUE)


ll <- lapply(ff, function(x){
  load(x)

  aa <- LTfitBest$data %>%
    filter(Yr==2018) %>%
    select(Long, Lat, starts_with("SP"), -SppMax) %>%
    mutate(Long = -Long) %>%
    pivot_longer(starts_with("SP"),
                 names_to="SpCode", values_to="somenumber") %>%
    filter(somenumber>0)

})

codes <- read.csv("./Data/MASTER_CCE_1991_2018_GS.csv") %>%
  mutate(SpCode = paste0("SP", str_pad(SpCode, 3, "left", "0"))) %>%
  select(Species, SpCode)

dat <- do.call(rbind, ll) %>%
  left_join(codes)


library(ggplot2)
library(sf)

dat_sf <- st_as_sf(dat, coords=c("Long", "Lat"), crs=4326) %>%
  filter(!is.na(Species))

ggplot() +
  geom_sf(data=dat_sf, aes(colour=Species)) +
  theme_minimal()

dat_sf$Type <- "Visual"

