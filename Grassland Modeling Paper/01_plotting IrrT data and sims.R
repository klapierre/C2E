### Compile and graph IrrT data and model output
###
### Author: Kevin Wilcox (kevin.wilcox@uwyo.edu)
### Created May 28, 2020; Last updated May 28, 2020

### Set up workspace
setwd("C:\\Users\\wilco\\Dropbox\\shared working groups\\C2E\\Products\\Modeling C2A\\")
library(tidyverse)
library(ggthemes)

###
### Read in and prep data
###

### Prepare model output for combining with annual empirical measurements

### teco ###
teco_output <- read.csv("Simulations\\TECO_output_long.csv") %>%
  rename(treatment=trt, year=yr) %>%
  mutate(anpp_simulated_gm2=Leaf.growth+Stem.growth)

teco_out_annual <- teco_output %>%
  group_by(year, treatment) %>%
  summarize(anpp_simulated_gm2 = sum(anpp_simulated_gm2)) %>%
  mutate(source="TECO")

### sdvgm ###

### BiomeE ###
biomee_out_annual <- read.csv("Simulations\\BiomeE-Konza-simulations-2018-11-27_long.csv") %>%
  dplyr::select(year, treatment, NPPlf1, NPPlf2, NPPwd1, NPPwd2) %>%
  mutate(anpp_simulated_gm2=(NPPlf1+NPPlf2+NPPwd1+NPPwd2)/0.45) %>%
  mutate(source="BiomeE")



### Prepare empirical data
anpp_empir <- read.csv("Data\\anpp_1991_2012.csv") %>%
  rename(anpp_observed_gm2=anpp_gm2) %>%
  mutate(source="Empirical")

### Combine empirical data and model output
teco_data_combined <- teco_out_annual %>%
  full_join(anpp_empir, by=c("year", "treatment"))

biomee_data_combined <- biomee_out_annual %>%
  full_join(anpp_empir, by=c("year", "treatment"))


###
### Plotting
###

### TECO ###
### Control plots versus model
ggplot(filter(teco_data_combined,treatment=="control"), aes(x=year, y=anpp_observed_gm2)) +
  geom_bar(stat="identity") +
  geom_point(data=filter(teco_data_combined,treatment=="control"), aes(x=year, y=anpp_simulated_gm2), inherit.aes=F, col="green", alpha=0.5, size=5) +
  theme_few()

### Irrigated plots vursus model
ggplot(filter(teco_data_combined,treatment=="irrigated"), aes(x=year, y=anpp_observed_gm2)) +
  geom_bar(stat="identity") +
  geom_point(data=filter(teco_data_combined,treatment=="irrigated"), aes(x=year, y=anpp_simulated_gm2), inherit.aes=F, col="blue", alpha=0.5, size=5)+
  theme_few()

### BiomeE ###
### Control plots
ggplot(filter(biomee_data_combined,treatment=="control"), aes(x=year, y=anpp_observed_gm2)) +
  geom_bar(stat="identity") +
  geom_point(data=filter(biomee_data_combined,treatment=="control"), aes(x=year, y=anpp_simulated_gm2), inherit.aes=F, col="green", alpha=0.5, size=5) +
  theme_few()

### Irrigated plots vursus model
ggplot(filter(biomee_data_combined,treatment=="irrigated"), aes(x=year, y=anpp_observed_gm2)) +
  geom_bar(stat="identity") +
  geom_point(data=filter(biomee_data_combined,treatment=="irrigated"), aes(x=year, y=anpp_simulated_gm2), inherit.aes=F, col="blue", alpha=0.5, size=5)+
  theme_few()


###
### Calcualte rmse and r2 for model-data comparisons
###

### TECO ###

### controls
teco_control_lm <- lm(anpp_simulated_gm2 ~ anpp_observed_gm2, data = filter(model_data_combined,treatment=="control"))
teco_control_r2 <- summary(teco_control_lm)$r.squared
teco_control_rmse <- modelr::rmse(model=teco_control_lm, data =filter(model_data_combined,treatment=="control"))

### irrigated
# all years
teco_irrigated_lm <- lm(anpp_simulated_gm2 ~ anpp_observed_gm2, data = filter(model_data_combined,treatment=="irrigated"))
teco_irrigated_r2 <- summary(teco_irrigated_lm)$r.squared
teco_irrigated_rmse <- modelr::rmse(model=teco_irrigated_lm, data =filter(model_data_combined,treatment=="irrigated"))

# pre-comm change (91-98)
teco_irrigated_prechange_lm <- lm(anpp_simulated_gm2 ~ anpp_observed_gm2, data = filter(model_data_combined,treatment=="irrigated"&year<1999))
teco_irrigated_prechange_r2 <- summary(teco_irrigated_prechange_lm)$r.squared
teco_irrigated_prechange_rmse <- modelr::rmse(model=teco_irrigated_prechange_lm, data =filter(model_data_combined,treatment=="irrigated"&year<1999))

# pre-comm change (99-12)
teco_irrigated_postchange_lm <- lm(anpp_simulated_gm2 ~ anpp_observed_gm2, data = filter(model_data_combined,treatment=="irrigated"&year>=1999))
teco_irrigated_postchange_r2 <- summary(teco_irrigated_postchange_lm)$r.squared
teco_irrigated_postchange_rmse <- modelr::rmse(model=teco_irrigated_postchange_lm, data =filter(model_data_combined,treatment=="irrigated"&year>=1999))













