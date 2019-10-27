### Simulating populations to assess mathematical dependencies of species loss on population asynchrony
###
### Author: Kevin Wilcox (kevin.wilcox@uwyo.edu)
### Created Oct 26, 2019; last updated Oct 26, 2019

library(reshape2)
library(tidyverse)
library(ggthemes)

###
### Stability metric function
###

var.partition <- function(metacomm_tsdata){   # metacomm_tsdata = arrayx
  
  ## The function "var.partition" performs the partitioning of variability 
  ## across hierarchical levesl within a metacommunity.
  ## The input array "metacomm_tsdata" is an N*T*M array. The first dimension represents N species, 
  ## the second represents time-series observations of length T, and the third represents M local communities. 
  ## The output includes four variability and four synchrony metrics as defined in the main text.
  ## Note that, to be able to handle large metacommunities, this code has avoided calculating all covariance.
  
  ts_metacom <- apply(metacomm_tsdata,2,sum)
  ts_patch <- apply(metacomm_tsdata,c(2,3),sum)
  ts_species <- apply(metacomm_tsdata,c(1,2),sum)
  
  sd_metacom <- sd(ts_metacom)
  sd_patch_k <- apply(ts_patch,2,sd)
  sd_species_i <- apply(ts_species,1,sd)
  sd_species_patch_ik <- apply(metacomm_tsdata,c(1,3),sd)
  
  mean_metacom <- mean(ts_metacom)
  
  CV_S_L <- sum(sd_species_patch_ik)/mean_metacom  # spp_var
  CV_C_L <- sum(sd_patch_k)/mean_metacom  # alpha_var
  CV_S_R <- sum(sd_species_i)/mean_metacom
  CV_C_R <- sd_metacom/mean_metacom  # gamma_var
  
  phi_S_L2R <- CV_S_R/CV_S_L # pop synch
  phi_C_L2R <- CV_C_R/CV_C_L # spatial synch
  phi_S2C_L <- CV_C_L/CV_S_L # spp_synch
  phi_S2C_R <- CV_C_R/CV_S_R
  
  partition_3level <- c(CV_S_L=CV_S_L, CV_C_L=CV_C_L, CV_S_R=CV_S_R, CV_C_R=CV_C_R, 
                        phi_S_L2R=phi_S_L2R, phi_C_L2R=phi_C_L2R, phi_S2C_L=phi_S2C_L, phi_S2C_R=phi_S2C_R)
  return(partition_3level)
}

###
### Simulate community datasets
###

###
### Loop to create species datasets with different levels of species differences
###
num_of_bootstraps = 10
num_of_years = 10
num_of_species = 10 ## Must be even number for now
num_of_plots = 10

spdiff_vector <- floor(1:num_of_species/2)

for(SPDIFF in 0:length(spdiff_vector)){
  bootstrap_vector <- 1:num_of_bootstraps
  
  for(BOOTSTRAP in 1:length(bootstrap_vector)){  
  
  df_temp_control <- data.frame(
    treatment = "control",
    species = sort(rep(1:num_of_species, num_of_years * num_of_plots)),
    year = rep(sort(rep(1:num_of_years, num_of_plots)), num_of_species),
    plot = paste0("c",rep(1:num_of_plots, num_of_species * num_of_years)),
    cover = c(rnorm(num_of_years * (num_of_species/2) * num_of_plots, mean = 0.50, sd = 0.15),
              rep(0, num_of_years * (num_of_species/2) * num_of_plots))
  ) 
  
  if(SPDIFF == 0){
    df_temp_treatment <- df_temp_control %>%
      mutate(plot = paste0("t",rep(1:num_of_plots, num_of_species * num_of_years))) %>%
      mutate(treatment = "treatment")
  }  
  
  if(SPDIFF>0){
    df_temp_treatment <- df_temp_control %>%
      mutate(plot = paste0("t",rep(1:num_of_plots, num_of_species * num_of_years))) %>%
      mutate(treatment = "treatment") %>%
      mutate(cover = replace(cover, species %in% 1:SPDIFF, 0)) %>% 
      mutate(cover = replace(cover, species %in% ((num_of_species/2)+1):(SPDIFF+5), 
                             rnorm(num_of_years * SPDIFF * num_of_plots, mean = 0.50, sd = 0.15)))  
  }
  
  df_temp_full <- df_temp_control %>%
    bind_rows(df_temp_treatment) %>%
    group_by(year, species) %>%
    mutate(sum_cover = sum(cover)) %>%
    ungroup() %>%
    filter(sum_cover > 0)
  
  spdiff_temp <- RAC_difference(df_temp_full, time.var = "year", species.var = "species", 
                                abundance.var = "cover", replicate.var = "plot", treatment.var = "treatment")

  control_array_temp <- acast(filter(df_temp_full, treatment=="control"), year ~ species ~ plot) %>%
    replace_na(0)
  treatment_array_temp <- acast(filter(df_temp_full, treatment=="treatment"), year ~ species ~ plot) %>%
    replace_na(0)
  
  }

  filter(df_temp_full, treatment == "treatment" & plot == "t7" & year == 1)
  
}

ggplot(df_temp_treatment, aes(x=year, y=cover, col=as.factor(species))) +
  geom_line() +
  facet_wrap(~plot)

###
### No species differences
### 

### Create community
comm_data_no_spdiff <- data.frame(
  year = sort(rep(1:10,9)),
  plot = rep(sort(rep(1:3,3)),10),
  species = rep(c("a","b","c"), 30),
  cover = rnorm(90, mean = .50, sd = 0.15)
)

### Visualize species cover through time
ggplot(comm_data_no_spdiff, aes(x=year, y=cover,  col=species)) +
  geom_line() +
  facet_wrap(~plot) +
  theme_few()

### Put dataframe into format for calculating stability metrics
no_spdiff_array <- acast(comm_data_no_spdiff, year ~ species ~ plot) %>%
  replace_na(0)

### Calculate stability metrics
no_spdiff_vars_object <- var.partition(no_spdiff_array)
no_spdiff_vars_df <- data.frame(t(as.data.frame(no_spdiff_vars_object)))

synchrony_vars_df_no_spdiff <- no_spdiff_vars_df %>%
  mutate(method = "no spdiff",
         spp_stab = 1/CV_S_L,
         alpha_stab = 1/CV_C_L,
         gamma_stab = 1/CV_C_R,
         pop_synch = phi_S_L2R,
         spatial_synch = phi_C_L2R,
         spp_synch = phi_S2C_L)

###
### With species differences
###
comm_data_with_spdiff <- comm_data_no_spdiff %>%
  mutate(cover = replace(cover, plot == 1 & species == "a", 0)) %>%
  mutate(cover = replace(cover, plot == 2 & species == "b", 0))
  
### visualize
ggplot(comm_data_with_spdiff, aes(x=year, y=cover,  col=species)) +
  geom_line() +
  facet_wrap(~plot) +
  theme_few()

### Put dataframe into format for calculating stability metrics
###
with_spdiff_array <- acast(comm_data_with_spdiff, year ~ species ~ plot) %>%
  replace_na(0)

###
### Calculate stability metrics
###
with_spdiff_vars_object <- var.partition(with_spdiff_array)
with_spdiff_vars_df <- data.frame(t(as.data.frame(with_spdiff_vars_object)))

synchrony_vars_df_with_spdiff <- with_spdiff_vars_df %>%
  mutate(method = "with spdiff",
         spp_stab = 1/CV_S_L,
         alpha_stab = 1/CV_C_L,
         gamma_stab = 1/CV_C_R,
         pop_synch = phi_S_L2R,
         spatial_synch = phi_C_L2R,
         spp_synch = phi_S2C_L)

synchrony_both <- synchrony_vars_df_with_spdiff %>%
  bind_rows(synchrony_vars_df_no_spdiff)

ggplot(synchrony_both, aes(method, pop_synch)) +
  geom_point() +
  ylim(0,0.5)
