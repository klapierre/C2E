### Simulating species differences to assess  mathematical dependencies of species loss on stability/synchrony metrics
###
### Author: Kevin Wilcox (kevin.wilcox@uwyo.edu)
### Created Apr 07, 2020; last updated Apr 07, 2020

rm(list=ls())
library(reshape2)
library(tidyverse)
library(ggthemes)
library(codyn)

###
### Stability/synchrony metric function from Shaopeng
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
### Simulate communities having different levels of species differences (same species across all replicates versus species in some but not all replicates) 
###

num_of_bootstraps <- 1
num_of_years = 10
num_of_species = 10 ## Must be even number for now
num_of_plots = 5

### Generate random covers to be used in all species diff simulations
df_temp_control <- data.frame(
  species = sort(rep(1:num_of_species, num_of_years * num_of_plots)),
  year = rep(sort(rep(1:num_of_years, num_of_plots)), num_of_species),
  plot = rep(1:num_of_plots, num_of_species * num_of_years),
  cover = c(rnorm(num_of_years * (num_of_species) * num_of_plots, mean = 0.50, sd = 0.15))
) 

### Combine cover values generated above in different ways to simulate different species difference
df_no_spdiff <- df_temp_control %>%
  bind_rows(
    df_temp_control %>%
      mutate(plot=plot+num_of_plots)
    )

df_one_spdiff <- df_temp_control %>%
  bind_rows(
    df_temp_control %>%
      mutate(plot=plot+num_of_plots) %>%
      mutate(species=replace(species, species == 1, 11))      
    )
  

df_two_spdiff <- df_temp_control %>%
  bind_rows(
    df_temp_control %>%
      mutate(plot=plot+num_of_plots) %>%
      mutate(species=replace(species, species == 1, 11),
             species=replace(species, species == 2, 12)
      )
  )

df_three_spdiff <- df_temp_control %>%
  bind_rows(
    df_temp_control %>%
      mutate(plot=plot+num_of_plots) %>%
      mutate(species=replace(species, species == 1, 11),
             species=replace(species, species == 2, 12),
             species=replace(species, species == 3, 13)
      )
  )

df_four_spdiff <- df_temp_control %>%
  bind_rows(
    df_temp_control %>%
      mutate(plot=plot+num_of_plots) %>%
      mutate(species=replace(species, species == 1, 11),
             species=replace(species, species == 2, 12),
             species=replace(species, species == 3, 13),
             species=replace(species, species == 4, 14)
      )
  )

df_five_spdiff <- df_temp_control %>%
  bind_rows(
    df_temp_control %>%
      mutate(plot=plot+num_of_plots) %>%
      mutate(species=replace(species, species == 1, 11),
             species=replace(species, species == 2, 12),
             species=replace(species, species == 3, 13),
             species=replace(species, species == 4, 14),
             species=replace(species, species == 5, 15)
      )
  )

df_six_spdiff <- df_temp_control %>%
  bind_rows(
    df_temp_control %>%
      mutate(plot=plot+num_of_plots) %>%
      mutate(species=replace(species, species == 1, 11),
             species=replace(species, species == 2, 12),
             species=replace(species, species == 3, 13),
             species=replace(species, species == 4, 14),
             species=replace(species, species == 5, 15),
             species=replace(species, species == 6, 16)
      )
  )

df_seven_spdiff <- df_temp_control %>%
  bind_rows(
    df_temp_control %>%
      mutate(plot=plot+num_of_plots) %>%
      mutate(species=replace(species, species == 1, 11),
             species=replace(species, species == 2, 12),
             species=replace(species, species == 3, 13),
             species=replace(species, species == 4, 14),
             species=replace(species, species == 5, 15),
             species=replace(species, species == 6, 16),
             species=replace(species, species == 7, 17)
      )
  )

df_eight_spdiff <- df_temp_control %>%
  bind_rows(
    df_temp_control %>%
      mutate(plot=plot+num_of_plots) %>%
      mutate(species=replace(species, species == 1, 11),
             species=replace(species, species == 2, 12),
             species=replace(species, species == 3, 13),
             species=replace(species, species == 4, 14),
             species=replace(species, species == 5, 15),
             species=replace(species, species == 6, 16),
             species=replace(species, species == 7, 17),
             species=replace(species, species == 8, 18)
      )
  )

df_nine_spdiff <- df_temp_control %>%
  bind_rows(
    df_temp_control %>%
      mutate(plot=plot+num_of_plots) %>%
      mutate(species=replace(species, species == 1, 11),
             species=replace(species, species == 2, 12),
             species=replace(species, species == 3, 13),
             species=replace(species, species == 4, 14),
             species=replace(species, species == 5, 15),
             species=replace(species, species == 6, 16),
             species=replace(species, species == 7, 17),
             species=replace(species, species == 8, 18),
             species=replace(species, species == 9, 19)
      )
  )


df_ten_spdiff <- df_temp_control %>%
  bind_rows(
    df_temp_control %>%
      mutate(plot=plot+num_of_plots) %>%
      mutate(species=replace(species, species == 1, 11),
             species=replace(species, species == 2, 12),
             species=replace(species, species == 3, 13),
             species=replace(species, species == 4, 14),
             species=replace(species, species == 5, 15),
             species=replace(species, species == 6, 16),
             species=replace(species, species == 7, 17),
             species=replace(species, species == 8, 18),
             species=replace(species, species == 9, 19),
             species=replace(species, species == 10, 20)
             
      )
  )

###
### Cast into arrays for var.partition function
###

no_spdiff_array_temp <- acast(df_no_spdiff, year ~ species ~ plot, 
                            value.var="cover") %>%
  replace_na(0)

one_spdiff_array_temp <- acast(df_one_spdiff, year ~ species ~ plot,
                                value.var="cover") %>%
  replace_na(0)

two_spdiff_array_temp <- acast(df_two_spdiff, year ~ species ~ plot,
                              value.var="cover") %>%
  replace_na(0)

three_spdiff_array_temp <- acast(df_three_spdiff, year ~ species ~ plot,
                                value.var="cover") %>%
  replace_na(0)

four_spdiff_array_temp <- acast(df_four_spdiff, year ~ species ~ plot,
                                value.var="cover") %>%
  replace_na(0)

five_spdiff_array_temp <- acast(df_five_spdiff, year ~ species ~ plot,
                                value.var="cover") %>%
  replace_na(0)

six_spdiff_array_temp <- acast(df_six_spdiff, year ~ species ~ plot,
                                value.var="cover") %>%
  replace_na(0)

seven_spdiff_array_temp <- acast(df_seven_spdiff, year ~ species ~ plot,
                                value.var="cover") %>%
  replace_na(0)

eight_spdiff_array_temp <- acast(df_eight_spdiff, year ~ species ~ plot,
                                value.var="cover") %>%
  replace_na(0)

nine_spdiff_array_temp <- acast(df_nine_spdiff, year ~ species ~ plot,
                                value.var="cover") %>%
  replace_na(0)

ten_spdiff_array_temp <- acast(df_ten_spdiff, year ~ species ~ plot,
                                value.var="cover") %>%
  replace_na(0)

###
### Calculate stability metrics from arrays
###

### Create synch variable objects
synch_vars_object_no_spdiff <- var.partition(no_spdiff_array_temp)
synch_vars_df_no_spdiff <- data.frame(t(as.data.frame(synch_vars_object_no_spdiff)))

synch_vars_object_one_spdiff <- var.partition(one_spdiff_array_temp)
synch_vars_df_one_spdiff <- data.frame(t(as.data.frame(synch_vars_object_one_spdiff)))

synch_vars_object_two_spdiff <- var.partition(two_spdiff_array_temp)
synch_vars_df_two_spdiff <- data.frame(t(as.data.frame(synch_vars_object_two_spdiff)))

synch_vars_object_three_spdiff <- var.partition(three_spdiff_array_temp)
synch_vars_df_three_spdiff <- data.frame(t(as.data.frame(synch_vars_object_three_spdiff)))

synch_vars_object_four_spdiff <- var.partition(four_spdiff_array_temp)
synch_vars_df_four_spdiff <- data.frame(t(as.data.frame(synch_vars_object_four_spdiff)))

synch_vars_object_five_spdiff <- var.partition(five_spdiff_array_temp)
synch_vars_df_five_spdiff <- data.frame(t(as.data.frame(synch_vars_object_five_spdiff)))

synch_vars_object_six_spdiff <- var.partition(six_spdiff_array_temp)
synch_vars_df_six_spdiff <- data.frame(t(as.data.frame(synch_vars_object_six_spdiff)))

synch_vars_object_seven_spdiff <- var.partition(seven_spdiff_array_temp)
synch_vars_df_seven_spdiff <- data.frame(t(as.data.frame(synch_vars_object_seven_spdiff)))

synch_vars_object_eight_spdiff <- var.partition(eight_spdiff_array_temp)
synch_vars_df_eight_spdiff <- data.frame(t(as.data.frame(synch_vars_object_eight_spdiff)))

synch_vars_object_nine_spdiff <- var.partition(nine_spdiff_array_temp)
synch_vars_df_nine_spdiff <- data.frame(t(as.data.frame(synch_vars_object_nine_spdiff)))

synch_vars_object_ten_spdiff <- var.partition(ten_spdiff_array_temp)
synch_vars_df_ten_spdiff <- data.frame(t(as.data.frame(synch_vars_object_ten_spdiff)))

### convert to stabilty and synchrony metrics
synchrony_vars_df_no_spdiff <- synch_vars_df_no_spdiff %>%
  mutate(spdiff = 0) %>%
  mutate(spp_stab = 1/CV_S_L,
         alpha_stab = 1/CV_C_L,
         gamma_stab = 1/CV_C_R,
         pop_synch = phi_S_L2R,
         spatial_synch = phi_C_L2R,
         spp_synch = phi_S2C_L) %>%
  dplyr::select(spdiff, spp_stab:spp_synch)

synchrony_vars_df_one_spdiff <- synch_vars_df_one_spdiff %>%
  mutate(spdiff = 1) %>%
  mutate(spp_stab = 1/CV_S_L,
         alpha_stab = 1/CV_C_L,
         gamma_stab = 1/CV_C_R,
         pop_synch = phi_S_L2R,
         spatial_synch = phi_C_L2R,
         spp_synch = phi_S2C_L) %>%
  dplyr::select(spdiff, spp_stab:spp_synch)

synchrony_vars_df_two_spdiff <- synch_vars_df_two_spdiff %>%
  mutate(spdiff = 2) %>%
  mutate(spp_stab = 1/CV_S_L,
         alpha_stab = 1/CV_C_L,
         gamma_stab = 1/CV_C_R,
         pop_synch = phi_S_L2R,
         spatial_synch = phi_C_L2R,
         spp_synch = phi_S2C_L) %>%
  dplyr::select(spdiff, spp_stab:spp_synch)

synchrony_vars_df_three_spdiff <- synch_vars_df_three_spdiff %>%
  mutate(spdiff = 3) %>%
  mutate(spp_stab = 1/CV_S_L,
         alpha_stab = 1/CV_C_L,
         gamma_stab = 1/CV_C_R,
         pop_synch = phi_S_L2R,
         spatial_synch = phi_C_L2R,
         spp_synch = phi_S2C_L) %>%
  dplyr::select(spdiff, spp_stab:spp_synch)

synchrony_vars_df_four_spdiff <- synch_vars_df_four_spdiff %>%
  mutate(spdiff = 4) %>%
  mutate(spp_stab = 1/CV_S_L,
         alpha_stab = 1/CV_C_L,
         gamma_stab = 1/CV_C_R,
         pop_synch = phi_S_L2R,
         spatial_synch = phi_C_L2R,
         spp_synch = phi_S2C_L) %>%
  dplyr::select(spdiff, spp_stab:spp_synch)

synchrony_vars_df_five_spdiff <- synch_vars_df_five_spdiff %>%
  mutate(spdiff = 5) %>%
  mutate(spp_stab = 1/CV_S_L,
         alpha_stab = 1/CV_C_L,
         gamma_stab = 1/CV_C_R,
         pop_synch = phi_S_L2R,
         spatial_synch = phi_C_L2R,
         spp_synch = phi_S2C_L) %>%
  dplyr::select(spdiff, spp_stab:spp_synch)

synchrony_vars_df_six_spdiff <- synch_vars_df_six_spdiff %>%
  mutate(spdiff = 6) %>%
  mutate(spp_stab = 1/CV_S_L,
         alpha_stab = 1/CV_C_L,
         gamma_stab = 1/CV_C_R,
         pop_synch = phi_S_L2R,
         spatial_synch = phi_C_L2R,
         spp_synch = phi_S2C_L) %>%
  dplyr::select(spdiff, spp_stab:spp_synch)

synchrony_vars_df_seven_spdiff <- synch_vars_df_seven_spdiff %>%
  mutate(spdiff = 7) %>%
  mutate(spp_stab = 1/CV_S_L,
         alpha_stab = 1/CV_C_L,
         gamma_stab = 1/CV_C_R,
         pop_synch = phi_S_L2R,
         spatial_synch = phi_C_L2R,
         spp_synch = phi_S2C_L) %>%
  dplyr::select(spdiff, spp_stab:spp_synch)

synchrony_vars_df_eight_spdiff <- synch_vars_df_eight_spdiff %>%
  mutate(spdiff = 8) %>%
  mutate(spp_stab = 1/CV_S_L,
         alpha_stab = 1/CV_C_L,
         gamma_stab = 1/CV_C_R,
         pop_synch = phi_S_L2R,
         spatial_synch = phi_C_L2R,
         spp_synch = phi_S2C_L) %>%
  dplyr::select(spdiff, spp_stab:spp_synch)

synchrony_vars_df_nine_spdiff <- synch_vars_df_nine_spdiff %>%
  mutate(spdiff = 9) %>%
  mutate(spp_stab = 1/CV_S_L,
         alpha_stab = 1/CV_C_L,
         gamma_stab = 1/CV_C_R,
         pop_synch = phi_S_L2R,
         spatial_synch = phi_C_L2R,
         spp_synch = phi_S2C_L) %>%
  dplyr::select(spdiff, spp_stab:spp_synch)

synchrony_vars_df_ten_spdiff <- synch_vars_df_ten_spdiff %>%
  mutate(spdiff = 10) %>%
  mutate(spp_stab = 1/CV_S_L,
         alpha_stab = 1/CV_C_L,
         gamma_stab = 1/CV_C_R,
         pop_synch = phi_S_L2R,
         spatial_synch = phi_C_L2R,
         spp_synch = phi_S2C_L) %>%
  dplyr::select(spdiff, spp_stab:spp_synch)

###
### Bind synch/stability metrics together
###
synchrony_vars_all <- synchrony_vars_df_no_spdiff %>%
  bind_rows(synchrony_vars_df_one_spdiff) %>%
  bind_rows(synchrony_vars_df_two_spdiff) %>%
  bind_rows(synchrony_vars_df_three_spdiff) %>%
  bind_rows(synchrony_vars_df_four_spdiff) %>%
  bind_rows(synchrony_vars_df_five_spdiff) %>%
  bind_rows(synchrony_vars_df_six_spdiff) %>%
  bind_rows(synchrony_vars_df_seven_spdiff) %>%
  bind_rows(synchrony_vars_df_eight_spdiff) %>%
  bind_rows(synchrony_vars_df_nine_spdiff) %>%
  bind_rows(synchrony_vars_df_ten_spdiff) %>%
  gather(key=metric_name, value=metric_value, -spdiff)
  
###
### Visualize
###
ggplot(synchrony_vars_all, aes(x=spdiff, y=metric_value)) +
  geom_point() +
  geom_path() +
  theme_few() +
  facet_wrap(~metric_name, scales="free")





