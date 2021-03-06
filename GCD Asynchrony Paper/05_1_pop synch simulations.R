### Simulating populations to assess mathematical dependencies of species loss on population asynchrony
###
### Author: Kevin Wilcox (kevin.wilcox@uwyo.edu)
### Created Oct 26, 2019; last updated Oct 26, 2019
rm(list=ls())
library(reshape2)
library(tidyverse)
library(ggthemes)
library(codyn)

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
num_of_bootstraps <- 1
num_of_years = 10
num_of_species = 10 ## Must be even number for now
num_of_plots = 10

spdiff_vector <- floor(1:(num_of_species/2))
synch_RR_master <- {}

pb = txtProgressBar(min = 0, max = num_of_bootstraps, initial = 0, style=3) 

for(SPDIFF in 0:length(spdiff_vector)){
  bootstrap_vector <- 1:num_of_bootstraps
  print(SPDIFF)
  
  for(BOOTSTRAP in 1:length(bootstrap_vector)){  
  setTxtProgressBar(pb,BOOTSTRAP)
  
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
    filter(sum_cover > 0) %>%
    dplyr::select(-sum_cover)
  
  spdiff_temp <- RAC_difference(df_temp_full, time.var = "year", species.var = "species", 
                                abundance.var = "cover", replicate.var = "plot", treatment.var = "treatment")
  
  spdiff_value <- max(spdiff_temp$species_diff)

  control_array_temp <- acast(filter(df_temp_full, treatment=="control"), year ~ species ~ plot, 
                              value.var="cover") %>%
    replace_na(0)
  treatment_array_temp <- acast(filter(df_temp_full, treatment=="treatment"), year ~ species ~ plot,
                                value.var="cover") %>%
    replace_na(0)

  ### Calculate stability metrics
  synch_vars_object_control_temp <- var.partition(control_array_temp)
  synch_vars_df_control_temp <- data.frame(t(as.data.frame(synch_vars_object_control_temp)))
  
  synchrony_vars_df_control_temp <- synch_vars_df_control_temp %>%
    mutate(spp_stab = 1/CV_S_L,
           alpha_stab = 1/CV_C_L,
           gamma_stab = 1/CV_C_R,
           pop_synch = phi_S_L2R,
           spatial_synch = phi_C_L2R,
           spp_synch = phi_S2C_L) %>%
    dplyr::select(spp_stab:spp_synch)
  
  
  synch_vars_object_treatment_temp <- var.partition(treatment_array_temp)
  synch_vars_df_treatment_temp <- data.frame(t(as.data.frame(synch_vars_object_treatment_temp)))
  
  synchrony_vars_df_treatment_temp <- synch_vars_df_treatment_temp %>%
    mutate(spp_stab = 1/CV_S_L,
           alpha_stab = 1/CV_C_L,
           gamma_stab = 1/CV_C_R,
           pop_synch = phi_S_L2R,
           spatial_synch = phi_C_L2R,
           spp_synch = phi_S2C_L) %>%
    dplyr::select(spp_stab:spp_synch)
  
  synch_RR_temp <- data.frame(iteration = BOOTSTRAP,
                              species_diff = spdiff_value,
                              log(synchrony_vars_df_treatment_temp / synchrony_vars_df_control_temp))
  synch_RR_master <- rbind(synch_RR_master, synch_RR_temp)
  
  rm(synch_RR_temp, synchrony_vars_df_treatment_temp, synchrony_vars_df_control_temp, spdiff_value,
     synch_vars_object_treatment_temp, synch_vars_object_control_temp, control_array_temp,
     treatment_array_temp,  df_temp_full,  df_temp_treatment,  df_temp_control, spdiff_temp,
     synch_vars_df_control_temp, synch_vars_df_treatment_temp)
  } ### End bootstrap loop ###

}
#################### End species diff loop (and full loop) ###############################

write.csv(synch_RR_master, file="species difference effects on sychrony metrics_simulation 1000 bootstraps_2019Nov04.csv", row.names=F)

###
### manipulate simulation data set
###
synch_RR_long <- synch_RR_master %>%
  gather(key=metric, value=response_ratio, -iteration:-species_diff)


###
### Plot simulation data to see if there are systematic effects of species differences on any of the stability/synchrony metrics
###
simulation_boxplot <- ggplot(synch_RR_long, aes(x=species_diff, y=response_ratio)) +
  geom_jitter(alpha = 0.5) +
  geom_boxplot(aes(group=factor(species_diff)), col="dodgerblue") +
  facet_wrap(~metric, scales="free") +
  theme_few()


pdf(file= "C:\\Users\\kwilcox4\\Dropbox\\Shared working groups\\C2E\\GCD asynchrony\\figures\\species difference bootstrap boxplot.pdf", width=7, height=5, useDingbats = F)
print(simulation_boxplot)
dev.off()




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
ggplot(synch_RR_master, aes(x=species_diff, y=pop_synch)) +
  geom_point() +
#  facet_wrap(~plot) +
  theme_few()

### Put dataframe into format for calculating stability metrics
no_spdiff_array <- acast(synch_RR_master, year ~ species ~ plot) %>%
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


### 
### Compare community with all the same species across replicates, with community using different species across replicates
###
num_of_bootstraps <- 1
num_of_years = 10
num_of_species = 10 ## Must be even number for now
num_of_plots = 5

spdiff_vector <- floor(1:(num_of_species/2))

df_temp_control <- data.frame(
  species = sort(rep(1:num_of_species, num_of_years * num_of_plots)),
  year = rep(sort(rep(1:num_of_years, num_of_plots)), num_of_species),
  plot = rep(1:num_of_plots, num_of_species * num_of_years),
  cover = c(rnorm(num_of_years * (num_of_species) * num_of_plots, mean = 0.50, sd = 0.15))
) 

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
### Cast into arrays
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

synchrony_vars_df_ten_spdiff <- synch_vars_df_ten_spdiff %>%
  mutate(spdiff = 10) %>%
  mutate(spp_stab = 1/CV_S_L,
         alpha_stab = 1/CV_C_L,
         gamma_stab = 1/CV_C_R,
         pop_synch = phi_S_L2R,
         spatial_synch = phi_C_L2R,
         spp_synch = phi_S2C_L) %>%
  dplyr::select(spdiff, spp_stab:spp_synch)

### Bind synch/stability metrics together
synchrony_vars_all <- synchrony_vars_df_no_spdiff %>%
  bind_rows(synchrony_vars_df_one_spdiff) %>%
  bind_rows(synchrony_vars_df_two_spdiff) %>%
  bind_rows(synchrony_vars_df_three_spdiff) %>%
  bind_rows(synchrony_vars_df_four_spdiff) %>%
  bind_rows(synchrony_vars_df_five_spdiff) %>%
  bind_rows(synchrony_vars_df_ten_spdiff) %>%
  gather(key=metric_name, value=metric_value, -spdiff)
  
ggplot(synchrony_vars_all, aes(x=spdiff, y=metric_value)) +
  geom_point() +
  geom_path() +
  theme_few() +
  facet_wrap(~metric_name, scales="free")





