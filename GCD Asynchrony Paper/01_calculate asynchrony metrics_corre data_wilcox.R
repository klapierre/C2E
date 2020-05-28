##### Load packages
library(reshape2)
library(tidyverse) # data munging
library(codyn)     # community dynamics metrics
library(vegan)     # diversity metrics
#library(piecewiseSEM)
library(nlme)

#############################################################################
#############################################################################
######################### 1. Load and clean data ############################
#############################################################################
#############################################################################

# Load data
setwd("C:\\Users\\wilco\\Dropbox\\shared working groups\\C2E\\GCD asynchrony\\")
cover_df <- read.csv('data\\SpeciesRawAbundance_March2019.csv', header = T)
exp_info <- read.csv('data\\ExperimentInformation_March2019.csv', header = T) %>%
  mutate(site_proj_comm_trt = paste(site_code, project_name, community_type, treatment, sep='_'))
site_info <- read.csv('data\\SiteExperimentDetails_March2019.csv')
exp_info_short <- unique(dplyr::select(exp_info, -X, -calendar_year, -treatment_year))

full_df <- cover_df %>%
  full_join(exp_info, by=c("site_code", 
                           "project_name", 
                           "calendar_year",
                           "treatment_year",
                           "treatment",
                           "community_type"))

################################################################################################################
################################################################################################################
### 2. Shaopeng - Variance partitioning method
################################################################################################################
################################################################################################################

###### here is the function for calculating all metrics defined in Wang et al. (2019 Ecography)
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
### Loop through experiments to calculate synchrony/stability metrics
###
  
site_vector <- sort(unique(exp_info$site_proj_comm_trt))
master_partition_vars_df <- {}


for(SITE in 1:length(site_vector)){
  temp_df <- full_df %>%
    filter(site_proj_comm_trt == site_vector[SITE]) %>%
    dplyr::select(calendar_year, plot_id, genus_species, abundance)
  temp_exp_info <- exp_info_short %>%
    filter(site_proj_comm_trt == site_vector[SITE]) %>%
    slice(1) ### THIS TAKES THE FIRST OBSERVATION BECAUSE OF CAR AND BOWMAN HAVING DIFFERENT NUTRIENT AMOUNTS ACROSS YEARS... CANNOT TRUST N P AND K AMOUNTS FOR CAR AND BOWMAN IN THE RESULTING DATAFRAME
  
  temp_array <- acast(temp_df, genus_species ~ calendar_year ~ plot_id) %>%
    replace_na(0)
  
  temp_vars_object <- var.partition(temp_array)
  temp_partition_vars_df <- data.frame(temp_exp_info, t(as.data.frame(temp_vars_object)))
  master_partition_vars_df <- rbind(master_partition_vars_df, temp_partition_vars_df)
  rm(temp_df, temp_exp_info, temp_array, temp_vars_object, temp_partition_vars_df)
}  

synchrony_vars_df <- master_partition_vars_df %>%
  mutate(spp_stab = 1/CV_S_L,
         alpha_stab = 1/CV_C_L,
         gamma_stab = 1/CV_C_R,
         pop_synch = phi_S_L2R,
         spatial_synch = phi_C_L2R,
         spp_synch = phi_S2C_L)

write.csv(synchrony_vars_df, file="partition variables from Shaopeng's function_22Apr2020.csv", row.names=F)
synchrony_vars_df <- read.csv("partition variables from Shaopeng's function_22Apr2020.csv") %>%
  mutate(pop_synch       = 1/pop_asynch,
         spatial_synch   = 1/spatial_asynch,
         spp_synch       = 1/spp_asynch)

###
### Take a quick look
###
synchrony_vars_long <- synchrony_vars_df %>%
  dplyr::select(-c(CV_S_L:phi_S2C_R)) %>%
  gather(key=metric_name, value=metric_value, -(site_code:site_proj_comm_trt))

ggplot(filter(synchrony_vars_long, plot_mani<=1 & (plot_mani==0|nutrients==1)), aes(x=factor(nutrients), y=metric_value)) +
  geom_boxplot() +
  geom_jitter() +
  facet_wrap(~metric_name, scales="free") +
  theme_classic()

ggplot(filter(synchrony_vars_long, plot_mani<=1 & (plot_mani==0|light==1)), aes(x=factor(light), y=metric_value)) +
  geom_boxplot() +
  geom_jitter() +
  facet_wrap(~metric_name, scales="free") +
  theme_classic()

ggplot(filter(synchrony_vars_long, plot_mani<=1 & (plot_mani==0|carbon==1)), aes(x=factor(carbon), y=metric_value)) +
  geom_boxplot() +
  geom_jitter() +
  facet_wrap(~metric_name, scales="free") +
  theme_classic()

ggplot(filter(synchrony_vars_long, plot_mani<=1 & (plot_mani==0|water==1)), aes(x=factor(water), y=metric_value)) +
  geom_boxplot() +
  geom_jitter() +
  facet_wrap(~metric_name, scales="free") +
  theme_classic()

ggplot(filter(synchrony_vars_long, plot_mani<=1 & (plot_mani==0|carbon==1)), aes(x=factor(carbon), y=metric_value)) +
  geom_boxplot() +
  geom_jitter() +
  facet_wrap(~metric_name, scales="free") +
  theme_classic()


###
### Calculate response ratios
###
synchrony_controls <- synchrony_vars_long %>%
  filter(plot_mani == 0) %>%
  rename(metric_value_control = metric_value) %>%
  dplyr::select(site_code, project_name, community_type, metric_name, metric_value_control)
synchrony_rr_long <- synchrony_vars_long %>%
  filter(plot_mani > 0) %>%
  rename(metric_value_trt = metric_value) %>%
  full_join(synchrony_controls, by=c("site_code", "project_name", "community_type", "metric_name")) %>%
  mutate(lnRR = log(metric_value_trt/metric_value_control))
  
write.csv(synchrony_rr_long, file="Synchrony metrics response ratios_long form_22Apr2020.csv", row.names=F)

###
### Visualize response ratios
###
table(filter(synchrony_rr_long, plot_mani<=1 & metric_name=="alpha_stab")$trt_type)

ggplot(filter(synchrony_rr_long, plot_mani==1 & trt_type %in% c("CO2", "drought", "irr", "N", "P", "temp")), 
              aes(x=factor(trt_type), y=lnRR, col=factor(trt_type))) +
  geom_boxplot() +
  geom_jitter(shape=21) +
  facet_wrap(~metric_name, scales="free") +
  theme_bw() +
  geom_hline(yintercept=0)

###
### means and se's for response ratios
###
synchrony_rr_subset_summary <- synchrony_rr_long %>%
  filter(trt_type %in% c("CO2", "drought", "irr", "N", "P", "temp")) %>%
  group_by(trt_type, metric_name) %>%
  summarize(lnRR_mean = mean(lnRR, na.rm=T),
            lnRR_se = sd(lnRR)/sqrt(length(lnRR)))


ggplot(synchrony_rr_subset_summary, aes(x=trt_type, y=lnRR_mean, ymin=lnRR_mean-lnRR_se,ymax=lnRR_mean+lnRR_se, col=factor(trt_type))) +
  geom_point() +
  geom_errorbar(width=0) +
  facet_wrap(~metric_name, scales="free") +
  theme_bw() +
  geom_hline(yintercept=0)

###
### t tests
###
spp_stab_model_out <- synch_RR_diff_metrics %>% 
  filter(trt_type2 %in% c("CO2", "drought", "Irrigation", "Precip. Vari.", "Temperature", "N", "P", "Mult. Nuts.")) %>%
  filter(metric_name == "spp_stab") %>%
  group_by(trt_type2) %>%
  do(data.frame(tidy(t.test(.$lnRR))))%>%
  mutate(metric = "spp_stab")

spp_synch_model_out <- synch_RR_diff_metrics %>% 
  filter(trt_type2 %in% c("CO2", "drought", "Irrigation", "Precip. Vari.", "Temperature", "N", "P", "Mult. Nuts.")) %>%
  filter(metric_name == "spp_synch") %>%
  group_by(trt_type2) %>%
  do(data.frame(tidy(t.test(.$lnRR))))%>%
  mutate(metric = "spp_synch")

alpha_stab_model_out <- synch_RR_diff_metrics %>% 
  filter(trt_type2 %in% c("CO2", "drought", "Irrigation", "Precip. Vari.", "Temperature", "N", "P", "Mult. Nuts.")) %>%
  filter(metric_name == "alpha_stab") %>%
  group_by(trt_type2) %>%
  do(data.frame(tidy(t.test(.$lnRR))))%>%
  mutate(metric = "alpha_stab")

spatial_synch_model_out <- synch_RR_diff_metrics %>% 
  filter(trt_type2 %in% c("CO2", "drought", "Irrigation", "Precip. Vari.", "Temperature", "N", "P", "Mult. Nuts.")) %>%
  filter(metric_name == "spatial_synch") %>%
  group_by(trt_type2) %>%
  do(data.frame(tidy(t.test(.$lnRR))))%>%
  mutate(metric = "spatial_synch")

pop_synch_model_out <- synch_RR_diff_metrics %>% 
  filter(trt_type2 %in% c("CO2", "drought", "Irrigation", "Precip. Vari.", "Temperature", "N", "P", "Mult. Nuts.")) %>%
  filter(metric_name == "pop_synch") %>%
  group_by(trt_type2) %>%
  do(data.frame(tidy(t.test(.$lnRR)))) %>%
  mutate(metric = "pop_synch")

gamma_stab_model_out <- synch_RR_diff_metrics %>% 
  filter(trt_type2 %in% c("CO2", "drought", "Irrigation", "Precip. Vari.", "Temperature", "N", "P", "Mult. Nuts.")) %>%
  filter(metric_name == "spp_stab") %>%
  group_by(trt_type2) %>%
  do(data.frame(tidy(t.test(.$lnRR))))%>%
  mutate(metric = "gamma_stab")

model_out_all <- spp_stab_model_out %>%
  bind_rows(spp_synch_model_out,
            alpha_stab_model_out,
            spatial_synch_model_out,
            pop_synch_model_out,
            gamma_stab_model_out) %>%
  mutate(bonf_adjP = p.adjust(p.value, method = "bonferroni"))
?p.adjust
write.table(model_out_all, file="..//model out//t test output_response ratios of stability metric.csv", row.names=F, sep=",")


###
### Look at response ratios across environmental gradients
###
synchrony_rr_long_with_env <- synchrony_rr_long %>%
  full_join(site_info, by=c("site_code", "project_name", "community_type"))

### Visualize RR across MAP MAT richness anpp and exp length
### Nitrogen alone
ggplot(filter(synchrony_rr_long_with_env, trt_type=="N"), aes(x=MAP, y=lnRR)) +
  geom_hline(yintercept=0) +
  geom_point(size=2) +
  geom_smooth(method='lm') +
  facet_wrap(~metric_name) +
  theme_bw()

ggplot(filter(synchrony_rr_long_with_env, trt_type=="N"), aes(x=MAT, y=lnRR)) +
  geom_hline(yintercept=0) +
  geom_point(size=2) +
  geom_smooth(method='lm') +
  facet_wrap(~metric_name) +
  theme_bw()

ggplot(filter(synchrony_rr_long_with_env, trt_type=="N"), aes(x=anpp, y=lnRR)) +
  geom_hline(yintercept=0) +
  geom_point(size=2) +
  geom_smooth(method='lm') +
  facet_wrap(~metric_name) +
  theme_bw()

ggplot(filter(synchrony_rr_long_with_env, trt_type=="N"), aes(x=experiment_length, y=lnRR)) +
  geom_hline(yintercept=0) +
  geom_point(size=2) +
  geom_smooth(method='lm') +
  facet_wrap(~metric_name) +
  theme_bw()

### other variables
ggplot(filter(synchrony_rr_long_with_env, trt_type=="temp"), aes(x=MAP, y=lnRR)) +
  geom_hline(yintercept=0) +
  geom_point(size=2) +
  geom_smooth(method='lm') +
  facet_wrap(~metric_name) +
  theme_bw()

ggplot(filter(synchrony_rr_long_with_env, trt_type=="temp"), aes(x=MAT, y=lnRR)) +
  geom_hline(yintercept=0) +
  geom_point(size=2) +
  geom_smooth(method='lm') +
  facet_wrap(~metric_name) +
  theme_bw()

ggplot(filter(synchrony_rr_long_with_env, trt_type=="temp"), aes(x=anpp, y=lnRR)) +
  geom_hline(yintercept=0) +
  geom_point(size=2) +
  geom_smooth(method='lm') +
  facet_wrap(~metric_name) +
  theme_bw()

ggplot(filter(synchrony_rr_long_with_env, trt_type=="temp"), aes(x=rrich, y=lnRR)) +
  geom_hline(yintercept=0) +
  geom_point(size=2) +
  geom_smooth(method='lm') +
  facet_wrap(~metric_name) +
  theme_bw()

ggplot(filter(synchrony_rr_long_with_env, trt_type=="temp"), aes(x=experiment_length, y=lnRR)) +
  geom_hline(yintercept=0) +
  geom_point(size=2) +
  geom_smooth(method='lm') +
  facet_wrap(~metric_name) +
  theme_bw()

###
### regression models
###
synchrony_rr_long_with_env %>%
  filter(trt_type %in% c("CO2", "drought", "irr", "N", "P", "temp")) %>%
  filter(metric_name == "spp_asynch") %>%
  nest(-trt_type)  %>%
  mutate(
    test = map(data, ~ cor.test(.x$anpp, .x$lnRR)), # S3 list-col
    tidied = map(test, tidy)
  ) %>% 
  unnest(tidied, .drop = TRUE)

###### select experiments that have >5 years data and >3 blocks (we used plots from each block with a given treatment (e.g. NPK) to form a metacommunity)

TH.year <- 5      ### only experiments with >5 years will be used (for all experiments, only the first 5 years were used)
TH.block <- 3     ### only experiments with >3 blocks will be used (for all experiments, only the first 3 blocks were used, to control for the effect of number of patches)



allsite <- unique(data$site_code)
site.ind <- lapply(allsite,function(x){
  datax <- data[data$site_code==x,]
  datax.year <- unique(datax$year_trt)
  datax.year <- datax.year[datax.year>0]
  datax.block <- unique(datax$block)
  c(length(datax.year)>=TH.year & length(datax.block)>=TH.block)
})
site.used <- allsite[unlist(site.ind)]

################################################################################################################
################################################################################################################
### 2. Shaopeng - Variance partitioning method
################################################################################################################
################################################################################################################

###### here is the function for calculating all metrics defined in Wang et al. (2019 Ecography)
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

#var.partition_live <- function(metacomm_tsdata){
#
#  ts_metacom_live <- apply(metacomm_tsdata,2,sum, na.rm=TRUE)
#  ts_patch_live <- apply(metacomm_tsdata,c(2,3),sum, na.rm=TRUE)
#
#  sd_metacom_live <- sd(ts_metacom_live, na.rm=TRUE)
#  sd_patch_k_live <- apply(ts_patch_live,2,sd, na.rm=TRUE)
#
#  mean_metacom <- mean(ts_metacom_live)
#
#  CV_C_L <- sum(sd_patch_k_live)/mean_metacom # alpha_var
#  CV_C_R <- sd_metacom_live/mean_metacom  # gamma_var
#
#  phi_C_L2R <- CV_C_R/CV_C_L  # spatial_synch
#
#  partition_3level <- c(CV_C_L_live=CV_C_L, CV_C_R_live=CV_C_R,
#                        phi_C_L2R_live=phi_C_L2R)
#  return(partition_3level)
#}

################################################################################################################
### End of Variance partitioning method
################################################################################################################



################################################################################################################
################################################################################################################
### 4. derive the stability and synchrony metrics for each experimental treatment
################################################################################################################
################################################################################################################

head(data)[1:9]

site.used
x = 'lancaster.uk'
Treatment = 'Control'

# note: no data for bnch.us, control, block 2, plot 20 year_trt 2
# note: doane.us. it looks like there is data on cover over time only for block 1? There is only data for 2012 for the other blocks

### Control
site.treatment_Control <- lapply(site.used, function(x, Treatment = 'Control'){
  datax <- data[data$site_code==x & data$trt==Treatment,]
  arrayx <- array(NA,dim=c(ncol(datax)-7,TH.year,TH.block))

  blockind <- sort(unique(datax$block))
  yearind <- sort(unique(datax$year_trt))
  yearind <- yearind[yearind>0]

  for(i in 1:TH.block)
  for(j in 1:TH.year)
  {
    tmpdata <- datax[datax$block==blockind[i] & datax$year_trt==yearind[j],-c(1:7)]
    if(nrow(tmpdata)==0){tmpdata <- rep(NA,ncol(datax)-7)}
    arrayx[,j,i] <- unlist(tmpdata)
  }

  var.partition(arrayx)
  })


### N
site.treatment_N <- lapply(site.used, function(x, Treatment = 'N'){
  datax <- data[data$site_code==x & data$trt==Treatment,]
  arrayx <- array(NA,dim=c(ncol(datax)-7,TH.year,TH.block))

  blockind <- sort(unique(datax$block))
  yearind <- sort(unique(datax$year_trt))
  yearind <- yearind[yearind>0]

  for(i in 1:TH.block)
  for(j in 1:TH.year)
  {
    tmpdata <- datax[datax$block==blockind[i] & datax$year_trt==yearind[j],-c(1:7)]
    if(nrow(tmpdata)==0){tmpdata <- rep(NA,ncol(datax)-7)}
    arrayx[,j,i] <- unlist(tmpdata)
  }

  var.partition(arrayx)
  })

### P
site.treatment_P <- lapply(site.used, function(x, Treatment = 'P'){
  datax <- data[data$site_code==x & data$trt==Treatment,]
  arrayx <- array(NA,dim=c(ncol(datax)-7,TH.year,TH.block))

  blockind <- sort(unique(datax$block))
  yearind <- sort(unique(datax$year_trt))
  yearind <- yearind[yearind>0]

  for(i in 1:TH.block)
  for(j in 1:TH.year)
  {
    tmpdata <- datax[datax$block==blockind[i] & datax$year_trt==yearind[j],-c(1:7)]
    if(nrow(tmpdata)==0){tmpdata <- rep(NA,ncol(datax)-7)}
    arrayx[,j,i] <- unlist(tmpdata)
  }

  var.partition(arrayx)
  })

### K
site.treatment_K <- lapply(site.used, function(x, Treatment = 'K'){
  datax <- data[data$site_code==x & data$trt==Treatment,]
  arrayx <- array(NA,dim=c(ncol(datax)-7,TH.year,TH.block))

  blockind <- sort(unique(datax$block))
  yearind <- sort(unique(datax$year_trt))
  yearind <- yearind[yearind>0]

  for(i in 1:TH.block)
  for(j in 1:TH.year)
  {
    tmpdata <- datax[datax$block==blockind[i] & datax$year_trt==yearind[j],-c(1:7)]
    if(nrow(tmpdata)==0){tmpdata <- rep(NA,ncol(datax)-7)}
    arrayx[,j,i] <- unlist(tmpdata)
  }

  var.partition(arrayx)
  })

### NP
site.treatment_NP <- lapply(site.used, function(x, Treatment = 'NP'){
  datax <- data[data$site_code==x & data$trt==Treatment,]
  arrayx <- array(NA,dim=c(ncol(datax)-7,TH.year,TH.block))

  blockind <- sort(unique(datax$block))
  yearind <- sort(unique(datax$year_trt))
  yearind <- yearind[yearind>0]

  for(i in 1:TH.block)
  for(j in 1:TH.year)
  {
    tmpdata <- datax[datax$block==blockind[i] & datax$year_trt==yearind[j],-c(1:7)]
    if(nrow(tmpdata)==0){tmpdata <- rep(NA,ncol(datax)-7)}
    arrayx[,j,i] <- unlist(tmpdata)
  }

  var.partition(arrayx)
  })

### NK
site.treatment_NK <- lapply(site.used, function(x, Treatment = 'NK'){
  datax <- data[data$site_code==x & data$trt==Treatment,]
  arrayx <- array(NA,dim=c(ncol(datax)-7,TH.year,TH.block))

  blockind <- sort(unique(datax$block))
  yearind <- sort(unique(datax$year_trt))
  yearind <- yearind[yearind>0]

  for(i in 1:TH.block)
  for(j in 1:TH.year)
  {
    tmpdata <- datax[datax$block==blockind[i] & datax$year_trt==yearind[j],-c(1:7)]
    if(nrow(tmpdata)==0){tmpdata <- rep(NA,ncol(datax)-7)}
    arrayx[,j,i] <- unlist(tmpdata)
  }

  var.partition(arrayx)
  })

### PK
site.used
x = 'lancaster.uk'
Treatment = 'PK'

site.treatment_PK <- lapply(site.used, function(x, Treatment = 'PK'){
  datax <- data[data$site_code==x & data$trt==Treatment,]
  arrayx <- array(NA,dim=c(ncol(datax)-7,TH.year,TH.block))

  blockind <- sort(unique(datax$block))
  yearind <- sort(unique(datax$year_trt))
  yearind <- yearind[yearind>0]

  for(i in 1:TH.block)
  for(j in 1:TH.year)
  {
    tmpdata <- datax[datax$block==blockind[i] & datax$year_trt==yearind[j],-c(1:7)]
    if(nrow(tmpdata)==0){
    tmpdata <- rep(NA,ncol(datax)-7)
    } else if(nrow(tmpdata)>1){
    tmpdata <- rep(NA,ncol(datax)-7)
    } # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    arrayx[,j,i] <- unlist(tmpdata)
  }

  var.partition(arrayx)
  })

### NPK
site.treatment_NPK <- lapply(site.used, function(x, Treatment = 'NPK'){
  datax <- data[data$site_code==x & data$trt==Treatment,]
  arrayx <- array(NA,dim=c(ncol(datax)-7,TH.year,TH.block))

  blockind <- sort(unique(datax$block))
  yearind <- sort(unique(datax$year_trt))
  yearind <- yearind[yearind>0]

  for(i in 1:TH.block)
  for(j in 1:TH.year)
  {
    tmpdata <- datax[datax$block==blockind[i] & datax$year_trt==yearind[j],-c(1:7)]
    if(nrow(tmpdata)==0){tmpdata <- rep(NA,ncol(datax)-7)}
    arrayx[,j,i] <- unlist(tmpdata)
  }

  var.partition(arrayx)
  })

## for biomass
#site.used
#x = 'rook.uk'
#Treatment = 'Control'
#
#site.treatment_Control_live <- lapply(site.used, function(x, Treatment = 'Control'){
#  datax <- data[data$site_code==x & data$trt==Treatment,]
#  arrayy <- array(NA,dim=c(1,TH.year,TH.block))
#
#  blockind <- sort(unique(datax$block))
#  yearind <- sort(unique(datax$year_trt))
#  yearind <- yearind[yearind>0]
#
#  for(i in 1:TH.block)
#  for(j in 1:TH.year)
#  {
#    tmpdata <- datax$live[datax$block==blockind[i] & datax$year_trt==yearind[j]]
#    if(length(tmpdata)==0){tmpdata <- rep(NA,1)}
#    arrayy[,j,i] <- unlist(tmpdata)
#  }
#
#  var.partition_live(arrayy)
#  })

varpart_Control <- t(data.frame(site.treatment_Control)); row.names(varpart_Control) <- site.used
varpart_Control <- data.frame(varpart_Control) 
varpart_Control$site_code <- site.used
varpart_Control$trt <- 'Control'

varpart_N <- t(data.frame(site.treatment_N)); row.names(varpart_N) <- site.used
varpart_N <- data.frame(varpart_N) 
varpart_N$site_code <- site.used
varpart_N$trt <- 'N'

varpart_P <- t(data.frame(site.treatment_P)); row.names(varpart_P) <- site.used
varpart_P <- data.frame(varpart_P) 
varpart_P$site_code <- site.used
varpart_P$trt <- 'P'

varpart_K <- t(data.frame(site.treatment_K)); row.names(varpart_K) <- site.used
varpart_K <- data.frame(varpart_K) 
varpart_K$site_code <- site.used
varpart_K$trt <- 'K'

varpart_NP <- t(data.frame(site.treatment_NP)); row.names(varpart_NP) <- site.used
varpart_NP <- data.frame(varpart_NP) 
varpart_NP$site_code <- site.used
varpart_NP$trt <- 'NP'

varpart_PK <- t(data.frame(site.treatment_PK)); row.names(varpart_PK) <- site.used
varpart_PK <- data.frame(varpart_PK) 
varpart_PK$site_code <- site.used
varpart_PK$trt <- 'PK'

varpart_NK <- t(data.frame(site.treatment_NK)); row.names(varpart_NK) <- site.used
varpart_NK <- data.frame(varpart_NK) 
varpart_NK$site_code <- site.used
varpart_NK$trt <- 'NK'
                     
varpart_NPK <- t(data.frame(site.treatment_NPK)); row.names(varpart_NPK) <- site.used
varpart_NPK <- data.frame(varpart_NPK) 
varpart_NPK$site_code <- site.used
varpart_NPK$trt <- 'NPK'

#varpart_Control_live <- t(data.frame(site.treatment_Control_live)); row.names(varpart_Control_live) <- site.used
#varpart_Control_live <- data.frame(varpart_Control_live) 
#varpart_Control_live$site_code <- site.used
#varpart_Control_live$trt <- 'Control'

### remove NAs
##varpart_Control <- na.omit(varpart_Control)
##varpart_NPK <- na.omit(varpart_NPK)
#
####### compare metrics between control vs. NPK treatments.
#plot(varpart_Control[,2],varpart_NPK[,2])    ### compare the average community-level stability between control vs. NPK: no difference!
#abline(0,1)
#plot(varpart_Control[,4],varpart_NPK[,4])    ### compare the total metacommunity stability between control vs. NPK: no difference!
#abline(0,1)
#plot(varpart_Control[,7],varpart_NPK[,7])    ### compare the (local) species synchrony between control vs. NPK: no difference!
#abline(0,1)
#plot(varpart_Control[,6],varpart_NPK[,6])    ### compare the spatial community synchrony (i.e. 1/beta variability) between control vs. NPK: no difference!
#abline(0,1)
#
#  ###### compare spatial synchrony between species-level (i.e. within metapopulation) vs. community-level (i.e. within metacommunity)
#  plot(ll<-varpart_Control[,6]~varpart_Control[,5]); summary(re<-lm(ll))  ### significant positive relation for Control(the dark-red arrow in SEM of slide 5)
#    abline(re$coef); abline(0,1,lty=2)
#  plot(ll<-varpart_NPK[,6]~varpart_NPK[,5]); summary(re<-lm(ll))  ### insignificant relation for NPK (the red arrow in SEM of slide 5)
#    abline(re$coef); abline(0,1,lty=2)
#
##   ###### compare spatial synchrony between species-level (i.e. within metapopulation) vs. community-level (i.e. within metacommunity)
##   plot(ll<-log(varpart_Control[,4])~log(varpart_Control[,6])); summary(re<-lm(ll))  ### r2=0.55
##     abline(re$coef);
##   plot(ll<-log(varpart_NPK[,4])~log(varpart_NPK[,6])); summary(re<-lm(ll))  ### r2=0.46
##     abline(re$coef); abline(0,1,lty=2)
##
##   ###### compare spatial synchrony between species-level (i.e. within metapopulation) vs. community-level (i.e. within metacommunity)
##   plot(ll<-log(varpart_Control[,4])~log(varpart_Control[,2])); summary(re<-lm(ll))  ### r2=0.93
##     abline(re$coef);
##   plot(ll<-log(varpart_NPK[,4])~log(varpart_NPK[,2])); summary(re<-lm(ll))  ### r2=0.75
##     abline(re$coef); abline(0,1,lty=2)


data_Control <- varpart_Control %>%
  left_join(alpha_diversity, by = c('site_code', 'trt')) %>%
  left_join(beta_dist, by = c('site_code', 'trt'))

data_N <- varpart_N %>%
  left_join(alpha_diversity, by = c('site_code', 'trt')) %>%
  left_join(beta_dist, by = c('site_code', 'trt'))

data_P <- varpart_P %>%
  left_join(alpha_diversity, by = c('site_code', 'trt')) %>%
  left_join(beta_dist, by = c('site_code', 'trt'))

data_K <- varpart_K %>%
  left_join(alpha_diversity, by = c('site_code', 'trt')) %>%
  left_join(beta_dist, by = c('site_code', 'trt'))

data_NP <- varpart_NP %>%
  left_join(alpha_diversity, by = c('site_code', 'trt')) %>%
  left_join(beta_dist, by = c('site_code', 'trt'))

data_NK <- varpart_NK %>%
  left_join(alpha_diversity, by = c('site_code', 'trt')) %>%
  left_join(beta_dist, by = c('site_code', 'trt'))

data_PK <- varpart_PK %>%
  left_join(alpha_diversity, by = c('site_code', 'trt')) %>%
  left_join(beta_dist, by = c('site_code', 'trt'))

data_NPK <- varpart_NPK %>%
  left_join(alpha_diversity, by = c('site_code', 'trt')) %>%
  left_join(beta_dist, by = c('site_code', 'trt'))

data_Control$spp_stab <- with(data_Control, 1/CV_S_L)
data_Control$alpha_stab <- with(data_Control, 1/CV_C_L)
data_Control$gamma_stab <- with(data_Control, 1/CV_C_R)
data_Control$pop_asynch <- with(data_Control, 1/phi_S_L2R)
data_Control$spatial_asynch <- with(data_Control, 1/phi_C_L2R)
data_Control$spp_asynch <- with(data_Control, 1/phi_S2C_L)

data_N$spp_stab <- with(data_N, 1/CV_S_L)
data_N$alpha_stab <- with(data_N, 1/CV_C_L)
data_N$gamma_stab <- with(data_N, 1/CV_C_R)
data_N$pop_asynch <- with(data_N, 1/phi_S_L2R)
data_N$spatial_asynch <- with(data_N, 1/phi_C_L2R)
data_N$spp_asynch <- with(data_N, 1/phi_S2C_L)

data_P$spp_stab <- with(data_P, 1/CV_S_L)
data_P$alpha_stab <- with(data_P, 1/CV_C_L)
data_P$gamma_stab <- with(data_P, 1/CV_C_R)
data_P$pop_asynch <- with(data_P, 1/phi_S_L2R)
data_P$spatial_asynch <- with(data_P, 1/phi_C_L2R)
data_P$spp_asynch <- with(data_P, 1/phi_S2C_L)

data_K$spp_stab <- with(data_K, 1/CV_S_L)
data_K$alpha_stab <- with(data_K, 1/CV_C_L)
data_K$gamma_stab <- with(data_K, 1/CV_C_R)
data_K$pop_asynch <- with(data_K, 1/phi_S_L2R)
data_K$spatial_asynch <- with(data_K, 1/phi_C_L2R)
data_K$spp_asynch <- with(data_K, 1/phi_S2C_L)

data_NP$spp_stab <- with(data_NP, 1/CV_S_L)
data_NP$alpha_stab <- with(data_NP, 1/CV_C_L)
data_NP$gamma_stab <- with(data_NP, 1/CV_C_R)
data_NP$pop_asynch <- with(data_NP, 1/phi_S_L2R)
data_NP$spatial_asynch <- with(data_NP, 1/phi_C_L2R)
data_NP$spp_asynch <- with(data_NP, 1/phi_S2C_L)

data_NK$spp_stab <- with(data_NK, 1/CV_S_L)
data_NK$alpha_stab <- with(data_NK, 1/CV_C_L)
data_NK$gamma_stab <- with(data_NK, 1/CV_C_R)
data_NK$pop_asynch <- with(data_NK, 1/phi_S_L2R)
data_NK$spatial_asynch <- with(data_NK, 1/phi_C_L2R)
data_NK$spp_asynch <- with(data_NK, 1/phi_S2C_L)

data_PK$spp_stab <- with(data_PK, 1/CV_S_L)
data_PK$alpha_stab <- with(data_PK, 1/CV_C_L)
data_PK$gamma_stab <- with(data_PK, 1/CV_C_R)
data_PK$pop_asynch <- with(data_PK, 1/phi_S_L2R)
data_PK$spatial_asynch <- with(data_PK, 1/phi_C_L2R)
data_PK$spp_asynch <- with(data_PK, 1/phi_S2C_L)

data_NPK$spp_stab <- with(data_NPK, 1/CV_S_L)
data_NPK$alpha_stab <- with(data_NPK, 1/CV_C_L)
data_NPK$gamma_stab <- with(data_NPK, 1/CV_C_R)
data_NPK$pop_asynch <- with(data_NPK, 1/phi_S_L2R)
data_NPK$spatial_asynch <- with(data_NPK, 1/phi_C_L2R)
data_NPK$spp_asynch <- with(data_NPK, 1/phi_S2C_L)

# remove NAs
data_Control <- na.omit(data_Control)
data_N <- na.omit(data_N)
data_P <- na.omit(data_P)
data_K <- na.omit(data_K)
data_NP <- na.omit(data_NP)
data_NK <- na.omit(data_NK)
data_PK <- na.omit(data_PK)
data_NPK <- na.omit(data_NPK)
data_Control$site_code
data_NPK$site_code

data_All <- rbind(data_Control, data_N, data_P, data_K, data_NP, data_NK, data_PK, data_NPK)
data_All$nb_trt <- factor(data_All$trt)
levels(data_All$nb_trt) <- c('0','1','1','2','2','3','1','2')
data_All$nb_trt <- as.numeric(as.character(data_All$nb_trt))

################################################################################################################
################################################################################################################
### 5. SEM analyses
################################################################################################################
################################################################################################################

####
####  LOG TRANSFORM METRICS FOR NORMALITY ----
####
data_Control <- data_Control %>%
  mutate(log_spp_asynch       = log(spp_asynch),
         log_spp_stab         = log(spp_stab),
         log_alpha_stab       = log(alpha_stab),
         log_spatial_asynch   = log(spatial_asynch),
         log_gamma_stab       = log(gamma_stab))

data_NPK <- data_NPK %>%
  mutate(log_spp_asynch       = log(spp_asynch),
         log_spp_stab         = log(spp_stab),
         log_alpha_stab       = log(alpha_stab),
         log_spatial_asynch   = log(spatial_asynch),
         log_gamma_stab       = log(gamma_stab))

data_All <- data_All %>%
  mutate(log_spp_asynch       = log(spp_asynch),
         log_spp_stab         = log(spp_stab),
         log_alpha_stab       = log(alpha_stab),
         log_spatial_asynch   = log(spatial_asynch),
         log_gamma_stab       = log(gamma_stab))

#####################
#  Control
#####################
model = psem(
  lm(log_gamma_stab ~ richness + log_alpha_stab + log_spatial_asynch, data=data_Control),
  lm(log_alpha_stab ~ log_spatial_asynch + richness + log_spp_stab + log_spp_asynch, data=data_Control),
  lm(log_spp_stab ~ beta_Bray + richness, data=data_Control),
  lm(log_spp_asynch ~ richness, data=data_Control),
  lm(log_spatial_asynch ~ richness + beta_Bray, data=data_Control),
  richness %~~% beta_Bray,
  log_spp_stab %~~% log_spp_asynch,
  log_spp_asynch %~~% log_spatial_asynch
)

summary(model, .progressBar = F) # by default the standardize estimates are shown as well as the unstandardized.

qplot(richness, log_alpha_stab, data=data_Control) + geom_smooth(method=lm)
qplot(richness, alpha_stab, data=data_Control) + geom_smooth(method=lm)
qplot(richness, log_gamma_stab, data=data_Control) + geom_smooth(method=lm)
qplot(richness, gamma_stab, data=data_Control) + geom_smooth(method=lm)
qplot(richness, log_spp_asynch, data=data_Control) + geom_smooth(method=lm)
qplot(alpha_stab, log_spp_asynch, data=data_Control) + geom_smooth(method=lm)

#####################
#  NPK
#####################
model = psem(
  lm(log_gamma_stab ~ richness + log_alpha_stab + log_spatial_asynch, data=data_NPK),
  lm(log_alpha_stab ~ log_spatial_asynch + richness + log_spp_stab + log_spp_asynch, data=data_NPK),
  lm(log_spp_stab ~ beta_Bray + richness, data=data_NPK),
  lm(log_spp_asynch ~ richness, data=data_NPK),
  lm(log_spatial_asynch ~ richness + beta_Bray, data=data_NPK),
  richness %~~% beta_Bray,
  log_spp_stab %~~% log_spp_asynch,
  log_spp_asynch %~~% log_spatial_asynch
)

summary(model, .progressBar = F) # by default the standardize estimates are shown as well as the unstandardized.

#####################
#  All
#####################
model = psem(
  lm(log_gamma_stab ~ log_spp_asynch + richness + nb_trt + log_alpha_stab + log_spatial_asynch, data=data_All),
  lm(log_alpha_stab ~ beta_Bray + richness + nb_trt + log_spp_stab + log_spp_asynch, data=data_All),
  lm(log_spp_stab ~ beta_Bray + nb_trt + richness, data=data_All),
  lm(log_spp_asynch ~ nb_trt + richness, data=data_All),
  lm(log_spatial_asynch ~ nb_trt + beta_Bray, data=data_All),
  lm(richness ~ nb_trt, data=data_All),
  lm(beta_Bray ~ nb_trt, data=data_All),
  richness %~~% beta_Bray,
  log_spp_stab %~~% log_spp_asynch,
  log_spp_asynch %~~% log_spatial_asynch
)

summary(model, .progressBar = F) # by default the standardize estimates are shown as well as the unstandardized.

qplot(spp_asynch, gamma_stab, data=data_All) + geom_smooth(method=lm)
qplot(richness, gamma_stab, data=data_All) + geom_smooth(method=lm)
qplot(alpha_stab, gamma_stab, data=data_All) + geom_smooth(method=lm)
qplot(spatial_asynch, gamma_stab, data=data_All) + geom_smooth(method=lm)

qplot(beta_Bray, alpha_stab, data=data_All) + geom_smooth(method=lm)
qplot(richness, alpha_stab, data=data_All) + geom_smooth(method=lm)
qplot(spp_stab, alpha_stab, data=data_All) + geom_smooth(method=lm)
qplot(spp_asynch, alpha_stab, data=data_All) + geom_smooth(method=lm)

qplot(beta_Bray, spp_stab, data=data_All) + geom_smooth(method=lm)

qplot(spp_stab, spp_asynch, data=data_All) + geom_smooth(method=lm)


################################################################################################################
################################################################################################################
### 3. Diversity indices
################################################################################################################
################################################################################################################

####
####  ALPHA DIVERSITY (SIMPSON'S D and SHANNON'S H) ----
####

cover_comm_id   <- site.used
cover_comm_trt   <- c('Control','N','P','K','NP','NK','PK','NPK') # levels(data$trt) # c('Control', 'NPK')

do_site <- 'bldr.us'
do_trt <- 'Control'

avg_richness <- data.frame()
avg_simpsons <- data.frame()
avg_shannons <- data.frame()

library(plyr) # have to do that due to problem of using the pipes from dplyr when plyr is loaded
i <- 1
for(do_site in cover_comm_id){
  for(do_trt in cover_comm_trt){
    richness_dat <- data %>%
      filter(site_code == do_site, trt == do_trt)%>%
      select(starts_with("sp"))
    if(dim(richness_dat)[1]>0) {
      avg_richness_t <- data.frame(richness=mean(specnumber(richness_dat)))
      avg_shannons_t <- data.frame(shannons=mean(diversity(richness_dat, index="shannon")))
      avg_simpsons_t <- data.frame(simpsons=mean(diversity(richness_dat, index="simpson")))
      avg_richness_t$site_code <- do_site
      avg_shannons_t$site_code <- do_site
      avg_simpsons_t$site_code <- do_site
      avg_richness_t$trt <- do_trt
      avg_shannons_t$trt <- do_trt
      avg_simpsons_t$trt <- do_trt
      avg_richness <- rbind.fill(avg_richness, avg_richness_t)
      avg_simpsons <- rbind.fill(avg_simpsons, avg_simpsons_t)
      avg_shannons <- rbind.fill(avg_shannons, avg_shannons_t)
      i <- i+1
    }
  }
}
detach("package:plyr", unload=TRUE)
alpha_diversity <- avg_richness %>%
  left_join(avg_simpsons, by = c('site_code', 'trt')) %>%
  left_join(avg_shannons, by = c('site_code', 'trt'))
alpha_diversity$trt <- factor(alpha_diversity$trt)
alpha_diversity$site_code <- factor(alpha_diversity$site_code)

alpha_diversity

####
####  BETA DIVERSITY (EUCLIDEAN DISTANCE)
####
do_site <- 'bldr.us'
do_trt <- 'Control'

beta_dist <- data.frame()

library(plyr)
for(do_site in cover_comm_id){
  for(do_trt in cover_comm_trt){
    site_dat <- data %>%
      filter(site_code == do_site, trt == do_trt)
    
    site_mat <- site_dat %>%
      select(starts_with("sp"))
    if(dim(site_mat)[1]>1) {
      
      site_mat_rel <- site_mat/rowSums(site_mat)
      site_mat_rel[is.na(site_mat_rel)] <- 0  #  add one line code to remove NA in data "site_mat_rel".
      vdist <- vegdist(site_mat_rel)
      beta_Bray <- vegdist(vdist)
      beta_dispersion <- betadisper(vdist, site_dat$year, type="centroid")
      distances_t <- data.frame(beta_Bray=mean(beta_Bray), euclid_dist=mean(beta_dispersion$distances))
      distances_t$site_code <- do_site
      distances_t$trt <- do_trt
      beta_dist <- rbind.fill(beta_dist,distances_t)
    }
  } 
}
detach("package:plyr", unload=TRUE)
beta_dist$trt <- factor(beta_dist$trt)
beta_dist$site_code <- factor(beta_dist$site_code)

beta_dist


# Select columns
coverDFsub <- coverDF[c('site_code', 'year', 'year_trt', 'block', 'plot', 'Taxon', 'max_cover')]
head(coverDFsub)

levels(coverDFsub$Taxon)
# remove unwanted Taxons
coverDFsubTax <- subset(coverDFsub, !(Taxon %in% c('BRYOPHYTE','BRYOPHYTE ','BRYOPHYTE SP.', 'DEER', 'FORB', 'FORB SP.',
                                                   'FUNGI', 'FUNGI SP.', 'GROUND', 'LICHEN', 'LICHEN SP.', 'OTHER','OTHER ANIMAL DIGGING', 'OTHER ANIMAL DIGGINGS',
                                                   'OTHER ANIMAL DROPPINGS', 'OTHER ANT NEST', 'OTHER ARISTIDA CONTORTA (DEAD)', 'OTHER CRUST', 'OTHER DISTURBED SOIL',
                                                   'OTHER LIGNOTUBER', 'OTHER LITTER', 'OTHER ROCK', 'OTHER SALSOLA KALI (DEAD)', 'OTHER SOIL BIOCRUST',
                                                   'OTHER STANDING DEAD', 'OTHER TRIODIA BASEDOWII (DEAD)', 'OTHER UNKNOWN SOIL_CRUST', 'OTHER WOOD',
                                                   'OTHER WOODY OVERSTORY', 'UNKNOWN', 'UNKNOWN ', 'UNKNOWN ACANTHACEAE','UNKNOWN AMARANTHACEAE SP.',
                                                   'UNKNOWN APIACEAE', 'UNKNOWN ASTERACEAE', 'UNKNOWN ASTERACEAE ', 'UNKNOWN ASTERACEAE SP.', 'UNKNOWN BRASSICACEAE',
                                                   'UNKNOWN BRASSICACEAE SP.', 'UNKNOWN CARYOPHYLLACEAE', 'UNKNOWN CARYOPHYLLACEAE SP.', 'UNKNOWN CRYPTOGAM',
                                                   'UNKNOWN CUPRESSACEAE SP.', 'UNKNOWN CYPERACEAE', 'UNKNOWN CYPERACEAE SP.', 'UNKNOWN DICOT', 'UNKNOWN EUPHORBIACEAE',
                                                   'UNKNOWN FABACEAE', 'UNKNOWN FABACEAE ', 'UNKNOWN FABACEAE SP.', 'UNKNOWN FORB', 'UNKNOWN GRASS',
                                                   'UNKNOWN GRASS ', 'UNKNOWN GRASS SP.', 'UNKNOWN LAMIACEAE', 'UNKNOWN LAMIACEAE ', 'UNKNOWN LILIACEAE',
                                                   'UNKNOWN ONAGRACEAE SP.', 'UNKNOWN ORCHIDACEAE SP.', 'UNKNOWN POLEMONIACEAE SP.', 'UNKNOWN POTTIACEAE SP.',
                                                   'UNKNOWN PTERIDOPHYTA', 'UNKNOWN ROSACEAE ', 'UNKNOWN RUBIACEAE', 'UNKNOWN SHRUB', 'UNKNOWN SP.')))
# drop unused levels
coverDFsubTax[] <- lapply(coverDFsubTax, function(x) x[drop=T])
levels(coverDFsubTax$Taxon)
head(coverDFsubTax)

# convert cover between 0 and 1 to 1
coverDFsubTax$max_cover[coverDFsubTax$max_cover>0 & coverDFsubTax$max_cover<1] <- 1

# Long- to wide-format data
coverDFsubTax_w <- reshape(coverDFsubTax, idvar = c('site_code', 'year', 'year_trt', 'block', 'plot'), timevar = "Taxon", direction = "wide")
head(coverDFsubTax_w)[1:9]

# replace NA with 0
coverDFsubTax_w[is.na(coverDFsubTax_w)] <- 0
head(coverDFsubTax_w)[1:9]

# Handy function to strip prefixes from the column names (of cosmetic value only!) since the reshape process will add a prefix:
colnames_removing_prefix <- function(df, prefix) {
  names <- colnames(df)
  indices <- (substr(names,1,nchar(prefix))==prefix)
  names[indices] <- substr(names[indices], nchar(prefix)+1, nchar(names[indices]))
  return(names)
}

colnames_replace_prefix <- function(df, prefix, replace) {
  names <- colnames(df)
  indices <- (substr(names,1,nchar(prefix))==prefix)
  names[indices] <- str_replace(names[indices], prefix, replace)
  return(names)
}

# replace prefix max_cover. by sp
colnames(coverDFsubTax_w) <- colnames_replace_prefix(coverDFsubTax_w, "max_cover.", "sp.")
head(coverDFsubTax_w)[1:9]
coverDFsubTax_w[coverDFsubTax_w$site_code=='bldr.us',][1:9]

# make a shorter name
spp <- coverDFsubTax_w
head(spp)[1:9]

data <- 
  # Combine
  data <- merge(combbyplot, spp)
data[1:20,c(1:10,ncol(data))]
data[data$site_code=='koffler.ca',][c(1:9)]

### only experimental years
#data <- data[data$year_trt != 0,]

data$ind <- with(data, paste(site_code, trt, plot, sep=''))

#koffler.ca: 3 control plots in each block, remove extra control plots
data[data$site_code=='koffler.ca'&data$trt=='Control',][c(1:9)]
data <- data[data$ind!='koffler.caControl11',]
data <- data[data$ind!='koffler.caControl9',]
data <- data[data$ind!='koffler.caControl17',]
data <- data[data$ind!='koffler.caControl21',]
data <- data[data$ind!='koffler.caControl34',]
data <- data[data$ind!='koffler.caControl36',]

#marc.ar: 3 control plots in blocks 1 and 2. Verified with J. Alberti on 170324, remove extra control plots
data[data$site_code=='marc.ar'&data$trt=='Control',][c(1:9)]
data <- data[data$ind!='marc.arControl6',]
data <- data[data$ind!='marc.arControl8',]
data <- data[data$ind!='marc.arControl17',]
data <- data[data$ind!='marc.arControl19',]

#sedg.us: 2 control, 2 NPK (no fences) in each block
data[data$site_code=='sedg.us'&data$trt=='Control',][c(1:9)]
data[data$site_code=='sedg.us'&data$trt=='NPK',][c(1:9)]
data <- data[data$ind!='sedg.usControl7',]
data <- data[data$ind!='sedg.usControl17',]
data <- data[data$ind!='sedg.usControl28',]
data <- data[data$ind!='sedg.usNPK10',]
data <- data[data$ind!='sedg.usNPK18',]
data <- data[data$ind!='sedg.usNPK27',]

#temple.us: extra plots in block 2
data[data$site_code=='temple.us'&data$block=='2',][c(1:9)]
data <- data[data$ind!='temple.usControl19',]
data <- data[data$ind!='temple.usNPK20',]

#ukul.za: 2 control, 2 NPK (no fences) in each block
data[data$site_code=='ukul.za'&data$trt=='Control',][c(1:9)]
data[data$site_code=='ukul.za'&data$trt=='NPK',][c(1:9)]
data <- data[data$ind!='ukul.zaControl8',]
data <- data[data$ind!='ukul.zaControl20',]
data <- data[data$ind!='ukul.zaControl30',]
data <- data[data$ind!='ukul.zaNPK10',]
data <- data[data$ind!='ukul.zaNPK19',]
data <- data[data$ind!='ukul.zaNPK25',]

data <- select(data, -ncol(data))
# drop unused levels
data[] <- lapply(data, function(x) x[drop=T])

data <- data[data$site_code!='doane.us',]





###
### Take a look at old patches site list and current full list to see what we removed
###

file.choose()
old_metrics_df <- read.csv("C:\\Users\\wilco\\Desktop\\Working groups\\C2E_May2019\\asynchrony\\Patches\\Patches\\full_analysis\\data\\Site level synch and var metrics_Jan2017_final.csv")
new_metrics_df <- read.csv("C:\\Users\\wilco\\Desktop\\Working groups\\C2E_May2019\\asynchrony\\Synchrony metrics response ratios_long form_13May2019.csv") %>%
  mutate(site_proj_comm = paste(site_code, project_name, community_type, sep=".."))

new_sites <- new_metrics_df %>%
  dplyr::select(site_proj_comm) %>%
  mutate(new_present=1) %>%
  unique()

sites_merged <- old_metrics_df %>%
  dplyr::select(site_proj_comm) %>%
  unique() %>%
  mutate(old_present=1) %>%
  full_join(new_sites, by="site_proj_comm") 

select_site_vector <- sites_merged %>%
  filter(old_present == 1) %>%
  .$site_proj_comm

diff_site_vec <- sites_merged %>%
  filter(is.na(old_present)) %>%
  .$site_proj_comm


synchrony_df_diff_sites <- synchrony_vars_df %>%
  mutate(site_proj_comm = paste(site_code, project_name, community_type, sep="_")) %>%
  filter(site_proj_comm %in% diff_site_vec)

plot_df <- full_df %>%
  dplyr::select(site_code:plot_id, community_type) %>%
  unique()

# CDR_e001_C
with(filter(plot_df, site_code=="CDR" & project_name == "e001" & community_type=="A"), table(calendar_year, factor(treatment)))

data <- read.csv("C:\\Users\\wilco\\Desktop\\Working groups\\C2E_May2019\\asynchrony\\data\\SpeciesRawAbundance_March2019.csv")

data %>%
  filter(site_code == "CAR" &
           project_name == "salt marsh" &)
