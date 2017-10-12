################################################################################
##
##  bootstrap_rac_metrics.R: Script to bootstrap within-year replicates of
##  community change metrics. Basic idea is to use the means and covariances
##  of the data to simulate larger sample sizes, with which we can conduct
##  variance partitioning. The end goal is to create a time series showing the
##  relative importance of different metrics of community change.
##
##  This script produces a simulated dataset where each site_proj_comm of the 
##  CoRRE database has 20 simulated within-year replicates. The result is a
##  long dataframe, saved as an RDS file named "bootstrapped_rac_metrics.RDS".
##
##  Authors: Andrew Tredennick (atredenn@gmail.com)
##           Kevin Wilcox (wilcoxkr@gmail.com)
##           Chris Lortie (pitched original idea)
##
##  Date created: October 10, 2017
##
################################################################################


##  Clear everything
rm(list = ls(all.names = TRUE))

##  Set number of desired replicates within a year
nsims <- 20


####
####  SET DIRECTORY INFORMATION ----
####
work_dir <- dirname(rstudioapi::getActiveDocumentContext()$path) # to file location
data_dir <- "/Users/atredenn/Google_Drive/C2E/data/"
out_dir  <- "/Users/atredenn/Google_Drive/C2E/"
setwd(work_dir)



####
####  LOAD LIBRARIES ----
####
library(tidyverse)
library(MASS)



####
####  LOAD AND FILTER METRICS DATA ----
####
##  Bring in the metrics data
fname    <- "CORRE_RAC_Metrics_Oct2017_allReplicates.csv"
rac_data <- read.csv(paste0(data_dir,fname), row.names = 1)

##  Find the site_project_comms with 10 or more years of data
long_sites <- rac_data %>%
  group_by(site_project_comm) %>%
  summarise(first_year = min(calendar_year),
            last_year  = max(calendar_year)) %>%
  mutate(year_diff = last_year - first_year) %>%
  dplyr::filter(year_diff >= 10)

##  Filter the metrics data to only long time series (10 years or more)
### TODO: DISCUSS REMOVING THIS SERC DATASET !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
long_data <- rac_data %>%
  # dplyr::filter(site_project_comm %in% long_sites$site_project_comm) %>%
  dplyr::filter(site_project_comm != "SERC_TMECE_SP" & treatment != "E") # nasty data



####
####  LOOP THROUGH SITE_PROJECT_COMMs AND BOOTSTRAP ----
####
groups       <- unique(long_data$site_project_comm)
boots_rac_df <- {} # empty object for storage

##  Progress bar
pb <- txtProgressBar(min=1, max=length(groups), char="+", style=3, width=65) 
counter <- 1

for(do_group in groups){
  proj_data <- long_data %>%
    dplyr::filter(site_project_comm == do_group)
  
  treatments <- unique(proj_data$treatment)
  
  for(do_treat in treatments){
    treat_data <- proj_data %>%
      dplyr::filter(treatment == do_treat)
    
    years <- unique(treat_data$calendar_year)
    
    for(do_year in years){
      year_data <- treat_data %>%
        dplyr::filter(calendar_year == do_year) %>%
        dplyr::select(appearance, disappearance, S, E_q, MRSc, bc_dissim) %>%
        mutate(S = as.numeric(S)) # coerce S to numeric for matrix conversion
      metric_names <- colnames(year_data)
      
      ##  Calculate empirical means and covariances
      ### TODO: DISCUSS WHAT TO DO ABOUT NA'S IN E_q !!!!!!!!!!!!!!!!!!!!!!!!!!
      ### TODO: IS 'S' REALLY 1 IN SOME CASES??????? !!!!!!!!!!!!!!!!!!!!!!!!!!
      year_means <- as.numeric(colMeans(year_data, na.rm = TRUE))
      bad_metric <- which(is.nan(year_means))
      bad_metric_name <- metric_names[bad_metric]
      
      if(length(bad_metric) > 0){
        bad_name   <- colnames(year_data)[bad_metric]
        year_data  <- year_data[-bad_metric]
        year_means <- as.numeric(colMeans(year_data, na.rm = TRUE))
        year_covar <- cov(as.matrix(year_data), use = "complete.obs")
      }else{
        year_covar <- cov(as.matrix(year_data), use = "complete.obs")
      }
      
      if(!is.na(sum(year_covar))){
        ##  Simulate data
        tmp_sim <- MASS::mvrnorm(n = nsims, mu = year_means, Sigma = year_covar)
        
        ##  Store as data.frame
        tmp_df <- as.data.frame(tmp_sim) %>%
          mutate(site_project_comm = do_group,
                 treatment      = do_treat,
                 calendar_year  = do_year) 
        
        if(length(bad_metric) > 0){
          tmp_df[,bad_metric_name] <- NA
        }
        
        tmp_df <- tmp_df %>%
          dplyr::select(site_project_comm, treatment, calendar_year, appearance,
                        disappearance, S, E_q, MRSc, bc_dissim)
        
        ##  rbind() the temporary data frame to the storage data frame
        boots_rac_df <- rbind(boots_rac_df, tmp_df)
      }
      
    } # end year loop, within treatment and site_project_comm
    
  } # end treatment loop, within site_project_comm
  
  ##  Update progess
  setTxtProgressBar(pb, counter)
  counter <- counter+1
  
} # end site_project_comm loop

##  Save the emprically-based simulated data
saveRDS(boots_rac_df, paste0(out_dir,"bootstrapped_rac_metrics.RDS"))
