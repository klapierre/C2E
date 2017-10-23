################################################################################
##
##  peak_detection.R: just trying this out for now...
##
##  Authors: Andrew Tredennick (atredenn@gmail.com)
##
##  Date created: October 22, 2017
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
data_dir <- "/Users/atredenn/Google_Drive/C2E/"
out_dir  <- "/Users/atredenn/Google_Drive/C2E/"
setwd(work_dir)



####
####  LOAD LIBRARIES ----
####
library(tidyverse)
library(peakPick)



####
####  LOAD AND FILTER METRICS DATA ----
####
##  Bring in the metrics data
fname       <- "varpart_output_with MRSc_15Oct2017.csv"
varpar_data <- read.csv(paste0(data_dir,fname))

##  Find the site_project_comms with 10 or more years of data
long_sites <- varpar_data %>%
  group_by(site_project_comm) %>%
  summarise(first_year = min(calendar_year),
            last_year  = max(calendar_year)) %>%
  mutate(year_diff = last_year - first_year) %>%
  dplyr::filter(year_diff >= 20)

##  Filter the metrics data to only long time series (10 years or more)
long_data <- varpar_data %>%
  dplyr::filter(site_project_comm %in% long_sites$site_project_comm) %>%
  dplyr::mutate(comm_trt = paste(site_project_comm, treatment, sep = "::"))



####
####  LOOK AT SINGLE COMMUNITY::TREATMENT ----
####
all_groups <- unique(long_data$comm_trt)
do_group <- all_groups[1]
group_dat <- dplyr::filter(long_data, comm_trt == do_group) %>%
  dplyr::select(-appearance_all_var_expl,
                -disappearance_all_var_expl,
                -MRSc_all_var_expl,
                -total_var_expl) %>%
  gather(key = metric, value = value,
         c(appearance_alone_var_expl,
         disappearance_alone_var_expl,
         MRSc_alone_var_expl))

ggplot(group_dat, aes(x=treatment_year, y=value, color=metric))+
  geom_line()



####
####  TRY PEAK DETECTION ----
####
ts_dat <- dplyr::filter(long_data, comm_trt == do_group) %>%
  dplyr::select(-appearance_all_var_expl,
                -disappearance_all_var_expl,
                -MRSc_all_var_expl,
                -total_var_expl)
ts_mat <- as.matrix(ts_dat[,5:7])
peaks <- peakPick::detect.spikes(ts_mat,roi = c(2,nrow(ts_mat)-1), winlen = 1, spike.min.sd=1.8)
which(peaks[, 1])
