# run_cumulative_gams.R
#  Script to run different versions of the cumulative GAM scripts.
#
# Author:
#  Andrew Tredennick (atredenn@gmail.com)


# Clear the workspace -----------------------------------------------------

rm(list = ls(all.names = TRUE))


# Load libraries ----------------------------------------------------------

library(tidyverse)
library(ggthemes)
library(mgcv)


# Set working directory info ----------------------------------------------

work_dir  <- "~/Repos/C2E/Community Paper/" # change as needed
data_dir  <- "~/Dropbox/C2E/Products/CommunityChange/March2018 WG/"
results_dir <- "~/Dropbox/C2E/Products/CommunityChange/Summer2018_Results/"
data_file <- "MetricsTrts_July2018.csv"
setwd(work_dir)


# Run GAM analysis for last year diffs ------------------------------------

diff_type <- "last_year"
source("./July2018/fit_cumulative_gam_interactions.R")


# Run GAM analysis for mid-year diffs -------------------------------------

diff_type <- "mid_year"
source("./July2018/fit_cumulative_gam_interactions.R")

