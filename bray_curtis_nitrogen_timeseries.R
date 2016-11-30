##  bray_curtis_nitrogen_timeseries.R: script to calculate Bray-Curtis 
##  distances between treatment and control plots. Then, uses GAMs to 
##  fit temporal trends to the time series. We only use datasets that
##  are at least 7 years long.
##
##  Author: Andrew Tredennick (atredenn@gmail.com)
##  Date created: November 30, 2016
##

rm(list=ls(all.names = TRUE))

diffwd <- "/Users/atredenn/Google Drive/C2E/data/"

####
####  LOAD LIBRARIES
####
library(plyr)
library(reshape2)
library(mgcv)
library(ggplot2)
library(ggthemes)



####
####  SOURCE DATA SUBSETTING SCRIPT
####
source("subsetting 7 yr N data.R")