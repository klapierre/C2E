##  bray_curtis_nitrogen_timeseries.R: script to calculate Bray-Curtis 
##  distances between treatment and control plots. Then, uses GAMs to 
##  fit temporal trends to the time series. We only use datasets that
##  are at least 7 years long.
##
##  Author: Andrew Tredennick (atredenn@gmail.com)
##  Date created: November 30, 2016
##

rm(list=ls(all.names = TRUE))



####
####  LOAD LIBRARIES
####
