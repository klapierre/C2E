################################################################################
##  plot_varpart.R: script to look at results from bootstrapped variance
##  partitioning.
##
##  Author: Andrew Tredennick (atredenn@gmail.com)
##          Kevin Wilcox
##  Date created: 10-24-2017
################################################################################
##  Clear everything...
rm(list=ls())


####
####  LOAD LIBRARIES ----
####
library(tidyverse)
library(vegan)
library(ggplot2)
library(ggthemes)



####
####  SETUP DIRECTORIES ----
##  KW
#setwd("C:\\Users\\Kevin.Wilcox\\Desktop\\C2E\\BC drivers\\")
#setwd("C:\\Users\\WilcoxKR.WILCOXKR-LAPTOP\\Desktop\\Working groups\\C2E\\BC drivers\\")
#data_dir <- paste0(getwd(),"/data/")

##  AT
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
data_dir <- "../../Google_Drive/C2E/data/"



####
####  LOAD OUTPUT AND MAYBE PLOT ----
####
varpart_results <- read.csv(vpart_out, file=paste0(data_dir,"varpart_output_",Sys.Date(),".csv"))

ggplot(subset(vpart_out, metric!="total_var_expl"&site_project_comm=="KNZ_pplots_0"),
       aes(x=treatment_year, y=var_expl_together, colour=metric)) +
  geom_line() + facet_wrap(~treatment) + theme_few()

vpart_means <- vpart_out %>%
  group_by(treatment_year, metric) %>%
  summarise(var_expl_alone=mean(var_expl_alone,na.rm=T))

ggplot(subset(vpart_out, metric!="total_var_expl"),
       aes(x=treatment_year, y=var_expl_together, colour=metric)) +
  geom_line() + theme_few()

ggplot(subset(vpart_means, metric!="total_var_expl"),
       aes(x=treatment_year, y=var_expl_alone, fill=metric)) +
  geom_col(position=position_dodge(0.9))

