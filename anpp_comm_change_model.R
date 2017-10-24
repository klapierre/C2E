################################################################################
##  anpp_comm_change_model.R: script to construct a linear model between
##  mean community change and ANPP change.
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
library(cowplot)



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
####  BRING IN DATA ----
####
all_dat <- read.csv(paste0(data_dir,"SEM_allyr.csv"))

mod <- lm(anpp_PC_transform ~ mean_change_transform*treatment_year, data = all_dat)
predic(mod)
summary(mod)

ggplot(all_dat, aes(x=mean_change_transform, y=anpp_PC_transform, color=as.character(treatment_year)))+
  geom_point()+
  stat_smooth(method = "lm", se=FALSE)


