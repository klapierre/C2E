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
####  LOAD OUTPUT AND MAYBE PLOT ----
####
##  Get variance partitioning results
infile <- list.files(data_dir)[grep("varpart_output_",list.files(data_dir))]
vpart_out <- read.csv(file=paste0(data_dir,infile))

##  Get treatment key
trt_key <- read.csv(paste0(data_dir, "ExperimentInformation_May2017.csv")) %>%
  mutate(site_project_comm = paste(site_code,project_name,community_type,sep="_")) %>%
  mutate(other_trt = as.numeric(as.character(other_trt)),
         other_trt = ifelse(is.na(other_trt),1,0)) %>%
  mutate(all_other_trt = temp+mow_clip+burn+herb_removal+management+other_trt) %>%
  mutate(all_other_trt = ifelse(all_other_trt > 0,1,0)) %>%
  dplyr::select(site_project_comm, treatment, n:precip, all_other_trt) %>%
  unique()

##  Get raw metrics data for Bray-Curtis values
metrics_df <-  readRDS(paste0(data_dir,"bootstrapped_rac_metrics.RDS")) %>%
  dplyr::select(site_project_comm, treatment, calendar_year, bc_dissim)

##  Add in treatment info
vpart_results <- vpart_out %>%
  filter(site_project_comm != "KNZ_BGP_0") %>%
  left_join(trt_key, by = c("site_project_comm", "treatment")) %>%
  left_join(metrics_df, by = c("site_project_comm", "treatment", "calendar_year"))



####
####  SOME PRELIM PLOTS ----
####
### Bray-Curtis Time Series Plots

##  Nitrogen
bc_nitrogen <- vpart_results %>%
  filter(n > 0) %>%
  mutate(non_n = p+k+CO2+precip+all_other_trt) %>%
  filter(non_n == 0) %>%
  dplyr::select(site_project_comm, treatment_year,bc_dissim) %>%
  unique() %>%
  group_by(treatment_year) %>%
  summarise(avg_bc = mean(bc_dissim))

bc_n <- ggplot(bc_nitrogen, aes(x=treatment_year, y=avg_bc))+
  geom_line()+
  geom_point(size=3)+
  labs(x="Treatment Year", y="Bray-Curtis Change")+
  scale_y_continuous(limits=c(0.1,0.7))

##  Precipitation
bc_ppt <- vpart_results %>%
  filter(precip > 0) %>%
  mutate(non_n = p+k+CO2+n+all_other_trt) %>%
  filter(non_n == 0) %>%
  dplyr::select(site_project_comm, treatment_year,bc_dissim) %>%
  unique() %>%
  group_by(treatment_year) %>%
  summarise(avg_bc = mean(bc_dissim))

bc_precip <- ggplot(bc_ppt, aes(x=treatment_year, y=avg_bc))+
  geom_line()+
  geom_point(size=3)+
  labs(x="Treatment Year", y="Bray-Curtis Change")+
  scale_y_continuous(limits=c(0.1,0.7))

##  Phosphorous
bc_p <- vpart_results %>%
  filter(p > 0) %>%
  mutate(non_n = precip+k+CO2+n+all_other_trt) %>%
  filter(non_n == 0) %>%
  dplyr::select(site_project_comm, treatment_year,bc_dissim) %>%
  unique() %>%
  group_by(treatment_year) %>%
  summarise(avg_bc = mean(bc_dissim))

bc_phos <- ggplot(bc_p, aes(x=treatment_year, y=avg_bc))+
  geom_line()+
  geom_point(size=3)+
  labs(x="Treatment Year", y="Bray-Curtis Change")+
  scale_y_continuous(limits=c(0.1,0.7))

##  CO2
bc_co2 <- vpart_results %>%
  filter(CO2 > 0) %>%
  mutate(non_n = precip+k+p+n+all_other_trt) %>%
  filter(non_n == 0) %>%
  dplyr::select(site_project_comm, treatment_year,bc_dissim) %>%
  unique() %>%
  group_by(treatment_year) %>%
  summarise(avg_bc = mean(bc_dissim))

bc_carbon <- ggplot(bc_co2, aes(x=treatment_year, y=avg_bc))+
  geom_line()+
  geom_point(size=3)+
  labs(x="Treatment Year", y="Bray-Curtis Change")+
  scale_y_continuous(limits=c(0.1,0.7))

### Variance Partitioning Plots
##  Nitrogen only
vpart_means <- vpart_results %>%
  filter(n > 0) %>%
  mutate(non_n = p+k+CO2+precip+all_other_trt) %>%
  filter(non_n == 0) %>%
  group_by(treatment_year, metric) %>%
  summarise(var_expl_alone=mean(var_expl_alone,na.rm=T))

nitrogen <- ggplot(subset(vpart_means, metric!="total_var_expl"),
       aes(x=treatment_year, y=var_expl_alone, color=metric)) +
  # geom_col(position=position_dodge(0.9))+
  geom_line()+
  geom_point()+
  labs(x="Treatment Year", y="Variance Explained")+
  ggtitle("Nitrogen Additions Only")+
  scale_color_brewer(palette = "Set1", name = "")+
  theme(legend.position = c(0.3, 0.8))

##  Phosporous only
vpart_means <- vpart_results %>%
  filter(p > 0) %>%
  mutate(non_p = n+k+CO2+precip+all_other_trt) %>%
  filter(non_p == 0) %>%
  group_by(treatment_year, metric) %>%
  summarise(var_expl_alone=mean(var_expl_alone,na.rm=T))

phosphorous <- ggplot(subset(vpart_means, metric!="total_var_expl"),
       aes(x=treatment_year, y=var_expl_alone, color=metric)) +
  # geom_col(position=position_dodge(0.9))+
  geom_line()+
  geom_point()+
  labs(x="Treatment Year", y="Variance Explained")+
  ggtitle("Phosporous Additions Only")+
  scale_color_brewer(palette = "Set1", name = "")+
  guides(color=FALSE)

##  CO2 only
vpart_means <- vpart_results %>%
  filter(CO2 > 0) %>%
  mutate(non_c = n+p+k+precip+all_other_trt) %>%
  filter(non_c == 0) %>%
  group_by(treatment_year, metric) %>%
  summarise(var_expl_alone=mean(var_expl_alone,na.rm=T))

co2 <- ggplot(subset(vpart_means, metric!="total_var_expl"),
       aes(x=treatment_year, y=var_expl_alone, color=metric)) +
  # geom_col(position=position_dodge(0.9))+
  geom_line()+
  geom_point()+
  labs(x="Treatment Year", y="Variance Explained")+
  ggtitle("CO2 Treatment Only")+
  scale_color_brewer(palette = "Set1", name = "")+
  guides(color=FALSE)

##  Precip only
vpart_means <- vpart_results %>%
  filter(precip > 0) %>%
  mutate(non_c = n+p+k+CO2+all_other_trt) %>%
  filter(non_c == 0) %>%
  group_by(treatment_year, metric) %>%
  summarise(var_expl_alone=mean(var_expl_alone,na.rm=T))

precip <- ggplot(subset(vpart_means, metric!="total_var_expl"),
       aes(x=treatment_year, y=var_expl_alone, color=metric)) +
  # geom_col(position=position_dodge(0.9))+
  geom_line()+
  geom_point()+
  labs(x="Treatment Year", y="Variance Explained")+
  ggtitle("Precipitation Treatment Only")+
  scale_color_brewer(palette = "Set1", name = "")+
  guides(color=FALSE)

# plot_grid(nitrogen, bc_n, phosphorous, bc_phos, co2, bc_carbon, precip, bc_precip, nrow = 4, ncol = 2)

bc_dynamics <- plot_grid(nitrogen, phosphorous, co2, precip, 
                         bc_n, bc_phos, bc_carbon, bc_precip, 
                         nrow = 2, ncol = 4)

outfile <- paste0(data_dir,"bc_dynamics.pdf")
ggsave(filename = outfile, plot = bc_dynamics, width = 14, height = 6, units = "in")


# 
# ggplot(subset(vpart_out, metric!="total_var_expl"&site_project_comm=="KNZ_pplots_0"),
#        aes(x=treatment_year, y=var_expl_together, colour=metric)) +
#   geom_line() + facet_wrap(~treatment) + theme_few()
# 
