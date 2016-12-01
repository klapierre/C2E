##  anpp_nitrogen_patterns.R: script to calculate ANPP differences among control and 
##  nitrogen treatments, then look at temporal patterns with GAMs.
##
##


rm(list=ls(all.names = TRUE))


####
####  LOAD LIBRARIES
####
library(plyr)
library(reshape2)
library(mgcv)
library(vegan)
library(ggplot2)
library(ggthemes)



####
####  MY PLOTTING THEME --------------------------------------------------------
####
my_theme <- theme_bw()+
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_line(color="white"),
        panel.background   = element_rect(fill = "#EFEFEF"),
        axis.text          = element_text(size=10, color="grey35", family = "Arial Narrow"),
        axis.title         = element_text(size=12, family = "Arial Narrow", face = "bold"),
        panel.border       = element_blank(),
        axis.line.x        = element_line(color="black"),
        axis.line.y        = element_line(color="black"),
        strip.background   = element_blank(),
        strip.text         = element_text(size=10, color="grey35", family = "Arial Narrow"))



####
####  SOURCE DATA SUBSETTING SCRIPT
####
diffwd <- "/Users/atredenn/Google Drive/C2E/data/"
source("subsetting 7 yr N data_anpp.R")
setwd("~/Repos/C2E/")
relAbundN$site_proj_comm <- with(relAbundN, paste0(site_code, project_name, community_type))


### Aggregate over plots
anppAgg <- relAbundN%>%
  ungroup()%>%
  select(site_code, project_name, community_type, treatment, treatment_year, calendar_year, plot_id, anpp, n, n_treatment)%>%
  group_by(site_code, project_name, community_type, treatment, treatment_year, calendar_year, n, n_treatment)%>%
  summarise(anpp_agg=mean(anpp))

anpp_sub <- anppAgg[,c(-4,-7)]

#calculate effect size (as ln response ratio)
anppLRR <- anpp_sub%>%
  spread(key=n_treatment, value=anpp_agg)%>%
  mutate(lrr=log(`1`/`0`))


anppLRR <- subset(anppLRR, community_type!="WetBowman")
anppLRR <- subset(anppLRR, project_name!="RMAPC" & project_name!="BGP")

saveRDS(anppLRR, "~/Desktop/anpp_timeseries.RDS")

ggplot(anppLRR, aes(x=calendar_year, y=lrr))+
  geom_point()+
  geom_line()+
  facet_wrap("project_name")
