##  plot_community_dynamics_nitrogen.R


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
library(viridis)



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




bray_curtis_gams <- readRDS("bray_curtis_gams.RDS")
mean_rank_shifts <- readRDS("rankdf.RDS")
comm_metrics <- read.csv("~/Desktop/community_LRR.csv")
torms <- which(comm_metrics$metric %in% c("appearance", "disappearance"))
comm_metrics <- comm_metrics[-torms,]
comm_metrics$site_proj_comm <- with(comm_metrics, paste0(site_code, project_name, community_type))

trt_mrs <- subset(mean_rank_shifts, trt == 1)
cntl_mrs <- subset(mean_rank_shifts, trt == 0)
cntl_mrs$MRS.corr.control <- cntl_mrs$MRS.corr
back_mrs <- cbind(trt_mrs, cntl_mrs$MRS.corr.control)
colnames(back_mrs)[ncol(back_mrs)] <- "MRS.corr.control"
back_mrs$lrr <- with(back_mrs, log(MRS.corr/MRS.corr.control))


ggplot(data=subset(bray_curtis_gams), aes(x=treat_year))+
  geom_point(aes(y=bc_dist), alpha=0.5)+
  geom_line(aes(y=yhat))+
  xlab("Treatment year")+
  ylab("Bray-Curtis Dissimilarity")+
  facet_wrap("site_proj_comm", scales = "free", nrow=1)+
  guides(color=FALSE)+
  my_theme
ggsave(filename = "~/Desktop/bray_curtis_ts_gams.png", width = 10, height = 2, units="in", dpi = 120)


ggplot(back_mrs, aes(x=as.numeric(year_pair), y=lrr))+
  geom_line(aes(group=site_proj_comm))+
  geom_hline(aes(yintercept=0), linetype="dotted", color="grey25")+
  xlab("Year")+
  ylab("Mean Rank Shift (LRR)")+
  facet_wrap("site_proj_comm", nrow=1, scales = "free_x")+
  my_theme+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=8))
ggsave(filename = "~/Desktop/mrs_ts.png", width = 10, height = 2, units="in", dpi = 120)


ttt <- back_mrs[,c("site_proj_comm", "year_pair", "lrr")]
ttt$metric <- "MRS"
ttt$treatment_year <- as.numeric(ttt$year_pair)-1

comm_metrics <- comm_metrics[,c("site_proj_comm", "treatment_year", "lrr", "metric")]
comm_metrics <- rbind(comm_metrics, ttt[,c("site_proj_comm", "treatment_year", "lrr", "metric")])

ggplot(comm_metrics, aes(x=treatment_year, y=lrr))+
  geom_line(aes(group=site_proj_comm))+
  geom_hline(aes(yintercept=0), linetype="dotted", color="grey25")+
  xlab("Year")+
  ylab("Log Response Ratio")+
  facet_grid(metric~site_proj_comm, scales = "free")+
  my_theme
ggsave(filename = "~/Desktop/comm_metrics_ts.png", width = 10, height = 6, units="in", dpi = 120)


ggplot(comm_metrics, aes(x=treatment_year, y=lrr))+
  geom_point(aes(group=site_proj_comm), alpha=0.5)+
  geom_smooth(method="loess", se=FALSE, color="black")+
  geom_hline(aes(yintercept=0), linetype="dotted", color="grey25")+
  xlab("Year")+
  ylab("Log Response Ratio")+
  facet_grid(metric~site_proj_comm, scales = "free")+
  my_theme
ggsave(filename = "~/Desktop/comm_metrics_ts_smooths.png", width = 10, height = 10, units="in", dpi = 120)
