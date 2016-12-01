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




anpp_lrr <- readRDS("~/Desktop/anpp_timeseries.RDS")
bray_curtis_gams <- readRDS("~/Desktop/bray_curtis_gams.RDS")
bray_curtis_obs_only <- bray_curtis_gams[which(is.na(bray_curtis_gams$bc_dist)!=TRUE),]
colnames(bray_curtis_obs_only)[which(colnames(bray_curtis_obs_only)=="bc_dist")] <- "lrr"
bray_curtis_obs_only$metric <- "bray_curtis"
comm_metrics <- read.csv("~/Desktop/community_LRR.csv")
torms <- which(comm_metrics$metric %in% c("appearance", "disappearance","bp_dominance","richness","turnover"))
comm_metrics <- comm_metrics[-torms,]
comm_metrics$site_proj_comm <- with(comm_metrics, paste0(site_code, project_name, community_type))
comm_metrics <- comm_metrics[,c("site_proj_comm", "calendar_year","metric","lrr")]
comm_metrics <- rbind(comm_metrics, bray_curtis_obs_only[,c("site_proj_comm", "calendar_year","metric","lrr")])

ggplot(comm_metrics, aes(x=treatment_year, y=lrr))+
  geom_point(aes(group=site_proj_comm), alpha=0.5)+
  # geom_smooth(method="loess", se=FALSE, color="black", size=0.5)+
  geom_hline(aes(yintercept=0), linetype="dotted", color="grey25")+
  xlab("Year")+
  ylab("Log Response Ratio")+
  facet_grid(site_proj_comm~metric, scales = "free")+
  my_theme


plot_list <- list()
for(do_site in unique(comm_metrics$site_proj_comm)){
  temp_dat <- subset(comm_metrics, site_proj_comm==do_site)
  minyr <- min(temp_dat$treatment_year)
  maxyr <- max(temp_dat$treatment_year)
  ggplot(temp_dat, aes(x=treatment_year, y=lrr))+
    geom_point(aes(group=site_proj_comm), alpha=0.5)+
    # geom_smooth(method="loess", se=FALSE, color="black", size=0.5)+
    geom_hline(aes(yintercept=0), linetype="dotted", color="grey25")+
    xlab("")+
    ylab("")+
    scale_x_continuous(breaks = seq(minyr, maxyr,2))+
    facet_wrap("metric", scales = "free", nrow=1)+
    my_theme
}

test <- temp_dat$metric
test2 <- factor(test, levels=c("invD/S", "MRS", "loss", "gain"), labels = c("Evenness", "Rank Corr.", "Loss", "Gain"))




####
####  FIT GAMs TO EACH METRIC BY SITE
####
comm_metrics_sub <- subset(comm_metrics, site_proj_comm!="PIETIDE0" & site_proj_comm!="SERCCXN0")
gam_df <- list()
for(do_unit in unique(comm_metrics_sub$site_proj_comm)){
  site_data <- subset(comm_metrics_sub, site_proj_comm==do_unit)
  for(do_metric in unique(site_data$metric)){
    todo_data <- subset(site_data, metric==do_metric)
    nyrs      <- length(unique(todo_data$treatment_year))
    gam_fit   <- gam(lrr ~ s(treatment_year, k=nyrs-1), method="REML", 
                     data = todo_data,
                     select = TRUE) 
    years_to_predict   <- seq(min(unique(todo_data$treatment_year)),max(unique(todo_data$treatment_year)), by=0.1)
    pred_df            <- data.frame(treatment_year = years_to_predict)
    pred_df$yhat       <- predict(object = gam_fit, 
                                  newdata = pred_df, 
                                  type = "response")
    pred_df            <- merge(pred_df, todo_data, all.x = TRUE)
    pred_df$site_proj_comm <- do_unit
    pred_df$metric     <- do_metric
    gam_df             <- rbind(gam_df, pred_df)
  }
}


## Plot the fits by site_proj_com
ggplot(data=gam_df, aes(x=treatment_year))+
  geom_point(aes(y=lrr), alpha=0.5)+
  geom_line(aes(y=yhat))+
  xlab("Treatment year")+
  ylab("Value")+
  facet_grid(metric~site_proj_comm, scales = "free")+
  guides(color=FALSE)+
  my_theme
