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
library(gridExtra)



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
        strip.text         = element_text(size=10, color="grey35", family = "Arial Narrow"),
        plot.title         = element_text(size=12, family = "Arial Narrow", face = "bold"))




anpp_lrr <- readRDS("~/Desktop/anpp_timeseries.RDS")
anpp_lrr$site_proj_comm <- with(anpp_lrr, paste0(site_code, project_name, community_type))
anpp_lrr$metric <- "ANPP"
sev_yrs <- length(anpp_lrr[which(anpp_lrr$site_proj_comm=="SEVNfert0"), "treatment_year"])
anpp_lrr[which(anpp_lrr$site_proj_comm=="SEVNfert0"), "treatment_year"] <- c(10:(10+sev_yrs-1))

bray_curtis_gams <- readRDS("~/Desktop/bray_curtis_gams.RDS")
bray_curtis_obs_only <- bray_curtis_gams[which(is.na(bray_curtis_gams$bc_dist)!=TRUE),]
colnames(bray_curtis_obs_only)[which(colnames(bray_curtis_obs_only)=="bc_dist")] <- "lrr"
colnames(bray_curtis_obs_only)[which(colnames(bray_curtis_obs_only)=="treat_year")] <- "treatment_year"
bray_curtis_obs_only$metric <- "bray_curtis"

comm_metrics <- read.csv("~/Desktop/community_metrics.csv")
# torms <- which(comm_metrics$metric %in% c("appearance", "disappearance","bp_dominance","richness","turnover"))
comm_metrics <- comm_metrics[-torms,]
comm_metrics$site_proj_comm <- with(comm_metrics, paste0(site_code, project_name, community_type))

clim_data <- read.csv("/Users/atredenn/Google Drive/C2E/data/clim_data_yearly.csv")
clim_data$calendar_year <- clim_data$year
clim_data <- ddply(clim_data, .(site_code, calendar_year), summarise,
                   ann_ppt = mean(ann_ppt))
ttt <- merge(comm_metrics, clim_data[,c("ann_ppt","site_code","calendar_year")], all.x = TRUE)
ttt$metric <- "ppt"
ttt$lrr <- ttt$ann_ppt+1

comm_metrics <- comm_metrics[,c("site_proj_comm", "treatment_year","metric","lrr")]
comm_metrics <- rbind(comm_metrics, bray_curtis_obs_only[,c("site_proj_comm", "treatment_year","metric","lrr")])
comm_metrics <- rbind(comm_metrics, anpp_lrr[,c("site_proj_comm", "treatment_year","metric","lrr")])
comm_metrics <- rbind(comm_metrics, ttt[,c("site_proj_comm", "treatment_year","metric","lrr")])
# 
# comm_metrics$metric <- factor(comm_metrics$metric, 
#                               levels=c("ANPP","bray_curtis","invD/S", "RankCor", "loss", "gain", "ppt"), 
#                               labels = c("ANPP","Bray-Curtis","Evenness", "Rank Corr.", "Loss", "Gain", "Precip."))
# 

sites_to_plot <- intersect(unique(bray_curtis_obs_only$site_proj_comm), unique(anpp_lrr$site_proj_comm))
comm_metrics <- comm_metrics[which(comm_metrics$site_proj_comm %in% sites_to_plot),]
# comm_metrics[which(comm_metrics$metric=="Precip."&comm_metrics$lrr==0),"lrr"]





####
####  FIT GAMs TO EACH METRIC BY SITE
####
gam_df <- list()
for(do_unit in unique(comm_metrics$site_proj_comm)){
  site_data <- subset(comm_metrics, site_proj_comm==do_unit)
  for(do_metric in unique(site_data$metric)){
    todo_data <- subset(site_data, metric==do_metric)
    nyrs      <- length(unique(todo_data$treatment_year))
    gam_fit   <- gam(lrr ~ s(treatment_year, k=nyrs-1), method="REML", 
                     data = todo_data,
                     select = TRUE) 
    if(do_metric=="Precip."){
      gam_fit   <- gam(lrr ~ s(treatment_year, k=nyrs-1), method="REML", 
                       data = todo_data,
                       select = TRUE,family = Gamma) 
    }
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



gam_sub <- subset(gam_df, site_proj_comm=="KNZpplots0" & metric%in%c("invD/S","RankCor","bp_dominance","appearance","disappearance", "bray_curtis"))
ggplot(data = gam_sub, aes(x=treatment_year))+
  geom_line(aes(y=yhat, color=metric), size=1)


