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




bray_curtis_gams <- readRDS("~/Desktop/bray_curtis_gams.RDS")
mean_rank_shifts <- readRDS("~/Desktop/rankdf.RDS")
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
x <- strsplit(as.character(ttt$year_pair), "-")
ttt$treatment_year <- unlist(lapply(x, `[[`, 1))

comm_metrics <- comm_metrics[,c("site_proj_comm", "treatment_year", "lrr", "metric")]
comm_metrics <- rbind(comm_metrics, ttt[,c("site_proj_comm", "treatment_year", "lrr", "metric")])
comm_metrics$treatment_year <- as.numeric(comm_metrics$treatment_year)

# ggplot(comm_metrics, aes(x=treatment_year, y=lrr))+
#   geom_line(aes(group=site_proj_comm))+
#   geom_hline(aes(yintercept=0), linetype="dotted", color="grey25")+
#   xlab("Year")+
#   ylab("Log Response Ratio")+
#   facet_grid(metric~site_proj_comm, scales = "free")+
#   my_theme
# ggsave(filename = "~/Desktop/comm_metrics_ts.png", width = 10, height = 6, units="in", dpi = 120)


ggplot(comm_metrics, aes(x=treatment_year, y=lrr))+
  geom_point(aes(group=site_proj_comm), alpha=0.5)+
  geom_smooth(method="loess", se=FALSE, color="black", size=0.5)+
  geom_hline(aes(yintercept=0), linetype="dotted", color="grey25")+
  xlab("Year")+
  ylab("Log Response Ratio")+
  facet_grid(metric~site_proj_comm, scales = "free")+
  my_theme
ggsave(filename = "~/Desktop/comm_metrics_ts_smooths.png", width = 10, height = 6, units="in", dpi = 120)




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
