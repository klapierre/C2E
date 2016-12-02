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
sev_yrs <- nrow(anpp_lrr[which(anpp_lrr$site_proj_comm=="SEVNfert0"), "treatment_year"])
anpp_lrr[which(anpp_lrr$site_proj_comm=="SEVNfert0"), "treatment_year"] <- c(10:(10+sev_yrs-1))

bray_curtis_gams <- readRDS("~/Desktop/bray_curtis_gams.RDS")
bray_curtis_obs_only <- bray_curtis_gams[which(is.na(bray_curtis_gams$bc_dist)!=TRUE),]
colnames(bray_curtis_obs_only)[which(colnames(bray_curtis_obs_only)=="bc_dist")] <- "lrr"
colnames(bray_curtis_obs_only)[which(colnames(bray_curtis_obs_only)=="treat_year")] <- "treatment_year"
bray_curtis_obs_only$metric <- "bray_curtis"

comm_metrics <- read.csv("~/Desktop/community_metrics.csv")
torms <- which(comm_metrics$metric %in% c("appearance", "disappearance","bp_dominance","richness","turnover"))
comm_metrics <- comm_metrics[-torms,]
comm_metrics$site_proj_comm <- with(comm_metrics, paste0(site_code, project_name, community_type))
comm_metrics <- comm_metrics[,c("site_proj_comm", "treatment_year","metric","lrr")]
comm_metrics <- rbind(comm_metrics, bray_curtis_obs_only[,c("site_proj_comm", "treatment_year","metric","lrr")])
comm_metrics <- rbind(comm_metrics, anpp_lrr[,c("site_proj_comm", "treatment_year","metric","lrr")])

comm_metrics$metric <- factor(comm_metrics$metric, 
                              levels=c("ANPP","bray_curtis","invD/S", "RankCor", "loss", "gain"), 
                              labels = c("ANPP","Bray-Curtis","Evenness", "Rank Corr.", "Loss", "Gain"))


sites_to_plot <- intersect(unique(bray_curtis_obs_only$site_proj_comm), unique(anpp_lrr$site_proj_comm))
comm_metrics <- comm_metrics[which(comm_metrics$site_proj_comm %in% sites_to_plot),]



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


mycols <- viridis(5, end=0.7)
plot_list <- list()
count <- 1
for(do_site in unique(gam_df$site_proj_comm)){
  temp_dat <- subset(gam_df, site_proj_comm==do_site)
  minyr <- min(temp_dat$treatment_year)
  maxyr <- max(temp_dat$treatment_year)
  p <- ggplot(temp_dat, aes(x=treatment_year))+
    geom_point(aes(y=lrr, group=site_proj_comm), alpha=0.5, color=mycols[count])+
    geom_line(aes(y=yhat), color=mycols[count])+
    xlab("Treatment year")+
    ylab("")+
    scale_x_continuous(breaks = seq(minyr, maxyr,2))+
    facet_wrap("metric", scales = "free", nrow=1)+
    my_theme
  plot_list[[do_site]] <- p
  count <- count+1
}


g <- arrangeGrob(plot_list[[1]], plot_list[[2]],plot_list[[3]],
                 plot_list[[4]],plot_list[[5]], ncol=1) #generates g
ggsave(file="~/Desktop/all_metrics_timeseries.png", g, width = 10, height=8, units="in", dpi=90) #saves g




####
####  PLOT BIVARIATE RELATIONSHIPS
####
serc_data <- subset(comm_metrics, site_proj_comm=="SEVNfert0")
serc_wide <- dcast(serc_data, treatment_year~metric, value.var = "lrr")
plot(serc_wide$Evenness, serc_wide$`Bray-Curtis`)
abline(lm(`Bray-Curtis`~Evenness, data=serc_wide))
summary(lm(`Bray-Curtis`~Evenness+Gain+Loss+`Rank Corr.`, data=serc_wide))
summary(lm(ANPP~Evenness+Gain+Loss+`Rank Corr.`, data=serc_wide))
summary(lm(ANPP~`Bray-Curtis`, data=serc_wide))



