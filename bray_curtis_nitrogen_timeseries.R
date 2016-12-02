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
source("subsetting 7 yr N data.R")
setwd("~/Repos/C2E/")
relAbundN$site_proj_comm <- with(relAbundN, paste0(site_code, project_name, community_type))


####
####  CALCULATE BRAY-CURTIS FOR TREATMENTS Vs. CONTROLS
####
all_sites_bc_time <- list()
for(do_unit in unique(relAbundN$site_proj_comm)){
  tmp_dat <- subset(relAbundN, site_proj_comm==do_unit)
  comm_matrix <- dcast(tmp_dat, formula = treatment_year+plot_id+n_treatment~genus_species,
                       value.var = "relcov", fill=0)
  id_cols <- which(colnames(comm_matrix)%in%c("treatment_year","plot_id","n_treatment"))
  out_bc_abund  <- list()
  for(do_year in unique(comm_matrix$treatment_year)){
    year_data <- subset(comm_matrix, treatment_year==do_year)
    bc_all    <- vegdist(year_data[,-id_cols], method="bray")
    disper    <- betadisper(bc_all, year_data$n_treatment, type="centroid")
    distances <- as.matrix(vegdist(disper$centroids, method = "euclidean"))
    tmp_df    <- data.frame(treat_year = rep(do_year, nrow(distances)),
                            bc_dist    = distances[,1],
                            treat_name = rownames(distances))
    out_bc_abund    <- rbind(out_bc_abund, tmp_df)
  }
  out_bc_abund <- subset(out_bc_abund, treat_name==1)
  out_bc_abund$site_proj_comm <- do_unit
  all_sites_bc_time <- rbind(all_sites_bc_time, out_bc_abund)
}



####
####  FIT GAMs FOR EACH SITE-PROJECT-COMMUNITY
####
gam_df <- list()
for(do_unit in unique(all_sites_bc_time$site_proj_comm)){
  todo_data <- subset(all_sites_bc_time, site_proj_comm==do_unit)
  nyrs      <- length(unique(todo_data$treat_year))
  gam_fit   <- gam(bc_dist ~ s(treat_year, k=nyrs-1), method="REML", 
                   data = todo_data, 
                   family = betar(link="logit"),
                   select = TRUE) 
  years_to_predict   <- seq(min(unique(todo_data$treat_year)),max(unique(todo_data$treat_year)), by=0.1)
  pred_df            <- data.frame(treat_year = years_to_predict)
  pred_df$yhat       <- predict(object = gam_fit, 
                                newdata = pred_df, 
                                type = "response")
  pred_df            <- merge(pred_df, todo_data, all.x = TRUE)
  pred_df$site_proj_comm <- do_unit
  gam_df             <- rbind(gam_df, pred_df)
}

saveRDS(gam_df, "bray_curtis_gams.RDS")
## Plot the fits by site_proj_com
ggplot(data=gam_df, aes(x=treat_year))+
  geom_point(aes(y=bc_dist), alpha=0.5)+
  geom_line(aes(y=yhat))+
  xlab("Treatment year")+
  ylab("Bray-Curtis Dissimilarity")+
  facet_wrap("site_proj_comm", scales = "free", nrow=1)+
  guides(color=FALSE)+
  my_theme
ggsave(filename = "~/Desktop/bray_curtis_ts_gams.png", width = 10, height = 2, units="in", dpi = 120)



