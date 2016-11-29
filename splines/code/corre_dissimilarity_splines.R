##  corre_dissimilarity_splines.R: script to fit GAMs to Cedar Creek and 
##  irrigation plot data from the CORRE dataset to test for species
##  reordering and periods of high rate of change.
##
##  Author: Andrew Tredennick (atredenn@gmail.com)
##  Date created: November 29, 2016


rm(list=ls(all.names = TRUE))

#########
## USER!! Reset this data path to match your machine...
########
data_path <- "/Users/atredenn/Google Drive/C2E/splines/data/"


####
####  LIBRARIES ----------------------------------------------------------------
####
library(mgcv)
library(vegan)
library(ggplot2)
library(ggthemes)
library(reshape2)
library(plyr)



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


#########
######### CEDAR CREEK
#########

####
####  LOAD AND SUBSET DATA -----------------------------------------------------
####
all_data    <- read.csv(paste0(data_path,"CORRE_irrigationL_e001D_subset.csv"))
cdr_data    <- subset(all_data, site_code=="CDR")
comm_matrix <- dcast(cdr_data, formula = treatment_year+plot_id+treatment~genus_species,
                     value.var = "abundance", fill=0)



####
####  CALCULATE BRAY-CURTIS ----------------------------------------------------
####
id_cols <- which(colnames(comm_matrix)%in%c("treatment_year","plot_id","treatment"))
out_df_abund  <- list()
out_df_pa <- list()
for(do_year in unique(comm_matrix$treatment_year)){
  ### Abundance
  year_data <- subset(comm_matrix, treatment_year==do_year)
  bc_all    <- vegdist(year_data[,-id_cols], method="bray")
  disper    <- betadisper(bc_all, year_data$treatment, type="centroid")
  distances <- as.matrix(vegdist(disper$centroids, method = "euclidean"))
  tmp_df    <- data.frame(treat_year = rep(do_year, nrow(distances)),
                          bc_dist    = distances[,1],
                          treat_name = rownames(distances))
  out_df_abund    <- rbind(out_df_abund, tmp_df)
  
  ### Presence/Absence (Jaccard)
  bc_all_pa    <- vegdist(year_data[,-id_cols], method="bray", binary = TRUE)
  disper_pa    <- betadisper(bc_all_pa, year_data$treatment, type="centroid")
  distances_pa <- as.matrix(vegdist(disper_pa$centroids, method = "euclidean"))
  tmp_df_pa    <- data.frame(treat_year = rep(do_year, nrow(distances)),
                          bc_dist    = distances_pa[,1],
                          treat_name = rownames(distances_pa))
  out_df_pa    <- rbind(out_df_pa, tmp_df_pa)
}



####
####  FIT THE GAMs -------------------------------------------------------------
####

### ABUNDANCE BRAY-CURTIS
gam_df <- list()
out_df_abund <- subset(out_df_abund, treat_name != "1_y_n")
for(do_treat in unique(out_df_abund$treat_name)){
  todo_data <- subset(out_df_abund, treat_name==do_treat)
  nyrs <- length(unique(todo_data$treat_year))
  gam_fit  <- gam(bc_dist ~ s(treat_year, k=nyrs-1), method="REML", 
                  data = todo_data, 
                  family = betar(link="logit"),
                  select = TRUE) 
  
  ##  Make predictions from best model 
  ##  (here we can just show the all separate model for demonstration)
  years_to_predict   <- 1:max(unique(todo_data$treat_year))
  pred_df            <- data.frame(treat_year = years_to_predict)
  pred_df$yhat       <- predict(object = gam_fit, 
                                newdata = pred_df, 
                                type = "response")
  pred_df            <- merge(pred_df, todo_data, all.x = TRUE)
  pred_df$treat_name <- do_treat
  gam_df             <- rbind(gam_df, pred_df)
}

### PRESENCE/ABSENCE JACCARD
gam_df_pa <- list()
out_df_pa <- subset(out_df_pa, treat_name != "1_y_n")
for(do_treat in unique(out_df_pa$treat_name)){
  todo_data <- subset(out_df_pa, treat_name==do_treat)
  nyrs <- length(unique(todo_data$treat_year))
  gam_fit  <- gam(bc_dist ~ s(treat_year, k=nyrs-1), method="REML", 
                  data = todo_data, 
                  family = betar(link="logit"),
                  select = TRUE) 
  
  ##  Make predictions from best model 
  ##  (here we can just show the all separate model for demonstration)
  years_to_predict   <- 1:max(unique(todo_data$treat_year))
  pred_df            <- data.frame(treat_year = years_to_predict)
  pred_df$yhat       <- predict(object = gam_fit, 
                                newdata = pred_df, 
                                type = "response")
  pred_df            <- merge(pred_df, todo_data, all.x = TRUE)
  pred_df$treat_name <- do_treat
  gam_df_pa         <- rbind(gam_df_pa, pred_df)
}

gam_df$type="Abundance"
gam_df_pa$type="Presence/Absence"
gam_df_all <- rbind(gam_df, gam_df_pa)


## Plot the fits by treatment
ggplot(data=subset(gam_df_all, treat_name!="9_y_n"), aes(x=treat_year, color=treat_name))+
  geom_point(aes(y=bc_dist), alpha=0.5)+
  geom_line(aes(y=yhat))+
  scale_color_viridis(discrete = TRUE, end=0.85)+
  xlab("Treatment year")+
  ylab("Bray-Curtis Dissimilarity")+
  facet_grid(type~treat_name)+
  guides(color=FALSE)+
  my_theme
ggsave(filename = "../figures/bray_curtis_gams_cdr.png", width = 8.5, height = 3, units="in", dpi = 120)





#########
######### IRRIGATION PLOTS
#########

####
####  LOAD AND SUBSET DATA -----------------------------------------------------
####
all_data    <- read.csv(paste0(data_path,"CORRE_irrigationL_e001D_subset.csv"))
irrig_data    <- subset(all_data, site_code=="KNZ")
comm_matrix <- dcast(irrig_data, formula = treatment_year+plot_id+treatment~genus_species,
                     value.var = "abundance", fill=0)



####
####  CALCULATE BRAY-CURTIS ----------------------------------------------------
####
id_cols <- which(colnames(comm_matrix)%in%c("treatment_year","plot_id","treatment"))
out_df_abund  <- list()
out_df_pa <- list()
for(do_year in unique(comm_matrix$treatment_year)){
  ### Abundance
  year_data <- subset(comm_matrix, treatment_year==do_year)
  bc_all    <- vegdist(year_data[,-id_cols], method="bray")
  disper    <- betadisper(bc_all, year_data$treatment, type="centroid")
  distances <- as.matrix(vegdist(disper$centroids, method = "euclidean"))
  tmp_df    <- data.frame(treat_year = rep(do_year, nrow(distances)),
                          bc_dist    = distances[,1],
                          treat_name = rownames(distances))
  out_df_abund    <- rbind(out_df_abund, tmp_df)
  
  ### Presence/Absence (Jaccard)
  bc_all_pa    <- vegdist(year_data[,-id_cols], method="bray", binary = TRUE)
  disper_pa    <- betadisper(bc_all_pa, year_data$treatment, type="centroid")
  distances_pa <- as.matrix(vegdist(disper_pa$centroids, method = "euclidean"))
  tmp_df_pa    <- data.frame(treat_year = rep(do_year, nrow(distances)),
                             bc_dist    = distances_pa[,1],
                             treat_name = rownames(distances_pa))
  out_df_pa    <- rbind(out_df_pa, tmp_df_pa)
}



####
####  FIT THE GAMs -------------------------------------------------------------
####

### ABUNDANCE BRAY-CURTIS
out_df_abund <- subset(out_df_abund, treat_name != "c")
todo_data <- out_df_abund
nyrs <- length(unique(todo_data$treat_year))
gam_fit  <- gam(bc_dist ~ s(treat_year, k=nyrs-1), method="REML", 
                data = todo_data, 
                family = betar(link="logit"),
                select = TRUE) 
coefs_abund <- coef(gam_fit)

##  Make predictions from best model 
##  (here we can just show the all separate model for demonstration)
years_to_predict   <- 1:max(unique(todo_data$treat_year))
pred_df            <- data.frame(treat_year = years_to_predict)
pred_df$yhat       <- predict(object = gam_fit, 
                              newdata = pred_df, 
                              type = "response")
pred_df            <- merge(pred_df, todo_data, all.x = TRUE)
pred_df$treat_name <- do_treat
gam_df             <- pred_df


### PRESENCE/ABSENCE JACCARD
out_df_pa <- subset(out_df_pa, treat_name != "c")
todo_data <-out_df_pa
nyrs <- length(unique(todo_data$treat_year))
gam_fit  <- gam(bc_dist ~ s(treat_year, k=nyrs-1), method="REML", 
                data = todo_data, 
                family = betar(link="logit"),
                select = TRUE) 

##  Make predictions from best model 
##  (here we can just show the all separate model for demonstration)
years_to_predict   <- 1:max(unique(todo_data$treat_year))
pred_df            <- data.frame(treat_year = years_to_predict)
pred_df$yhat       <- predict(object = gam_fit, 
                              newdata = pred_df, 
                              type = "response")
pred_df            <- merge(pred_df, todo_data, all.x = TRUE)
pred_df$treat_name <- do_treat
gam_df_pa         <- pred_df


gam_df$type="Abundance"
gam_df_pa$type="Presence/Absence"
gam_df_all <- rbind(gam_df, gam_df_pa)


## Plot the fits by treatment
ggplot(data=subset(gam_df_all), aes(x=treat_year, color=treat_name))+
  geom_point(aes(y=bc_dist), alpha=0.5)+
  geom_line(aes(y=yhat))+
  scale_color_viridis(discrete = TRUE, end=0.85)+
  xlab("Treatment year")+
  ylab("Bray-Curtis Dissimilarity")+
  facet_wrap("type", ncol=1)+
  guides(color=FALSE)+
  my_theme
ggsave(filename = "../figures/bray_curtis_gams_irrig.png", width = 2.5, height = 4, units="in", dpi = 120)



