################################################################################
##  fit_cumulative_gams.R: Script that fits GAMs to the cumulative RAC metrics
##  data. The goal is to compare models with and without a treatment effect.
##
##  Author: Andrew Tredennick (atredenn@gmail.com)
##  Date created: March 21, 2018
################################################################################

##  Clear the workspace
rm(list = ls(all.names = TRUE))



####
####  LOAD LIBRARIES -----------------------------------------------------------
####
library(tidyverse)
library(ggthemes)
library(mgcv)



####
####  SET WORKING DIRECTORIES AND FILENAMES ------------------------------------
####
work_dir  <- "~/Repos/C2E/Community Paper/" # change as needed
data_dir  <- "~/Dropbox/C2E/Products/CommunityChange/March2018 WG/"
data_file <- "CORRE_RACS_Subset_Perm.csv" # change as needed
setwd(work_dir)



####
####  READ IN DATA AND CALCULATE CUMULATIVE CHANGE -----------------------------
####
change_metrics <- read.csv(paste0(data_dir,data_file))

##  Calculate cumulative sums of each metric (from Kevin)
change_cumsum <- change_metrics %>%
  group_by(site_project_comm, treatment, plot_id) %>%
  mutate(richness_change_abs = abs(richness_change)) %>%
  mutate(evenness_change_abs = abs(evenness_change)) %>%
  mutate_at(vars(richness_change, 
                 richness_change_abs, 
                 evenness_change,
                 evenness_change_abs, 
                 rank_change, 
                 gains, 
                 losses), 
            funs(cumsum)) %>%
  mutate(control = ifelse(plot_mani==0,"control","treatment"))



####
####  LOOP OVER SITE_PROJECT_COMMS AND COMPARE CONTROLS VS. TREATMENTS ---------
####
all_sites <- unique(change_cumsum$site_project_comm)
all_delta_aics <- {} # empty object for storage

for(do_site in all_sites){
  site_data <- filter(change_cumsum, site_project_comm == do_site)
  site_controls <- filter(site_data, plot_mani == 0)
  site_treatments <- filter(site_data, plot_mani != 0)
  all_treatments <- unique(site_treatments$treatment)
  
  for(do_treatment in all_treatments){
    treatment_data <- filter(site_treatments, treatment == do_treatment)
    model_data <- rbind(site_controls, treatment_data)
    num_years <- length(unique(model_data$treatment_year))
    
    if(num_years < 4){
      rich_delta_aic <- NA
      even_delta_aic <- NA
      rank_delta_aic <- NA
      gain_delta_aic <- NA
      loss_delta_aic <- NA
    }
    
    if(num_years > 3){
      ##  Richness
      gam_test <- gam(richness_change_abs ~ s(treatment_year, treatment, 
                                              bs = "fs", k = (num_years-1)) + 
                        s(plot_id, bs="re"), 
                      data = model_data)
      gam_null <- gam(richness_change_abs ~ s(treatment_year, k = (num_years-1)) + 
                        s(plot_id, bs="re"),
                      data = model_data)
      rich_aics <- AIC(gam_test, gam_null)$AIC
      rich_delta_aic <- rich_aics[1] - rich_aics[2] # full - null
      rm(gam_test)
      rm(gam_null)
      
      ##  Evenness
      gam_test <- gam(evenness_change_abs ~ s(treatment_year, treatment, 
                                              bs = "fs", k = (num_years-1)) + 
                        s(plot_id, bs="re"), 
                      data = model_data)
      gam_null <- gam(evenness_change_abs ~ s(treatment_year, k = (num_years-1)) + 
                        s(plot_id, bs="re"),
                      data = model_data)
      even_aics <- AIC(gam_test, gam_null)$AIC
      even_delta_aic <- even_aics[1] - even_aics[2] # full - null 
      rm(gam_test)
      rm(gam_null)
      
      ##  Rank change
      gam_test <- gam(rank_change ~ s(treatment_year, treatment, 
                                              bs = "fs", k = (num_years-1)) + 
                        s(plot_id, bs="re"), 
                      data = model_data)
      gam_null <- gam(rank_change ~ s(treatment_year, k = (num_years-1)) + 
                        s(plot_id, bs="re"),
                      data = model_data)
      rank_aics <- AIC(gam_test, gam_null)$AIC
      rank_delta_aic <- rank_aics[1] - rank_aics[2] # full - null 
      rm(gam_test)
      rm(gam_null)
      
      ##  Gains
      gam_test <- gam(gains ~ s(treatment_year, treatment, 
                                      bs = "fs", k = (num_years-1)) + 
                        s(plot_id, bs="re"), 
                      data = model_data)
      gam_null <- gam(gains ~ s(treatment_year, k = (num_years-1)) + 
                        s(plot_id, bs="re"),
                      data = model_data)
      gain_aics <- AIC(gam_test, gam_null)$AIC
      gain_delta_aic <- gain_aics[1] - gain_aics[2] # full - null 
      rm(gam_test)
      rm(gam_null)
      
      ##  Losses
      gam_test <- gam(losses ~ s(treatment_year, treatment, 
                                      bs = "fs", k = (num_years-1)) + 
                        s(plot_id, bs="re"), 
                      data = model_data)
      gam_null <- gam(losses ~ s(treatment_year, k = (num_years-1)) + 
                        s(plot_id, bs="re"),
                      data = model_data)
      loss_aics <- AIC(gam_test, gam_null)$AIC
      loss_delta_aic <- loss_aics[1] - loss_aics[2] # full - null 
      rm(gam_test)
      rm(gam_null)
      
    } # end if/then for number of years
    
    tmp_out <- data.frame(site_project_comm = do_site,
                          treatment = do_treatment,
                          rich_delta_aic = rich_delta_aic,
                          even_delta_aic = even_delta_aic,
                          rank_delta_aic = rank_delta_aic,
                          gain_delta_aic = gain_delta_aic,
                          loss_delta_aic = loss_delta_aic)
    
    all_delta_aics <- rbind(all_delta_aics, tmp_out)
    
    ##  Remove AIC objects, suppressing warnings if they don't exist
    suppressWarnings(rm(rich_aics))
    suppressWarnings(rm(rich_delta_aic))
    suppressWarnings(rm(even_aics))
    suppressWarnings(rm(even_delta_aic))
    suppressWarnings(rm(rank_aics))
    suppressWarnings(rm(rank_delta_aic))
    suppressWarnings(rm(gain_aics))
    suppressWarnings(rm(gain_delta_aic))
    suppressWarnings(rm(loss_aics))
    suppressWarnings(rm(loss_delta_aic))
    
  } # end treatment loop
  
  print(paste("Done with site:", do_site))
  
} # end site loop



####
####  SAVE DELTA_AIC TABLE -----------------------------------------------------
####
write.csv(x = all_delta_aics, 
          file = paste0(data_dir,"gam_delta_aic_table.csv"))



####
####  VISUALIZE THE DLETA AIC TABLE --------------------------------------------
####
delta_aics <- read.csv(paste0(data_dir,"gam_delta_aic_table.csv"), 
                       row.names = 1) %>%
  gather(key = metric, value = delta_aic, rich_delta_aic:loss_delta_aic) %>%
  mutate(site_treatment = paste(site_project_comm, treatment, sep = "::"),
         different = ifelse(delta_aic < -2, "yes", "no"))

ggplot(delta_aics, aes(y = site_treatment, x = metric))+
  geom_tile(aes(fill = different))+
  scale_fill_brewer(type = "qual", 
                    labels = c("C and T not different", "C and T different", "NA"), 
                    name = NULL)+
  scale_x_discrete(labels = c("Evenness","Gains","Losses","Rank Change", "Richness"))+
  xlab("Metric")+
  ylab("Site and Treatment")
ggsave(filename = paste0(data_dir,"figures/delta_aic_figure.pdf"), height = 14, width = 7, units = "in")


####
####  TEST FIT USING JUST PPLOTS -----------------------------------------------
####
##  Subset controls and one treatment for pplots
# do_site <- "CDR_e001_D"
# do_treatment <- 6
# pplots_cumsum  <- filter(change_cumsum, site_project_comm == do_site)
# pplot_controls <- filter(pplots_cumsum, treatment == 9 | treatment == do_treatment)
# 
# ##  Fit nested GAMs
# ### TODO: CHECK THESE MODEL STRUCTURES!!! Email Gavin?
# gam_test <- gam(evenness_change_abs ~ s(treatment_year, treatment, bs = "fs") + 
#                   s(plot_id, bs="re"), data = pplot_controls)
# gam_null <- gam(evenness_change_abs ~ s(treatment_year) + s(plot_id, bs="re"),
#                 data = pplot_controls)
# 
# ##  Compare nested models
# AIC(gam_test, gam_null) # less conservative
# BIC(gam_test, gam_null) # more conservative
# deltaAIC <- diff(AIC(gam_test, gam_null)$AIC)
# 
# ##  Make example plot of predictions from the complex model
# newdata <- data.frame(treatment_year = rep(0:11,2),
#                       treatment = rep(c("N1P0","N2P0"), each = 12),
#                       plot_id = rep(factor(27),times = 12*2))
# gam_preds <- predict(gam_test, 
#                      newdata = newdata, 
#                      type = "response", 
#                      exclude = "s(plot_id)")
# plot_preds <- data.frame(predictions = gam_preds,
#                          treatment_year = rep(0:11,2),
#                          treatment = rep(c("N1P0","N2P0"), each = 12))
# plot_data <- pplot_controls %>%
#   ungroup() %>%
#   select(treatment_year, rank_change, treatment)
# 
# ggplot()+
#   geom_point(data = plot_data, 
#              aes(x = treatment_year, y = rank_change, color = treatment))+
#   geom_line(data = plot_preds, 
#             aes(x = treatment_year, y = predictions, color = treatment),
#             size = 1.2)+
#   scale_color_brewer(type = "qual", name = "Treatment")+
#   xlab("Treatment year")+
#   ylab("Rank change")+
#   ggtitle("KNZ PPlots Example", subtitle = "Control (N1P0) vs. N2P0")
# ggsave(filename = paste0(data_dir,"figures/pplot_gam_example.pdf"),
#        height = 4,
#        width = 5,
#        units = "in")
# 
