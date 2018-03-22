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
####  TEST FIT USING JUST PPLOTS -----------------------------------------------
####
##  Subset controls and one treatment for pplots
pplots_cumsum  <- filter(change_cumsum, site_project_comm == "KNZ_pplots_0")
pplot_controls <- filter(pplots_cumsum, plot_mani == 0 | plot_mani == 1)

##  Fit nested GAMs
### TODO: CHECK THESE MODEL STRUCTURES!!! Email Gavin?
gam_test <- gam(rank_change ~ s(treatment_year, treatment, bs = "fs") + 
                              s(plot_id, bs="re"), data = pplot_controls)
gam_null <- gam(rank_change ~ s(treatment_year) + s(plot_id, bs="re"),
                data = pplot_controls)

##  Compare nested models
AIC(gam_test, gam_null) # less conservative
BIC(gam_test, gam_null) # more conservative
deltaAIC <- diff(AIC(gam_test, gam_null)$AIC)

##  Make example plot of predictions from the complex model
newdata <- data.frame(treatment_year = rep(0:11,2),
                      treatment = rep(c("N1P0","N2P0"), each = 12),
                      plot_id = rep(factor(27),times = 12*2))
gam_preds <- predict(gam_test, 
                     newdata = newdata, 
                     type = "response", 
                     exclude = "s(plot_id)")
plot_preds <- data.frame(predictions = gam_preds,
                        treatment_year = rep(0:11,2),
                        treatment = rep(c("N1P0","N2P0"), each = 12))
plot_data <- pplot_controls %>%
  ungroup() %>%
  select(treatment_year, rank_change, treatment)

ggplot()+
  geom_point(data = plot_data, 
             aes(x = treatment_year, y = rank_change, color = treatment))+
  geom_line(data = plot_preds, 
            aes(x = treatment_year, y = predictions, color = treatment),
            size = 1.2)+
  scale_color_brewer(type = "qual", name = "Treatment")+
  xlab("Treatment year")+
  ylab("Rank change")+
  ggtitle("KNZ PPlots Example", subtitle = "Control (N1P0) vs. N2P0")
ggsave(filename = paste0(data_dir,"figures/pplot_gam_example.pdf"),
       height = 4,
       width = 5,
       units = "in")


