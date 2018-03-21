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
  mutate_at(vars(richness_change, richness_change_abs, 
                 evenness_change,evenness_change_abs, 
                 rank_change, gains, losses), 
            funs(cumsum) ) %>%
  mutate(control = ifelse(plot_mani==0,"control","treatment"))



####
####  TEST FIT USING JUST PPLOTS -----------------------------------------------
####
pplots_cumsum <- filter(change_cumsum, site_project_comm == "KNZ_pplots_0")

pplot_controls <- filter(pplots_cumsum, plot_mani == 0 | plot_mani == 1)

gamtest <- gam(rank_change ~ s(treatment_year, treatment, bs = "fs") + s(plot_id, bs="re"), data = pplot_controls)
gamsim <- gam(rank_change ~ s(treatment_year) + s(plot_id, bs="re"), data = pplot_controls)
AIC(gamtest, gamsim)
summary(gamtest)
plot(pplot_controls$treatment_year,predict(gamtest, type = "response", exclude = "s(plot_id)"), type = "l")
points(pplot_controls$treatment_year, pplot_controls$rank_change)

ggplot(pplot_controls, aes(x = treatment_year, y = rank_change))+
  geom_point()



