# plot_gam_diff_results.R
#  Script to plot the results from estimating the difference between
#  control and treatment cumulative change for each change metric.
#  "Difference" is based on the difference between the two fitted GAMs
#  for each treatment-control pair. The "difference" itself is a smooth
#  curve over the course of the time series, but here we investigate the
#  mean difference over all years, the difference at year 5, and the
#  difference at the end of the time series for each site-experiment.

# Author:
#  Andrew Tredennick (atredenn@gmail.com)


# Clear the workspace -----------------------------------------------------

rm(list = ls(all.names = TRUE))


# Load libraries ----------------------------------------------------------

library(tidyverse)
library(ggthemes)


# Set working directory info ----------------------------------------------

work_dir  <- "~/Repos/C2E/Community Paper/July2018/" # change as needed
data_dir  <- "~/Dropbox/C2E/Products/CommunityChange/March2018 WG/"
results_dir <- "~/Dropbox/C2E/Products/CommunityChange/Summer2018_Results/"
data_file <- "MetricsTrts_July2018.csv"
setwd(work_dir)


# Load and collate results ------------------------------------------------

all_files <- list.files(results_dir)
gam_ids <- grep("gam_comparison_table_", all_files)
gam_files <- all_files[gam_ids]

gam_diffs <- tibble()  # empty storage

for(do_file in gam_files){
  if(length(grep("all_years", do_file)) > 0) diff_type = "all_years"
  if(length(grep("last_year", do_file)) > 0) diff_type = "clast_year"
  if(length(grep("mid_year", do_file)) > 0) diff_type = "bmid_year"
  
  tmp <- read_csv(paste0(results_dir, do_file)) %>%
    filter(sig_diff_cntrl_trt == "yes") %>%
    dplyr::select(site_proj_comm, treatment, response_var, diff, 
                  diff_lower, diff_upper) %>%
    mutate(diff_type = diff_type)
  
  gam_diffs <- bind_rows(gam_diffs, tmp)
}

gam_exp <- gam_diffs %>%
  select(site_proj_comm)%>%
  unique()

# Filter the press treatments for the experiments with enough data
trt_file <- paste0(data_dir, "treatment interactions_July2018.csv")
trt_touse <- read.csv(trt_file)%>%
  select(site_proj_comm, treatment)%>%
  unique()%>%
  right_join(gam_exp)

gam_touse <- gam_diffs %>%
  right_join(trt_touse)

gam_trt <- gam_touse%>%
  select(site_proj_comm, treatment)%>%
  unique()%>%
  mutate(site_project_comm=site_proj_comm)

trts_interactions<-read.csv(paste0(data_dir, "treatment interactions_July2018.csv"))%>%
  right_join(gam_trt)%>%
  filter(use==1)

gam_all_info <- gam_touse %>%
  left_join(trts_interactions)

plot_diffs <- gam_all_info %>%
  mutate(sig_diff = ifelse(sign(diff_lower*diff_upper) == 1, TRUE, FALSE)) %>%
  drop_na()


# Plot the average diffs --------------------------------------------------

# avg_diffs <- gam_diffs %>%
#   group_by(diff_type, response_var) %>%
#   summarise(mean_diff = mean(diff, na.rm = TRUE),
#             sd_diff = sd(diff, na.rm = TRUE))

ggplot(plot_diffs, aes(x = trt_type2, y = diff)) +
  geom_boxplot(aes(fill = sig_diff), outlier.size = 0.5) +
  geom_hline(aes(yintercept = 0), color = "dodgerblue4")+
  facet_grid(response_var~diff_type) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
