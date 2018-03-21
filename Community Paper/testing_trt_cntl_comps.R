################################################################################
##  testing_trt_cntl_comps.R: Exploratory script looking at different ways to
##  compare treatment and control RAC metrics. This script only works with
##  RAC metrics from communities identified to be different from controls as
##  per La Pierre et al. 201x (aka, the Bayes approach).
##
##  Author: Andrew Tredennick (atredenn@gmail.com)
##  Date: March 21, 2018
################################################################################

##  Clear the workspace
rm(list = ls(all.names = TRUE))



####
####  LOAD LIBRARIES -----------------------------------------------------------
####
library(tidyverse)
library(ggthemes)
library(cowplot)



####
####  SET WORKING DIRECTORIES AND FILENAMES ------------------------------------
####
work_dir  <- "~/Repos/C2E/Community Paper/" # change as needed
data_dir  <- "~/Dropbox/C2E/Products/CommunityChange/March2018 WG/"
data_file <- "CORRE_RACS_Subset_Bayes.csv" # change as needed
setwd(work_dir)



####
####  READ IN DATA AND REFORMAT ------------------------------------------------
####
metric_dat <- read.csv(paste0(data_dir,data_file))

mean_metrics <- metric_dat %>%
  mutate_at(vars(richness_change, evenness_change, rank_change, gains, losses), abs) %>%
  group_by(treatment_year, treatment_year2, site_project_comm, treatment, plot_mani) %>%
  summarize_at(vars(richness_change, evenness_change, rank_change, gains, losses), funs(mean), na.rm = T) %>%
  mutate(trt=ifelse(plot_mani==0, "control","treatment"))

sd_control_metrics <- metric_dat %>%
  mutate_at(vars(richness_change, evenness_change, rank_change, gains, losses), abs) %>%
  group_by(treatment_year, treatment_year2, site_project_comm, treatment, plot_mani) %>%
  summarize_at(vars(richness_change, evenness_change, rank_change, gains, losses), funs(sd), na.rm = T) %>%
  mutate(trt=ifelse(plot_mani==0, "control","treatment")) %>%
  filter(trt == "control") %>%
  ungroup() %>%
  mutate(sd_richness = richness_change,
         sd_even = evenness_change,
         sd_gains = gains,
         sd_losses = losses,
         sd_ranks = rank_change) %>%
  select(site_project_comm, treatment_year, treatment_year2,
         sd_richness, sd_even, sd_gains, sd_losses, sd_ranks)

control_means <- mean_metrics %>%
  filter(plot_mani==0) %>%
  ungroup() %>%
  mutate(C_S=richness_change, 
         C_even=evenness_change, 
         C_gain=gains, 
         C_loss=losses, 
         C_rank=rank_change) %>%
  select(site_project_comm, treatment_year, treatment_year2,
         C_S, C_even, C_gain, C_loss, C_rank)



####
####  CALCULATE RESPONSE RATIOS ------------------------------------------------
####
response_ratios <- merge(mean_metrics, control_means, 
                         by=c("site_project_comm",
                              "treatment_year", 
                              "treatment_year2"))%>%
  left_join(sd_control_metrics, by = c("site_project_comm",
                                       "treatment_year", 
                                       "treatment_year2")) %>%
  ungroup() %>%
  filter(plot_mani!=0) %>%
  mutate(Richness = ((richness_change-C_S)/sd_richness),
         Evenness = ((evenness_change-C_even)/sd_even),
         `Rank Change` = ((rank_change-C_rank)/sd_ranks),
         Losses = ((losses-C_loss)/sd_losses),
         Gains = ((gains-C_gain)/sd_gains))%>%
  select(site_project_comm, treatment_year, treatment_year2, treatment,
         Richness, Evenness, `Rank Change`, Losses, Gains) %>%
  gather(key = metric, value = ratio, Richness:Gains) %>%
  mutate(id = paste(site_project_comm, treatment, sep = "::"))

ggplot(response_ratios, aes(x = treatment_year, y = ratio))+
  geom_point(aes(group = id))+
  geom_smooth(se = FALSE, method = "loess")+
  guides(color = FALSE)+
  facet_wrap(~metric, scales = "free")+
  ylab(expression(paste("Glass's ", Delta)))+
  xlab("Treatment year")+
  theme_few()
# ggsave(filename = paste0(data_dir,"figures/glass_by_metric.pdf"), height = 5, width = 8.5, units = "in")

ggplot(response_ratios, aes(x = treatment_year, y = ratio, color = metric))+
  # geom_point(aes(group = metric))+
  geom_smooth(se = FALSE, method = "loess")+
  facet_wrap(~id, scales = "free") +
  theme_few()
# ggsave(filename = paste0(data_dir,"figures/glass_loess_large.pdf"), height = 15, width = 15, units = "in")




####
####  GENERATE FAKE DATA TO GAIN INTUITION -------------------------------------
####
n <- 1000
sd_c <- 1
sd_t <- 1
mu_c <- 10
mu_t <- 14
treatments <- rnorm(n,mu_t,sd_t)
controls <- rnorm(n,mu_c,sd_c)
glass_delta <- (mean(treatments) - mean(controls)) / sd_c

plot(density(controls), xlim = c(8,20), main = paste("Glass's = ",round(glass_delta,2)))
lines(density(treatments), lty=2)

