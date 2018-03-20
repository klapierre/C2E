################################################################################
##  comm_change_concept_fig.R: This script generates fake data for a conceptual
##  figure showing the information about metrics that we are pulling out and
##  summarizing from the time series.
##
##  Author: Andrew Tredennick (atredenn@gmail.com)
##  Date: March 19, 2018
################################################################################

rm(list = ls(all.names = TRUE))



####
####  LOAD LIBRARIES -----------------------------------------------------------
####
library(tidyverse)
library(dplyr)
library(ggthemes)
library(cowplot)



####
####  MAKE A BUNCH OF FAKE DELTA TIME SERIES -----------------------------------
####
n_sites <- 50
ts_means <- c(0.2,0.6,0.5,0.2,0.1,0.1,0.1,0.2,0.1,0.1) # arithmetic scale
ts_sd <- 0.3 # log scale
ts_length <- length(ts_means)
out_ts <- {} # empty object for storage

for(i in 1:ts_length){
  tmpout <- rlnorm(n_sites, log(ts_means[i]), ts_sd)
  tmpout[tmpout>1] <- 1
  tmpdf <- data.frame(site = as.character(1:n_sites),
                      year = i,
                      delta_metric = tmpout)
  out_ts <- rbind(out_ts, tmpdf)
}



####
####  PLOT TRAJECTORIES OF ALL SITES -------------------------------------------
####
all_sites_delta <- ggplot(out_ts, aes(x = year, y = delta_metric, group = site))+
  geom_line(alpha = 0.8)+
  guides(color = FALSE)+
  scale_x_continuous(breaks = seq(1,10,1))+
  xlab("Year of Experiment")+
  ylab(expression(paste(Delta,"Mean Rank Shift")))+
  theme_classic()+
  ggtitle("All sites")

cumulative_change <- out_ts %>%
  group_by(site) %>%
  mutate(cumulative_sum = cumsum(delta_metric))

all_sites_cumul_delta <- ggplot(cumulative_change, aes(x = year, y = cumulative_sum, group = site))+
  geom_line(alpha = 0.8)+
  guides(color = FALSE)+
  scale_x_continuous(breaks = seq(1,10,1))+
  xlab("Year of Experiment")+
  ylab(expression(paste("Cumulative ",Delta,"Mean Rank Shift")))+
  theme_classic()+
  ggtitle("All sites")



####
####  PLOT TRAJECTORY OF SINGLE EXEMPLAR SITE ----------------------------------
####
one_site_delta <- ggplot(filter(out_ts, site == "1"), aes(x = year, y = delta_metric))+
  geom_line()+
  geom_point(size = 2)+
  guides(color = FALSE)+
  scale_x_continuous(breaks = seq(1,10,1))+
  xlab("Year of Experiment")+
  ylab(expression(paste(Delta,"Mean Rank Shift")))+
  theme_classic()+
  ggtitle("One site")

one_site_cumul_delta <- ggplot(filter(cumulative_change, site == "1"), aes(x = year, y = cumulative_sum))+
  geom_line()+
  geom_point(size =  2)+
  guides(color = FALSE)+
  scale_x_continuous(breaks = seq(1,10,1))+
  xlab("Year of Experiment")+
  ylab(expression(paste("Cumulative ",Delta,"Mean Rank Shift")))+
  theme_classic()+
  ggtitle("One site")



####
####  GET YEAR OF MAXIMUM CHANGE -----------------------------------------------
####
max_change <- out_ts %>%
  group_by(site) %>%
  filter(delta_metric == max(delta_metric)) %>%
  ungroup() %>%
  group_by(year) %>%
  summarise(count = n_distinct(site)/n_sites)

ref_df <- data.frame(year = 1:10, count = 0)
rm_ids <- which(ref_df$year %in% max_change$year)

max_change <- max_change %>%
  full_join(ref_df[-rm_ids,])

max_freq <- ggplot(max_change, aes(x = year, y = count))+
  geom_line()+
  geom_point(size = 2)+
  ylab("Proportion of sites at maximum change")+
  scale_x_continuous(breaks = seq(1,10,1))+
  xlab("Year of Experiment")+
  theme_classic()+
  ggtitle("Community change summary")


####
####  COMBINE PLOTS AND SAVE ---------------------------------------------------
####
outplot <- plot_grid(one_site_delta, one_site_cumul_delta,
          all_sites_delta, all_sites_cumul_delta, 
          max_freq,NULL,
          ncol = 2)

ggsave(filename = "comm_change_concept_graphs.pdf", plot = outplot, height = 10, width = 7, units = "in")
