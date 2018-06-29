#emily's working directory
setwd("/Users/egrman/Dropbox/C2E/Products/CommunityChange/March2018 WG")

#kevin's working directory
#setwd("C:\\Users\\wilco\\Dropbox\\C2E\\Products\\CommunityChange\\March2018 WG\\")

library(tidyverse)
library(ggplot2)
library(ggthemes)
library(grid)
library(vegan)
library(pls)

theme_set(theme_bw())
theme_update(axis.title.x=element_text(size=20, vjust=-0.35), axis.text.x=element_text(size=16),
             axis.title.y=element_text(size=20, angle=90, vjust=0.5), axis.text.y=element_text(size=16),
             plot.title = element_text(size=24, vjust=2),
             panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
             legend.title=element_blank(), legend.text=element_text(size=20))

### stealing Kevin's code for creating Glass's delta to compare T vs C at each timestep

### Read in data 
change_metrics_perm <- read.csv("CORRE_RACS_Subset_Perm.csv") %>%
  mutate(abs_richness_change = abs(richness_change),
         abs_evenness_change = abs(evenness_change))

### Control data
change_control <- change_metrics_perm %>%
  filter(plot_mani==0) %>%
  dplyr::select(treatment_year, treatment_year2, abs_richness_change, abs_evenness_change, 
                rank_change, gains, losses, site_project_comm, treatment, plot_mani) %>%
  rename(abs_richness_change_ctrl = abs_richness_change,
         abs_evenness_change_ctrl = abs_evenness_change,
         rank_change_ctrl = rank_change,
         gains_ctrl = gains,
         losses_ctrl = losses
  ) %>%
  group_by(site_project_comm, treatment, treatment_year2) %>%
  summarise_at(vars(abs_richness_change_ctrl:losses_ctrl), funs(mean, sd), na.rm=T)

change_glass_d <- change_metrics_perm %>%
  filter(plot_mani != 0) %>%
  group_by(site_project_comm, treatment, treatment_year2, plot_mani) %>%
  summarise(abs_richness_change = mean(abs_richness_change,na.rm=T),
            abs_evenness_change = mean(abs_evenness_change, na.rm=T),
            rank_change = mean(rank_change, na.rm=T),
            gains = mean(gains, na.rm=T),
            losses = mean(losses, na.rm=T)) %>%
  left_join(change_control, by=c("site_project_comm","treatment_year2")) %>%
  mutate(abs_richness_glass = (abs_richness_change-abs_richness_change_ctrl_mean)/abs_richness_change_ctrl_sd,
         abs_evenness_glass = (abs_evenness_change-abs_evenness_change_ctrl_mean)/abs_evenness_change_ctrl_sd,
         rank_glass = (rank_change-rank_change_ctrl_mean)/rank_change_ctrl_sd,
         gains_glass = (gains-gains_ctrl_mean)/gains_ctrl_sd,
         losses_glass = (losses-losses_ctrl_mean)/losses_ctrl_sd
  ) %>%
  dplyr::select(site_project_comm:plot_mani, abs_richness_glass:losses_glass) %>%
  ungroup()

#change_glass_d is the thing that we want

## replace Inf with NAs in change_glass_d
change_glass_d <- change_glass_d %>%
  mutate(gains_glass=replace(gains_glass, gains_glass=="Inf", NA)) %>%
  mutate(losses_glass=replace(losses_glass, losses_glass=="Inf", NA))
  

# reading in predictor variables for PLSR
info.spc=read.csv("SiteExperimentDetails_Dec2016.csv") %>%
  mutate(site_project_comm = paste(site_code, project_name, community_type, sep="_"))

info.trt=read.csv("ExperimentInformation_Nov2017.csv") %>%
  mutate(site_project_comm = paste(site_code, project_name, community_type, sep="_")) %>%
  group_by(site_project_comm, treatment) %>%
  summarise_at(vars(n, p, k, CO2, precip, temp), funs(mean))

### calculate mean change through time and combine with predictor variables
change_glass_d_mean <- change_glass_d %>%
  group_by(site_project_comm, treatment.x, plot_mani) %>%
  summarise_at(vars(abs_richness_glass, abs_evenness_glass, rank_glass, gains_glass, losses_glass), funs(mean), na.rm=T) %>%
  rename(treatment=treatment.x) %>%
  left_join(info.spc, by=c("site_project_comm")) %>%
  left_join(info.trt, by=c("site_project_comm","treatment"))

#subsetting out predictor variables
pred=as.matrix(change_glass_d_mean[, c("MAP", "MAT", "rrich", "anpp", "n", "p", "k", "CO2", "precip", "temp")])
#note that some response var (evenness, losses, gains) have NA (for every year there was no variation among the controls; sd=0 and glass's delta was undefined for every year)


#subsetting out predictor variables
pred=as.matrix(change_glass_d_mean[, c("MAP", "MAT", "rrich", "anpp", "n", "p", "k", "CO2", "precip", "temp")])

#attempting PLSR
rich=plsr(abs_richness_glass~pred, data=change_glass_d_mean, validation="LOO")
plot(RMSEP(rich)) #none of the components have any predictive value
summary(rich)

even=plsr(abs_evenness_glass~pred, data=change_glass_d_mean, validation="LOO")
plot(RMSEP(even)) #none of the components have any predictive value
summary(even)

rank=plsr(rank_glass~pred, data=change_glass_d_mean, validation="LOO")
plot(RMSEP(rank)) #none of the components have any predictive value
summary(rank)

gains=plsr(gains_glass~pred, data=change_glass_d_mean, validation="LOO")
plot(RMSEP(gains))
summary(gains) #first component explains 5% of the variation in gains
plot(gains, "loadings", comps = 1:2, legendpos = "topright")

losses=plsr(losses_glass~pred, data=change_glass_d_mean, validation="LOO")
plot(RMSEP(losses))
summary(losses) #first component explains 8% of variation in losses
plot(losses, ncomp=2, asp=1, line=T) #terrible!
plot(losses, "loadings", comps = 1:2, legendpos = "topright")

#-----Multiple Regression
#we have 125 datapoints and 10 predictors--why can't we just do multiple regression?
library(car)

cor(pred)

#partial R2=(SSE(all other terms in model) - SSE(all terms in model))/SSE(all other terms in model)

rich=lm(abs_richness_glass ~ MAP + MAT + rrich + anpp + n + p + k + CO2 + precip + temp, data=change_glass_d_mean)
vif(rich)
Anova(rich)

even=lm(abs_evenness_glass~MAP + MAT + rrich + anpp + n + p + k + CO2 + precip + temp, data=change_glass_d_mean)
Anova(even)

rank=lm(rank_glass~MAP + MAT + rrich + anpp + n + p + k + CO2 + precip + temp, data=change_glass_d_mean)
Anova(rank)

gains=lm(gains_glass~MAP + MAT + rrich + anpp + n + p + k + CO2 + precip + temp, data=change_glass_d_mean); Anova(gains)
gains.noMAP=lm(gains_glass~MAT + rrich + anpp + n + p + k + CO2 + precip + temp, data=change_glass_d_mean); Anova(gains.noMAP)
(40.752-36.936)/40.752 #R2=0.09, which is higher than predictive power of first compnent in PLSR

losses=lm(losses_glass~MAP + MAT + rrich + anpp + n + p + k + CO2 + precip + temp, data=change_glass_d_mean); Anova(losses)
losses.noN=lm(losses_glass~MAP + MAT + rrich + anpp + p + k + CO2 + precip + temp, data=change_glass_d_mean); Anova(losses.noN)
(40.190-36.936)/40.190 #R2=0.08, which matches predictive power of first component in PLSR
