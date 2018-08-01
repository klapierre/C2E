#emily's working directory
setwd("/Users/egrman/Dropbox/C2E/Products/CommunityChange/March2018 WG")

#kevin's working directory
#setwd("C:\\Users\\wilco\\Dropbox\\C2E\\Products\\CommunityChange\\March2018 WG\\")

library(tidyverse)
library(ggplot2)
library(ggthemes)
library(grid)
library(vegan)
library(car)
library(rsq)
library(lme4)

theme_set(theme_bw())
theme_update(axis.title.x=element_text(size=20, vjust=-0.35), axis.text.x=element_text(size=16),
             axis.title.y=element_text(size=20, angle=90, vjust=0.5), axis.text.y=element_text(size=16),
             plot.title = element_text(size=24, vjust=2),
             panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
             legend.title=element_blank(), legend.text=element_text(size=20))

### stealing Kevin's code for creating Glass's delta to compare T vs C at each timestep

### Read in data 
change_metrics_perm <- read.csv("MetricsTrts_July2018.csv") %>%
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
  
# read in site level predictor variables
info.spc=read.csv("SiteExperimentDetails_Dec2016.csv") %>%
  mutate(site_project_comm = paste(site_code, project_name, community_type, sep="_"))

### calculate mean change through time and combine with predictor variables
change_glass_d_mean <- change_glass_d %>%
  group_by(site_project_comm, treatment.x, plot_mani) %>%
  summarise_at(vars(abs_richness_glass, abs_evenness_glass, rank_glass, gains_glass, losses_glass), funs(mean), na.rm=T) %>%
  rename(treatment=treatment.x) %>%
  left_join(info.spc, by=c("site_project_comm")) 

#looking for correlations among predictor variables
pred=as.matrix(change_glass_d_mean[, c("MAP", "MAT", "rrich", "anpp")])
cor(pred)
pairs(pred)
png(paste0("MR predictor variables SITE LEVEL pairs plot.png"), width=11, height=8, units="in", res=600)
print(pairs(pred))
dev.off()



#note that some response var (evenness, losses, gains) have NA (for every year there was no variation among the controls; sd=0 and glass's delta was undefined for every year)



#-------------Multiple regression with only site predictors

#with only some of the treatments (so we don't have too many treatments/experiments at a single site)




#partial R2=(SSE(all other terms in model) - SSE(all terms in model))/SSE(all other terms in model)

rich=lm(abs_richness_glass ~ MAP + MAT + rrich + anpp + n + p + k + CO2 + precip + temp, data=change_glass_d_mean)
vif(rich)
rsq.partial(rich)

rich=lmer(abs_richness_glass ~ MAP + MAT + rrich + anpp + n + p + k + CO2 + precip + temp + (1|site_code), data=change_glass_d_mean)
Anova(rich)
summary(rich)

even=lm(abs_evenness_glass~MAP + MAT + rrich + anpp + n + p + k + CO2 + precip + temp, data=change_glass_d_mean)
Anova(even)
rsq.partial(even)

even=lmer(abs_evenness_glass~MAP + MAT + rrich + anpp + n + p + k + CO2 + precip + temp + (1|site_code), data=change_glass_d_mean)
Anova(even)
summary(even)

rank=lm(rank_glass~MAP + MAT + rrich + anpp + n + p + k + CO2 + precip + temp, data=change_glass_d_mean)
Anova(rank)
rsq.partial(rank)

rank=lmer(rank_glass~MAP + MAT + rrich + anpp + n + p + k + CO2 + precip + temp + (1|site_code), data=change_glass_d_mean)
Anova(rank)
summary(rank)

gains=lm(gains_glass~MAP + MAT + rrich + anpp + n + p + k + CO2 + precip + temp, data=change_glass_d_mean); Anova(gains)
rsq.partial(gains)

gains=lmer(gains_glass~MAP + MAT + rrich + anpp + n + p + k + CO2 + precip + temp + (1|site_code), data=change_glass_d_mean); Anova(gains); summary(gains)

losses=lm(losses_glass~MAP + MAT + rrich + anpp + n + p + k + CO2 + precip + temp, data=change_glass_d_mean); Anova(losses)
rsq.partial(losses)

losses=lmer(losses_glass~MAP + MAT + rrich + anpp + n + p + k + CO2 + precip + temp + (1|site_code), data=change_glass_d_mean); Anova(losses); summary(losses)