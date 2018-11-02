#emily's working directory
setwd("/Users/egrman/Dropbox/C2E/Products/CommunityChange/March2018 WG")

#meghan's working directory
setwd("/Users/megha/Dropbox/C2E/Products/CommunityChange/March2018 WG")
setwd("~/Dropbox/C2E/Products/CommunityChange/March2018 WG")

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
change_metrics <- read.csv("MetricsTrts_July2018.csv") %>%
  mutate(abs_richness_change = abs(richness_change),
         abs_evenness_change = abs(evenness_change))

### Control data
change_control <- change_metrics %>%
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

change_glass_d <- change_metrics %>%
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
  mutate(abs_richness_glass=replace(abs_richness_glass, abs_richness_glass=="Inf", NA)) %>%
  mutate(rank_glass=replace(rank_glass, rank_glass=="Inf", NA)) %>%
  mutate(gains_glass=replace(gains_glass, gains_glass=="Inf", NA)) %>%
  mutate(losses_glass=replace(losses_glass, losses_glass=="Inf", NA))
  
# read in site level predictor variables
info.spc=read.csv("SiteExperimentDetails_Dec2016.csv") %>%
  mutate(site_project_comm = paste(site_code, project_name, community_type, sep="_"))

# read in treatment variables for subsetting later
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

#looking for correlations among predictor variables
pred=as.matrix(change_glass_d_mean[, c("MAP", "MAT", "rrich", "anpp")])
cor(pred)
pairs(pred)
png(paste0("MR predictor variables SITE LEVEL pairs plot.png"), width=11, height=8, units="in", res=600)
print(pairs(pred))
dev.off()

#-------------Multiple regression with only site predictors

#with only some of the treatments (so we don't have too many treatments/experiments at a single site)
#meghan picked out which ones we want (resource manipulations only?)

##for this first analysis I think we should use all the data - meghan, thus I am skipping this step and running on change_glass_d_mean
usethese=change_metrics[change_metrics$use==1, c("site_project_comm", "treatment", "use")]
use_change_glass_d_mean=merge(change_glass_d_mean, unique(usethese), by=c("site_project_comm", "treatment"))
#write.csv(use_change_glass_d_mean, "use for site predictors of SERGL.csv")

#note that some response var (evenness, losses, gains) have NA (for every year there was no variation among the controls; sd=0 and glass's delta was undefined for every year)

#-----1a) treating all experiments in this subset as independent data points:

rich=lm(abs_richness_glass ~ MAP + MAT + rrich + anpp, data=change_glass_d_mean)
vif(rich)
summary(rich)
rsq.partial(rich)
ggplot(data=change_glass_d_mean, aes(x = MAP, y = losses_glass))+
  geom_point()+
  geom_smooth(method = 'lm')

even=lm(abs_evenness_glass~MAP + MAT + rrich + anpp, data=change_glass_d_mean)
summary(even)
rsq.partial(even)

rank=lm(rank_glass~MAP + MAT + rrich + anpp, data=change_glass_d_mean)
summary(rank)
rsq.partial(rank)

gains=lm(gains_glass~MAP + MAT + rrich + anpp, data=change_glass_d_mean)
summary(gains)
rsq.partial(gains)

losses=lm(losses_glass~MAP + MAT + rrich + anpp, data=change_glass_d_mean)
summary(losses)
rsq.partial(losses)

#-------1b) Only N addition studies

rich=lm(abs_richness_glass ~ MAP + MAT + rrich + anpp, data=use_change_glass_d_mean[use_change_glass_d_mean$n>0,])
vif(rich)
summary(rich)
as.data.frame(rsq.partial(rich)[2:3])

even=lm(abs_evenness_glass~MAP + MAT + rrich + anpp, data=use_change_glass_d_mean[use_change_glass_d_mean$n>0,])
summary(even)
as.data.frame(rsq.partial(even)[2:3])

rank=lm(rank_glass~MAP + MAT + rrich + anpp, data=use_change_glass_d_mean[use_change_glass_d_mean$n>0,])
summary(rank)
as.data.frame(rsq.partial(rank)[2:3])

gains=lm(gains_glass~MAP + MAT + rrich + anpp, data=use_change_glass_d_mean[use_change_glass_d_mean$n>0,])
summary(gains)
as.data.frame(rsq.partial(gains)[2:3])

losses=lm(losses_glass~MAP + MAT + rrich + anpp, data=use_change_glass_d_mean[use_change_glass_d_mean$n>0,])
summary(losses)
as.data.frame(rsq.partial(losses)[2:3])



#------2) including site as a random factor to group treatments occurring at the same site

rich=lmer(abs_richness_glass ~ MAP + MAT + rrich + anpp + (1|site_code), data=use_change_glass_d_mean)
Anova(rich)
summary(rich)

even=lmer(abs_evenness_glass~MAP + MAT + rrich + anpp + (1|site_code), data=use_change_glass_d_mean)
Anova(even)
summary(even)

rank=lmer(rank_glass~MAP + MAT + rrich + anpp + (1|site_code), data=use_change_glass_d_mean)
Anova(rank)
summary(rank)

gains=lmer(gains_glass~MAP + MAT + rrich + anpp + (1|site_code), data=use_change_glass_d_mean); Anova(gains); summary(gains)

losses=lmer(losses_glass~MAP + MAT + rrich + anpp + (1|site_code), data=use_change_glass_d_mean); Anova(losses); summary(losses)