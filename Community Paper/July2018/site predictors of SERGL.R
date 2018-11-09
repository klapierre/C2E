#emily's working directory
setwd("/Users/egrman/Dropbox/C2E/Products/CommunityChange/March2018 WG")

#meghan's working directory
setwd("C:\\Users\\megha\\Dropbox\\C2E\\Products\\CommunityChange\\March2018 WG")
setwd("~/Dropbox/C2E/Products/CommunityChange/March2018 WG")


#kevin's working directory
#setwd("C:\\Users\\wilco\\Dropbox\\C2E\\Products\\CommunityChange\\March 2018 WG")

library(tidyverse)
library(ggthemes)
library(grid)
library(vegan)
library(rsq)
library(car)


theme_set(theme_bw())
theme_update(axis.title.x=element_text(size=12, vjust=-0.35), axis.text.x=element_text(size=12),
             axis.title.y=element_text(size=12, angle=90, vjust=0.5), axis.text.y=element_text(size=12),
             plot.title = element_text(size=24, vjust=2),
             panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
             legend.title=element_text(size=12), legend.text=element_text(size=12))

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
png(paste0("Summer2018_Results/site predictors of SERGL/MR predictor variables SITE LEVEL pairs plot.png"), width=11, height=8, units="in", res=600)
print(pairs(pred))
dev.off()

#-------------Standardized multiple regression with only site predictors

usethese=change_metrics[change_metrics$use==1, c("site_project_comm", "treatment", "use")]
use_change_glass_d_mean=merge(change_glass_d_mean, unique(usethese), by=c("site_project_comm", "treatment"))
#write.csv(use_change_glass_d_mean, "use for site predictors of SERGL.csv")

##for this first analysis I think we should use all the data 
length(use_change_glass_d_mean$site_project_comm) #202 studies

#note that some response var (evenness, losses, gains) have NA (for every year there was no variation among the controls; sd=0 and glass's delta was undefined for every year)

#-----1a) treating all experiments in this subset as independent data points:

change_glass_d_mean$sMAP<-scale(change_glass_d_mean$MAP)
change_glass_d_mean$sMAT<-scale(change_glass_d_mean$MAT)
change_glass_d_mean$srrich<-scale(change_glass_d_mean$rrich)
change_glass_d_mean$sanpp<-scale(change_glass_d_mean$anpp)

rich=lm(abs_richness_glass ~ sMAP + sMAT + srrich + sanpp, data=change_glass_d_mean)
#vif(rich) #need the car package for this.
summary(rich)
rsq.partial(rich)

#making object to contain results
richresults=data.frame(response="rich", predictor=names(rich$coefficients), slope=as.numeric(rich$coefficients), pval=as.numeric(summary(rich)$coef[,4]), rsq=c(NA, rsq.partial(rich)$partial.rsq))


even=lm(abs_evenness_glass~sMAP + sMAT + srrich + sanpp, data=change_glass_d_mean)
summary(even)
rsq.partial(even)
evenresults=data.frame(response="even", predictor=names(even$coefficients), slope=as.numeric(even$coefficients), pval=as.numeric(summary(even)$coef[,4]), rsq=c(NA, rsq.partial(even)$partial.rsq))

rank=lm(rank_glass~sMAP + sMAT + srrich + sanpp, data=change_glass_d_mean)
summary(rank)
rsq.partial(rank)
rankresults=data.frame(response="rank", predictor=names(rank$coefficients), slope=as.numeric(rank$coefficients), pval=as.numeric(summary(rank)$coef[,4]), rsq=c(NA, rsq.partial(rank)$partial.rsq))

gains=lm(gains_glass~sMAP + sMAT + srrich + sanpp, data=change_glass_d_mean)
summary(gains)
rsq.partial(gains)
gainsresults=data.frame(response="gains", predictor=names(gains$coefficients), slope=as.numeric(gains$coefficients), pval=as.numeric(summary(gains)$coef[,4]), rsq=c(NA, rsq.partial(gains)$partial.rsq))

losses=lm(losses_glass~sMAP + sMAT + srrich + sanpp, data=change_glass_d_mean)
summary(losses)
rsq.partial(losses)
lossesresults=data.frame(response="losses", predictor=names(losses$coefficients), slope=as.numeric(losses$coefficients), pval=as.numeric(summary(losses)$coef[,4]), rsq=c(NA, rsq.partial(losses)$partial.rsq))

fulldataset=rbind(richresults, evenresults, rankresults, gainsresults, lossesresults)
#write.csv(fulldataset, "Summer2018_Results/site predictors of SERGL/site predictors of SERGL, all studies.csv", row.names=F)


#-------1b) Only a subset of studies

summary(use_change_glass_d_mean) #so we can look at relationships separately for experiments manipulating each of these factors: N, P, K, CO2, precip, temp?
length(use_change_glass_d_mean[use_change_glass_d_mean$n>0,]$site_project_comm) #131 studies
length(use_change_glass_d_mean[use_change_glass_d_mean$p>0,]$site_project_comm) #73 studies
length(use_change_glass_d_mean[use_change_glass_d_mean$k>0,]$site_project_comm) #only 21 studies
length(use_change_glass_d_mean[use_change_glass_d_mean$CO2>0,]$site_project_comm) #only 10 studies
length(use_change_glass_d_mean[use_change_glass_d_mean$precip>0,]$site_project_comm) #38 studies
length(use_change_glass_d_mean[use_change_glass_d_mean$temp>0,]$site_project_comm) #only 2 studies

#----------N addition studies:
use_change_glass_d_mean$sMAP<-scale(use_change_glass_d_mean$MAP)
use_change_glass_d_mean$sMAT<-scale(use_change_glass_d_mean$MAT)
use_change_glass_d_mean$srrich<-scale(use_change_glass_d_mean$rrich)
use_change_glass_d_mean$sanpp<-scale(use_change_glass_d_mean$anpp)


rich=lm(abs_richness_glass ~ sMAP + sMAT + srrich + sanpp, data=use_change_glass_d_mean[use_change_glass_d_mean$n>0,])
vif(rich)
richresults=data.frame(response="rich", predictor=names(rich$coefficients), slope=as.numeric(rich$coefficients), pval=as.numeric(summary(rich)$coef[,4]), rsq=c(NA, rsq.partial(rich)$partial.rsq))

even=lm(abs_evenness_glass~sMAP + sMAT + srrich + sanpp, data=use_change_glass_d_mean[use_change_glass_d_mean$n>0,])
evenresults=data.frame(response="even", predictor=names(even$coefficients), slope=as.numeric(even$coefficients), pval=as.numeric(summary(even)$coef[,4]), rsq=c(NA, rsq.partial(even)$partial.rsq))

rank=lm(rank_glass~sMAP + sMAT + srrich + sanpp, data=use_change_glass_d_mean[use_change_glass_d_mean$n>0,])
rankresults=data.frame(response="rank", predictor=names(rank$coefficients), slope=as.numeric(rank$coefficients), pval=as.numeric(summary(rank)$coef[,4]), rsq=c(NA, rsq.partial(rank)$partial.rsq))

gains=lm(gains_glass~sMAP + sMAT + srrich + sanpp,data=use_change_glass_d_mean[use_change_glass_d_mean$n>0,])
gainsresults=data.frame(response="gains", predictor=names(gains$coefficients), slope=as.numeric(gains$coefficients), pval=as.numeric(summary(gains)$coef[,4]), rsq=c(NA, rsq.partial(gains)$partial.rsq))

losses=lm(losses_glass~sMAP + sMAT + srrich + sanpp, data=use_change_glass_d_mean[use_change_glass_d_mean$n>0,])
lossesresults=data.frame(response="losses", predictor=names(losses$coefficients), slope=as.numeric(losses$coefficients), pval=as.numeric(summary(losses)$coef[,4]), rsq=c(NA, rsq.partial(losses)$partial.rsq))

onlyNadditions=rbind(richresults, evenresults, rankresults, gainsresults, lossesresults)
#write.csv(onlyNadditions, "Summer2018_Results/site predictors of SERGL/site predictors of SERGL, only N additions.csv", row.names=F)

#----------P addition studies:

rich=lm(abs_richness_glass ~ sMAP + sMAT + srrich + sanpp,data=use_change_glass_d_mean[use_change_glass_d_mean$p>0,])
vif(rich)
richresults=data.frame(response="rich", predictor=names(rich$coefficients), slope=as.numeric(rich$coefficients), pval=as.numeric(summary(rich)$coef[,4]), rsq=c(NA, rsq.partial(rich)$partial.rsq))

even=lm(abs_evenness_glass~sMAP + sMAT + srrich + sanpp, data=use_change_glass_d_mean[use_change_glass_d_mean$p>0,])
evenresults=data.frame(response="even", predictor=names(even$coefficients), slope=as.numeric(even$coefficients), pval=as.numeric(summary(even)$coef[,4]), rsq=c(NA, rsq.partial(even)$partial.rsq))

rank=lm(rank_glass~sMAP + sMAT + srrich + sanpp, data=use_change_glass_d_mean[use_change_glass_d_mean$p>0,])
rankresults=data.frame(response="rank", predictor=names(rank$coefficients), slope=as.numeric(rank$coefficients), pval=as.numeric(summary(rank)$coef[,4]), rsq=c(NA, rsq.partial(rank)$partial.rsq))

gains=lm(gains_glass~sMAP + sMAT + srrich + sanpp, data=use_change_glass_d_mean[use_change_glass_d_mean$p>0,])
gainsresults=data.frame(response="gains", predictor=names(gains$coefficients), slope=as.numeric(gains$coefficients), pval=as.numeric(summary(gains)$coef[,4]), rsq=c(NA, rsq.partial(gains)$partial.rsq))

losses=lm(losses_glasss~MAP + sMAT + srrich + sanpp,data=use_change_glass_d_mean[use_change_glass_d_mean$p>0,])
lossesresults=data.frame(response="losses", predictor=names(losses$coefficients), slope=as.numeric(losses$coefficients), pval=as.numeric(summary(losses)$coef[,4]), rsq=c(NA, rsq.partial(losses)$partial.rsq))

onlyPadditions=rbind(richresults, evenresults, rankresults, gainsresults, lossesresults)

#write.csv(onlyPadditions, "Summer2018_Results/site predictors of SERGL/site predictors of SERGL, only P additions.csv", row.names=F)

#----------K addition studies:

rich=lm(abs_richness_glass ~ sMAP + sMAT + srrich + sanpp, data=use_change_glass_d_mean[use_change_glass_d_mean$k>0,])
vif(rich)
richresults=data.frame(response="rich", predictor=names(rich$coefficients), slope=as.numeric(rich$coefficients), pval=as.numeric(summary(rich)$coef[,4]), rsq=c(NA, rsq.partial(rich)$partial.rsq))

even=lm(abs_evenness_glass~sMAP + sMAT + srrich + sanpp, data=use_change_glass_d_mean[use_change_glass_d_mean$k>0,])
evenresults=data.frame(response="even", predictor=names(even$coefficients), slope=as.numeric(even$coefficients), pval=as.numeric(summary(even)$coef[,4]), rsq=c(NA, rsq.partial(even)$partial.rsq))

rank=lm(rank_glass~sMAP + sMAT + srrich + sanpp, data=use_change_glass_d_mean[use_change_glass_d_mean$k>0,])
rankresults=data.frame(response="rank", predictor=names(rank$coefficients), slope=as.numeric(rank$coefficients), pval=as.numeric(summary(rank)$coef[,4]), rsq=c(NA, rsq.partial(rank)$partial.rsq))

gains=lm(gains_glass~sMAP + sMAT + srrich + sanpp,data=use_change_glass_d_mean[use_change_glass_d_mean$k>0,])
gainsresults=data.frame(response="gains", predictor=names(gains$coefficients), slope=as.numeric(gains$coefficients), pval=as.numeric(summary(gains)$coef[,4]), rsq=c(NA, rsq.partial(gains)$partial.rsq))

losses=lm(losses_glass~sMAP + sMAT + srrich + sanpp, data=use_change_glass_d_mean[use_change_glass_d_mean$k>0,])
lossesresults=data.frame(response="losses", predictor=names(losses$coefficients), slope=as.numeric(losses$coefficients), pval=as.numeric(summary(losses)$coef[,4]), rsq=c(NA, rsq.partial(losses)$partial.rsq))

onlyKadditions=rbind(richresults, evenresults, rankresults, gainsresults, lossesresults)

#write.csv(onlyKadditions, "Summer2018_Results/site predictors of SERGL/site predictors of SERGL, only K additions.csv", row.names=F)

#----------CO2 manipulation studies:

rich=lm(abs_richness_glass ~ sMAP + sMAT + srrich + sanpp, data=use_change_glass_d_mean[use_change_glass_d_mean$CO2>0,])
vif(rich)
richresults=data.frame(response="rich", predictor=names(rich$coefficients), slope=as.numeric(rich$coefficients), pval=as.numeric(summary(rich)$coef[,4]), rsq=c(NA, rsq.partial(rich)$partial.rsq))

even=lm(abs_evenness_glass~sMAP + sMAT + srrich + sanpp,data=use_change_glass_d_mean[use_change_glass_d_mean$CO2>0,])
evenresults=data.frame(response="even", predictor=names(even$coefficients), slope=as.numeric(even$coefficients), pval=as.numeric(summary(even)$coef[,4]), rsq=c(NA, rsq.partial(even)$partial.rsq))

rank=lm(rank_glass~sMAP + sMAT + srrich + sanpp, data=use_change_glass_d_mean[use_change_glass_d_mean$CO2>0,])
rankresults=data.frame(response="rank", predictor=names(rank$coefficients), slope=as.numeric(rank$coefficients), pval=as.numeric(summary(rank)$coef[,4]), rsq=c(NA, rsq.partial(rank)$partial.rsq))

gains=lm(gains_glass~sMAP + sMAT + srrich + sanpp, data=use_change_glass_d_mean[use_change_glass_d_mean$CO2>0,])
gainsresults=data.frame(response="gains", predictor=names(gains$coefficients), slope=as.numeric(gains$coefficients), pval=as.numeric(summary(gains)$coef[,4]), rsq=c(NA, rsq.partial(gains)$partial.rsq))

losses=lm(losses_glass~sMAP + sMAT + srrich + sanpp,data=use_change_glass_d_mean[use_change_glass_d_mean$CO2>0,])
lossesresults=data.frame(response="losses", predictor=names(losses$coefficients), slope=as.numeric(losses$coefficients), pval=as.numeric(summary(losses)$coef[,4]), rsq=c(NA, rsq.partial(losses)$partial.rsq))

onlyCO2manipulations=rbind(richresults, evenresults, rankresults, gainsresults, lossesresults)

#write.csv(onlyCO2manipulations, "Summer2018_Results/site predictors of SERGL/site predictors of SERGL, only CO2 manip.csv", row.names=F)

#----------precip studies:

rich=lm(abs_richness_glass ~ sMAP + sMAT + srrich + sanpp, data=use_change_glass_d_mean[use_change_glass_d_mean$precip>0,])
vif(rich)
richresults=data.frame(response="rich", predictor=names(rich$coefficients), slope=as.numeric(rich$coefficients), pval=as.numeric(summary(rich)$coef[,4]), rsq=c(NA, rsq.partial(rich)$partial.rsq))

even=lm(abs_evenness_glass~sMAP + sMAT + srrich + sanpp, data=use_change_glass_d_mean[use_change_glass_d_mean$precip>0,])
evenresults=data.frame(response="even", predictor=names(even$coefficients), slope=as.numeric(even$coefficients), pval=as.numeric(summary(even)$coef[,4]), rsq=c(NA, rsq.partial(even)$partial.rsq))

rank=lm(rank_glass~sMAP + sMAT + srrich + sanpp, data=use_change_glass_d_mean[use_change_glass_d_mean$precip>0,])
rankresults=data.frame(response="rank", predictor=names(rank$coefficients), slope=as.numeric(rank$coefficients), pval=as.numeric(summary(rank)$coef[,4]), rsq=c(NA, rsq.partial(rank)$partial.rsq))

gains=lm(gains_glass~sMAP + sMAT + srrich + sanpp, data=use_change_glass_d_mean[use_change_glass_d_mean$precip>0,])
gainsresults=data.frame(response="gains", predictor=names(gains$coefficients), slope=as.numeric(gains$coefficients), pval=as.numeric(summary(gains)$coef[,4]), rsq=c(NA, rsq.partial(gains)$partial.rsq))

losses=lm(losses_glass~sMAP + sMAT + srrich + sanpp, data=use_change_glass_d_mean[use_change_glass_d_mean$precip>0,])
lossesresults=data.frame(response="losses", predictor=names(losses$coefficients), slope=as.numeric(losses$coefficients), pval=as.numeric(summary(losses)$coef[,4]), rsq=c(NA, rsq.partial(losses)$partial.rsq))

onlyprecipmanipulations=rbind(richresults, evenresults, rankresults, gainsresults, lossesresults)
#write.csv(onlyprecipmanipulations, "Summer2018_Results/site predictors of SERGL/site predictors of SERGL, only precip manip.csv", row.names=F)

#-------make a figure

fulldataset$studies="All manipulations"
onlyNadditions$studies="N addition"
onlyPadditions$studies="P addition"
onlyKadditions$studies="K addition"
onlyCO2manipulations$studies="CO2 manipulation"
onlyprecipmanipulations$studies="Water manipulation"
forbigfig=rbind(fulldataset, onlyNadditions, onlyPadditions, onlyKadditions, onlyCO2manipulations, onlyprecipmanipulations)
levels(forbigfig$response)=c("Richness change", "Evenness change", "Rank change", "Species gains", "Species losses")
levels(forbigfig$predictor)=c("(Intercept)", "ANPP", "MAP", "MAT", "Regional SR")
forbigfig$significant=as.factor(1*(forbigfig$pval<0.05))
forbigfig$star.location=ifelse(forbigfig$slope>0, forbigfig$slope+0.1, forbigfig$slope-0.1)

ggplot(aes(predictor, slope, fill=rsq), data=forbigfig[!forbigfig$predictor=="(Intercept)",]) + geom_col() + geom_point(aes(predictor, star.location, shape=significant)) + facet_grid(studies~response) + theme(axis.text.x=element_text(angle = 90, vjust = 0.4, hjust=1)) + xlab("Aspect of community change") + ylab("Effect on aspect of community change\n(slope from standardized multiple regression)") + guides(fill = guide_colorbar(title = "Partial R2")) + scale_shape_manual(values=c(NA, 8), guide=FALSE)
#ggsave("Summer2018_Results/site predictors of SERGL/site predictors of SERGL.pdf", width=7, height=9.5)
#ggsave("Summer2018_Results/site predictors of SERGL/site predictors of SERGL.png", width=7, height=9.5)

ggplot(aes(response, slope, fill=rsq), data=forbigfig[!forbigfig$predictor=="(Intercept)" & !forbigfig$pval>0.05,]) + geom_col() + facet_grid(studies~predictor) + theme(axis.text.x=element_text(angle = 90, vjust = 0.4)) + xlab("") + ylab("Effect on aspect of community change\n(slope from standardized multiple regression)") + guides(fill = guide_colorbar(title = "Partial R2"))
ggsave("Summer2018_Results/site predictors of SERGL/significant site predictors of SERGL.pdf", width=8, height=9.5)

