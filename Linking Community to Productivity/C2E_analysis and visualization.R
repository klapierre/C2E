################################################################################
##  C2E_analysis and visualization.R: Data anaylsis and visualization.
##
##  Author: Kimberly Komatsu
##  Date created: July 12, 2021
################################################################################

# #for boosted regression trees
# library(mgcv)
# library(gratia)
# library(gbm)
# library(dismo)
# library(magrittr)
# library(rsq)

#for causal modeling
library(data.table)
library(fixest)
library(lme4)
library(tidyverse)

# #other
# library(emmeans)
# library(MuMIn)
# library(lmerTest)
# library(performance)
# library(tidyverse)
# library(ggpubr)

#options, themes, and functions

#kim's laptop
setwd('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\CoRRE\\C2E\\Products\\testing HRF\\data')

theme_set(theme_bw())
theme_update(axis.title.x=element_text(size=20, vjust=-0.35), axis.text.x=element_text(size=20),
             axis.title.y=element_text(size=20, angle=90, vjust=0.5), axis.text.y=element_text(size=20),
             plot.title = element_text(size=24, vjust=2),
             panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
             legend.title=element_text(size=16), legend.text=element_text(size=16))

###bar graph summary statistics function
#barGraphStats(data=, variable="", byFactorNames=c(""))
barGraphStats <- function(data, variable, byFactorNames) {
  count <- length(byFactorNames)
  N <- aggregate(data[[variable]], data[byFactorNames], FUN=length)
  names(N)[1:count] <- byFactorNames
  names(N) <- sub("^x$", "N", names(N))
  mean <- aggregate(data[[variable]], data[byFactorNames], FUN=mean)
  names(mean)[1:count] <- byFactorNames
  names(mean) <- sub("^x$", "mean", names(mean))
  sd <- aggregate(data[[variable]], data[byFactorNames], FUN=sd)
  names(sd)[1:count] <- byFactorNames
  names(sd) <- sub("^x$", "sd", names(sd))
  preSummaryStats <- merge(N, mean, by=byFactorNames)
  finalSummaryStats <- merge(preSummaryStats, sd, by=byFactorNames)
  finalSummaryStats$se <- finalSummaryStats$sd / sqrt(finalSummaryStats$N)
  return(finalSummaryStats)
}


options(contrasts=c('contr.sum','contr.poly'))

#not in function
`%!in%` = Negate(`%in%`)


##### import data #####
correAllDataTrt <- read.csv('CoRRE_comm and anpp diff_07122021.csv')%>%
  select(site_project_comm, site_project_comm_trt, site_code, project_name, community_type, treatment_year, calendar_year, treatment, richness_diff, evenness_diff, rank_diff, species_diff, composition_diff, abs_dispersion_diff, richness_change, evenness_change, rank_change, gains, losses, composition_change, dispersion_change, anpp_pdiff, plot_mani)%>%
  # drop_na()%>%
  filter(anpp_pdiff<4)%>%
  mutate(abs_anpp_pdiff=abs(anpp_pdiff))

correDiff <- correAllDataTrt%>%
  select(-richness_change, -evenness_change, -rank_change, -gains, -losses, -composition_change, -dispersion_change)%>%
  drop_na()

correChange <- correAllDataTrt%>%
  select(-richness_diff, -evenness_diff, -rank_diff, -species_diff, -composition_diff, -abs_dispersion_diff)%>%
  drop_na()

# ***marginal R2 estimate is based on only fixed factors, while conditional R2 estimate is based on full model (fixed and random effects)
# AIC greater than 6 to prefer polynomical model, due to sample size of 196 (trt within site_project_comm)

##### ANPP differences through time #####
summary(anppTemporal <- lmer(abs(anpp_pdiff)~poly(treatment_year,2) + (1+treatment_year|site_project_comm/treatment), data=correDiff))  #significant quadratic effect
check_model(anppTemporal)
r.squaredGLMM(anppTemporal) #marginal R2=0.01678326, conditional R2=0.7479032 
AIC(anppTemporal) #polynomial model (AIC=898.1287) better than liner model (AIC=935.4736)

ggplot(data=correDiff, aes(x=treatment_year, y=abs(anpp_pdiff))) +
  geom_point(color='grey') +
  geom_smooth(method='lm', formula=y~x+I(x^2), se=F, color='grey',aes(group=interaction(site_project_comm, treatment))) +
  geom_smooth(method='lm', formula=y~x+I(x^2), se=F, size=3, color='black') +
  xlab('Treatment Year') + ylab('|ANPP % Difference|') +
  theme(legend.position='none') +
  coord_cartesian(ylim=c(0,4))

##### ANPP difference related to community difference #####
###overall community change
summary(anppCommDiff <- lmer(abs(anpp_pdiff)~composition_diff + (1+treatment_year|site_project_comm/treatment), data=correDiff)) #significant effect
check_model(anppCommDiff)
r.squaredGLMM(anppCommDiff) #marginal R2=0.0875811, conditional R2=0.7393942 
AIC(anppCommDiff) #polynomial model (AIC=823.8902) NOT better than liner model (AIC=828.4461)

ggplot(data=correDiff, aes(x=composition_diff, y=abs(anpp_pdiff))) +
    geom_point(color='grey') +
    geom_smooth(method='lm', se=F, color='grey',aes(group=interaction(site_project_comm, treatment))) +
    geom_smooth(method='lm', se=F, size=3, color='black') +
    xlab('Composition Difference') + ylab('|ANPP % Difference|') +
    theme(legend.position='none')

###richness only
summary(anppRichDiff <- lmer(abs(anpp_pdiff)~poly(richness_diff,2) + (1+treatment_year|site_project_comm/treatment), data=correDiff)) #significant effect
check_model(anppRichDiff)
r.squaredGLMM(anppRichDiff) #marginal R2=0.06009032, conditional R2=0.6841572 
AIC(anppRichDiff) #polynomial model (AIC=848.7159) better than liner model (AIC=878.3956)

ggplot(data=correDiff, aes(x=richness_diff, y=abs(anpp_pdiff))) +
  geom_point(color='grey') +
  geom_smooth(method='lm', formula=y~x+I(x^2), se=F, color='grey',aes(group=interaction(site_project_comm, treatment))) +
  geom_smooth(method='lm', formula=y~x+I(x^2), se=F, size=3, color='black') +
  xlab('Richness Difference') + ylab('|ANPP % Difference|') +
  theme(legend.position='none') +
  coord_cartesian(ylim=c(0,4))

###multiple change metrics
#individually to determine quadratic vs linear
AIC(anppEvenDiff <- lmer(abs(anpp_pdiff)~evenness_diff + (1+treatment_year|site_project_comm/treatment), data=correDiff)) #polynomial model (AIC=929.7081) NOT better than linear model (AIC=931.1781)
AIC(anppRankDiff <- lmer(abs(anpp_pdiff)~rank_diff + (1+treatment_year|site_project_comm/treatment), data=correDiff)) #polynomial model (AIC=909.5364) NOT better than linear model (AIC=910.5558)
AIC(anppGainDiff <- lmer(abs(anpp_pdiff)~species_diff + (1+treatment_year|site_project_comm/treatment), data=correDiff)) #polynomial model (AIC=929.8562) better than linear model (AIC=931.3553)

#all together
summary(anppCommMetricDiff <- lmer(abs(anpp_pdiff)~evenness_diff+rank_diff+species_diff + (1+treatment_year|site_project_comm/treatment), data=correDiff)) #significant effect of rank diff
check_model(anppCommMetricDiff)
r.squaredGLMM(anppCommMetricDiff) #marginal R2=0.01977212, conditional R2=0.7415795

ggplot(data=correDiff, aes(x=evenness_diff, y=abs(anpp_pdiff))) +
  geom_point(color='grey') +
  geom_smooth(method='lm', formula=y~x+I(x^2), se=F, color='grey',aes(group=interaction(site_project_comm, treatment))) +
  geom_smooth(method='lm', formula=y~x+I(x^2), se=F, size=3, color='black') +
  xlab('Evenness Difference') + ylab('|ANPP % Difference|') +
  theme(legend.position='none') +
  coord_cartesian(ylim=c(0,4))

ggplot(data=correDiff, aes(x=rank_diff, y=abs(anpp_pdiff))) +
  geom_point(color='grey') +
  geom_smooth(method='lm', formula=y~x+I(x^2), se=F, color='grey',aes(group=interaction(site_project_comm, treatment))) +
  geom_smooth(method='lm', formula=y~x+I(x^2), se=F, size=3, color='black') +
  xlab('Rank Difference') + ylab('|ANPP % Difference|') +
  theme(legend.position='none') +
  coord_cartesian(ylim=c(0,4))

ggplot(data=correDiff, aes(x=species_diff, y=abs(anpp_pdiff))) +
  geom_point(color='grey') +
  geom_smooth(method='lm', formula=y~x+I(x^2), se=F, color='grey',aes(group=interaction(site_project_comm, treatment))) +
  geom_smooth(method='lm', formula=y~x+I(x^2), se=F, size=3, color='black') +
  xlab('Species Difference') + ylab('|ANPP % Difference|') +
  theme(legend.position='none') +
  coord_cartesian(ylim=c(0,4))


##### site drivers of ANPP difference #####
siteInfo <- read.csv('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\CoRRE\\converge_diverge\\datasets\\LongForm\\SiteExperimentDetails_March2019.csv')%>%
  select(-X)

correDiffSite <- correDiff%>%
  left_join(siteInfo)

summary(anppCommMetricDiff <- lmer(abs(anpp_pdiff) ~ evenness_diff + rank_diff + species_diff + rrich + (1+treatment_year|site_project_comm/treatment), data=correDiffSite)) #significant effects of rank diff, species diff, and rrich
check_model(anppCommMetricDiff)
r.squaredGLMM(anppCommMetricDiff)

summary(anppSiteDrivers <- lmer(abs(anpp_pdiff)~MAP+MAT+rrich+anpp+(1+treatment_year|site_project_comm/treatment), data=correDiffSite)) #significant positive effect of rrich, negative effect of anpp

ggplot(data=correDiffSite, aes(x=rrich, y=abs(anpp_pdiff))) +
  geom_point(color='grey') +
  geom_smooth(method='lm', se=F, color='grey',aes(group=interaction(site_project_comm, treatment))) +
  geom_smooth(method='lm', se=F, size=3, color='black') +
  xlab('Relative Richness') + ylab('|ANPP % Difference|') +
  theme(legend.position='none') +
  coord_cartesian(ylim=c(0,4))

ggplot(data=correDiffSite, aes(x=anpp, y=abs(anpp_pdiff))) +
  geom_point(color='grey') +
  geom_smooth(method='lm', se=F, color='grey',aes(group=interaction(site_project_comm, treatment))) +
  geom_smooth(method='lm', se=F, size=3, color='black') +
  xlab('Site Productivity') + ylab('|ANPP % Difference|') +
  theme(legend.position='none') +
  coord_cartesian(ylim=c(0,4))

ggplot(data=correDiffSite, aes(x=rank_diff, y=abs(anpp_pdiff))) +
  geom_point() +
  geom_smooth(method='lm', formula=y~x+I(x^2), se=F, aes(group=interaction(site_project_comm, treatment), color=rrich)) +
  geom_smooth(method='lm', formula=y~x+I(x^2), se=F, size=3, color='black') +
  xlab('Rank Difference') + ylab('|ANPP % Difference|') +
  theme(legend.position='none') +
  coord_cartesian(ylim=c(0,4))

ggplot(data=correDiffSite, aes(x=species_diff, y=abs(anpp_pdiff))) +
  geom_point(color='grey') +
  geom_smooth(method='lm', formula=y~x+I(x^2), se=F, aes(group=interaction(site_project_comm, treatment), color=rrich)) +
  geom_smooth(method='lm', formula=y~x+I(x^2), se=F, size=3, color='black') +
  xlab('Species Difference') + ylab('|ANPP % Difference|') +
  theme(legend.position='none') +
  coord_cartesian(ylim=c(0,4))


##### boosted regression trees #####
#bidirectional response
brtModel <-  gbm.step(data=correDiffSite, gbm.x = c('MAP', 'MAT', 'rrich', 'anpp', 'richness_diff', 'evenness_diff', 'rank_diff', 'species_diff', 'plot_mani'),
                      gbm.y = "anpp_pdiff", family = "gaussian", #can consider other distributions
                      tree.complexity = 5, learning.rate = 0.001, #these values punish the models for over fitting
                      bag.fraction = 0.5, step.size=100, max.trees = 50000) #learning rate (0.5 is high) and step size are important

#get some metrics from the model:
brtModel$gbm.call$best.trees #14650 trees
hist(brtModel$residuals)
contributions <- brtModel$contributions
cor(correDiffSite$anpp_pdiff, brtModel$fitted)^2 #pseudo-R2 = 0.6414528
mean(brtModel$residuals * brtModel$residuals) #MSE = 0.1152091
brtModel$cv.statistics$deviance.mean #minimum cv deviance = 0.1979288
brtModel$cv.statistics$deviance.se #cv deviance se = 0.01132706

# contributions plot
ggplot(data=contributions, aes(x=var, y=rel.inf)) + 
  geom_bar(stat='identity') +
  coord_flip() + 
  scale_x_discrete(limits = c('MAP', 'MAT', 'rrich', 'anpp', 'richness_diff', 'evenness_diff', 'rank_diff', 'species_diff', 'plot_mani')) +
  xlab('Variable') + ylab('Relative Influence')

# partial dependence plots
gbm.plot(brtModel, smooth=T, write.title=F, y.label="ANPP % Difference")

find.int <- gbm.interactions(brtModel) #this is the interactions
find.int$rank.list
find.int$interactions
# The returned object is a list. The first 2 components summarise the results, first as a ranked list of the 5 most important pairwise interactions, and the second tabulating all pairwise interactions. The variable index numbers in $rank.list can be used for plotting.
# You can plot pairwise interactions like this:
gbm.perspec(brtModel,9,6) #plot_mani vs evenness_diff
gbm.perspec(brtModel,7,6) #rank_diff vs evenness_diff
gbm.perspec(brtModel,8,4) #species_diff vs anpp - uninformative
gbm.perspec(brtModel,5,1) #richness_diff vs MAP - uninformative

ggplot(data=correDiffSite, aes(x=richness_diff, y=anpp_pdiff)) + geom_point() + ylab('ANPP % Difference')


###abs value (magnitude) of response
brtModel <-  gbm.step(data=correDiffSite, gbm.x = c('MAP', 'MAT', 'rrich', 'anpp', 'richness_diff', 'evenness_diff', 'rank_diff', 'species_diff', 'plot_mani'),
                      gbm.y = "abs_anpp_pdiff", family = "gaussian", #can consider other distributions
                      tree.complexity = 5, learning.rate = 0.001, #these values punish the models for over fitting
                      bag.fraction = 0.5, step.size=100, max.trees = 50000) #learning rate (0.5 is high) and step size are important

#get some metrics from the model:
brtModel$gbm.call$best.trees #14450 trees
hist(brtModel$residuals)
contributions <- brtModel$contributions
cor(correDiffSite$anpp_pdiff, brtModel$fitted)^2 #pseudo-R2 = 0.4582868
mean(brtModel$residuals * brtModel$residuals) #MSE = 0.0664859
brtModel$cv.statistics$deviance.mean #minimum cv deviance = 0.1162645
brtModel$cv.statistics$deviance.se #cv deviance se = 0.008857221

# contributions plot
ggplot(data=contributions, aes(x=var, y=rel.inf)) + 
  geom_bar(stat='identity') +
  coord_flip() + 
  scale_x_discrete(limits = c('MAP', 'MAT', 'rrich', 'anpp', 'richness_diff', 'evenness_diff', 'rank_diff', 'species_diff', 'plot_mani')) +
  xlab('Variable') + ylab('Relative Influence')

# partial dependence plots
gbm.plot(brtModel, smooth=T, write.title=F, y.label="ANPP % Difference")

find.int <- gbm.interactions(brtModel) #this is the interactions
find.int$rank.list
find.int$interactions
# The returned object is a list. The first 2 components summarise the results, first as a ranked list of the 5 most important pairwise interactions, and the second tabulating all pairwise interactions. The variable index numbers in $rank.list can be used for plotting.
# You can plot pairwise interactions like this:
gbm.perspec(brtModel,7,5) #rank_diff vs richness_diff
gbm.perspec(brtModel,8,5) #species_diff vs richness_diff
gbm.perspec(brtModel,6,3) #evenness_diff vs rrich -- uninformative
gbm.perspec(brtModel,7,6) #rank_diff vs evenness_diff

ggplot(data=correDiffSite, aes(x=richness_diff, y=abs_anpp_pdiff)) + geom_point() + ylab('ANPP % Difference')



##### causal modeling - difference #####
correDiffFixed <- correDiff%>%
  mutate(site_year=paste(site_project_comm, treatment_year, sep="::"))%>% #creates dummy site-by-year variable (time variant site effects)
  mutate(site_treatment=paste(site_project_comm, treatment, sep='::')) #creates dummy site-by-plot variable (time invariant plot effects)

diffCausalModel <- feols(anpp_pdiff ~ richness_diff + evenness_diff + rank_diff + species_diff | site_treatment + site_year, data=correDiffFixed)
etable(diffCausalModel,
       cluster="site_treatment")


##### causal modeling - change #####
#to do: try doing change on residuals from weather data? or just change overall for each plot to anpp because the time variant site effects should account for rainfall?
correChangeFixed <- correChange%>%
  mutate(site_year=paste(site_project_comm, treatment_year, sep="::"))%>% #creates dummy site-by-year variable (time variant site effects)
  mutate(site_treatment=paste(site_project_comm, treatment, sep='::')) #creates dummy site-by-plot variable (time invariant plot effects)

changeCausalModel <- feols(anpp_pdiff ~ richness_change + evenness_change + rank_change + losses | site_treatment + site_year, data=correChangeFixed)
etable(changeCausalModel,
       cluster="site_treatment")



##### causal modeling - change by plot #####
correCommChangePlotANPP <- read.csv('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\CoRRE\\C2E\\Products\\testing HRF\\data\\CoRRE_comm and anpp change_by plot_07202021.csv')%>%
  mutate(site_year=paste(site_project_comm, treatment_year, sep="::"))%>% #creates dummy site-by-year variable (time variant site effects)
  mutate(site_plot=paste(site_project_comm, plot_id, sep='::')) #creates dummy site-by-plot variable (time invariant plot effects)

#need to calculate change in anpp from year 0 as well to include in the model instead of raw anpp values
#also try including the treatment variables as interactors with the change variables? like N effects alone vs interating with richness, evenness etc to get at the physiological responses vs the community driven responses?
changePlotCausalModel <- feols(anpp ~ richness_change + evenness_change + rank_change + gains + losses | site_plot + site_year, data=subset(correCommChangePlotANPP, anpp<4000))
etable(changePlotCausalModel,
       cluster="site_plot")

# significant effects of evenness change (negative) and rank change (positive) on ANPP
with(subset(correCommChangePlotANPP, anpp<4000), plot(anpp, evenness_change))
with(subset(correCommChangePlotANPP, anpp<4000), plot(anpp, rank_change))

#losses removed from model because collinearity

with(correCommChangePlotANPP, plot(gains, losses))
with(correCommChangePlotANPP, plot(richness_change, losses))
with(correCommChangePlotANPP, plot(evenness_change, losses))
with(correCommChangePlotANPP, plot(rank_change, losses))
with(correCommChangePlotANPP, plot(anpp, losses))

etable(feols(losses ~ gains | site_plot + site_year, data=correCommChangePlotANPP), cluster="site_plot") #sig related
etable(feols(losses ~ richness_change | site_plot + site_year, data=correCommChangePlotANPP), cluster="site_plot") #sig related
etable(feols(losses ~ evenness_change | site_plot + site_year, data=correCommChangePlotANPP), cluster="site_plot") #sig related
etable(feols(losses ~ rank_change | site_plot + site_year, data=correCommChangePlotANPP), cluster="site_plot") #sig related






