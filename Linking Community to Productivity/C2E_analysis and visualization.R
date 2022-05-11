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

#for data formatting and figures
library(tidyverse)
library(grid)

#for causal modeling
library(data.table)
library(fixest)
library(lme4)

# #other
library(emmeans)
library(MuMIn)
library(lmerTest)
# library(performance)
# library(tidyverse)
# library(ggpubr)


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
#site details
siteData <- read.csv('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\CoRRE\\CoRRE_database\\Data\\CompiledData\\siteBiotic.csv')%>%
  left_join(read.csv('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\CoRRE\\CoRRE_database\\Data\\CompiledData\\siteLocationClimate.csv'))

correAllDataTrt <- read.csv('CoRRE_comm and anpp diff_05022022.csv')%>%
  select(site_project_comm, site_project_comm_trt, site_code, project_name, community_type, treatment_year, calendar_year, treatment, trt_type, trt_type_2, richness_diff, evenness_diff, rank_diff, species_diff, composition_diff, abs_dispersion_diff, richness_change, evenness_change, rank_change, gains, losses, composition_change, dispersion_change, anpp_pdiff, plot_mani)%>%
  # drop_na()%>%
  filter(anpp_pdiff<4)%>%
  mutate(abs_anpp_pdiff=abs(anpp_pdiff))%>%
  filter(!(project_name %in% c('GFP', 'TIDE', 'NPKDNet', 'gb')))%>% #filter out datasets with less than 3 datapoints
  filter(!(project_name %in% c('TON', 'RMAPC') & !(site_code %in% c('ORNL')))) #filter out datasets with a biiiiiig gap before final year of data, which then drives all trends

correDiff <- correAllDataTrt%>%
  select(-richness_change, -evenness_change, -rank_change, -gains, -losses, -composition_change, -dispersion_change)%>%
  drop_na()

#dataset details
correDiffDetails <- correDiff%>%
  group_by(site_project_comm, site_project_comm_trt, site_code, project_name, community_type, treatment, trt_type, trt_type_2, plot_mani)%>%
  summarise(data_points=length(anpp_pdiff), exp_length=max(treatment_year), first_yr=min(treatment_year))%>%
  ungroup()%>%
  left_join(siteData)

#filter to just the projects with enough data points
correDiffLong <- correDiff%>%
  left_join(correDiffDetails)%>%
  filter(data_points>4)

#filter to just the treatments with large replication
correDiffMain <- correDiffLong%>%
  filter(trt_type %in% c('CO2', 'herb_removal', 'irr', 'K', 'mult_nutrient', 'N', 'N*irr', 'N*P', 'P', 'precip_vari', 'temp'))

# correChange <- correAllDataTrt%>%
#   select(-richness_diff, -evenness_diff, -rank_diff, -species_diff, -composition_diff, -abs_dispersion_diff)%>%
#   drop_na()

# ***marginal R2 estimate is based on only fixed factors, while conditional R2 estimate is based on full model (fixed and random effects)
# AIC greater than 6 to prefer polynomial model, due to sample size of 196 (trt within site_project_comm)

##### ANPP differences through time #####
summary(anppTemporal <- lmer(abs(anpp_pdiff)~poly(treatment_year,2) + (1+treatment_year|site_project_comm/treatment), data=correDiffLong))  #no effect (marginally significant negative quadratic term)
r.squaredGLMM(anppTemporal)
AIC(anppTemporal)

ggplot(data=correDiffLong, aes(x=treatment_year, y=anpp_pdiff)) +
  geom_point(color='grey') +
  geom_smooth(method='lm', formula=y~x+I(x^2), se=F, color='grey',aes(group=interaction(site_project_comm, treatment))) +
  geom_smooth(method='lm', formula=y~x+I(x^2), se=F, size=3, color='black') +
  xlab('Treatment Year') + ylab('ANPP % Difference') +
  theme(legend.position='none')

#look at directionality by trt type
ggplot(data=subset(correDiffLong, trt_type %in% c('CO2', 'drought', 'herb_removal', 'irr', 'K', 'mult_nutrient', 'N', 'N*irr', 'N*P', 'P', 'precip_vari', 'temp')), aes(x=treatment_year, y=anpp_pdiff)) +
  geom_point(color='grey') +
  geom_smooth(method='lm', formula=y~x+I(x^2), se=F, color='grey',aes(group=interaction(site_project_comm, treatment))) +
  geom_smooth(method='lm', formula=y~x+I(x^2), se=F, size=3, aes(color=trt_type)) +
  xlab('Treatment Year') + ylab('ANPP % Difference') +
  facet_wrap(~trt_type)

##### checking for patterns #####
site_project_comm_vector <- unique(correDiffLong$site_project_comm) 
site_project_comm_trt_vector <- unique(correDiffLong$site_project_comm_trt) 

##### ANPP over time #####
anppTimeModels=data.frame(row.names=1) 
for(PROJ in 1:length(site_project_comm_trt_vector)){ 
  subset=filter(correDiffLong, site_project_comm_trt == site_project_comm_trt_vector[PROJ])
  model=lm(anpp_pdiff~poly(treatment_year,2), data=subset)
  cf=as.data.frame(coef(summary(model)))%>%
    mutate(component=c('intercept','linear','quadratic'),
           r_sq=summary(model)$r.squared,
           site_project_comm_trt=site_project_comm_trt_vector[PROJ])
  anppTimeModels=rbind(cf, anppTimeModels) 
}
anppTimeModels <- anppTimeModels%>%
  rename(estimate=1,
         std_error=2,
         t=3,
         p=4)%>%
  select(site_project_comm_trt, r_sq, component, estimate, std_error, p)%>%
  mutate(weight=1/std_error)%>%
  mutate(significant=ifelse(p<0.05, 'yes', 'no'))%>%
  pivot_wider(names_from=component, values_from = c(estimate, std_error, weight, p, significant))
  # %>%
  # mutate(pattern=ifelse(significant_intercept=='no'&significant_slope=='yes'&estimate_slope>0, 'increase through time',
  #                       ifelse(significant_intercept=='no'&significant_slope=='yes'&estimate_slope<0, 'decrease through time',
  #                              ifelse(significant_intercept=='yes'&significant_slope=='no'&estimate_intercept>0, 'consistently higher',
  #                                     ifelse(significant_intercept=='yes'&significant_slope=='no'&estimate_intercept<0, 'consistently lower',
  #                                            ifelse(significant_intercept=='yes'&significant_slope=='yes'&estimate_intercept>0&estimate_slope>0, 'initially higher and increasing',
  #                                                   ifelse(significant_intercept=='yes'&significant_slope=='yes'&estimate_intercept>0&estimate_slope<0, 'initially higher and decreasing',
  #                                                          ifelse(significant_intercept=='yes'&significant_slope=='yes'&estimate_intercept<0&estimate_slope<0, 'initially lower and decreasing',
  #                                                                 ifelse(significant_intercept=='yes'&significant_slope=='yes'&estimate_intercept<0&estimate_slope>0, 'initially lower and increasing', 'no pattern')))))))))
colnames(anppTimeModels) <-paste(colnames(anppTimeModels), 'anpp', sep='_')
anppTimeModels <- anppTimeModels%>%
  rename(site_project_comm_trt=site_project_comm_trt_anpp)

# #plots of all ANPP through time
# for(PROJ in 1:length(site_project_comm_trt_vector)){
#   ggplot(data=filter(correDiffLong, site_project_comm == site_project_comm_vector[PROJ]),
#          aes(x=treatment_year, y=anpp_pdiff, color=treatment)) +
#     geom_point() +
#     geom_smooth(method='lm', se=F, formula=y~poly(x,2)) +
#     ggtitle(site_project_comm_vector[PROJ]) +
#     theme_bw() +
#     facet_wrap(~treatment)
#   ggsave(filename=paste0("C:\\Users\\lapie\\Desktop\\anpp figs\\",
#                          site_project_comm_vector[PROJ], "_anpp.png"))
# }

##### community difference over time #####
commTimeModels=data.frame(row.names=1) 
for(PROJ in 1:length(site_project_comm_trt_vector)){ 
  subset=filter(correDiffLong, site_project_comm_trt == site_project_comm_trt_vector[PROJ])
  model=lm(composition_diff~poly(treatment_year,2), data=subset)
  cf=as.data.frame(coef(summary(model)))%>%
    mutate(component=c('intercept','linear','quadratic'),
           r_sq=summary(model)$r.squared,
           site_project_comm_trt=site_project_comm_trt_vector[PROJ])
  commTimeModels=rbind(cf, commTimeModels) 
}
commTimeModels <- commTimeModels%>%
  rename(estimate=1,
         std_error=2,
         t=3,
         p=4)%>%
  select(site_project_comm_trt, r_sq, component, estimate, std_error, p)%>%
  mutate(weight=1/std_error)%>%
  mutate(significant=ifelse(p<0.05, 'yes', 'no'))%>%
  pivot_wider(names_from=component, values_from = c(estimate, std_error, weight, p, significant))
  # %>%
  # mutate(pattern=ifelse(significant_intercept=='no'&significant_slope=='yes'&estimate_slope>0, 'increase through time',
  #                       ifelse(significant_intercept=='no'&significant_slope=='yes'&estimate_slope<0, 'decrease through time',
  #                              ifelse(significant_intercept=='yes'&significant_slope=='no'&estimate_intercept>0, 'consistently higher',
  #                                     ifelse(significant_intercept=='yes'&significant_slope=='no'&estimate_intercept<0, 'consistently lower',
  #                                            ifelse(significant_intercept=='yes'&significant_slope=='yes'&estimate_intercept>0&estimate_slope>0, 'initially higher and increasing',
  #                                                   ifelse(significant_intercept=='yes'&significant_slope=='yes'&estimate_intercept>0&estimate_slope<0, 'initially higher and decreasing',
  #                                                          ifelse(significant_intercept=='yes'&significant_slope=='yes'&estimate_intercept<0&estimate_slope<0, 'initially lower and decreasing',
  #                                                                 ifelse(significant_intercept=='yes'&significant_slope=='yes'&estimate_intercept<0&estimate_slope>0, 'initially lower and increasing', 'no pattern')))))))))
colnames(commTimeModels) <-paste(colnames(commTimeModels), 'comm', sep='_')
commTimeModels <- commTimeModels%>%
  rename(site_project_comm_trt=site_project_comm_trt_comm)

# #graphs for each site
# for(PROJ in 1:length(site_project_comm_vector)){
#   ggplot(data=filter(correDiffLong, site_project_comm == site_project_comm_vector[PROJ]),
#          aes(x=treatment_year, y=composition_diff, color=treatment)) +
#     geom_point() +
#     geom_smooth(method='lm', se=F, formula=y~poly(x,2)) +
#     ggtitle(site_project_comm_vector[PROJ]) +
#     theme_bw() +
#     facet_wrap(~treatment)
#   ggsave(filename=paste0("C:\\Users\\lapie\\Desktop\\comm figs\\",
#                          site_project_comm_vector[PROJ], "_comm.png"))
# }

##### richness difference over time #####
richTimeModels=data.frame(row.names=1) 
for(PROJ in 1:length(site_project_comm_trt_vector)){ 
  subset=filter(correDiffLong, site_project_comm_trt == site_project_comm_trt_vector[PROJ])
  model=lm(richness_diff~poly(treatment_year,2), data=subset)
  cf=as.data.frame(coef(summary(model)))%>%
    mutate(component=c('intercept','linear','quadratic'),
           r_sq=summary(model)$r.squared,
           site_project_comm_trt=site_project_comm_trt_vector[PROJ])
  richTimeModels=rbind(cf, richTimeModels) 
}
richTimeModels <- richTimeModels%>%
  rename(estimate=1,
         std_error=2,
         t=3,
         p=4)%>%
  select(site_project_comm_trt, r_sq, component, estimate, std_error, p)%>%
  mutate(weight=1/std_error)%>%
  mutate(significant=ifelse(p<0.05, 'yes', 'no'))%>%
  pivot_wider(names_from=component, values_from = c(estimate, std_error, weight, p, significant))
  # %>%
  # mutate(pattern=ifelse(significant_intercept=='no'&significant_slope=='yes'&estimate_slope>0, 'increase through time',
  #                       ifelse(significant_intercept=='no'&significant_slope=='yes'&estimate_slope<0, 'decrease through time',
  #                              ifelse(significant_intercept=='yes'&significant_slope=='no'&estimate_intercept>0, 'consistently higher',
  #                                     ifelse(significant_intercept=='yes'&significant_slope=='no'&estimate_intercept<0, 'consistently lower',
  #                                            ifelse(significant_intercept=='yes'&significant_slope=='yes'&estimate_intercept>0&estimate_slope>0, 'initially higher and increasing',
  #                                                   ifelse(significant_intercept=='yes'&significant_slope=='yes'&estimate_intercept>0&estimate_slope<0, 'initially higher and decreasing',
  #                                                          ifelse(significant_intercept=='yes'&significant_slope=='yes'&estimate_intercept<0&estimate_slope<0, 'initially lower and decreasing',
  #                                                                 ifelse(significant_intercept=='yes'&significant_slope=='yes'&estimate_intercept<0&estimate_slope>0, 'initially lower and increasing', 'no pattern')))))))))
colnames(richTimeModels) <-paste(colnames(richTimeModels), 'rich', sep='_')
richTimeModels <- richTimeModels%>%
  rename(site_project_comm_trt=site_project_comm_trt_rich)

# #plots for each site
# for(PROJ in 1:length(site_project_comm_vector)){
#   ggplot(data=filter(correDiffLong, site_project_comm == site_project_comm_vector[PROJ]),
#          aes(x=treatment_year, y=richness_diff, color=treatment)) +
#     geom_point() +
#     geom_smooth(method='lm', se=F, formula=y~poly(x,2)) +
#     ggtitle(site_project_comm_vector[PROJ]) +
#     theme_bw()
#   ggsave(filename=paste0("C:\\Users\\lapie\\Desktop\\richness figs\\",
#                          site_project_comm_vector[PROJ], "_rich.png"))
# }

##### evenness difference over time #####
evenTimeModels=data.frame(row.names=1) 
for(PROJ in 1:length(site_project_comm_trt_vector)){ 
  subset=filter(correDiffLong, site_project_comm_trt == site_project_comm_trt_vector[PROJ]) #breaks at SERC because evenness=0 for all plots in SERC_TMECE_SP_E
  model=lm(evenness_diff~poly(treatment_year,2), data=subset)
  cf=as.data.frame(coef(summary(model)))%>%
    mutate(component=c('intercept','linear','quadratic'),
           r_sq=summary(model)$r.squared,
           site_project_comm_trt=site_project_comm_trt_vector[PROJ])
  evenTimeModels=rbind(cf, evenTimeModels) 
}
evenTimeModels <- evenTimeModels%>%
  rename(estimate=1,
         std_error=2,
         t=3,
         p=4)%>%
  select(site_project_comm_trt, r_sq, component, estimate, std_error, p)%>%
  mutate(weight=1/std_error)%>%
  mutate(significant=ifelse(p<0.05, 'yes', 'no'))%>%
  pivot_wider(names_from=component, values_from = c(estimate, std_error, weight, p, significant))
  # %>%
  # mutate(pattern=ifelse(significant_intercept=='no'&significant_slope=='yes'&estimate_slope>0, 'increase through time',
  #                       ifelse(significant_intercept=='no'&significant_slope=='yes'&estimate_slope<0, 'decrease through time',
  #                              ifelse(significant_intercept=='yes'&significant_slope=='no'&estimate_intercept>0, 'consistently higher',
  #                                     ifelse(significant_intercept=='yes'&significant_slope=='no'&estimate_intercept<0, 'consistently lower',
  #                                            ifelse(significant_intercept=='yes'&significant_slope=='yes'&estimate_intercept>0&estimate_slope>0, 'initially higher and increasing',
  #                                                   ifelse(significant_intercept=='yes'&significant_slope=='yes'&estimate_intercept>0&estimate_slope<0, 'initially higher and decreasing',
  #                                                          ifelse(significant_intercept=='yes'&significant_slope=='yes'&estimate_intercept<0&estimate_slope<0, 'initially lower and decreasing',
  #                                                                 ifelse(significant_intercept=='yes'&significant_slope=='yes'&estimate_intercept<0&estimate_slope>0, 'initially lower and increasing', 'no pattern')))))))))
colnames(evenTimeModels) <-paste(colnames(evenTimeModels), 'even', sep='_')
evenTimeModels <- evenTimeModels%>%
  rename(site_project_comm_trt=site_project_comm_trt_even)

# #plotting all sites
# for(PROJ in 1:length(site_project_comm_vector)){
#   ggplot(data=filter(correDiffLong, site_project_comm == site_project_comm_vector[PROJ]),
#          aes(x=treatment_year, y=evenness_diff, color=treatment)) +
#     geom_point() +
#     geom_smooth(method='lm', se=F, formula=y~poly(x,2)) +
#     ggtitle(site_project_comm_vector[PROJ]) +
#     theme_bw()
#   ggsave(filename=paste0("C:\\Users\\lapie\\Desktop\\evenness figs\\",
#                          site_project_comm_vector[PROJ], "_even.png"))
# }

##### rank difference over time #####
rankTimeModels=data.frame(row.names=1) 
for(PROJ in 1:length(site_project_comm_trt_vector)){ 
  subset=filter(correDiffLong, site_project_comm_trt == site_project_comm_trt_vector[PROJ])
  model=lm(rank_diff~poly(treatment_year,2), data=subset)
  cf=as.data.frame(coef(summary(model)))%>%
    mutate(component=c('intercept','linear','quadratic'),
           r_sq=summary(model)$r.squared,
           site_project_comm_trt=site_project_comm_trt_vector[PROJ])
  rankTimeModels=rbind(cf, rankTimeModels) 
}
rankTimeModels <- rankTimeModels%>%
  rename(estimate=1,
         std_error=2,
         t=3,
         p=4)%>%
  select(site_project_comm_trt, r_sq, component, estimate, std_error, p)%>%
  mutate(weight=1/std_error)%>%
  mutate(significant=ifelse(p<0.05, 'yes', 'no'))%>%
  pivot_wider(names_from=component, values_from = c(estimate, std_error, weight, p, significant))
  # %>%
  # mutate(pattern=ifelse(significant_intercept=='no'&significant_slope=='yes'&estimate_slope>0, 'increase through time',
  #                       ifelse(significant_intercept=='no'&significant_slope=='yes'&estimate_slope<0, 'decrease through time',
  #                              ifelse(significant_intercept=='yes'&significant_slope=='no'&estimate_intercept>0, 'consistently higher',
  #                                     ifelse(significant_intercept=='yes'&significant_slope=='no'&estimate_intercept<0, 'consistently lower',
  #                                            ifelse(significant_intercept=='yes'&significant_slope=='yes'&estimate_intercept>0&estimate_slope>0, 'initially higher and increasing',
  #                                                   ifelse(significant_intercept=='yes'&significant_slope=='yes'&estimate_intercept>0&estimate_slope<0, 'initially higher and decreasing',
  #                                                          ifelse(significant_intercept=='yes'&significant_slope=='yes'&estimate_intercept<0&estimate_slope<0, 'initially lower and decreasing',
  #                                                                 ifelse(significant_intercept=='yes'&significant_slope=='yes'&estimate_intercept<0&estimate_slope>0, 'initially lower and increasing', 'no pattern')))))))))
colnames(rankTimeModels) <-paste(colnames(rankTimeModels), 'rank', sep='_')
rankTimeModels <- rankTimeModels%>%
  rename(site_project_comm_trt=site_project_comm_trt_rank)

# #plotting all sites
# for(PROJ in 1:length(site_project_comm_vector)){
#   ggplot(data=filter(correDiffLong, site_project_comm == site_project_comm_vector[PROJ]),
#          aes(x=treatment_year, y=rank_diff, color=treatment)) +
#     geom_point() +
#     geom_smooth(method='lm', se=F, formula=y~poly(x,2)) +
#     ggtitle(site_project_comm_vector[PROJ]) +
#     theme_bw()
#   ggsave(filename=paste0("C:\\Users\\lapie\\Desktop\\rank figs\\",
#                          site_project_comm_vector[PROJ], "_rank.png"))
# }

##### spp difference over time #####
sppDiffTimeModels=data.frame(row.names=1) 
for(PROJ in 1:length(site_project_comm_trt_vector)){ 
  subset=filter(correDiffLong, site_project_comm_trt == site_project_comm_trt_vector[PROJ])
  model=lm(species_diff~poly(treatment_year,2), data=subset)
  cf=as.data.frame(coef(summary(model)))%>%
    mutate(component=c('intercept','linear','quadratic'),
           r_sq=summary(model)$r.squared,
           site_project_comm_trt=site_project_comm_trt_vector[PROJ])
  sppDiffTimeModels=rbind(cf, sppDiffTimeModels) 
}
sppDiffTimeModels <- sppDiffTimeModels%>%
  rename(estimate=1,
         std_error=2,
         t=3,
         p=4)%>%
  select(site_project_comm_trt, r_sq, component, estimate, std_error, p)%>%
  mutate(weight=1/std_error)%>%
  mutate(significant=ifelse(p<0.05, 'yes', 'no'))%>%
  pivot_wider(names_from=component, values_from = c(estimate, std_error, weight, p, significant))
  # %>%
  # mutate(pattern=ifelse(significant_intercept=='no'&significant_slope=='yes'&estimate_slope>0, 'increase through time',
  #                       ifelse(significant_intercept=='no'&significant_slope=='yes'&estimate_slope<0, 'decrease through time',
  #                              ifelse(significant_intercept=='yes'&significant_slope=='no'&estimate_intercept>0, 'consistently higher',
  #                                     ifelse(significant_intercept=='yes'&significant_slope=='no'&estimate_intercept<0, 'consistently lower',
  #                                            ifelse(significant_intercept=='yes'&significant_slope=='yes'&estimate_intercept>0&estimate_slope>0, 'initially higher and increasing',
  #                                                   ifelse(significant_intercept=='yes'&significant_slope=='yes'&estimate_intercept>0&estimate_slope<0, 'initially higher and decreasing',
  #                                                          ifelse(significant_intercept=='yes'&significant_slope=='yes'&estimate_intercept<0&estimate_slope<0, 'initially lower and decreasing',
  #                                                                 ifelse(significant_intercept=='yes'&significant_slope=='yes'&estimate_intercept<0&estimate_slope>0, 'initially lower and increasing', 'no pattern')))))))))
colnames(sppDiffTimeModels) <-paste(colnames(sppDiffTimeModels), 'spp_diff', sep='_')
sppDiffTimeModels <- sppDiffTimeModels%>%
  rename(site_project_comm_trt=site_project_comm_trt_spp_diff)

# #plotting all sites
# for(PROJ in 1:length(site_project_comm_vector)){
#   ggplot(data=filter(correDiffLong, site_project_comm == site_project_comm_vector[PROJ]),
#          aes(x=treatment_year, y=species_diff, color=treatment)) +
#     geom_point() +
#     geom_smooth(method='lm', se=F, formula=y~poly(x,2)) +
#     ggtitle(site_project_comm_vector[PROJ]) +
#     theme_bw()
#   ggsave(filename=paste0("C:\\Users\\lapie\\Desktop\\spp diff figs\\",
#                          site_project_comm_vector[PROJ], "_sppDiff.png"))
# }


##### putting together trends through time #####
allTimeModels <- anppTimeModels%>%
  left_join(commTimeModels)%>%
  left_join(richTimeModels)%>%
  left_join(evenTimeModels)%>%
  left_join(rankTimeModels)%>%
  left_join(sppDiffTimeModels)%>%
  left_join(correDiffDetails)%>%
  mutate(linear_anpp_modifier=ifelse(site_project_comm_trt %in% c('Bt_DroughtNet_0_drought', 'Bt_EVENT2_0_D1-N1', 'Bt_EVENT2_0_D2-N1'), -1, 1),
         estimate_linear_anpp_alt=estimate_linear_anpp*linear_anpp_modifier)


##### comparing anpp diff across years 1-3, 4-6, and 7-9 of exp #####
diffByYear <- correDiffLong%>%
  mutate(year_set=ifelse(treatment_year %in% c(1,2,3), 'yrs_1to3',
                         ifelse(treatment_year %in% c(4,5,6), 'yrs_4to6',
                                ifelse(treatment_year %in% c(7,8,9), 'yrs_7to9', 'yrs10plus'))))

ggplot(data=diffByYear, aes(x=year_set, y=anpp_pdiff)) +
  geom_boxplot()
ggplot(data=diffByYear, aes(x=as.factor(treatment_year), y=anpp_pdiff)) +
  geom_boxplot()


#try looking at max diff vs first and last yr diff
do this

##### comparing anpp slope to community differences -- some go up and some go down with increasing community change #####
ggplot(data=allTimeModels, aes(x=estimate_linear_comm, y=estimate_linear_anpp_alt, color=MAP)) +
  geom_point(size=4) +
  geom_smooth(method='lm', formula=y~poly(x,2)) +
  # geom_text(hjust=0, vjust=0, size=4) +
  geom_hline(yintercept=0) + geom_vline(xintercept=0)
summary(lm(estimate_linear_anpp_alt~poly(estimate_linear_comm,2), weights=weight_linear_anpp, data=allTimeModels))


##### anpp trends with site predictors -- more positive anpp responses at sites with greater MAP #####
ggplot(data=subset(allTimeModels, data_points>4), aes(x=MAP, y=estimate_slope_anpp_alt)) + geom_point() + geom_smooth(method='lm', formula=y~poly(x,2))
summary(model <- lm(estimate_slope_anpp_alt~poly(MAP,2), weights=weight_slope_anpp, data=subset(allTimeModels, data_points>4))) #negative quadratic trend
ggplot(data=subset(allTimeModels, data_points>4), aes(x=MAT, y=estimate_slope_anpp_alt)) + geom_point() + geom_smooth(method='lm', formula=y~poly(x,2))
summary(model <- lm(estimate_slope_anpp_alt~poly(MAT,2), weights=weight_slope_anpp, data=subset(allTimeModels, data_points>4))) #no trend
ggplot(data=subset(allTimeModels, data_points>4), aes(x=rrich, y=estimate_slope_anpp_alt)) + geom_point() + geom_smooth(method='lm', formula=y~poly(x,2))
summary(model <- lm(estimate_slope_anpp_alt~poly(rrich,2), weights=weight_slope_anpp, data=subset(allTimeModels, data_points>4))) #positive trend
ggplot(data=subset(allTimeModels, data_points>4), aes(x=anpp, y=estimate_slope_anpp_alt)) + geom_point() + geom_smooth(method='lm', formula=y~poly(x,2))
summary(model <- lm(estimate_slope_anpp_alt~poly(anpp,2), weights=weight_slope_anpp, data=subset(allTimeModels, data_points>4))) #negative quadratic trend

#community trends with site predictors -- more positive community responses at sites with higher MAP/anpp
ggplot(data=subset(allTimeModels, data_points>4), aes(x=MAP, y=estimate_slope_comm)) + geom_point() + geom_smooth(method='lm', formula=y~poly(x,2))
summary(model <- lm(estimate_slope_comm~poly(MAP,2), weights=weight_slope_comm, data=subset(allTimeModels, data_points>4))) #no trend
ggplot(data=subset(allTimeModels, data_points>4), aes(x=MAT, y=estimate_slope_comm)) + geom_point() + geom_smooth(method='lm', formula=y~poly(x,2))
summary(model <- lm(estimate_slope_comm~poly(MAT,2), weights=weight_slope_comm, data=subset(allTimeModels, data_points>4))) #positive quadratic trend
ggplot(data=subset(allTimeModels, data_points>4), aes(x=rrich, y=estimate_slope_comm)) + geom_point() + geom_smooth(method='lm', formula=y~poly(x,2))
summary(model <- lm(estimate_slope_comm~poly(rrich,2), weights=weight_slope_comm, data=subset(allTimeModels, data_points>4))) #negative quadratic trend
ggplot(data=subset(allTimeModels, data_points>4), aes(x=anpp, y=estimate_slope_comm)) + geom_point() + geom_smooth(method='lm', formula=y~poly(x,2))
summary(model <- lm(estimate_slope_comm~poly(anpp,2), weights=weight_slope_comm, data=subset(allTimeModels, data_points>4))) #no trend

 
#comparing anpp slope to different metrics of community difference
rich <- ggplot(data=allTimeModels, aes(x=estimate_slope_rich, y=estimate_slope_anpp_alt, label=site_project_comm_trt)) +
  geom_point() +
  geom_hline(yintercept=0) + geom_vline(xintercept=0)
even <- ggplot(data=allTimeModels, aes(x=estimate_slope_even, y=estimate_slope_anpp_alt, label=site_project_comm_trt)) +
  geom_point() +
  geom_hline(yintercept=0) + geom_vline(xintercept=0)
rank <- ggplot(data=allTimeModels, aes(x=estimate_slope_rank, y=estimate_slope_anpp_alt, label=site_project_comm_trt)) +
  geom_point() +
  geom_hline(yintercept=0) + geom_vline(xintercept=0)
sppDiff <- ggplot(data=allTimeModels, aes(x=estimate_slope_spp_diff, y=estimate_slope_anpp_alt, label=site_project_comm_trt)) +
  geom_point() +
  geom_hline(yintercept=0) + geom_vline(xintercept=0)

pushViewport(viewport(layout=grid.layout(2,2)))
print(rich, vp=viewport(layout.pos.row=1, layout.pos.col=1))
print(even, vp=viewport(layout.pos.row=2, layout.pos.col=1))
print(rank, vp=viewport(layout.pos.row=1, layout.pos.col=2))
print(sppDiff, vp=viewport(layout.pos.row=2, layout.pos.col=2))




###### anpp vs composition diff #####
anppCompModels=data.frame(row.names=1) 
for(PROJ in 1:length(site_project_comm_trt_vector)){ 
  subset=filter(correDiff, site_project_comm_trt == site_project_comm_trt_vector[PROJ])
  model=lm(anpp_pdiff~composition_diff, data=subset)
  cf=as.data.frame(coef(summary(model)))%>%
    mutate(component=c('intercept','slope'),
           r_sq=summary(model)$r.squared,
           site_project_comm_trt=site_project_comm_trt_vector[PROJ])
  anppCompModels=rbind(cf, anppCompModels) 
}

for(PROJ in 1:length(site_project_comm_vector)){
  ggplot(data=filter(correDiff, site_project_comm == site_project_comm_vector[PROJ]),
         aes(x=composition_diff, y=anpp_pdiff, color=treatment)) +
    geom_point() +
    geom_smooth(method='lm', se=F) +
    ggtitle(site_project_comm_vector[PROJ]) +
    theme_bw()
  ggsave(filename=paste0("C:\\Users\\lapie\\Desktop\\anpp v comp figs\\",
                         site_project_comm_vector[PROJ], "_anppComp.png"))
}

#anpp vs richness diff
anppRichModels=data.frame(row.names=1) 
for(PROJ in 1:length(site_project_comm_trt_vector)){ 
  subset=filter(correDiff, site_project_comm_trt == site_project_comm_trt_vector[PROJ])
  model=lm(anpp_pdiff~richness_diff, data=subset)
  cf=as.data.frame(coef(summary(model)))%>%
    mutate(component=c('intercept','slope'),
           r_sq=summary(model)$r.squared,
           site_project_comm_trt=site_project_comm_trt_vector[PROJ])
  anppRichModels=rbind(cf, anppRichModels) 
}

for(PROJ in 1:length(site_project_comm_vector)){
  ggplot(data=filter(correDiff, site_project_comm == site_project_comm_vector[PROJ]),
         aes(x=richness_diff, y=anpp_pdiff, color=treatment)) +
    geom_point() +
    geom_smooth(method='lm', se=F) +
    ggtitle(site_project_comm_vector[PROJ]) +
    theme_bw()
  ggsave(filename=paste0("C:\\Users\\lapie\\Desktop\\anpp v rich figs\\",
                         site_project_comm_vector[PROJ], "_anppRich.png"))
}

#anpp vs evenness diff
anppEvenModels=data.frame(row.names=1) 
for(PROJ in 1:length(site_project_comm_trt_vector)){ 
  subset=filter(correDiff, site_project_comm_trt == site_project_comm_trt_vector[PROJ]) #fails at SERC TMECE because evenness=0
  model=lm(anpp_pdiff~evenness_diff, data=subset)
  cf=as.data.frame(coef(summary(model)))%>%
    mutate(component=c('intercept','slope'),
           r_sq=summary(model)$r.squared,
           site_project_comm_trt=site_project_comm_trt_vector[PROJ])
  anppEvenModels=rbind(cf, anppEvenModels) 
}

for(PROJ in 1:length(site_project_comm_vector)){
  ggplot(data=filter(correDiff, site_project_comm == site_project_comm_vector[PROJ]),
         aes(x=evenness_diff, y=anpp_pdiff, color=treatment)) +
    geom_point() +
    geom_smooth(method='lm', se=F) +
    ggtitle(site_project_comm_vector[PROJ]) +
    theme_bw()
  ggsave(filename=paste0("C:\\Users\\lapie\\Desktop\\anpp v even figs\\",
                         site_project_comm_vector[PROJ], "_anppEven.png"))
}

#anpp vs rank diff
anppRankModels=data.frame(row.names=1) 
for(PROJ in 1:length(site_project_comm_trt_vector)){ 
  subset=filter(correDiff, site_project_comm_trt == site_project_comm_trt_vector[PROJ])
  model=lm(anpp_pdiff~rank_diff, data=subset)
  cf=as.data.frame(coef(summary(model)))%>%
    mutate(component=c('intercept','slope'),
           r_sq=summary(model)$r.squared,
           site_project_comm_trt=site_project_comm_trt_vector[PROJ])
  anppRankModels=rbind(cf, anppRankModels) 
}

for(PROJ in 1:length(site_project_comm_vector)){
  ggplot(data=filter(correDiff, site_project_comm == site_project_comm_vector[PROJ]),
         aes(x=rank_diff, y=anpp_pdiff, color=treatment)) +
    geom_point() +
    geom_smooth(method='lm', se=F) +
    ggtitle(site_project_comm_vector[PROJ]) +
    theme_bw()
  ggsave(filename=paste0("C:\\Users\\lapie\\Desktop\\anpp v rank figs\\",
                         site_project_comm_vector[PROJ], "_anppRank.png"))
}

#anpp vs spp diff
anppSppDiffModels=data.frame(row.names=1) 
for(PROJ in 1:length(site_project_comm_trt_vector)){ 
  subset=filter(correDiff, site_project_comm_trt == site_project_comm_trt_vector[PROJ])
  model=lm(anpp_pdiff~species_diff, data=subset)
  cf=as.data.frame(coef(summary(model)))%>%
    mutate(component=c('intercept','slope'),
           r_sq=summary(model)$r.squared,
           site_project_comm_trt=site_project_comm_trt_vector[PROJ])
  anppSppDiffModels=rbind(cf, anppSppDiffModels) 
}

for(PROJ in 1:length(site_project_comm_vector)){
  ggplot(data=filter(correDiff, site_project_comm == site_project_comm_vector[PROJ]),
         aes(x=species_diff, y=anpp_pdiff, color=treatment)) +
    geom_point() +
    geom_smooth(method='lm', se=F) +
    ggtitle(site_project_comm_vector[PROJ]) +
    theme_bw()
  ggsave(filename=paste0("C:\\Users\\lapie\\Desktop\\anpp v spp diff figs\\",
                         site_project_comm_vector[PROJ], "_anppSpp.png"))
}

##### ANPP difference related to community difference #####
###overall community change
summary(anppCommDiff <- lmer(abs(anpp_pdiff)~composition_diff + (1+treatment_year|site_project_comm/treatment), data=correDiffMain)) #significant effect
r.squaredGLMM(anppCommDiff) #marginal R2=0.1634425 , conditional R2=0.6189774 
AIC(anppCommDiff) #polynomial model (AIC=960.5969) NOT better than linear model (AIC=962.7837)

ggplot(data=correDiffMain, aes(x=composition_diff, y=abs(anpp_pdiff))) +
    geom_point(color='grey') +
    geom_smooth(method='lm', se=F, color='grey',aes(group=interaction(site_project_comm, treatment))) +
    geom_smooth(method='lm', se=F, size=3, color='black') +
    xlab('Composition Difference') + ylab('|ANPP % Difference|') +
    theme(legend.position='none')

###richness only
summary(anppRichDiff <- lmer(abs(anpp_pdiff)~poly(richness_diff,2) + (1+treatment_year|site_project_comm/treatment), data=correDiffMain)) #significant effect
r.squaredGLMM(anppRichDiff) #marginal R2=0.1173921 , conditional R2=0.5799932 
AIC(anppRichDiff) #polynomial model (AIC=971.4683) better than linear model (AIC=975.0343)

ggplot(data=correDiffMain, aes(x=richness_diff, y=abs(anpp_pdiff))) +
  geom_point(color='grey') +
  geom_smooth(method='lm', formula=y~x+I(x^2), se=F, color='grey',aes(group=interaction(site_project_comm, treatment))) +
  geom_smooth(method='lm', formula=y~x+I(x^2), se=F, size=3, color='black') +
  xlab('Richness Difference') + ylab('|ANPP % Difference|') +
  theme(legend.position='none') +
  coord_cartesian(ylim=c(0,4))

###multiple change metrics
#individually to determine quadratic vs linear
AIC(anppEvenDiff <- lmer(abs(anpp_pdiff)~evenness_diff + (1+treatment_year|site_project_comm/treatment), data=correDiffMain)) #polynomial model (AIC=1385.726) NOT better than linear model (AIC=1386.233)
AIC(anppRankDiff <- lmer(abs(anpp_pdiff)~rank_diff + (1+treatment_year|site_project_comm/treatment), data=correDiffMain)) #polynomial model (AIC=1378.816) NOT better than linear model (AIC=1384.194)
AIC(anppGainDiff <- lmer(abs(anpp_pdiff)~species_diff + (1+treatment_year|site_project_comm/treatment), data=correDiffMain)) #polynomial model (AIC=1392.324) NOT better than linear model (AIC=1395.222)

#all together
summary(anppCommMetricDiff <- lmer(anpp_pdiff~evenness_diff+rank_diff+species_diff + (1+treatment_year|site_project_comm/treatment), data=correDiffMain)) #significant effect of rank diff
r.squaredGLMM(anppCommMetricDiff) #marginal R2=0.02157451, conditional R2=0.517346

ggplot(data=correDiffMain, aes(x=evenness_diff, y=anpp_pdiff)) +
  geom_point(color='grey') +
  geom_smooth(method='lm', se=F, color='grey',aes(group=interaction(site_project_comm, treatment))) +
  geom_smooth(method='lm', se=F, size=3, color='black') +
  xlab('Evenness Difference') + ylab('ANPP % Difference') +
  theme(legend.position='none') +
  coord_cartesian(ylim=c(0,4))

ggplot(data=correDiffMain, aes(x=rank_diff, y=anpp_pdiff)) +
  geom_point(color='grey') +
  geom_smooth(method='lm', se=F, color='grey',aes(group=interaction(site_project_comm, treatment))) +
  geom_smooth(method='lm', se=F, size=3, color='black') +
  xlab('Rank Difference') + ylab('ANPP % Difference') +
  theme(legend.position='none')

ggplot(data=correDiffMain, aes(x=species_diff, y=anpp_pdiff)) +
  geom_point(color='grey') +
  geom_smooth(method='lm', se=F, color='grey',aes(group=interaction(site_project_comm, treatment))) +
  geom_smooth(method='lm', se=F, size=3, color='black') +
  xlab('Species Difference') + ylab('ANPP % Difference') +
  theme(legend.position='none')


##### site drivers of ANPP difference #####
siteInfo <- read.csv('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\CoRRE\\converge_diverge\\datasets\\LongForm\\SiteExperimentDetails_March2019.csv')%>%
  select(-X)

correDiffSite <- correDiffMain%>%
  left_join(siteInfo)

summary(anppCommMetricDiff <- lmer(abs(anpp_pdiff) ~ evenness_diff + rank_diff + species_diff + rrich + (1+treatment_year|site_project_comm/treatment), data=correDiffSite)) #significant effects of rank diff and species diff
r.squaredGLMM(anppCommMetricDiff) #marginal R2=0.03223563  , conditional R2=0.7324602

summary(anppSiteDrivers <- lmer(abs(anpp_pdiff)~MAP+MAT+rrich+anpp+(1+treatment_year|site_project_comm/treatment), data=correDiffSite)) #no significant drivers

# ggplot(data=correDiffSite, aes(x=rrich, y=abs(anpp_pdiff))) +
#   geom_point(color='grey') +
#   geom_smooth(method='lm', se=F, color='grey',aes(group=interaction(site_project_comm, treatment))) +
#   geom_smooth(method='lm', se=F, size=3, color='black') +
#   xlab('Relative Richness') + ylab('|ANPP % Difference|') +
#   theme(legend.position='none') +
#   coord_cartesian(ylim=c(0,4))

# ggplot(data=correDiffSite, aes(x=anpp, y=abs(anpp_pdiff))) +
#   geom_point(color='grey') +
#   geom_smooth(method='lm', se=F, color='grey',aes(group=interaction(site_project_comm, treatment))) +
#   geom_smooth(method='lm', se=F, size=3, color='black') +
#   xlab('Site Productivity') + ylab('|ANPP % Difference|') +
#   theme(legend.position='none') +
#   coord_cartesian(ylim=c(0,4))

ggplot(data=correDiffSite, aes(x=rank_diff, y=abs(anpp_pdiff))) +
  # geom_point(color='grey') +
  geom_smooth(method='lm', formula=y~x, se=F, color='grey', aes(group=interaction(site_project_comm, treatment))) +
  geom_smooth(method='lm', se=F, size=3, color='black') +
  xlab('Rank Difference') + ylab('|ANPP % Difference|') +
  theme(legend.position='none') +
  coord_cartesian(ylim=c(0,2))


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

correANPPctl <- correCommChangePlotANPP%>%
  filter(trt_type=='control')%>%
  group_by(site_year)%>%
  summarise(anpp_ctl=mean(anpp))%>%
  ungroup()

correCommChangePlotANPPDiff <- correCommChangePlotANPP%>%
  filter(trt_type!='control')%>%
  left_join(correANPPctl)%>%
  mutate(anpp_diff=(anpp-anpp_ctl)/anpp_ctl)

#need to calculate change in anpp from year 0 as well to include in the model instead of raw anpp values
#also try including the treatment variables as interactors with the change variables? like N effects alone vs interating with richness, evenness etc to get at the physiological responses vs the community driven responses?
changePlotCausalModel <- feols(abs(anpp_diff) ~ richness_change + evenness_change + rank_change + gains + losses | site_plot + site_year, data=subset(correCommChangePlotANPPDiff, anpp<4000))
etable(changePlotCausalModel,
       cluster="site_plot")

# # significant effects of evenness change (negative) and rank change (positive) on ANPP
# with(subset(correCommChangePlotANPPDiff, anpp<4000), plot(anpp_diff, evenness_change))
# with(subset(correCommChangePlotANPPDiff, anpp<4000), plot(anpp_diff, rank_change))

ggplot(data=subset(correCommChangePlotANPPDiff, anpp<4000), aes(x=evenness_change, y=abs(anpp_diff))) +
  # geom_point(color='grey') +
  geom_smooth(method='lm', formula=y~x, color='grey', se=F, aes(group=interaction(site_project_comm, treatment))) +
  geom_smooth(method='lm', formula=y~x, se=F, size=3, color='black') +
  xlab('Evenness Change') + ylab('|ANPP % Difference|') +
  theme(legend.position='none') +
  coord_cartesian(ylim=c(0,4))

ggplot(data=subset(correCommChangePlotANPPDiff, anpp<4000), aes(x=rank_change, y=abs(anpp_diff))) +
  # geom_point(color='grey') +
  geom_smooth(method='lm', formula=y~x, color='grey', se=F, aes(group=interaction(site_project_comm, treatment))) +
  geom_smooth(method='lm', formula=y~x, se=F, size=3, color='black') +
  xlab('Rank Change') + ylab('|ANPP % Difference|') +
  theme(legend.position='none') +
  coord_cartesian(ylim=c(0,4))

ggplot(data=subset(correCommChangePlotANPPDiff, anpp<4000), aes(x=gains, y=abs(anpp_diff))) +
  # geom_point(color='grey') +
  geom_smooth(method='lm', formula=y~x, color='grey', se=F, aes(group=interaction(site_project_comm, treatment))) +
  geom_smooth(method='lm', formula=y~x, se=F, size=3, color='black') +
  xlab('Species Gains') + ylab('|ANPP % Difference|') +
  theme(legend.position='none') +
  coord_cartesian(ylim=c(0,2.5))

ggplot(data=subset(correCommChangePlotANPPDiff, anpp<4000), aes(x=losses, y=abs(anpp_diff))) +
  # geom_point(color='grey') +
  geom_smooth(method='lm', formula=y~x, color='grey', se=F, aes(group=interaction(site_project_comm, treatment))) +
  geom_smooth(method='lm', formula=y~x, se=F, size=3, color='black') +
  xlab('Species Losses') + ylab('|ANPP % Difference|') +
  theme(legend.position='none') +
  coord_cartesian(ylim=c(0,2.5))


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






