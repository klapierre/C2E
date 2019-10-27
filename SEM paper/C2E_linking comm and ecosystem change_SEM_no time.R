library(grid)
library(PerformanceAnalytics)
library(piecewiseSEM)
library(tidyverse)


#kim's desktop
setwd('C:\\Users\\la pierrek\\Dropbox (Smithsonian)\\working groups\\CoRRE\\converge_diverge\\datasets\\LongForm')

#kim's laptop
setwd('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\CoRRE\\converge_diverge\\datasets\\LongForm')


theme_set(theme_bw())
theme_update(axis.title.x=element_text(size=40, vjust=-0.35, margin=margin(t=15)), axis.text.x=element_text(size=34, color='black'),
             axis.title.y=element_text(size=40, angle=90, vjust=0.5, margin=margin(r=15)), axis.text.y=element_text(size=34, color='black'),
             plot.title = element_text(size=40, vjust=2),
             panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
             legend.title=element_blank(), legend.text=element_text(size=20))


###merge CoRRE and NutNet datasets
correSEMdataTrt <- read.csv('CoRRE_comm and anpp diff_07160219.csv')%>%
  filter(trt_type %in% c('N','P','N*P','mult_nutrient', 'drought', 'CO2', 'irr', 'N*CO2', 'N*irr'))%>%
  select(site_project_comm, treatment_year, treatment, n, p, k, drought, irr, CO2, anpp_pdiff, composition_diff, richness_difference, evenness_diff, rank_difference, species_difference)%>%
  mutate(dataset='CoRRE')

semData <- read.csv('NutNet_comm and anpp diff_07160219.csv')%>%
  filter(! treatment %in% c('Fence', 'NPK+Fence'))%>%
  rename(site_project_comm=site_code)%>%
  select(site_project_comm, treatment_year, treatment, n, p, k, anpp_pdiff, composition_diff, richness_difference, evenness_diff, rank_difference, species_difference)%>%
  mutate(n=ifelse(treatment %in% c('N', 'NP', 'NK', 'NPK'), 10, 0), p=ifelse(treatment %in% c('P', 'NP', 'PK', 'NPK'), 10, 0), k=ifelse(treatment %in% c('K', 'NK', 'PK', 'NPK'), 10, 0))%>%
  mutate(drought=0, irr=0, CO2=0)%>%
  mutate(dataset='NutNet')%>%
  rbind(correSEMdataTrt)%>%
  separate(col=site_project_comm, into=c('site_code', 'project_name', 'community_type'), sep='::', remove=F)%>%
  mutate(anpp_pdiff_transform=log(anpp_pdiff+(1-min(anpp_pdiff))), composition_diff_transform=log(composition_diff))

dataVis <- (semData)%>%
  select(anpp_pdiff_transform, composition_diff_transform, richness_difference, evenness_diff, rank_difference, species_difference) #make visualization dataframe
chart.Correlation(dataVis, histogram=T, pch=19)





#######MODEL STRUCTURE: use differences between treatment and control plots in each year (horizontal arrows)-------------------

###all data, all years-----------

###difference metrics not through composition model
#all years
summary(compositionModel <- psem(
  lm(anpp_pdiff_transform ~ n + p + k + drought + irr + CO2 + evenness_diff + rank_difference + species_difference, data=semData),
  lm(evenness_diff ~ n + p + k + drought + irr + CO2, data=semData),
  lm(rank_difference ~ n + p + k + drought + irr + CO2, data=semData),
  lm(species_difference ~ n + p + k + drought + irr + CO2, data=semData),
  evenness_diff %~~% rank_difference,
  evenness_diff %~~% species_difference,
  rank_difference %~~% species_difference,
  data=semData
))
coefs1 <- coefs(compositionModel, standardize = "scale", standardize.type = "latent.linear", intercepts = FALSE)%>%
  select(Response, Predictor, Estimate, Std.Error, DF, Crit.Value, P.Value, Std.Estimate)


#all years - no richness
summary(compositionModel <- psem(
  lm(anpp_pdiff_transform ~ n + p + k + drought + irr + CO2 + evenness_diff + rank_difference + species_difference, data=semData),
  lm(richness_difference ~ n + p + k + drought + irr + CO2, data=semData),
  lm(evenness_diff ~ n + p + k + drought + irr + CO2, data=semData),
  lm(rank_difference ~ n + p + k + drought + irr + CO2, data=semData),
  lm(species_difference ~ n + p + k + drought + irr + CO2, data=semData),
  richness_difference %~~% evenness_diff,
  richness_difference %~~% rank_difference,
  richness_difference %~~% species_difference,
  evenness_diff %~~% rank_difference,
  evenness_diff %~~% species_difference,
  rank_difference %~~% species_difference,
  data=semData
))
coefs1 <- coefs(compositionModel, standardize = "scale", standardize.type = "latent.linear", intercepts = FALSE)%>%
  select(Response, Predictor, Estimate, Std.Error, DF, Crit.Value, P.Value, Std.Estimate)
