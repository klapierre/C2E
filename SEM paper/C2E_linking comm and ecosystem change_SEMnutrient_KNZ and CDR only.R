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
  select(site_project_comm, treatment_year, treatment, n, p, k, anpp_pdiff, composition_diff, richness_difference, evenness_diff, rank_difference, species_difference)%>%
  mutate(dataset='CoRRE')

CDRsem <- read.csv('NutNet_comm and anpp diff_07160219.csv')%>%
  filter(! treatment %in% c('Fence', 'NPK+Fence'))%>%
  rename(site_project_comm=site_code)%>%
  select(site_project_comm, treatment_year, treatment, n, p, k, anpp_pdiff, composition_diff, richness_difference, evenness_diff, rank_difference, species_difference)%>%
  mutate(n=ifelse(treatment %in% c('N', 'NP', 'NK', 'NPK'), 10, 0), p=ifelse(treatment %in% c('P', 'NP', 'PK', 'NPK'), 10, 0), k=ifelse(treatment %in% c('K', 'NK', 'PK', 'NPK'), 10, 0))%>%
  mutate(dataset='NutNet')%>%
  rbind(correSEMdataTrt)%>%
  separate(col=site_project_comm, into=c('site_code', 'project_name', 'community_type'), sep='::', remove=F)%>%
  mutate(anpp_pdiff_transform=log(anpp_pdiff+(1-min(anpp_pdiff))), composition_diff_transform=log(composition_diff))%>%
  filter(site_code %in% c('CDR', 'cdcr.us'))

dataVisCDR <- (CDRsem)%>%
  select(anpp_pdiff_transform, composition_diff_transform, richness_difference, evenness_diff, rank_difference, species_difference) #make visualization dataframe
chart.Correlation(dataVisCDR, histogram=T, pch=19)


KNZsem <- read.csv('NutNet_comm and anpp diff_07160219.csv')%>%
  filter(! treatment %in% c('Fence', 'NPK+Fence'))%>%
  rename(site_project_comm=site_code)%>%
  select(site_project_comm, treatment_year, treatment, n, p, k, anpp_pdiff, composition_diff, richness_difference, evenness_diff, rank_difference, species_difference)%>%
  mutate(n=ifelse(treatment %in% c('N', 'NP', 'NK', 'NPK'), 10, 0), p=ifelse(treatment %in% c('P', 'NP', 'PK', 'NPK'), 10, 0), k=ifelse(treatment %in% c('K', 'NK', 'PK', 'NPK'), 10, 0))%>%
  mutate(dataset='NutNet')%>%
  rbind(correSEMdataTrt)%>%
  separate(col=site_project_comm, into=c('site_code', 'project_name', 'community_type'), sep='::', remove=F)%>%
  mutate(anpp_pdiff_transform=log(anpp_pdiff+(1-min(anpp_pdiff))), composition_diff_transform=log(composition_diff))%>%
  filter(site_code %in% c('KNZ', 'konz.us'))

dataVisKNZ <- (allSEMdata)%>%
  select(anpp_pdiff_transform, composition_diff_transform, richness_difference, evenness_diff, rank_difference, species_difference) #make visualization dataframe
chart.Correlation(dataVisKNZ, histogram=T, pch=19)
  


correSEMdataTrt <- read.csv('CoRRE_comm and anpp diff_07160219.csv')%>%
  filter(trt_type %in% c('N','P','N*P','mult_nutrient'))%>%
  select(site_project_comm, treatment_year, treatment, n, p, k, anpp_pdiff, composition_diff, richness_difference, evenness_diff, rank_difference, species_difference)%>%
  mutate(dataset='CoRRE')

tallgrassSEM <- read.csv('NutNet_comm and anpp diff_07160219.csv')%>%
  filter(! treatment %in% c('Fence', 'NPK+Fence'))%>%
  rename(site_project_comm=site_code)%>%
  select(site_project_comm, treatment_year, treatment, n, p, k, anpp_pdiff, composition_diff, richness_difference, evenness_diff, rank_difference, species_difference)%>%
  mutate(n=ifelse(treatment %in% c('N', 'NP', 'NK', 'NPK'), 10, 0), p=ifelse(treatment %in% c('P', 'NP', 'PK', 'NPK'), 10, 0), k=ifelse(treatment %in% c('K', 'NK', 'PK', 'NPK'), 10, 0))%>%
  mutate(dataset='NutNet')%>%
  rbind(correSEMdataTrt)%>%
  separate(col=site_project_comm, into=c('site_code', 'project_name', 'community_type'), sep='::', remove=F)%>%
  mutate(anpp_pdiff_transform=log(anpp_pdiff+(1-min(anpp_pdiff))), composition_diff_transform=log(composition_diff))%>%
  filter(site_code %in% c('CDR', 'cdcr.us', 'KNZ', 'konz.us'))

#removing a subset of CDR e001 and e002 treatments to prevent the majority of data being from CDR; keeping lowest, highest, and 10 gm-2 N additions (levels most comparable to other studies)
divCDRe001<-tallgrassSEM%>%
  filter(site_code=="CDR"&treatment %in% c(1,6,8,9))
divCDRe002<-tallgrassSEM%>%
  filter(site_code=="CDR"&treatment %in% c('1_f_u_n','6_f_u_n','8_f_u_n','9_f_u_n'))

tallgrassSEMsubset <- tallgrassSEM%>%
  filter(site_code!='CDR')%>%
  rbind(divCDRe001, divCDRe002)

#combine the two CDR experiments
divTrt<-rbind(divTrt1, divCDRe002, divCDRe001)

#######MODEL STRUCTURE: use differences between treatment and control plots in each year (horizontal arrows)-------------------

###all data, all years-----------

###difference metrics not through composition model
#all years
summary(compositionModel1 <- psem(
  lm(anpp_pdiff_transform ~ n + p + k + richness_difference + evenness_diff + rank_difference + species_difference, data=allSEMdata),
  lm(richness_difference ~ n + p + k, data=allSEMdata),
  lm(evenness_diff ~ n + p + k, data=allSEMdata),
  lm(rank_difference ~ n + p + k, data=allSEMdata),
  lm(species_difference ~ n + p + k, data=allSEMdata),
  richness_difference %~~% evenness_diff,
  richness_difference %~~% rank_difference,
  richness_difference %~~% species_difference,
  evenness_diff %~~% rank_difference,
  evenness_diff %~~% species_difference,
  rank_difference %~~% species_difference,
  data=allSEMdata
))
coefs1 <- coefs(compositionModel1, standardize = "scale", standardize.type = "latent.linear", intercepts = FALSE)%>%
  select(Response, Predictor, Estimate, Std.Error, DF, Crit.Value, P.Value, Std.Estimate)%>%
  mutate(treatment_year=1)


####################
# Konza and Cedar Creek
####################
#year 1
summary(compositionModel1 <- psem(
  lm(anpp_pdiff_transform ~ richness_difference + evenness_diff + rank_difference + species_difference + n + p + k, data=subset(tallgrassSEMsubset, treatment_year==1)),
  lm(richness_difference ~ n + p + k, data=subset(tallgrassSEMsubset, treatment_year==1)),
  lm(evenness_diff ~ n + p + k, data=subset(tallgrassSEMsubset, treatment_year==1)),
  lm(rank_difference ~ n + p + k, data=subset(tallgrassSEMsubset, treatment_year==1)),
  lm(species_difference ~ n + p + k, data=subset(tallgrassSEMsubset, treatment_year==1)),
  richness_difference %~~% evenness_diff,
  richness_difference %~~% rank_difference,
  richness_difference %~~% species_difference,
  evenness_diff %~~% rank_difference,
  evenness_diff %~~% species_difference,
  rank_difference %~~% species_difference,
  data=subset(tallgrassSEMsubset, treatment_year==1)
))
coefs1 <- coefs(compositionModel1, standardize = "scale", standardize.type = "latent.linear", intercepts = FALSE)%>%
  select(Response, Predictor, Estimate, Std.Error, DF, Crit.Value, P.Value, Std.Estimate)%>%
  mutate(treatment_year=1)

#year 2
summary(compositionModel2 <- psem(
  lm(anpp_pdiff_transform ~ richness_difference + evenness_diff + rank_difference + species_difference + n + p + k, data=subset(tallgrassSEMsubset, treatment_year==2)),
  lm(richness_difference ~ n + p + k, data=subset(tallgrassSEMsubset, treatment_year==2)),
  lm(evenness_diff ~ n + p + k, data=subset(tallgrassSEMsubset, treatment_year==2)),
  lm(rank_difference ~ n + p + k, data=subset(tallgrassSEMsubset, treatment_year==2)),
  lm(species_difference ~ n + p + k, data=subset(tallgrassSEMsubset, treatment_year==2)),
  richness_difference %~~% evenness_diff,
  richness_difference %~~% rank_difference,
  richness_difference %~~% species_difference,
  evenness_diff %~~% rank_difference,
  evenness_diff %~~% species_difference,
  rank_difference %~~% species_difference,
  data=subset(tallgrassSEMsubset, treatment_year==2)
))
coefs2 <- coefs(compositionModel2, standardize = "scale", standardize.type = "latent.linear", intercepts = FALSE)%>%
  select(Response, Predictor, Estimate, Std.Error, DF, Crit.Value, P.Value, Std.Estimate)%>%
  mutate(treatment_year=2)

#year 3
summary(compositionModel3 <- psem(
  lm(anpp_pdiff_transform ~ richness_difference + evenness_diff + rank_difference + species_difference + n + p + k, data=subset(tallgrassSEMsubset, treatment_year==3)),
  lm(richness_difference ~ n + p + k, data=subset(tallgrassSEMsubset, treatment_year==3)),
  lm(evenness_diff ~ n + p + k, data=subset(tallgrassSEMsubset, treatment_year==3)),
  lm(rank_difference ~ n + p + k, data=subset(tallgrassSEMsubset, treatment_year==3)),
  lm(species_difference ~ n + p + k, data=subset(tallgrassSEMsubset, treatment_year==3)),
  richness_difference %~~% evenness_diff,
  richness_difference %~~% rank_difference,
  richness_difference %~~% species_difference,
  evenness_diff %~~% rank_difference,
  evenness_diff %~~% species_difference,
  rank_difference %~~% species_difference,
  data=subset(tallgrassSEMsubset, treatment_year==3)
))
coefs3 <- coefs(compositionModel3, standardize = "scale", standardize.type = "latent.linear", intercepts = FALSE)%>%
  select(Response, Predictor, Estimate, Std.Error, DF, Crit.Value, P.Value, Std.Estimate)%>%
  mutate(treatment_year=3)

#year 4
summary(compositionModel4 <- psem(
  lm(anpp_pdiff_transform ~ richness_difference + evenness_diff + rank_difference + species_difference + n + p + k, data=subset(tallgrassSEMsubset, treatment_year==4)),
  lm(richness_difference ~ n + p + k, data=subset(tallgrassSEMsubset, treatment_year==4)),
  lm(evenness_diff ~ n + p + k, data=subset(tallgrassSEMsubset, treatment_year==4)),
  lm(rank_difference ~ n + p + k, data=subset(tallgrassSEMsubset, treatment_year==4)),
  lm(species_difference ~ n + p + k, data=subset(tallgrassSEMsubset, treatment_year==4)),
  richness_difference %~~% evenness_diff,
  richness_difference %~~% rank_difference,
  richness_difference %~~% species_difference,
  evenness_diff %~~% rank_difference,
  evenness_diff %~~% species_difference,
  rank_difference %~~% species_difference,
  data=subset(tallgrassSEMsubset, treatment_year==4)
))
coefs4 <- coefs(compositionModel4, standardize = "scale", standardize.type = "latent.linear", intercepts = FALSE)%>%
  select(Response, Predictor, Estimate, Std.Error, DF, Crit.Value, P.Value, Std.Estimate)%>%
  mutate(treatment_year=4)

#year 5
summary(compositionModel5 <- psem(
  lm(anpp_pdiff_transform ~ richness_difference + evenness_diff + rank_difference + species_difference + n + p + k, data=subset(tallgrassSEMsubset, treatment_year==5)),
  lm(richness_difference ~ n + p + k, data=subset(tallgrassSEMsubset, treatment_year==5)),
  lm(evenness_diff ~ n + p + k, data=subset(tallgrassSEMsubset, treatment_year==5)),
  lm(rank_difference ~ n + p + k, data=subset(tallgrassSEMsubset, treatment_year==5)),
  lm(species_difference ~ n + p + k, data=subset(tallgrassSEMsubset, treatment_year==5)),
  richness_difference %~~% evenness_diff,
  richness_difference %~~% rank_difference,
  richness_difference %~~% species_difference,
  evenness_diff %~~% rank_difference,
  evenness_diff %~~% species_difference,
  rank_difference %~~% species_difference,
  data=subset(tallgrassSEMsubset, treatment_year==5)
))
coefs5 <- coefs(compositionModel5, standardize = "scale", standardize.type = "latent.linear", intercepts = FALSE)%>%
  select(Response, Predictor, Estimate, Std.Error, DF, Crit.Value, P.Value, Std.Estimate)%>%
  mutate(treatment_year=5)

#year 6
summary(compositionModel6 <- psem(
  lm(anpp_pdiff_transform ~ richness_difference + evenness_diff + rank_difference + species_difference + n + p + k, data=subset(tallgrassSEMsubset, treatment_year==6)),
  lm(richness_difference ~ n + p + k, data=subset(tallgrassSEMsubset, treatment_year==6)),
  lm(evenness_diff ~ n + p + k, data=subset(tallgrassSEMsubset, treatment_year==6)),
  lm(rank_difference ~ n + p + k, data=subset(tallgrassSEMsubset, treatment_year==6)),
  lm(species_difference ~ n + p + k, data=subset(tallgrassSEMsubset, treatment_year==6)),
  richness_difference %~~% evenness_diff,
  richness_difference %~~% rank_difference,
  richness_difference %~~% species_difference,
  evenness_diff %~~% rank_difference,
  evenness_diff %~~% species_difference,
  rank_difference %~~% species_difference,
  data=subset(tallgrassSEMsubset, treatment_year==6)
))
coefs6 <- coefs(compositionModel6, standardize = "scale", standardize.type = "latent.linear", intercepts = FALSE)%>%
  select(Response, Predictor, Estimate, Std.Error, DF, Crit.Value, P.Value, Std.Estimate)%>%
  mutate(treatment_year=6)

#year 7
summary(compositionModel7 <- psem(
  lm(anpp_pdiff_transform ~ richness_difference + evenness_diff + rank_difference + species_difference + n + p + k, data=subset(tallgrassSEMsubset, treatment_year==7)),
  lm(richness_difference ~ n + p + k, data=subset(tallgrassSEMsubset, treatment_year==7)),
  lm(evenness_diff ~ n + p + k, data=subset(tallgrassSEMsubset, treatment_year==7)),
  lm(rank_difference ~ n + p + k, data=subset(tallgrassSEMsubset, treatment_year==7)),
  lm(species_difference ~ n + p + k, data=subset(tallgrassSEMsubset, treatment_year==7)),
  richness_difference %~~% evenness_diff,
  richness_difference %~~% rank_difference,
  richness_difference %~~% species_difference,
  evenness_diff %~~% rank_difference,
  evenness_diff %~~% species_difference,
  rank_difference %~~% species_difference,
  data=subset(tallgrassSEMsubset, treatment_year==7)
))
coefs7 <- coefs(compositionModel7, standardize = "scale", standardize.type = "latent.linear", intercepts = FALSE)%>%
  select(Response, Predictor, Estimate, Std.Error, DF, Crit.Value, P.Value, Std.Estimate)%>%
  mutate(treatment_year=7)

#year 8
summary(compositionModel8 <- psem(
  lm(anpp_pdiff_transform ~ richness_difference + evenness_diff + rank_difference + species_difference + n + p + k, data=subset(tallgrassSEMsubset, treatment_year==8)),
  lm(richness_difference ~ n + p + k, data=subset(tallgrassSEMsubset, treatment_year==8)),
  lm(evenness_diff ~ n + p + k, data=subset(tallgrassSEMsubset, treatment_year==8)),
  lm(rank_difference ~ n + p + k, data=subset(tallgrassSEMsubset, treatment_year==8)),
  lm(species_difference ~ n + p + k, data=subset(tallgrassSEMsubset, treatment_year==8)),
  richness_difference %~~% evenness_diff,
  richness_difference %~~% rank_difference,
  richness_difference %~~% species_difference,
  evenness_diff %~~% rank_difference,
  evenness_diff %~~% species_difference,
  rank_difference %~~% species_difference,
  data=subset(tallgrassSEMsubset, treatment_year==8)
))
coefs8 <- coefs(compositionModel8, standardize = "scale", standardize.type = "latent.linear", intercepts = FALSE)%>%
  select(Response, Predictor, Estimate, Std.Error, DF, Crit.Value, P.Value, Std.Estimate)%>%
  mutate(treatment_year=8)

#year 9
summary(compositionModel9 <- psem(
  lm(anpp_pdiff_transform ~ richness_difference + evenness_diff + rank_difference + species_difference + n + p + k, data=subset(tallgrassSEMsubset, treatment_year==9)),
  lm(richness_difference ~ n + p + k, data=subset(tallgrassSEMsubset, treatment_year==9)),
  lm(evenness_diff ~ n + p + k, data=subset(tallgrassSEMsubset, treatment_year==9)),
  lm(rank_difference ~ n + p + k, data=subset(tallgrassSEMsubset, treatment_year==9)),
  lm(species_difference ~ n + p + k, data=subset(tallgrassSEMsubset, treatment_year==9)),
  richness_difference %~~% evenness_diff,
  richness_difference %~~% rank_difference,
  richness_difference %~~% species_difference,
  evenness_diff %~~% rank_difference,
  evenness_diff %~~% species_difference,
  rank_difference %~~% species_difference,
  data=subset(tallgrassSEMsubset, treatment_year==9)
))
coefs9 <- coefs(compositionModel9, standardize = "scale", standardize.type = "latent.linear", intercepts = FALSE)%>%
  select(Response, Predictor, Estimate, Std.Error, DF, Crit.Value, P.Value, Std.Estimate)%>%
  mutate(treatment_year=9)

#year 10
summary(compositionModel10 <- psem(
  lm(anpp_pdiff_transform ~ richness_difference + evenness_diff + rank_difference + species_difference + n + p + k, data=subset(tallgrassSEMsubset, treatment_year==10)),
  lm(richness_difference ~ n + p + k, data=subset(tallgrassSEMsubset, treatment_year==10)),
  lm(evenness_diff ~ n + p + k, data=subset(tallgrassSEMsubset, treatment_year==10)),
  lm(rank_difference ~ n + p + k, data=subset(tallgrassSEMsubset, treatment_year==10)),
  lm(species_difference ~ n + p + k, data=subset(tallgrassSEMsubset, treatment_year==10)),
  richness_difference %~~% evenness_diff,
  richness_difference %~~% rank_difference,
  richness_difference %~~% species_difference,
  evenness_diff %~~% rank_difference,
  evenness_diff %~~% species_difference,
  rank_difference %~~% species_difference,
  data=subset(tallgrassSEMsubset, treatment_year==10)
))
coefs10 <- coefs(compositionModel10, standardize = "scale", standardize.type = "latent.linear", intercepts = FALSE)%>%
  select(Response, Predictor, Estimate, Std.Error, DF, Crit.Value, P.Value, Std.Estimate)%>%
  mutate(treatment_year=10)

#combine the coeffecients
metricsModelCoef <- rbind(coefs1, coefs2, coefs3, coefs4, coefs5, coefs6, coefs7, coefs8, coefs9, coefs10)
metricsSum <- metricsModelCoef%>%
  group_by(treatment_year)%>%
  summarise(total_std_estimate=sum(abs(Std.Estimate)))%>%
  ungroup()

###figure of drivers of each response - export at 1000x600
ggplot(data=subset(metricsModelCoef, Response=='anpp_pdiff_transform'), aes(x=treatment_year, y=Std.Estimate)) +
  geom_point(size=5) +
  geom_smooth(method='lm') +
  ylab('ANPP Response') +
  facet_wrap(~Predictor, scales='free')

ggplot(data=subset(metricsModelCoef, Response=='richness_difference'), aes(x=treatment_year, y=Std.Estimate)) +
  geom_point(size=5) +
  geom_smooth(method='lm') +
  ylab('Richness Response') +
  facet_wrap(~Predictor)

ggplot(data=subset(metricsModelCoef, Response=='evenness_diff'), aes(x=treatment_year, y=Std.Estimate)) +
  geom_point(size=5) +
  geom_smooth(method='lm') +
  ylab('Evenness Response') +
  facet_wrap(~Predictor)

ggplot(data=subset(metricsModelCoef, Response=='rank_difference'), aes(x=treatment_year, y=Std.Estimate)) +
  geom_point(size=5) +
  geom_smooth(method='lm') +
  ylab('Rank Diff Response') +
  facet_wrap(~Predictor)

ggplot(data=subset(metricsModelCoef, Response=='species_difference'), aes(x=treatment_year, y=Std.Estimate)) +
  geom_point(size=5) +
  geom_smooth(method='lm') +
  ylab('Species Diff Response') +
  facet_wrap(~Predictor, scales='free')

# ###figure of drivers of community difference
# ggplot(data=subset(metricsModelCoef, Response=='composition_diff_transform'), aes(x=treatment_year, y=Std.Estimate)) +
#   geom_point(size=5) +
#   geom_smooth(method='lm') +
#   facet_wrap(~Predictor)

####################
# Konza only
####################
#year 1
summary(compositionModel1 <- psem(
  lm(anpp_pdiff_transform ~ richness_difference + evenness_diff + rank_difference + species_difference + n + p + k, data=subset(KNZsem, treatment_year==1)),
  lm(richness_difference ~ n + p + k, data=subset(KNZsem, treatment_year==1)),
  lm(evenness_diff ~ n + p + k, data=subset(KNZsem, treatment_year==1)),
  lm(rank_difference ~ n + p + k, data=subset(KNZsem, treatment_year==1)),
  lm(species_difference ~ n + p + k, data=subset(KNZsem, treatment_year==1)),
  richness_difference %~~% evenness_diff,
  richness_difference %~~% rank_difference,
  richness_difference %~~% species_difference,
  evenness_diff %~~% rank_difference,
  evenness_diff %~~% species_difference,
  rank_difference %~~% species_difference,
  data=subset(KNZsem, treatment_year==1)
  ))
coefs1 <- coefs(compositionModel1, standardize = "scale", standardize.type = "latent.linear", intercepts = FALSE)%>%
  select(Response, Predictor, Estimate, Std.Error, DF, Crit.Value, P.Value, Std.Estimate)%>%
  mutate(treatment_year=1)

#year 2
summary(compositionModel2 <- psem(
  lm(anpp_pdiff_transform ~ richness_difference + evenness_diff + rank_difference + species_difference + n + p + k, data=subset(KNZsem, treatment_year==2)),
  lm(richness_difference ~ n + p + k, data=subset(KNZsem, treatment_year==2)),
  lm(evenness_diff ~ n + p + k, data=subset(KNZsem, treatment_year==2)),
  lm(rank_difference ~ n + p + k, data=subset(KNZsem, treatment_year==2)),
  lm(species_difference ~ n + p + k, data=subset(KNZsem, treatment_year==2)),
  richness_difference %~~% evenness_diff,
  richness_difference %~~% rank_difference,
  richness_difference %~~% species_difference,
  evenness_diff %~~% rank_difference,
  evenness_diff %~~% species_difference,
  rank_difference %~~% species_difference,
  data=subset(KNZsem, treatment_year==2)
))
coefs2 <- coefs(compositionModel2, standardize = "scale", standardize.type = "latent.linear", intercepts = FALSE)%>%
  select(Response, Predictor, Estimate, Std.Error, DF, Crit.Value, P.Value, Std.Estimate)%>%
  mutate(treatment_year=2)

#year 3
summary(compositionModel3 <- psem(
  lm(anpp_pdiff_transform ~ richness_difference + evenness_diff + rank_difference + species_difference + n + p + k, data=subset(KNZsem, treatment_year==3)),
  lm(richness_difference ~ n + p + k, data=subset(KNZsem, treatment_year==3)),
  lm(evenness_diff ~ n + p + k, data=subset(KNZsem, treatment_year==3)),
  lm(rank_difference ~ n + p + k, data=subset(KNZsem, treatment_year==3)),
  lm(species_difference ~ n + p + k, data=subset(KNZsem, treatment_year==3)),
  richness_difference %~~% evenness_diff,
  richness_difference %~~% rank_difference,
  richness_difference %~~% species_difference,
  evenness_diff %~~% rank_difference,
  evenness_diff %~~% species_difference,
  rank_difference %~~% species_difference,
  data=subset(KNZsem, treatment_year==3)
))
coefs3 <- coefs(compositionModel3, standardize = "scale", standardize.type = "latent.linear", intercepts = FALSE)%>%
  select(Response, Predictor, Estimate, Std.Error, DF, Crit.Value, P.Value, Std.Estimate)%>%
  mutate(treatment_year=3)

#year 4
summary(compositionModel4 <- psem(
  lm(anpp_pdiff_transform ~ richness_difference + evenness_diff + rank_difference + species_difference + n + p + k, data=subset(KNZsem, treatment_year==4)),
  lm(richness_difference ~ n + p + k, data=subset(KNZsem, treatment_year==4)),
  lm(evenness_diff ~ n + p + k, data=subset(KNZsem, treatment_year==4)),
  lm(rank_difference ~ n + p + k, data=subset(KNZsem, treatment_year==4)),
  lm(species_difference ~ n + p + k, data=subset(KNZsem, treatment_year==4)),
  richness_difference %~~% evenness_diff,
  richness_difference %~~% rank_difference,
  richness_difference %~~% species_difference,
  evenness_diff %~~% rank_difference,
  evenness_diff %~~% species_difference,
  rank_difference %~~% species_difference,
  data=subset(KNZsem, treatment_year==4)
))
coefs4 <- coefs(compositionModel4, standardize = "scale", standardize.type = "latent.linear", intercepts = FALSE)%>%
  select(Response, Predictor, Estimate, Std.Error, DF, Crit.Value, P.Value, Std.Estimate)%>%
  mutate(treatment_year=4)

#year 5
summary(compositionModel5 <- psem(
  lm(anpp_pdiff_transform ~ richness_difference + evenness_diff + rank_difference + species_difference + n + p + k, data=subset(KNZsem, treatment_year==5)),
  lm(richness_difference ~ n + p + k, data=subset(KNZsem, treatment_year==5)),
  lm(evenness_diff ~ n + p + k, data=subset(KNZsem, treatment_year==5)),
  lm(rank_difference ~ n + p + k, data=subset(KNZsem, treatment_year==5)),
  lm(species_difference ~ n + p + k, data=subset(KNZsem, treatment_year==5)),
  richness_difference %~~% evenness_diff,
  richness_difference %~~% rank_difference,
  richness_difference %~~% species_difference,
  evenness_diff %~~% rank_difference,
  evenness_diff %~~% species_difference,
  rank_difference %~~% species_difference,
  data=subset(KNZsem, treatment_year==5)
))
coefs5 <- coefs(compositionModel5, standardize = "scale", standardize.type = "latent.linear", intercepts = FALSE)%>%
  select(Response, Predictor, Estimate, Std.Error, DF, Crit.Value, P.Value, Std.Estimate)%>%
  mutate(treatment_year=5)

#year 6
summary(compositionModel6 <- psem(
  lm(anpp_pdiff_transform ~ richness_difference + evenness_diff + rank_difference + species_difference + n + p + k, data=subset(KNZsem, treatment_year==6)),
  lm(richness_difference ~ n + p + k, data=subset(KNZsem, treatment_year==6)),
  lm(evenness_diff ~ n + p + k, data=subset(KNZsem, treatment_year==6)),
  lm(rank_difference ~ n + p + k, data=subset(KNZsem, treatment_year==6)),
  lm(species_difference ~ n + p + k, data=subset(KNZsem, treatment_year==6)),
  richness_difference %~~% evenness_diff,
  richness_difference %~~% rank_difference,
  richness_difference %~~% species_difference,
  evenness_diff %~~% rank_difference,
  evenness_diff %~~% species_difference,
  rank_difference %~~% species_difference,
  data=subset(KNZsem, treatment_year==6)
))
coefs6 <- coefs(compositionModel6, standardize = "scale", standardize.type = "latent.linear", intercepts = FALSE)%>%
  select(Response, Predictor, Estimate, Std.Error, DF, Crit.Value, P.Value, Std.Estimate)%>%
  mutate(treatment_year=6)

#year 7
summary(compositionModel7 <- psem(
  lm(anpp_pdiff_transform ~ richness_difference + evenness_diff + rank_difference + species_difference + n + p + k, data=subset(KNZsem, treatment_year==7)),
  lm(richness_difference ~ n + p + k, data=subset(KNZsem, treatment_year==7)),
  lm(evenness_diff ~ n + p + k, data=subset(KNZsem, treatment_year==7)),
  lm(rank_difference ~ n + p + k, data=subset(KNZsem, treatment_year==7)),
  lm(species_difference ~ n + p + k, data=subset(KNZsem, treatment_year==7)),
  richness_difference %~~% evenness_diff,
  richness_difference %~~% rank_difference,
  richness_difference %~~% species_difference,
  evenness_diff %~~% rank_difference,
  evenness_diff %~~% species_difference,
  rank_difference %~~% species_difference,
  data=subset(KNZsem, treatment_year==7)
))
coefs7 <- coefs(compositionModel7, standardize = "scale", standardize.type = "latent.linear", intercepts = FALSE)%>%
  select(Response, Predictor, Estimate, Std.Error, DF, Crit.Value, P.Value, Std.Estimate)%>%
  mutate(treatment_year=7)

#year 8
summary(compositionModel8 <- psem(
  lm(anpp_pdiff_transform ~ richness_difference + evenness_diff + rank_difference + species_difference + n + p + k, data=subset(KNZsem, treatment_year==8)),
  lm(richness_difference ~ n + p + k, data=subset(KNZsem, treatment_year==8)),
  lm(evenness_diff ~ n + p + k, data=subset(KNZsem, treatment_year==8)),
  lm(rank_difference ~ n + p + k, data=subset(KNZsem, treatment_year==8)),
  lm(species_difference ~ n + p + k, data=subset(KNZsem, treatment_year==8)),
  richness_difference %~~% evenness_diff,
  richness_difference %~~% rank_difference,
  richness_difference %~~% species_difference,
  evenness_diff %~~% rank_difference,
  evenness_diff %~~% species_difference,
  rank_difference %~~% species_difference,
  data=subset(KNZsem, treatment_year==8)
))
coefs8 <- coefs(compositionModel8, standardize = "scale", standardize.type = "latent.linear", intercepts = FALSE)%>%
  select(Response, Predictor, Estimate, Std.Error, DF, Crit.Value, P.Value, Std.Estimate)%>%
  mutate(treatment_year=8)

#year 9
summary(compositionModel9 <- psem(
  lm(anpp_pdiff_transform ~ richness_difference + evenness_diff + rank_difference + species_difference + n + p + k, data=subset(KNZsem, treatment_year==9)),
  lm(richness_difference ~ n + p + k, data=subset(KNZsem, treatment_year==9)),
  lm(evenness_diff ~ n + p + k, data=subset(KNZsem, treatment_year==9)),
  lm(rank_difference ~ n + p + k, data=subset(KNZsem, treatment_year==9)),
  lm(species_difference ~ n + p + k, data=subset(KNZsem, treatment_year==9)),
  richness_difference %~~% evenness_diff,
  richness_difference %~~% rank_difference,
  richness_difference %~~% species_difference,
  evenness_diff %~~% rank_difference,
  evenness_diff %~~% species_difference,
  rank_difference %~~% species_difference,
  data=subset(KNZsem, treatment_year==9)
))
coefs9 <- coefs(compositionModel9, standardize = "scale", standardize.type = "latent.linear", intercepts = FALSE)%>%
  select(Response, Predictor, Estimate, Std.Error, DF, Crit.Value, P.Value, Std.Estimate)%>%
  mutate(treatment_year=9)

#year 10
summary(compositionModel10 <- psem(
  lm(anpp_pdiff_transform ~ richness_difference + evenness_diff + rank_difference + species_difference + n + p + k, data=subset(KNZsem, treatment_year==10)),
  lm(richness_difference ~ n + p + k, data=subset(KNZsem, treatment_year==10)),
  lm(evenness_diff ~ n + p + k, data=subset(KNZsem, treatment_year==10)),
  lm(rank_difference ~ n + p + k, data=subset(KNZsem, treatment_year==10)),
  lm(species_difference ~ n + p + k, data=subset(KNZsem, treatment_year==10)),
  richness_difference %~~% evenness_diff,
  richness_difference %~~% rank_difference,
  richness_difference %~~% species_difference,
  evenness_diff %~~% rank_difference,
  evenness_diff %~~% species_difference,
  rank_difference %~~% species_difference,
  data=subset(KNZsem, treatment_year==10)
))
coefs10 <- coefs(compositionModel10, standardize = "scale", standardize.type = "latent.linear", intercepts = FALSE)%>%
  select(Response, Predictor, Estimate, Std.Error, DF, Crit.Value, P.Value, Std.Estimate)%>%
  mutate(treatment_year=10)

#combine the coeffecients
metricsModelCoef <- rbind(coefs1, coefs2, coefs3, coefs4, coefs5, coefs6, coefs7, coefs8, coefs9, coefs10)
metricsSum <- metricsModelCoef%>%
  group_by(treatment_year)%>%
  summarise(total_std_estimate=sum(abs(Std.Estimate)))%>%
  ungroup()

###figure of drivers of each response - export at 1000x600
ggplot(data=subset(metricsModelCoef, Response=='anpp_pdiff_transform'), aes(x=treatment_year, y=Std.Estimate)) +
  geom_point(size=5) +
  geom_smooth(method='lm') +
  ylab('ANPP Response') +
  facet_wrap(~Predictor, scales='free')

ggplot(data=subset(metricsModelCoef, Response=='richness_difference'), aes(x=treatment_year, y=Std.Estimate)) +
  geom_point(size=5) +
  geom_smooth(method='lm') +
  ylab('Richness Response') +
  facet_wrap(~Predictor)

ggplot(data=subset(metricsModelCoef, Response=='evenness_diff'), aes(x=treatment_year, y=Std.Estimate)) +
  geom_point(size=5) +
  geom_smooth(method='lm') +
  ylab('Evenness Response') +
  facet_wrap(~Predictor)

ggplot(data=subset(metricsModelCoef, Response=='rank_difference'), aes(x=treatment_year, y=Std.Estimate)) +
  geom_point(size=5) +
  geom_smooth(method='lm') +
  ylab('Rank Diff Response') +
  facet_wrap(~Predictor)

ggplot(data=subset(metricsModelCoef, Response=='species_difference'), aes(x=treatment_year, y=Std.Estimate)) +
  geom_point(size=5) +
  geom_smooth(method='lm') +
  ylab('Species Diff Response') +
  facet_wrap(~Predictor, scales='free')

# ###figure of drivers of community difference
# ggplot(data=subset(metricsModelCoef, Response=='composition_diff_transform'), aes(x=treatment_year, y=Std.Estimate)) +
#   geom_point(size=5) +
#   geom_smooth(method='lm') +
#   facet_wrap(~Predictor)



####################
# Cedar Creek only
####################
#year 1
summary(compositionModel1 <- psem(
  lm(anpp_pdiff_transform ~ richness_difference + evenness_diff + rank_difference + species_difference + n + p + k, data=subset(CDRsem, treatment_year==1)),
  lm(richness_difference ~ n + p + k, data=subset(CDRsem, treatment_year==1)),
  lm(evenness_diff ~ n + p + k, data=subset(CDRsem, treatment_year==1)),
  lm(rank_difference ~ n + p + k, data=subset(CDRsem, treatment_year==1)),
  lm(species_difference ~ n + p + k, data=subset(CDRsem, treatment_year==1)),
  richness_difference %~~% evenness_diff,
  richness_difference %~~% rank_difference,
  richness_difference %~~% species_difference,
  evenness_diff %~~% rank_difference,
  evenness_diff %~~% species_difference,
  rank_difference %~~% species_difference,
  data=subset(CDRsem, treatment_year==1)
))
coefs1 <- coefs(compositionModel1, standardize = "scale", standardize.type = "latent.linear", intercepts = FALSE)%>%
  select(Response, Predictor, Estimate, Std.Error, DF, Crit.Value, P.Value, Std.Estimate)%>%
  mutate(treatment_year=1)

#year 2
summary(compositionModel2 <- psem(
  lm(anpp_pdiff_transform ~ richness_difference + evenness_diff + rank_difference + species_difference + n + p + k, data=subset(CDRsem, treatment_year==2)),
  lm(richness_difference ~ n + p + k, data=subset(CDRsem, treatment_year==2)),
  lm(evenness_diff ~ n + p + k, data=subset(CDRsem, treatment_year==2)),
  lm(rank_difference ~ n + p + k, data=subset(CDRsem, treatment_year==2)),
  lm(species_difference ~ n + p + k, data=subset(CDRsem, treatment_year==2)),
  richness_difference %~~% evenness_diff,
  richness_difference %~~% rank_difference,
  richness_difference %~~% species_difference,
  evenness_diff %~~% rank_difference,
  evenness_diff %~~% species_difference,
  rank_difference %~~% species_difference,
  data=subset(CDRsem, treatment_year==2)
))
coefs2 <- coefs(compositionModel2, standardize = "scale", standardize.type = "latent.linear", intercepts = FALSE)%>%
  select(Response, Predictor, Estimate, Std.Error, DF, Crit.Value, P.Value, Std.Estimate)%>%
  mutate(treatment_year=2)

#year 3
summary(compositionModel3 <- psem(
  lm(anpp_pdiff_transform ~ richness_difference + evenness_diff + rank_difference + species_difference + n + p + k, data=subset(CDRsem, treatment_year==3)),
  lm(richness_difference ~ n + p + k, data=subset(CDRsem, treatment_year==3)),
  lm(evenness_diff ~ n + p + k, data=subset(CDRsem, treatment_year==3)),
  lm(rank_difference ~ n + p + k, data=subset(CDRsem, treatment_year==3)),
  lm(species_difference ~ n + p + k, data=subset(CDRsem, treatment_year==3)),
  richness_difference %~~% evenness_diff,
  richness_difference %~~% rank_difference,
  richness_difference %~~% species_difference,
  evenness_diff %~~% rank_difference,
  evenness_diff %~~% species_difference,
  rank_difference %~~% species_difference,
  data=subset(CDRsem, treatment_year==3)
))
coefs3 <- coefs(compositionModel3, standardize = "scale", standardize.type = "latent.linear", intercepts = FALSE)%>%
  select(Response, Predictor, Estimate, Std.Error, DF, Crit.Value, P.Value, Std.Estimate)%>%
  mutate(treatment_year=3)

#year 4
summary(compositionModel4 <- psem(
  lm(anpp_pdiff_transform ~ richness_difference + evenness_diff + rank_difference + species_difference + n + p + k, data=subset(CDRsem, treatment_year==4)),
  lm(richness_difference ~ n + p + k, data=subset(CDRsem, treatment_year==4)),
  lm(evenness_diff ~ n + p + k, data=subset(CDRsem, treatment_year==4)),
  lm(rank_difference ~ n + p + k, data=subset(CDRsem, treatment_year==4)),
  lm(species_difference ~ n + p + k, data=subset(CDRsem, treatment_year==4)),
  richness_difference %~~% evenness_diff,
  richness_difference %~~% rank_difference,
  richness_difference %~~% species_difference,
  evenness_diff %~~% rank_difference,
  evenness_diff %~~% species_difference,
  rank_difference %~~% species_difference,
  data=subset(CDRsem, treatment_year==4)
))
coefs4 <- coefs(compositionModel4, standardize = "scale", standardize.type = "latent.linear", intercepts = FALSE)%>%
  select(Response, Predictor, Estimate, Std.Error, DF, Crit.Value, P.Value, Std.Estimate)%>%
  mutate(treatment_year=4)

#year 5
summary(compositionModel5 <- psem(
  lm(anpp_pdiff_transform ~ richness_difference + evenness_diff + rank_difference + species_difference + n + p + k, data=subset(CDRsem, treatment_year==5)),
  lm(richness_difference ~ n + p + k, data=subset(CDRsem, treatment_year==5)),
  lm(evenness_diff ~ n + p + k, data=subset(CDRsem, treatment_year==5)),
  lm(rank_difference ~ n + p + k, data=subset(CDRsem, treatment_year==5)),
  lm(species_difference ~ n + p + k, data=subset(CDRsem, treatment_year==5)),
  richness_difference %~~% evenness_diff,
  richness_difference %~~% rank_difference,
  richness_difference %~~% species_difference,
  evenness_diff %~~% rank_difference,
  evenness_diff %~~% species_difference,
  rank_difference %~~% species_difference,
  data=subset(CDRsem, treatment_year==5)
))
coefs5 <- coefs(compositionModel5, standardize = "scale", standardize.type = "latent.linear", intercepts = FALSE)%>%
  select(Response, Predictor, Estimate, Std.Error, DF, Crit.Value, P.Value, Std.Estimate)%>%
  mutate(treatment_year=5)

#year 6
summary(compositionModel6 <- psem(
  lm(anpp_pdiff_transform ~ richness_difference + evenness_diff + rank_difference + species_difference + n + p + k, data=subset(CDRsem, treatment_year==6)),
  lm(richness_difference ~ n + p + k, data=subset(CDRsem, treatment_year==6)),
  lm(evenness_diff ~ n + p + k, data=subset(CDRsem, treatment_year==6)),
  lm(rank_difference ~ n + p + k, data=subset(CDRsem, treatment_year==6)),
  lm(species_difference ~ n + p + k, data=subset(CDRsem, treatment_year==6)),
  richness_difference %~~% evenness_diff,
  richness_difference %~~% rank_difference,
  richness_difference %~~% species_difference,
  evenness_diff %~~% rank_difference,
  evenness_diff %~~% species_difference,
  rank_difference %~~% species_difference,
  data=subset(CDRsem, treatment_year==6)
))
coefs6 <- coefs(compositionModel6, standardize = "scale", standardize.type = "latent.linear", intercepts = FALSE)%>%
  select(Response, Predictor, Estimate, Std.Error, DF, Crit.Value, P.Value, Std.Estimate)%>%
  mutate(treatment_year=6)

#year 7
summary(compositionModel7 <- psem(
  lm(anpp_pdiff_transform ~ richness_difference + evenness_diff + rank_difference + species_difference + n + p + k, data=subset(CDRsem, treatment_year==7)),
  lm(richness_difference ~ n + p + k, data=subset(CDRsem, treatment_year==7)),
  lm(evenness_diff ~ n + p + k, data=subset(CDRsem, treatment_year==7)),
  lm(rank_difference ~ n + p + k, data=subset(CDRsem, treatment_year==7)),
  lm(species_difference ~ n + p + k, data=subset(CDRsem, treatment_year==7)),
  richness_difference %~~% evenness_diff,
  richness_difference %~~% rank_difference,
  richness_difference %~~% species_difference,
  evenness_diff %~~% rank_difference,
  evenness_diff %~~% species_difference,
  rank_difference %~~% species_difference,
  data=subset(CDRsem, treatment_year==7)
))
coefs7 <- coefs(compositionModel7, standardize = "scale", standardize.type = "latent.linear", intercepts = FALSE)%>%
  select(Response, Predictor, Estimate, Std.Error, DF, Crit.Value, P.Value, Std.Estimate)%>%
  mutate(treatment_year=7)

#year 8
summary(compositionModel8 <- psem(
  lm(anpp_pdiff_transform ~ richness_difference + evenness_diff + rank_difference + species_difference + n + p + k, data=subset(CDRsem, treatment_year==8)),
  lm(richness_difference ~ n + p + k, data=subset(CDRsem, treatment_year==8)),
  lm(evenness_diff ~ n + p + k, data=subset(CDRsem, treatment_year==8)),
  lm(rank_difference ~ n + p + k, data=subset(CDRsem, treatment_year==8)),
  lm(species_difference ~ n + p + k, data=subset(CDRsem, treatment_year==8)),
  richness_difference %~~% evenness_diff,
  richness_difference %~~% rank_difference,
  richness_difference %~~% species_difference,
  evenness_diff %~~% rank_difference,
  evenness_diff %~~% species_difference,
  rank_difference %~~% species_difference,
  data=subset(CDRsem, treatment_year==8)
))
coefs8 <- coefs(compositionModel8, standardize = "scale", standardize.type = "latent.linear", intercepts = FALSE)%>%
  select(Response, Predictor, Estimate, Std.Error, DF, Crit.Value, P.Value, Std.Estimate)%>%
  mutate(treatment_year=8)

#year 9
summary(compositionModel9 <- psem(
  lm(anpp_pdiff_transform ~ richness_difference + evenness_diff + rank_difference + species_difference + n + p + k, data=subset(CDRsem, treatment_year==9)),
  lm(richness_difference ~ n + p + k, data=subset(CDRsem, treatment_year==9)),
  lm(evenness_diff ~ n + p + k, data=subset(CDRsem, treatment_year==9)),
  lm(rank_difference ~ n + p + k, data=subset(CDRsem, treatment_year==9)),
  lm(species_difference ~ n + p + k, data=subset(CDRsem, treatment_year==9)),
  richness_difference %~~% evenness_diff,
  richness_difference %~~% rank_difference,
  richness_difference %~~% species_difference,
  evenness_diff %~~% rank_difference,
  evenness_diff %~~% species_difference,
  rank_difference %~~% species_difference,
  data=subset(CDRsem, treatment_year==9)
))
coefs9 <- coefs(compositionModel9, standardize = "scale", standardize.type = "latent.linear", intercepts = FALSE)%>%
  select(Response, Predictor, Estimate, Std.Error, DF, Crit.Value, P.Value, Std.Estimate)%>%
  mutate(treatment_year=9)

#year 10
summary(compositionModel10 <- psem(
  lm(anpp_pdiff_transform ~ richness_difference + evenness_diff + rank_difference + species_difference + n + p + k, data=subset(CDRsem, treatment_year==10)),
  lm(richness_difference ~ n + p + k, data=subset(CDRsem, treatment_year==10)),
  lm(evenness_diff ~ n + p + k, data=subset(CDRsem, treatment_year==10)),
  lm(rank_difference ~ n + p + k, data=subset(CDRsem, treatment_year==10)),
  lm(species_difference ~ n + p + k, data=subset(CDRsem, treatment_year==10)),
  richness_difference %~~% evenness_diff,
  richness_difference %~~% rank_difference,
  richness_difference %~~% species_difference,
  evenness_diff %~~% rank_difference,
  evenness_diff %~~% species_difference,
  rank_difference %~~% species_difference,
  data=subset(CDRsem, treatment_year==10)
))
coefs10 <- coefs(compositionModel10, standardize = "scale", standardize.type = "latent.linear", intercepts = FALSE)%>%
  select(Response, Predictor, Estimate, Std.Error, DF, Crit.Value, P.Value, Std.Estimate)%>%
  mutate(treatment_year=10)

#combine the coeffecients
metricsModelCoef <- rbind(coefs1, coefs2, coefs3, coefs4, coefs5, coefs6, coefs7, coefs8, coefs9, coefs10)
metricsSum <- metricsModelCoef%>%
  group_by(treatment_year)%>%
  summarise(total_std_estimate=sum(abs(Std.Estimate)))%>%
  ungroup()

###figure of drivers of each response - export at 1000x600
ggplot(data=subset(metricsModelCoef, Response=='anpp_pdiff_transform'), aes(x=treatment_year, y=Std.Estimate)) +
  geom_point(size=5) +
  geom_smooth(method='lm') +
  ylab('ANPP Response') +
  facet_wrap(~Predictor, scales='free')

ggplot(data=subset(metricsModelCoef, Response=='richness_difference'), aes(x=treatment_year, y=Std.Estimate)) +
  geom_point(size=5) +
  geom_smooth(method='lm') +
  ylab('Richness Response') +
  facet_wrap(~Predictor)

ggplot(data=subset(metricsModelCoef, Response=='evenness_diff'), aes(x=treatment_year, y=Std.Estimate)) +
  geom_point(size=5) +
  geom_smooth(method='lm') +
  ylab('Evenness Response') +
  facet_wrap(~Predictor)

ggplot(data=subset(metricsModelCoef, Response=='rank_difference'), aes(x=treatment_year, y=Std.Estimate)) +
  geom_point(size=5) +
  geom_smooth(method='lm') +
  ylab('Rank Diff Response') +
  facet_wrap(~Predictor)

ggplot(data=subset(metricsModelCoef, Response=='species_difference'), aes(x=treatment_year, y=Std.Estimate)) +
  geom_point(size=5) +
  geom_smooth(method='lm') +
  ylab('Species Diff Response') +
  facet_wrap(~Predictor, scales='free')
