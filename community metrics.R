library(vegan)
library(codyn)
library(ggplot2)
library(grid)
library(plyr)
library(dplyr)
library(tidyr)

#1212 this is just a test

#kim's wd
setwd('C:\\Users\\Kim\\Dropbox\\working groups\\converge diverge working group\\converge_diverge\\datasets\\LongForm\\to yang')

source('C:\\Users\\Kim\\Dropbox\\working groups\\converge diverge working group\\C2E\\C2E\\aggregate across replicates.R')

#calculate dominance: berger-parker method
dominance <- relAbundAgg%>%
  group_by(site_code, project_name, community_type, treatment, treatment_year, calendar_year, n, n_treatment)%>%
  summarise(bp_dominance=max(relcov_agg))

#calculate richness
richness <- relAbundAgg%>%
  filter(relcov_agg>0)%>%
  group_by(site_code, project_name, community_type, treatment, treatment_year, calendar_year, n, n_treatment)%>%
  summarise(richness=length(relcov_agg))

#calculate evenness
relAbundAggWide <- relAbundAgg%>%
  group_by(site_code, project_name, community_type, treatment, treatment_year, calendar_year, n, n_treatment)%>%
  spread(key=genus_species, value=relcov_agg, fill=0)

S <- specnumber(relAbundAggWide[,9:ncol(relAbundAggWide)])
invD <- diversity(relAbundAggWide[,9:ncol(relAbundAggWide)],"inv")
evenness <- as.data.frame(invD/S)%>%
  cbind(relAbundAggWide)%>%
  select(site_code, project_name, community_type, treatment, treatment_year, calendar_year, n, n_treatment, `invD/S`)


#calculate spp gain and loss through space (trt compared to ctl)
turnoverSpace <- relAbundAgg[,c(-4,-8)]
turnoverSpace <- turnoverSpace%>%
  mutate(relcov_agg=ifelse(relcov_agg>0, 1, 0))%>%
  spread(key=n_treatment, value=relcov_agg, fill=0)%>%
  mutate(diff=`1`-`0`)%>%
  mutate(diff2=diff)

gain <- turnoverSpace%>%
  group_by(site_code, project_name, community_type, treatment_year, calendar_year)%>%
#need to count up -1 as loss and +1 as gain
  summarise(gain=sum(diff[diff2==1]))

loss <- turnoverSpace%>%
  group_by(site_code, project_name, community_type, treatment_year, calendar_year)%>%
  #need to count up -1 as loss and +1 as gain
  summarise(loss=sum(diff[diff2==-1]))

turnoverSpace2 <- gain%>%
  left_join(loss)

#calculate spp loss and gain through time
#make a new dataframe with just the label
expTrt <- relAbundAgg
expTrt$exp_trt <- with(relAbundAgg, (paste(site_code, project_name, community_type, treatment, sep=':')))
expTrt2 <- unique(expTrt$exp_trt)

#make a new dataframe to collect the turnover metrics
turnoverAll=data.frame(row.names=1)

for(i in 1:length(expTrt2)) {
  
  #creates a dataset for each unique year, trt, exp combo
  subset=expTrt[expTrt$exp_trt==as.character(expTrt2[i]),]%>%
    select(exp_trt, treatment, n, n_treatment, genus_species, relcov_agg)
  
  #need this to keep track of n treatment
  labels=as.data.frame(unique(subset[c('exp_trt', 'calendar_year')]))
  
  #calculate disappearance
  disappearance=turnover(df=subset, time.var='calendar_year', species.var='genus_species', abundance.var='relcov_agg', replicate.var=NA, metric='disappearance')
  
  #calculate appearance
  appearance=turnover(df=subset, time.var='calendar_year', species.var='genus_species', abundance.var='relcov_agg', replicate.var=NA, metric='appearance')
  
  #calculate turnover
  total=turnover(df=subset, time.var='calendar_year', species.var='genus_species', abundance.var='relcov_agg', replicate.var=NA, metric='total')
  
  #merging back with labels to get back plot_mani
  turnover=cbind(labels[-1,], disappearance$disappearance, appearance$appearance, total$total)
  
  #pasting variables into the dataframe made for this analysis
  turnoverAll=rbind(turnoverAll, turnover)
}

names(turnoverAll)[names(turnoverAll)=='appearance$appearance'] <- 'appearance'
names(turnoverAll)[names(turnoverAll)=='disappearance$disappearance'] <- 'disappearance'
names(turnoverAll)[names(turnoverAll)=='total$total'] <- 'turnover'

turnoverAll <- turnoverAll%>%
  separate(exp_trt, c('site_code', 'project_name', 'community_type', 'treatment'), sep=":")%>%
  select(site_code, project_name, community_type, treatment, disappearance, appearance, turnover, calendar_year)

#rank correlations
rankCorr <- read.csv('kendall_rank.csv')%>%
  select(-X)

#merge all community metrics
communityMetrics <- dominance%>%
  left_join(richness)%>%
  left_join(evenness)%>%
  left_join(turnoverAll, by=c('site_code', 'project_name', 'community_type', 'treatment', 'calendar_year'))%>%
  left_join(turnoverSpace2)%>%
  mutate(siteprojcom=paste(site_code, project_name, community_type, sep=''))%>%
  left_join(rankCorr, by=c('siteprojcom', 'treatment_year'))%>%
  select(-siteprojcom)%>%
  gather(key=metric, value=value, bp_dominance:RankCor)

communityMetrics <- communityMetrics[,c(-4,-7)]

#calculate effect size (as ln response ratio)
communityLRR <- communityMetrics%>%
  spread(key=n_treatment, value=value)%>%
  mutate(lrr=log(`1`/`0`), lrr=ifelse(metric=='gain'|metric=='loss', `0`, lrr))

