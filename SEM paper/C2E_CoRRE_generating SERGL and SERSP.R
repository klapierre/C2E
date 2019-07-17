library(codyn)
library(Kendall)
library(gridExtra)
library(reldist)
library(grid)
library(gtable)
library(tidyverse)


#kim's desktop
setwd('C:\\Users\\la pierrek\\Dropbox (Smithsonian)\\working groups\\CoRRE\\converge_diverge\\datasets\\LongForm')

#kim's laptop
setwd('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\CoRRE\\converge_diverge\\datasets\\LongForm')


trt <- read.csv('ExperimentInformation_March2019.csv')%>%
  select(-X)

siteInfo <- read.csv('SiteExperimentDetails_March2019.csv')%>%
  select(-X)

correRelCover <- read.csv('SpeciesRelativeAbundance_March2019.csv')%>%
  select(-X)%>%
  left_join(trt)%>%
  left_join(siteInfo)%>%
  filter(treatment_year>0)%>%
  rename(time=treatment_year, replicate=plot_id, species=genus_species, abundance=relcov)%>%
  mutate(C_T=as.factor(ifelse(plot_mani==0, 'Control', 'Treatment')), site_project_comm=paste(site_code, project_name, community_type, sep='::'))%>%
  select(site_project_comm, time, treatment, replicate, species, abundance, C_T)%>%
  group_by(site_project_comm, treatment)%>%
  mutate(min_year=min(time))%>%
  ungroup()%>%
  #drop datasets where first year isn't year 1
  filter(min_year==1)


#make a new dataframe with just the label;
site=correRelCover%>%
  select(site_project_comm)%>%
  unique()

#makes an empty dataframe
correMultDiff=data.frame(row.names=1)


###first get SERSp
###second get multivariate difference
for(i in 1:length(site$site_project_comm)) {
  #creates a dataset for each unique year, trt, exp combo
  subset=correRelCover[correRelCover$site_project_comm==as.character(site$site_project_comm[i]),]%>%
    mutate(treatment=ifelse(C_T=='Control', 'Control', as.character(treatment)))
  
  #get site_project_comm label to paste back on after functions are run
  repLabels=subset%>%
    select(site_project_comm, treatment, replicate)%>%
    unique()
  labels=subset%>%
    select(site_project_comm, treatment)%>%
    unique()
  
  subset2=subset%>%
    select(-site_project_comm)
  
  #calculate SERGL metrics
  sersp=RAC_difference(subset2, time.var='time', species.var='species', abundance.var='abundance', replicate.var='replicate', treatment.var='treatment', reference.treatment='Control')%>%
    left_join(labels)%>%
  #calculate mean across replicates
    group_by(time, treatment, treatment2)%>%
    summarise(richness_difference=mean(richness_diff), evenness_diff=mean(evenness_diff), rank_difference=mean(rank_diff), species_difference=mean(species_diff))
    
  #calculate multivariate change
  mult_diff=multivariate_difference(subset2, time.var='time', species.var='species', abundance.var='abundance', replicate.var='replicate', treatment.var='treatment', reference.treatment='Control')%>%
    left_join(sersp)%>%
    left_join(labels)
  
  #pasting dispersions into the dataframe made for this analysis
  correMultDiff=rbind(mult_diff, correMultDiff)  
}

# write.csv(correMultDiff, 'corre_community_difference_July2019.csv', row.names=F)

