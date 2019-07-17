library(tidyverse)
library(codyn)
library(vegan)
library(Kendall)
library(gridExtra)
library(reldist)
library(grid)
library(gtable)

setwd('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\CoRRE\\converge_diverge\\datasets\\LongForm')

nutnetAbsCover <- read.csv('C:\\Users\\lapie\\Dropbox (Smithsonian)\\nutrient network\\NutNet data\\full-cover-16-July-2019.csv')%>%
  filter(year_trt!=0, live!=0)

nutnetTotalCover <- nutnetAbsCover%>%
  group_by(year, site_name, plot)%>%
  summarise(total_cover=sum(max_cover))

numYears <- nutnetAbsCover%>%
  select(year, site_code)%>%
  unique()%>%
  group_by(site_code)%>%
  summarise(num_years=length(year))

nutnetRelCover <- nutnetAbsCover%>%
  left_join(nutnetTotalCover)%>%
  mutate(abundance=max_cover/total_cover)%>%
  left_join(numYears)%>%
  filter(num_years>2)%>%
  mutate(time=year_trt, treatment=trt, replicate=plot, species=Taxon, C_T=as.factor(ifelse(treatment=='Control', 'Control', 'Treatment')))%>%
  select(site_code, time, treatment, replicate, species, abundance, C_T)


#make a new dataframe with just the label;
site=nutnetRelCover%>%
  select(site_code)%>%
  unique()

#makes an empty dataframe
nutnetMultDiff=data.frame(row.names=1)


###first get SERSp
###second get multivariate differences
for(i in 1:length(site$site_code)) {
  #creates a dataset for each unique year, trt, exp combo
  subset=nutnetRelCover[nutnetRelCover$site_code==as.character(site$site_code[i]),]
  
  #get site_project_comm label to paste back on after functions are run
  labels=subset%>%
    select(site_code, treatment)%>%
    unique()
  
  subset2=subset%>%
    select(-site_code)
  
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
  nutnetMultDiff=rbind(mult_diff, nutnetMultDiff)  
}


# write.csv(nutnetMultDiff, 'NutNet_community differences_07162019.csv')

