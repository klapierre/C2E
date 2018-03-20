library(tidyverse)
library(codyn)
library(vegan)
library(Kendall)
library(gridExtra)
library(reldist)
library(grid)
library(gtable)

setwd('C:\\Users\\la pierrek\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\datasets\\LongForm')

nutnetAbsCover <- read.csv('C:\\Users\\la pierrek\\Dropbox (Smithsonian)\\NutNet data\\full-cover-04-December-2017.csv')%>%
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
nutnetSERSp=data.frame(row.names=1) 
nutnetMultDiff=data.frame(row.names=1)


###first get SERSp
###second get multivariate differences
for(i in 1:length(site$site_code)) {
  #creates a dataset for each unique year, trt, exp combo
  subset=nutnetRelCover[nutnetRelCover$site_code==as.character(site$site_code[i]),]
  
  #get site_code label to paste back on after functions are run
  labels=subset%>%
    select(site_code, treatment)%>%
    unique()
  
  subset2=subset%>%
    select(-site_code)
  
  #calculate metrics
  sersp=SERSp_func(subset2)
  mult_diff=multivariate_diff_func(subset2)%>%
    left_join(sersp)%>%
    left_join(labels)%>%
    filter(treatment!='Control')
  
  #pasting dispersions into the dataframe made for this analysis
  nutnetMultDiff=rbind(mult_diff, nutnetMultDiff)  
}

write.csv(nutnetMultDiff, 'NutNet_community differences_12052017.csv')

