library(grid)
library(PerformanceAnalytics)
library(tidyverse)

#kim's desktop
# setwd('C:\\Users\\la pierrek\\Dropbox (Smithsonian)\\working groups\\CoRRE\\converge_diverge\\datasets\\LongForm')

#kim's laptop
setwd('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\CoRRE\\converge_diverge\\datasets\\LongForm')


theme_set(theme_bw())
theme_update(axis.title.x=element_text(size=40, vjust=-0.35, margin=margin(t=15)), axis.text.x=element_text(size=34, color='black'),
             axis.title.y=element_text(size=40, angle=90, vjust=0.5, margin=margin(r=15)), axis.text.y=element_text(size=34, color='black'),
             plot.title = element_text(size=40, vjust=2),
             panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
             legend.title=element_blank(), legend.text=element_text(size=20))

####data input

###community data
nutnetCommChange <- read.csv('NutNet_community differences_07162019.csv')%>%
  mutate(treatment_year=time)%>%
  select(-X, -time)%>%
  mutate(n=ifelse(treatment=='N'|treatment=='NP'|treatment=='NK'|treatment=='NPK'|treatment=='NPK+Fence', 10, 0), p=ifelse(treatment=='P'|treatment=='NP'|treatment=='PK'|treatment=='NPK'|treatment=='NPK+Fence', 10, 0), k=ifelse(treatment=='K'|treatment=='NK'|treatment=='PK'|treatment=='NPK'|treatment=='NPK+Fence', 10, 0), fence=ifelse(treatment=='Fence'|treatment=='NPK+Fence', 1, 0))%>%
  select(-treatment)%>%
  rename(treatment=treatment2)

# #checking community data for outliers
# ggplot(nutnetCommChange, aes(comp_diff)) + geom_histogram() + facet_wrap(~site_code, scales='free')

###anpp outliers
nutnetANPPoutliers <- read.csv('C:\\Users\\lapie\\Dropbox (Smithsonian)\\nutrient network\\NutNet data\\La Pierre_NutNet_anpp_potential outliers_12122017.csv')%>%
  filter(checked.with.PI=='incorrect')%>%
  select(-live, -notes)

###anpp data
#remove outliers
nutnetANPP <- read.csv('C:\\Users\\lapie\\Dropbox (Smithsonian)\\nutrient network\\NutNet data\\full-biomass-16-July-2019.csv')%>%
  filter(year_trt!=0, live==1)%>%
  merge(nutnetANPPoutliers, by=c("year", "year_trt", "trt", "site_name", "site_code", "block", "plot", "subplot", "mass", "category"), all.x=T)%>%
  filter(is.na(checked.with.PI))%>%
  group_by(site_code, plot, year_trt, trt)%>%
  summarise(anpp=sum(mass))%>%
  ungroup()


# #checking site level data for outliers in anpp
# ggplot(subset(nutnetANPP, category=='GRAMINOID'), aes(mass)) + geom_histogram() + facet_wrap(~site_code, scales='free')
# ggplot(subset(nutnetANPP, category=='FORB'), aes(mass)) + geom_histogram() + facet_wrap(~site_code, scales='free')
# ggplot(subset(nutnetANPP, category=='LEGUME'), aes(mass)) + geom_histogram() + facet_wrap(~site_code, scales='free')
# ggplot(subset(nutnetANPP, category=='BRYOPHYTE'), aes(mass)) + geom_histogram() + facet_wrap(~site_code, scales='free')
# ggplot(subset(nutnetANPP, category=='CACTUS'), aes(mass)) + geom_histogram() + facet_wrap(~site_code, scales='free')
# ggplot(subset(nutnetANPP, category=='LIVE'), aes(mass)) + geom_histogram() + facet_wrap(~site_code, scales='free')
# ggplot(subset(nutnetANPP, category=='PERENNIAL'), aes(mass)) + geom_histogram() + facet_wrap(~site_code, scales='free')
# ggplot(subset(nutnetANPP, category=='WOODY'), aes(mass)) + geom_histogram() + facet_wrap(~site_code, scales='free')
# ggplot(subset(nutnetANPP, category=='VASCULAR'), aes(mass)) + geom_histogram() + facet_wrap(~site_code, scales='free')

#calculating anpp difference
#anpp ctl data
nutnetANPPctl <- nutnetANPP%>%
  filter(trt=='Control')%>%
  select(site_code, year_trt, anpp)%>%
  group_by(site_code, year_trt)%>%
  summarise(anpp_ctl=mean(anpp))%>%
  ungroup()
#anpp change
nutnetANPPchange <- nutnetANPP%>%
  filter(trt!='Control')%>%
  group_by(site_code, year_trt, trt)%>%
  summarise(anpp=mean(anpp))%>%
  ungroup()%>%
  left_join(nutnetANPPctl)%>%
  #calculate anpp change as percent change from ctl in each year
  mutate(anpp_pdiff=(anpp-anpp_ctl)/anpp_ctl, treatment=trt, treatment_year=year_trt)%>%
  select(-anpp, -anpp_ctl, -trt, -year_trt)

#calculate dataset lengths
numYears <- nutnetANPPchange%>%
  select(treatment_year, site_code)%>%
  unique()%>%
  group_by(site_code)%>%
  summarise(num_years=length(treatment_year))

#merge
nutnetSEMdata <- nutnetANPPchange%>%
  left_join(nutnetCommChange)%>%
  na.omit()%>%
  #filter out experiments less than 3 years old
  left_join(numYears)%>%
  filter(num_years>2)%>%
  #drop fencing treatment
  filter(fence==0)
#all of these compare treatment to control!

###exploratory correlations and histograms (all variables compare treatment to control plots)
dataVis <- nutnetSEMdata%>%
  select(anpp_pdiff, composition_diff, richness_difference, evenness_diff, rank_difference, species_difference) #make visualization dataframe
chart.Correlation(dataVis, histogram=T, pch=19)


rm(list=setdiff(ls(), "nutnetSEMdata"))
# write.csv(nutnetSEMdata, 'NutNet_comm and anpp diff_07160219.csv')
