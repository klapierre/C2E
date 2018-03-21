library(tidyverse)
library(ggplot2)
library(grid)

#kim's working directory
setwd('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\C2E\\Products\\CommunityChange\\March2018 WG')

#emily's working directory
setwd("/Users/egrman/Dropbox/C2E/Products/CommunityChange/March2018 WG")

theme_set(theme_bw())
theme_update(axis.title.x=element_text(size=20, vjust=-0.35), axis.text.x=element_text(size=16),
             axis.title.y=element_text(size=20, angle=90, vjust=0.5), axis.text.y=element_text(size=16),
             plot.title = element_text(size=24, vjust=2),
             panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
             legend.title=element_blank(), legend.text=element_text(size=20))

lnRR <- read.csv('CORRE_RAC_LogRR_March2018_trtyr.csv')
metrics <- read.csv('CORRE_RAC_Metrics_March2018_trtyr.csv')

trt <- read.csv('ExperimentInformation_Nov2017.csv')%>%
  mutate(site_project_comm=as.factor(paste(site_code, project_name, community_type, sep='_')))%>%
  select(-X)

rawAbund <- read.csv('SpeciesRawAbundance_Oct2017.csv')%>%
  select(-X, -abundance, -genus_species, -block, -calendar_year, -treatment_year)%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep='_'))%>%
  unique()

lnRRtrt <- lnRR%>%
  select(-X)%>%
  left_join(trt)

metricsTrt <- metrics%>%
  select(-X)%>%
  left_join(rawAbund)%>%
  left_join(trt)

#lines of change for KBS across all treatements
metricsSFig <- ggplot(data=subset(metricsTrt, site_project_comm=='KBS_T7_0'), aes(x=treatment_year, y=abs(richness_change), color=treatment)) +
  # geom_point() +
  geom_smooth(method='loess', se=T, size=2) +
  geom_abline(slope=0, intercept=0)
metricsEFig <- ggplot(data=subset(metricsTrt, site_project_comm=='KBS_T7_0'), aes(x=treatment_year, y=evenness_change, color=treatment)) +
  # geom_point() +
  geom_smooth(method='loess', se=T, size=2) +
  geom_abline(slope=0, intercept=0)
metricsRFig <- ggplot(data=subset(metricsTrt, site_project_comm=='KBS_T7_0'), aes(x=treatment_year, y=rank_change, color=treatment)) +
  # geom_point() +
  geom_smooth(method='loess', se=T, size=2) +
  geom_abline(slope=0, intercept=0)
metricsLFig <- ggplot(data=subset(metricsTrt, site_project_comm=='KBS_T7_0'), aes(x=treatment_year, y=losses, color=treatment)) +
  # geom_point() +
  geom_smooth(method='loess', se=T, size=2) +
  geom_abline(slope=0, intercept=0)
metricsGFig <- ggplot(data=subset(metricsTrt, site_project_comm=='KBS_T7_0'), aes(x=treatment_year, y=gains, color=treatment)) +
  # geom_point() +
  geom_smooth(method='loess', se=T, size=2) +
  geom_abline(slope=0, intercept=0)
pushViewport(viewport(layout=grid.layout(5,1)))
print(metricsSFig, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
print(metricsEFig, vp=viewport(layout.pos.row = 2, layout.pos.col = 1))
print(metricsRFig, vp=viewport(layout.pos.row = 3, layout.pos.col = 1))
print(metricsLFig, vp=viewport(layout.pos.row = 4, layout.pos.col = 1))
print(metricsGFig, vp=viewport(layout.pos.row = 5, layout.pos.col = 1))


#lines of change for KUFS across all treatements
metricsSFig <- ggplot(data=subset(metricsTrt, site_project_comm=='KUFS_E6_type2'), aes(x=treatment_year, y=abs(richness_change), color=treatment)) +
  # geom_point() +
  geom_smooth(method='loess', se=T, size=2) +
  geom_abline(slope=0, intercept=0)
metricsEFig <- ggplot(data=subset(metricsTrt, site_project_comm=='KUFS_E6_type2'), aes(x=treatment_year, y=evenness_change, color=treatment)) +
  # geom_point() +
  geom_smooth(method='loess', se=T, size=2) +
  geom_abline(slope=0, intercept=0)
metricsRFig <- ggplot(data=subset(metricsTrt, site_project_comm=='KUFS_E6_type2'), aes(x=treatment_year, y=rank_change, color=treatment)) +
  # geom_point() +
  geom_smooth(method='loess', se=T, size=2) +
  geom_abline(slope=0, intercept=0)
metricsLFig <- ggplot(data=subset(metricsTrt, site_project_comm=='KUFS_E6_type2'), aes(x=treatment_year, y=losses, color=treatment)) +
  # geom_point() +
  geom_smooth(method='loess', se=T, size=2) +
  geom_abline(slope=0, intercept=0)
metricsGFig <- ggplot(data=subset(metricsTrt, site_project_comm=='KUFS_E6_type2'), aes(x=treatment_year, y=gains, color=treatment)) +
  # geom_point() +
  geom_smooth(method='loess', se=T, size=2) +
  geom_abline(slope=0, intercept=0)
pushViewport(viewport(layout=grid.layout(5,1)))
print(metricsSFig, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
print(metricsEFig, vp=viewport(layout.pos.row = 2, layout.pos.col = 1))
print(metricsRFig, vp=viewport(layout.pos.row = 3, layout.pos.col = 1))
print(metricsLFig, vp=viewport(layout.pos.row = 4, layout.pos.col = 1))
print(metricsGFig, vp=viewport(layout.pos.row = 5, layout.pos.col = 1))


#lines of change for NANT across all treatements
metricsSFig <- ggplot(data=subset(metricsTrt, site_project_comm=='NANT_wet_Nant-Mrsh_BSA_S'), aes(x=treatment_year, y=abs(richness_change), color=treatment)) +
  # geom_point() +
  geom_smooth(method='loess', se=T, size=2) +
  geom_abline(slope=0, intercept=0)
metricsEFig <- ggplot(data=subset(metricsTrt, site_project_comm=='NANT_wet_Nant-Mrsh_BSA_S'), aes(x=treatment_year, y=evenness_change, color=treatment)) +
  # geom_point() +
  geom_smooth(method='loess', se=T, size=2) +
  geom_abline(slope=0, intercept=0)
metricsRFig <- ggplot(data=subset(metricsTrt, site_project_comm=='NANT_wet_Nant-Mrsh_BSA_S'), aes(x=treatment_year, y=rank_change, color=treatment)) +
  # geom_point() +
  geom_smooth(method='loess', se=T, size=2) +
  geom_abline(slope=0, intercept=0)
metricsLFig <- ggplot(data=subset(metricsTrt, site_project_comm=='NANT_wet_Nant-Mrsh_BSA_S'), aes(x=treatment_year, y=losses, color=treatment)) +
  # geom_point() +
  geom_smooth(method='loess', se=T, size=2) +
  geom_abline(slope=0, intercept=0)
metricsGFig <- ggplot(data=subset(metricsTrt, site_project_comm=='NANT_wet_Nant-Mrsh_BSA_S'), aes(x=treatment_year, y=gains, color=treatment)) +
  # geom_point() +
  geom_smooth(method='loess', se=T, size=2) +
  geom_abline(slope=0, intercept=0)
pushViewport(viewport(layout=grid.layout(5,1)))
print(metricsSFig, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
print(metricsEFig, vp=viewport(layout.pos.row = 2, layout.pos.col = 1))
print(metricsRFig, vp=viewport(layout.pos.row = 3, layout.pos.col = 1))
print(metricsLFig, vp=viewport(layout.pos.row = 4, layout.pos.col = 1))
print(metricsGFig, vp=viewport(layout.pos.row = 5, layout.pos.col = 1))


pplots <- metricsTrt%>%
  filter(treatment=='N1P0'|treatment=='N2P0'|treatment=='N2P3')

pplots=metricsTrt[metricsTrt$treatment %in% c("N1P0", "N2P0", "N2P3") & metricsTrt$site_project_comm=="KNZ_pplots_0",]

#lines of change for pplots across all treatements
metricsSFig <- ggplot(data=pplots, aes(x=treatment_year, y=abs(richness_change), color=treatment)) +
  # geom_point() +
  geom_smooth(method='loess', se=T, size=2) +
  geom_abline(slope=0, intercept=0)
metricsEFig <- ggplot(data=pplots, aes(x=treatment_year, y=evenness_change, color=treatment)) +
  # geom_point() +
  geom_smooth(method='loess', se=T, size=2) +
  geom_abline(slope=0, intercept=0)
metricsRFig <- ggplot(data=pplots, aes(x=treatment_year, y=rank_change, color=treatment)) +
  # geom_point() +
  geom_smooth(method='loess', se=T, size=2) +
  geom_abline(slope=0, intercept=0)
metricsLFig <- ggplot(data=pplots, aes(x=treatment_year, y=losses, color=treatment)) +
  # geom_point() +
  geom_smooth(method='loess', se=T, size=2) +
  geom_abline(slope=0, intercept=0)
metricsGFig <- ggplot(data=pplots, aes(x=treatment_year, y=gains, color=treatment)) +
  # geom_point() +
  geom_smooth(method='loess', se=T, size=2) +
  geom_abline(slope=0, intercept=0)
pushViewport(viewport(layout=grid.layout(5,1)))
print(metricsSFig, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
print(metricsEFig, vp=viewport(layout.pos.row = 2, layout.pos.col = 1))
print(metricsRFig, vp=viewport(layout.pos.row = 3, layout.pos.col = 1))
print(metricsLFig, vp=viewport(layout.pos.row = 4, layout.pos.col = 1))
print(metricsGFig, vp=viewport(layout.pos.row = 5, layout.pos.col = 1))

#histograms for N trts only
lrSFig <- ggplot(data=subset(lnRRtrt, n>0), aes(x=as.factor(treatment_year), y=abs(lrS))) +
  geom_boxplot() +
  geom_abline(slope=0, intercept=0) +
  coord_cartesian(ylim=c(0,1))
lrEFig <- ggplot(data=subset(lnRRtrt, n>0), aes(x=as.factor(treatment_year), y=abs(lrE))) +
  geom_boxplot() +
  geom_abline(slope=0, intercept=0) +
  coord_cartesian(ylim=c(0,1))
lrRFig <- ggplot(data=subset(lnRRtrt, n>0), aes(x=as.factor(treatment_year), y=abs(lrR))) +
  geom_boxplot() +
  geom_abline(slope=0, intercept=0) +
  coord_cartesian(ylim=c(0,1))
lrLFig <- ggplot(data=subset(lnRRtrt, n>0), aes(x=as.factor(treatment_year), y=abs(lrL))) +
  geom_boxplot() +
  geom_abline(slope=0, intercept=0) +
  coord_cartesian(ylim=c(0,1))
lrGFig <- ggplot(data=subset(lnRRtrt, n>0), aes(x=as.factor(treatment_year), y=abs(lrG))) +
  geom_boxplot() +
  geom_abline(slope=0, intercept=0) +
  coord_cartesian(ylim=c(0,1))
pushViewport(viewport(layout=grid.layout(5,1)))
print(lrSFig, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
print(lrEFig, vp=viewport(layout.pos.row = 2, layout.pos.col = 1))
print(lrRFig, vp=viewport(layout.pos.row = 3, layout.pos.col = 1))
print(lrLFig, vp=viewport(layout.pos.row = 4, layout.pos.col = 1))
print(lrGFig, vp=viewport(layout.pos.row = 5, layout.pos.col = 1))




#N trts only - lines
lrSFig <- ggplot(data=subset(lnRRtrt, precip!=0), aes(x=treatment_year, y=lrS)) +
  geom_point() +
  geom_smooth(method='loess', se=F) +
  geom_abline(slope=0, intercept=0) 
lrEFig <- ggplot(data=subset(lnRRtrt, n>0), aes(x=as.factor(treatment_year), y=abs(lrE))) +
  geom_boxplot() +
  geom_abline(slope=0, intercept=0) +
  coord_cartesian(ylim=c(0,1))
lrRFig <- ggplot(data=subset(lnRRtrt, n>0), aes(x=as.factor(treatment_year), y=abs(lrR))) +
  geom_boxplot() +
  geom_abline(slope=0, intercept=0) +
  coord_cartesian(ylim=c(0,1))
lrLFig <- ggplot(data=subset(lnRRtrt, n>0), aes(x=as.factor(treatment_year), y=abs(lrL))) +
  geom_boxplot() +
  geom_abline(slope=0, intercept=0) +
  coord_cartesian(ylim=c(0,1))
lrGFig <- ggplot(data=subset(lnRRtrt, n>0), aes(x=as.factor(treatment_year), y=abs(lrG))) +
  geom_boxplot() +
  geom_abline(slope=0, intercept=0) +
  coord_cartesian(ylim=c(0,1))
pushViewport(viewport(layout=grid.layout(5,1)))
print(lrSFig, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
print(lrEFig, vp=viewport(layout.pos.row = 2, layout.pos.col = 1))
print(lrRFig, vp=viewport(layout.pos.row = 3, layout.pos.col = 1))
print(lrLFig, vp=viewport(layout.pos.row = 4, layout.pos.col = 1))
print(lrGFig, vp=viewport(layout.pos.row = 5, layout.pos.col = 1))

