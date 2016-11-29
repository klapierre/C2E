library(plyr)
library(dplyr)
library(tidyr)
library(codyn)
library(ggplot2)
library(grid)

#kim's wd
setwd('C:\\Users\\Kim\\Dropbox\\working groups\\converge diverge working group\\converge_diverge\\datasets\\LongForm\\to yang')

#ggplot theme set
theme_set(theme_bw())
theme_update(axis.title.x=element_text(size=40, vjust=-0.35, margin=margin(t=15)), axis.text.x=element_text(size=34),
             axis.title.y=element_text(size=40, angle=90, vjust=0.5, margin=margin(r=15)), axis.text.y=element_text(size=34),
             plot.title = element_text(size=24, vjust=2),
             panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
             legend.title=element_blank(), legend.text=element_text(size=20))

#load raw abundance data
rawAbundance <- read.csv('CORRE_irrigationL_e001D_subset.csv')%>%
  select(-X)


#generate presence/absence data
presenceAbsence <- rawAbundance%>%
  mutate(presence=1)%>%
  select(-abundance)%>%
  spread(key=genus_species, value=presence, fill=0)%>%
  gather(key=genus_species, value=presence, -site_code, -project_name, -calendar_year, -treatment_year, -block, -treatment, -plot_id, -community_type, -n, -p, -k, -precip, -other_trt, -plot_mani)

#generate treatment list
trt <- rawAbundance%>%
  select(site_code, project_name, treatment, plot_id, n, p, k, precip, other_trt)%>%
  unique()%>%
  mutate(plot_id=as.character(plot_id))


#subsetting datarames by project
#irrigation plots raw abundance subset
rawAbundanceIrr <- rawAbundance%>%
  filter(project_name=='IRG')

#irrigation plots presence/absence subset
presenceAbsenceIrr <- presenceAbsence%>%
  filter(project_name=='IRG')

#e001 raw abundance subset
rawAbundanceE001 <- rawAbundance%>%
  filter(project_name=='e001')%>%
  #get only the control and highest N treatments
  filter(treatment=='9_y_n'|treatment=='8_y_n')

#e001 raw abundance subset
presenceAbsenceE001 <- presenceAbsence%>%
  filter(project_name=='e001')%>%
  #get only the control and highest N treatments
  filter(treatment=='9_y_n'|treatment=='8_y_n')


###irrigation plots calculations
###raw abundance data
#rate change slopes
rateChangeIrrAbund <- rate_change(rawAbundanceIrr, time.var='treatment_year', species.var='genus_species', abundance.var='abundance', replicate.var='plot_id')

#rate change intervals for visualization
rateChangeIntervalIrrAbund <- rate_change_interval(rawAbundanceIrr, time.var='treatment_year', species.var='genus_species', abundance.var='abundance', replicate.var='plot_id')%>%
  mutate(project_name='IRG')%>%
  left_join(trt)%>%
  mutate(type='abundance')

#rate change figure
rateChangeIntervalIrrAbundFig <- ggplot(rateChangeIntervalIrrAbund, aes(interval, distance, linetype=plot_id)) +
  geom_point(aes(color=treatment)) + stat_smooth(method = "loess", se = F, size = 2, aes(color=treatment)) +
  xlab('Time Interval') +
  ylab('Euclidean Distance') +
  scale_linetype(guide='none') +
  scale_color_manual(values=c('#FF3333', '#3300FF'), name='Treatment', labels=c('Control', 'Irrigation')) +
  annotate('text', x=0, y=150, label='(a) Abundance-based', size=10, hjust='left') +
  theme(legend.position=c(0.1,0.8))

rateChangeIntervalIrrAbundFigB <- ggplot(rateChangeIntervalIrrAbund, aes(interval, distance)) +
  geom_point(aes(color=treatment)) + stat_smooth(method = "loess", se = F, size = 2, aes(color=treatment)) +
  xlab('Time Interval') +
  ylab('Euclidean Distance') +
  scale_linetype(guide='none') +
  scale_color_manual(values=c('#FF3333', '#3300FF'), name='Treatment', labels=c('Control', 'Irrigation')) +
  annotate('text', x=0, y=150, label='(c) Abundance-based', size=10, hjust='left') +
  theme(legend.position='none')

###presence/absence data
#rate change slopes
rateChangeIrrPres <- rate_change(presenceAbsenceIrr, time.var='treatment_year', species.var='genus_species', abundance.var='presence', replicate.var='plot_id')

#rate change intervals for visualization
rateChangeIntervalIrrPres <- rate_change_interval(presenceAbsenceIrr, time.var='treatment_year', species.var='genus_species', abundance.var='presence', replicate.var='plot_id')%>%
  mutate(project_name='IRG')%>%
  left_join(trt)%>%
  mutate(type='presence')

#rate change figure
rateChangeIntervalIrrPresFig <- ggplot(rateChangeIntervalIrrPres, aes(interval, distance, linetype=plot_id)) +
  geom_point(aes(color=treatment)) + stat_smooth(method = "loess", se = F, size = 2, aes(color=treatment)) +
  xlab('Time Interval') +
  ylab('Euclidean Distance') +
  scale_linetype(guide='none') +
  scale_color_manual(values=c('#FF3333', '#3300FF'), name='Treatment', labels=c('Control', 'Irrigation')) +
  annotate('text', x=0, y=5, label='(b) Presence-based', size=10, hjust='left') +
  theme(legend.position='none')

rateChangeIntervalIrrPresFigB <- ggplot(rateChangeIntervalIrrPres, aes(interval, distance)) +
  geom_point(aes(color=treatment)) + stat_smooth(method = "loess", se = F, size = 2, aes(color=treatment)) +
  xlab('Time Interval') +
  ylab('Euclidean Distance') +
  scale_linetype(guide='none') +
  scale_color_manual(values=c('#FF3333', '#3300FF'), name='Treatment', labels=c('Control', 'Irrigation')) +
  annotate('text', x=0, y=5, label='(d) Presence-based', size=10, hjust='left') +
  theme(legend.position='none')

pushViewport(viewport(layout=grid.layout(2,2)))
print(rateChangeIntervalIrrAbundFig, vp=viewport(layout.pos.row=1, layout.pos.col=1))
print(rateChangeIntervalIrrPresFig, vp=viewport(layout.pos.row=1, layout.pos.col=2))
print(rateChangeIntervalIrrAbundFigB, vp=viewport(layout.pos.row=2, layout.pos.col=1))
print(rateChangeIntervalIrrPresFigB, vp=viewport(layout.pos.row=2, layout.pos.col=2))
#export as 1500x1500


###e001
#e001 calculations
###raw abundance data
#rate change slopes
rateChangeE001Abund <- rate_change(rawAbundanceE001, time.var='treatment_year', species.var='genus_species', abundance.var='abundance', replicate.var='plot_id')

#rate change intervals for visualization
rateChangeIntervalE001Abund <- rate_change_interval(rawAbundanceE001, time.var='treatment_year', species.var='genus_species', abundance.var='abundance', replicate.var='plot_id')%>%
  mutate(project_name='e001')%>%
  left_join(trt)%>%
  mutate(type='abundance')

#rate change figure
rateChangeIntervalE001AbundFig <- ggplot(rateChangeIntervalE001Abund, aes(interval, distance, linetype=plot_id)) +
  geom_point(aes(color=treatment)) + stat_smooth(method = "loess", se = F, size = 2, aes(color=treatment)) +
  xlab('Time Interval') +
  ylab('Euclidean Distance') +
  scale_linetype(guide='none') +
  scale_color_manual(values=c('#3300FF', '#FF3333'), name='Treatment', labels=c('high NPK+', 'control')) +
  annotate('text', x=0, y=3000, label='(a) Abundance-based', size=10, hjust='left') +
  theme(legend.position=c(0.9,0.9))

rateChangeIntervalE001AbundFigB <- ggplot(rateChangeIntervalE001Abund, aes(interval, distance)) +
  geom_point(aes(color=treatment)) + stat_smooth(method = "loess", se = F, size = 2, aes(color=treatment)) +
  xlab('Time Interval') +
  ylab('Euclidean Distance') +
  scale_linetype(guide='none') +
  scale_color_manual(values=c('#3300FF', '#FF3333'), name='Treatment', labels=c('high NPK+', 'control')) +
  annotate('text', x=0, y=3000, label='(c) Abundance-based', size=10, hjust='left') +
  theme(legend.position='none')

###presence/absence data
#rate change slopes
rateChangeE001Pres <- rate_change(presenceAbsenceE001, time.var='treatment_year', species.var='genus_species', abundance.var='presence', replicate.var='plot_id')

#rate change intervals for visualization
rateChangeIntervalE001Pres <- rate_change_interval(presenceAbsenceE001, time.var='treatment_year', species.var='genus_species', abundance.var='presence', replicate.var='plot_id')%>%
  mutate(project_name='e001')%>%
  left_join(trt)%>%
  mutate(type='presence')

#rate change figure
rateChangeIntervalE001PresFig <- ggplot(rateChangeIntervalE001Pres, aes(interval, distance, linetype=plot_id)) +
  geom_point(aes(color=treatment)) + stat_smooth(method = "loess", se = F, size = 2, aes(color=treatment)) +
  xlab('Time Interval') +
  ylab('Euclidean Distance') +
  scale_linetype(guide='none') +
  scale_color_manual(values=c('#3300FF', '#FF3333'), name='Treatment', labels=c('high NPK+', 'control')) +
  annotate('text', x=0, y=5.2, label='(b) Presence-based', size=10, hjust='left') +
  theme(legend.position='none')

rateChangeIntervalE001PresFigB <- ggplot(rateChangeIntervalE001Pres, aes(interval, distance)) +
  geom_point(aes(color=treatment)) + stat_smooth(method = "loess", se = F, size = 2, aes(color=treatment)) +
  xlab('Time Interval') +
  ylab('Euclidean Distance') +
  scale_linetype(guide='none') +
  scale_color_manual(values=c('#3300FF', '#FF3333'), name='Treatment', labels=c('high NPK+', 'control')) +
  annotate('text', x=0, y=5.2, label='(d) Presence-based', size=10, hjust='left') +
  theme(legend.position='none')

pushViewport(viewport(layout=grid.layout(2,2)))
print(rateChangeIntervalE001AbundFig, vp=viewport(layout.pos.row=1, layout.pos.col=1))
print(rateChangeIntervalE001PresFig, vp=viewport(layout.pos.row=1, layout.pos.col=2))
print(rateChangeIntervalE001AbundFigB, vp=viewport(layout.pos.row=2, layout.pos.col=1))
print(rateChangeIntervalE001PresFigB, vp=viewport(layout.pos.row=2, layout.pos.col=2))
#export as 1500x1500





rawAbundanceFieldC <- read.csv('CORRE_raw_abundance.csv')%>%
  select(-X)%>%
  filter(site_code=='CDR'&project_name=='e001'&community_type=='C')


#generate presence/absence data
presenceAbsenceC <- rawAbundanceFieldC%>%
  mutate(presence=1)%>%
  select(-abundance)%>%
  spread(key=genus_species, value=presence, fill=0)%>%
  gather(key=genus_species, value=presence, -site_code, -project_name, -calendar_year, -treatment_year, -block, -treatment, -plot_id, -community_type)

#generate treatment list
trt <- rawAbundanceFieldC%>%
  select(site_code, project_name, treatment, plot_id)%>%
  unique()%>%
  mutate(plot_id=as.character(plot_id))

#e001 field C raw abundance subset
rawAbundanceE001C <- rawAbundanceFieldC%>%
  filter(project_name=='e001')%>%
  #get only the control and highest N treatments
  filter(treatment=='9_y_n'|treatment=='8_y_n')

#e001 raw abundance subset
presenceAbsenceE001C <- presenceAbsenceC%>%
  filter(project_name=='e001')%>%
  #get only the control and highest N treatments
  filter(treatment=='9_y_n'|treatment=='8_y_n')


###calculations
###raw abundance data
#rate change slopes
rateChangeCAbund <- rate_change(rawAbundanceE001C, time.var='treatment_year', species.var='genus_species', abundance.var='abundance', replicate.var='plot_id')

#rate change intervals for visualization
rateChangeIntervalCAbund <- rate_change_interval(rawAbundanceE001C, time.var='treatment_year', species.var='genus_species', abundance.var='abundance', replicate.var='plot_id')%>%
  mutate(project_name='e001')%>%
  left_join(trt)%>%
  mutate(type='abundance')

#rate change figure
rateChangeIntervalCAbundFig <- ggplot(rateChangeIntervalCAbund, aes(interval, distance, linetype=plot_id)) +
  geom_point(aes(color=treatment)) + stat_smooth(method = "loess", se = F, size = 2, aes(color=treatment)) +
  xlab('Time Interval') +
  ylab('Euclidean Distance') +
  scale_linetype(guide='none') +
  scale_color_manual(values=c('#3300FF', '#FF3333'), name='Treatment', labels=c('high NPK+', 'control')) +
  annotate('text', x=0, y=900, label='(a) Abundance-based', size=10, hjust='left') +
  theme(legend.position=c(0.1,0.8))

rateChangeIntervalCAbundFigB <- ggplot(rateChangeIntervalCAbund, aes(interval, distance)) +
  geom_point(aes(color=treatment)) + stat_smooth(method = "loess", se = F, size = 2, aes(color=treatment)) +
  xlab('Time Interval') +
  ylab('Euclidean Distance') +
  scale_linetype(guide='none') +
  scale_color_manual(values=c('#3300FF', '#FF3333'), name='Treatment', labels=c('high NPK+', 'control')) +
  annotate('text', x=0, y=900, label='(c) Abundance-based', size=10, hjust='left') +
  theme(legend.position='none')

###presence/absence data
#rate change slopes
rateChangeCPres <- rate_change(presenceAbsenceE001C, time.var='treatment_year', species.var='genus_species', abundance.var='presence', replicate.var='plot_id')

#rate change intervals for visualization
rateChangeIntervalCPres <- rate_change_interval(presenceAbsenceE001C, time.var='treatment_year', species.var='genus_species', abundance.var='presence', replicate.var='plot_id')%>%
  mutate(project_name='e001')%>%
  left_join(trt)%>%
  mutate(type='presence')

#rate change figure
rateChangeIntervalCPresFig <- ggplot(rateChangeIntervalCPres, aes(interval, distance, linetype=plot_id)) +
  geom_point(aes(color=treatment)) + stat_smooth(method = "loess", se = F, size = 2, aes(color=treatment)) +
  xlab('Time Interval') +
  ylab('Euclidean Distance') +
  scale_linetype(guide='none') +
  scale_color_manual(values=c('#3300FF', '#FF3333'), name='Treatment', labels=c('high NPK+', 'control')) +
  annotate('text', x=0, y=5, label='(b) Presence-based', size=10, hjust='left') +
  theme(legend.position='none')

rateChangeIntervalCPresFigB <- ggplot(rateChangeIntervalCPres, aes(interval, distance)) +
  geom_point(aes(color=treatment)) + stat_smooth(method = "loess", se = F, size = 2, aes(color=treatment)) +
  xlab('Time Interval') +
  ylab('Euclidean Distance') +
  scale_linetype(guide='none') +
  scale_color_manual(values=c('#3300FF', '#FF3333'), name='Treatment', labels=c('high NPK+', 'control')) +
  annotate('text', x=0, y=5, label='(d) Presence-based', size=10, hjust='left') +
  theme(legend.position='none')

pushViewport(viewport(layout=grid.layout(2,2)))
print(rateChangeIntervalCAbundFig, vp=viewport(layout.pos.row=1, layout.pos.col=1))
print(rateChangeIntervalCPresFig, vp=viewport(layout.pos.row=1, layout.pos.col=2))
print(rateChangeIntervalCAbundFigB, vp=viewport(layout.pos.row=2, layout.pos.col=1))
print(rateChangeIntervalCPresFigB, vp=viewport(layout.pos.row=2, layout.pos.col=2))
#export as 1500x1500







rawAbundancePplots <- read.csv('CORRE_raw_abundance.csv')%>%
  select(-X)%>%
  filter(project_name=='pplots')


#generate presence/absence data
presenceAbsencePplots <- rawAbundancePplots%>%
  mutate(presence=1)%>%
  select(-abundance)%>%
  spread(key=genus_species, value=presence, fill=0)%>%
  gather(key=genus_species, value=presence, -site_code, -project_name, -calendar_year, -treatment_year, -block, -treatment, -plot_id, -community_type)

#generate treatment list
trt <- rawAbundancePplots%>%
  select(site_code, project_name, treatment, plot_id)%>%
  unique()%>%
  mutate(plot_id=as.character(plot_id))

#pplots raw abundance subset
rawAbundancePplotsB <- rawAbundancePplots%>%
  #get only the control and highest N treatments
  filter(treatment=='N2P0'|treatment=='N1P0')

#pplots raw abundance subset
presenceAbsencePplotsB <- presenceAbsencePplots%>%
  #get only the control and highest N treatments
  filter(treatment=='N2P0'|treatment=='N1P0')


###calculations
###raw abundance data
#rate change slopes
rateChangePplotsAbund <- rate_change(rawAbundancePplotsB, time.var='treatment_year', species.var='genus_species', abundance.var='abundance', replicate.var='plot_id')

#rate change intervals for visualization
rateChangeIntervalPplotsAbund <- rate_change_interval(rawAbundancePplotsB, time.var='treatment_year', species.var='genus_species', abundance.var='abundance', replicate.var='plot_id')%>%
  mutate(project_name='pplots')%>%
  left_join(trt)%>%
  mutate(type='abundance')

#rate change figure
rateChangeIntervalPplotsAbundFig <- ggplot(rateChangeIntervalPplotsAbund, aes(interval, distance, linetype=plot_id)) +
  geom_point(aes(color=treatment)) + stat_smooth(method = "loess", se = F, size = 2, aes(color=treatment)) +
  xlab('Time Interval') +
  ylab('Euclidean Distance') +
  scale_linetype(guide='none') +
  scale_color_manual(values=c('#FF3333', '#3300FF'), name='Treatment', labels=c('control', '+N')) +
  annotate('text', x=0, y=150, label='(a) Abundance-based', size=10, hjust='left') +
  theme(legend.position=c(0.1,0.8))

rateChangeIntervalPplotsAbundFigB <- ggplot(rateChangeIntervalPplotsAbund, aes(interval, distance)) +
  geom_point(aes(color=treatment)) + stat_smooth(method = "loess", se = F, size = 2, aes(color=treatment)) +
  xlab('Time Interval') +
  ylab('Euclidean Distance') +
  scale_linetype(guide='none') +
  scale_color_manual(values=c('#FF3333', '#3300FF'), name='Treatment', labels=c('control', '+N')) +
  annotate('text', x=0, y=150, label='(c) Abundance-based', size=10, hjust='left') +
  theme(legend.position='none')

###presence/absence data
#rate change slopes
rateChangePplotsPres <- rate_change(presenceAbsencePplotsB, time.var='treatment_year', species.var='genus_species', abundance.var='presence', replicate.var='plot_id')

#rate change intervals for visualization
rateChangeIntervalPplotsPres <- rate_change_interval(presenceAbsencePplotsB, time.var='treatment_year', species.var='genus_species', abundance.var='presence', replicate.var='plot_id')%>%
  mutate(project_name='pplots')%>%
  left_join(trt)%>%
  mutate(type='presence')

#rate change figure
rateChangeIntervalPplotsPresFig <- ggplot(rateChangeIntervalPplotsPres, aes(interval, distance, linetype=plot_id)) +
  geom_point(aes(color=treatment)) + stat_smooth(method = "loess", se = F, size = 2, aes(color=treatment)) +
  xlab('Time Interval') +
  ylab('Euclidean Distance') +
  scale_linetype(guide='none') +
  scale_color_manual(values=c('#FF3333', '#3300FF'), name='Treatment', labels=c('control', '+N')) +
  annotate('text', x=0, y=5, label='(b) Presence-based', size=10, hjust='left') +
  theme(legend.position='none')

rateChangeIntervalPplotsPresFigB <- ggplot(rateChangeIntervalPplotsPres, aes(interval, distance)) +
  geom_point(aes(color=treatment)) + stat_smooth(method = "loess", se = F, size = 2, aes(color=treatment)) +
  xlab('Time Interval') +
  ylab('Euclidean Distance') +
  scale_linetype(guide='none') +
  scale_color_manual(values=c('#FF3333', '#3300FF'), name='Treatment', labels=c('control', '+N')) +
  annotate('text', x=0, y=5, label='(d) Presence-based', size=10, hjust='left') +
  theme(legend.position='none')

pushViewport(viewport(layout=grid.layout(2,2)))
print(rateChangeIntervalPplotsAbundFig, vp=viewport(layout.pos.row=1, layout.pos.col=1))
print(rateChangeIntervalPplotsPresFig, vp=viewport(layout.pos.row=1, layout.pos.col=2))
print(rateChangeIntervalPplotsAbundFigB, vp=viewport(layout.pos.row=2, layout.pos.col=1))
print(rateChangeIntervalPplotsPresFigB, vp=viewport(layout.pos.row=2, layout.pos.col=2))
#export as 1500x1500


