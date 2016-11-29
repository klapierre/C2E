library(plyr)
library(dplyr)
library(tidyr)
library(codyn)
library(ggplot2)

#kim's wd
setwd('C:\\Users\\Kim\\Dropbox\\working groups\\converge diverge working group\\converge_diverge\\datasets\\LongForm\\to yang')

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

#irrigation plots raw abundance subset
rawAbundanceE001 <- rawAbundance%>%
  filter(project_name=='e001')

#irrigation plots raw abundance subset
rawAbundanceE001 <- rawAbundance%>%
  filter(project_name=='e001')


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
ggplot(rateChangeIntervalIrrAbund, aes(interval, distance, linetype=plot_id)) +
  geom_point(aes(color=treatment)) + theme_bw() + stat_smooth(method = "lm", se = F, size = 2, aes(color=treatment)) +
  xlab('Time Interval') +
  ylab('Euclidean Distance') +
  scale_linetype(guide='none') +
  scale_color_manual(values=c('#FF3333', '#3300FF'), name='Treatment', labels=c('Control', 'Irrigation'))

###presence/absence data
#rate change slopes
rateChangeIrrPres <- rate_change(presenceAbsenceIrr, time.var='treatment_year', species.var='genus_species', abundance.var='presence', replicate.var='plot_id')

#rate change intervals for visualization
rateChangeIntervalIrrPres <- rate_change_interval(presenceAbsenceIrr, time.var='treatment_year', species.var='genus_species', abundance.var='presence', replicate.var='plot_id')%>%
  mutate(project_name='IRG')%>%
  left_join(trt)%>%
  mutate(type='presence')

#rate change figure
ggplot(rateChangeIntervalIrrPres, aes(interval, distance, linetype=plot_id)) +
  geom_point(aes(color=treatment)) + theme_bw() + stat_smooth(method = "lm", se = F, size = 2, aes(color=treatment)) +
  xlab('Time Interval') +
  ylab('Euclidean Distance') +
  scale_linetype(guide='none') +
  scale_color_manual(values=c('#FF3333', '#3300FF'), name='Treatment', labels=c('Control', 'Irrigation'))







