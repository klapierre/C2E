library(plyr)
library(dplyr)
library(tidyr)
library(codyn)
library(ggplot2)

#kim's wd
setwd('C:\\Users\\Kim\\Dropbox\\working groups\\converge diverge working group\\converge_diverge\\datasets\\LongForm\\to yang')

#load data
rawAbundance <- read.csv('CORRE_irrigationL_e001D_subset.csv')%>%
  select(-X)

presenceAbsence <- rawAbundance%>%
  mutate()

trt <- rawAbundance%>%
  select(site_code, project_name, treatment, plot_id, n, p, k, precip, other_trt)%>%
  unique()%>%
  mutate(plot_id=as.character(plot_id))

rawAbundanceIrr <- rawAbundance%>%
  filter(project_name=='IRG')

rateChange <- rate_change(rawAbundanceIrr, time.var='treatment_year', species.var='genus_species', abundance.var='abundance', replicate.var='plot_id')

rateChangeInterval <- rate_change_interval(rawAbundanceIrr, time.var='treatment_year', species.var='genus_species', abundance.var='abundance', replicate.var='plot_id')%>%
  mutate(project_name='IRG')%>%
  left_join(trt)

ggplot(rateChangeInterval, aes(interval, distance, linetype=plot_id)) +
  geom_point(aes(color=treatment)) + theme_bw() + stat_smooth(method = "lm", se = F, size = 2, aes(color=treatment)) +
  xlab('Time Interval') +
  ylab('Euclidean Distance') +
  scale_linetype(guide='none') +
  scale_color_manual(values=c('#FF3333', '#3300FF'), name='Treatment', labels=c('Control', 'Irrigation'))









