library(grid)
library(PerformanceAnalytics)
library(lme4)
library(lmerTest)
library(r2glmm)
library(MuMIn)
library(tidyverse)

#kim's desktop
setwd('C:\\Users\\komatsuk\\Dropbox (Smithsonian)\\working groups\\CoRRE\\converge_diverge\\datasets\\LongForm')
#kim's laptop
setwd('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\CoRRE\\converge_diverge\\datasets\\LongForm')

#import data
correSEMdataTrt <- read.csv('CoRRE_comm and anpp diff_04182021.csv')

###all treatments
###using absolute value of anpp percent difference
#temporal trends in ANPP -- tried polynomial equation (AIC=625.7441), which was better than linear (AIC=632.3642)
modelTime <- lmer(abs(anpp_pdiff) ~ poly(treatment_year,2) + (1|site_project_comm),
                  data=subset(correSEMdataTrt, treatment_year<11))
summary(modelTime)
anova(modelTime, type='I') #composition difference
r.squaredGLMM(modelTime)
AIC(modelTime)
#R2m (marginal)=0.02102201, R2c (conditional)=0.1909004; ***marginal R2 estimate is based on only fixed factors, while conditional R2 estimate is based on full model (fixed and random effects)
ggplot(data=subset(correSEMdataTrt, treatment_year<11), aes(x=treatment_year, y=abs(anpp_pdiff))) +
  # geom_point(color='grey') +
  geom_smooth(method='lm', formula=y~x+I(x^2), se=F, color='grey',aes(group=interaction(site_project_comm,treatment))) +
  geom_smooth(method='lm', formula=y~x+I(x^2), se=F, size=3, color='black') +
  xlab('Treatment Year') + ylab('|ANPP Difference|') +
  theme(legend.position='none')
#export 1000x1000

# #overall compositional difference -- tried polynomial equation (AIC=260.8949), which was better than linear (AIC=263.6193)
# modelComp <- lmer(abs(anpp_pdiff) ~ poly(composition_diff,2) + (1|site_project_comm),
#                  data=subset(correSEMdataTrt, treatment_year<11))
# summary(modelComp)
# anova(modelComp, type='I') #composition difference
# r.squaredGLMM(modelComp)
# AIC(modelComp)
# #R2m (marginal)=0.1090317, R2c (conditional)=0.352812; ***marginal R2 estimate is based on only fixed factors, while conditional R2 estimate is based on full model (fixed and random effects)
# 
# ggplot(data=subset(correSEMdataTrt, treatment_year<11), aes(x=composition_diff, y=abs(anpp_pdiff))) +
#   # geom_point(color='grey') +
#   geom_smooth(method='lm', formula=y~x+I(x^2), se=F, color='grey',aes(group=interaction(site_project_comm, treatment))) +
#   geom_smooth(method='lm', formula=y~x+I(x^2), se=F, size=3, color='black') +
#   xlab('Composition Difference') + ylab('|ANPP Difference|') +
#   theme(legend.position='none') +
#   coord_cartesian(ylim=c(0,2.0))
# 
# 
# #specific community difference metrics -- partial polynomial equation (AIC=288.838), which was better than linear (AIC=300.6363) or fully polynomial (AIC=292.117)
# modelAll <- lmer(abs(anpp_pdiff) ~ evenness_diff + rank_difference + poly(species_difference,2) + poly(richness_difference,2) + (1|site_project_comm),
#                   data=subset(correSEMdataTrt, treatment_year<11))
# summary(modelAll)
# anova(modelAll, type='I') #rank, richness, and spp diff
# r.squaredGLMM(modelAll)
# AIC(modelAll)
# r2 = r2beta(model=modelAll,partial=TRUE,method='sgv')
# #R2m (marginal)=0.1833687, R2c (conditional)=0.3427348; ***marginal R2 estimate is based on only fixed factors, while conditional R2 estimate is based on full model (fixed and random effects)
# 
# 
# 
# evennessPlot <- ggplot(data=subset(correSEMdataTrt, treatment_year<11), aes(x=evenness_diff, y=abs(anpp_pdiff))) +
#   # geom_point(color='grey') +
#   geom_smooth(method='lm', se=F, color='grey',aes(group=interaction(site_project_comm, treatment))) +
#   # geom_smooth(method='lm', se=F, size=3, color='black') +
#   xlab('Evenness Difference') + ylab('|ANPP Difference|') +
#   theme(legend.position='none') +
#   annotate('text', x=-0.3, y=2.2, label='ns', size=8, hjust='left')
# rankPlot <- ggplot(data=subset(correSEMdataTrt, treatment_year<11), aes(x=rank_difference, y=abs(anpp_pdiff))) +
#   # geom_point(color='grey') +
#   geom_smooth(method='lm', se=F, color='grey',aes(group=interaction(site_project_comm, treatment))) +
#   geom_smooth(method='lm', se=F, size=3, color='black') +
#   xlab('Rank Difference') + ylab('|ANPP Difference|') +
#   theme(legend.position='none') +
#   annotate('text', x=0.07, y=2.2, label=expression(paste(R^2,'=0.004',sep='')), size=8, hjust='left')
# sppdiffPlot <- ggplot(data=subset(correSEMdataTrt, treatment_year<11), aes(x=species_difference, y=abs(anpp_pdiff))) +
#   # geom_point(color='grey') +
#   geom_smooth(method='lm', formula=y~x+I(x^2), se=F, color='grey',aes(group=interaction(site_project_comm, treatment))) +
#   geom_smooth(method='lm', formula=y~x+I(x^2), se=F, size=3, color='black') +
#   xlab('Species Difference') + ylab('|ANPP Difference|') +
#   theme(legend.position='none') +
#   annotate('text', x=0.13, y=2.2, label=expression(paste(R^2,'=0.036',sep='')), size=8, hjust='left')
# richnessPlot <- ggplot(data=subset(correSEMdataTrt, treatment_year<11), aes(x=richness_difference, y=abs(anpp_pdiff))) +
#   # geom_point(color='grey') +
#   geom_smooth(method='lm', formula=y~x+I(x^2), se=F, color='grey',aes(group=interaction(site_project_comm, treatment))) +
#   geom_smooth(method='lm', formula=y~x+I(x^2), se=F, size=3, color='black') +
#   xlab('Richness Difference') + ylab('|ANPP Difference|') +
#   theme(legend.position='none') +
#   annotate('text', x=-0.6, y=2.2, label=expression(paste(R^2,'=0.108',sep='')), size=8, hjust='left')
# 
# pushViewport(viewport(layout=grid.layout(2,2)))
# print(evennessPlot, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
# print(rankPlot, vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
# print(sppdiffPlot, vp=viewport(layout.pos.row = 2, layout.pos.col = 1))
# print(richnessPlot, vp=viewport(layout.pos.row = 2, layout.pos.col = 2))
# #export at 2000x1000





###overall compositional difference, gains/losses included -- tried polynomial equation (AIC=474.742), which was better than linear (AIC=480.3047)
modelComp <- lmer(abs(anpp_pdiff) ~ poly(composition_diff,2) + (1|site_project_comm/treatment),
                  data=subset(correSEMdataTrt, treatment_year<11 & !is.na(composition_diff)))
summary(modelComp)
anova(modelComp, type='I') #composition difference
r.squaredGLMM(modelComp)
AIC(modelComp)
#R2m (marginal)=0.1482056, R2c (conditional)=0.3975152; ***marginal R2 estimate is based on only fixed factors, while conditional R2 estimate is based on full model (fixed and random effects)

ggplot(data=subset(correSEMdataTrt, treatment_year<11), aes(x=composition_diff, y=abs(anpp_pdiff))) +
  geom_smooth(method='lm', se=F, color='grey',aes(group=interaction(site_project_comm, treatment))) +
  geom_smooth(method='lm', formula=y~x+I(x^2), se=F, size=3, color='black') +
  xlab('Composition Difference') + ylab('|ANPP Difference|') +
  theme(legend.position='none') +
  coord_cartesian(ylim=c(0,2.0))
#export at 1000x1000

#specific community change metrics -- fully polynomial (AIC=410.2774), which was better than linear (AIC=436.9822)
modelAll <- lmer(abs(anpp_pdiff) ~ poly(evenness_diff,2) + poly(rank_diff,2) + poly(losses,2) + poly(gains,2) +poly(richness_diff,2) + (1|site_project_comm/treatment),
                 data=subset(correSEMdataTrt, treatment_year<11&!is.na(evenness_diff)&!is.na(rank_diff)&!is.na(gains)&!is.na(losses)))
summary(modelAll)
anova(modelAll, type='I') #evenness diff, rank diff, losses, richness
r.squaredGLMM(modelAll)
AIC(modelAll)
r2 = r2beta(model=modelAll,partial=TRUE,method='sgv')
#R2m (marginal)=0.2258467, R2c (conditional)=0.4996504; ***marginal R2 estimate is based on only fixed factors, while conditional R2 estimate is based on full model (fixed and random effects)
# partial R2 :
# evenness - 0.964 / (0.9864+2.5351+7.1190+0.0064+1.9795) = 0.076
# rank - 2.5351 / (0.9864+2.5351+7.1190+0.0064+1.9795) = 0.201
# losses - 7.1190 / (0.9864+2.5351+7.1190+0.0064+1.9795) = 0.564
# gains - 0.0064 / (0.9864+2.5351+7.1190+0.0064+1.9795) = 0.0005
# richness - 1.9795 / (0.9864+2.5351+7.1190+0.0064+1.9795) = 0.157

evennessPlot <- ggplot(data=subset(correSEMdataTrt, treatment_year<11), aes(x=evenness_diff, y=abs(anpp_pdiff))) +
  geom_smooth(method='lm', se=F, color='grey', aes(group=interaction(site_project_comm, treatment))) +
  geom_smooth(method='lm', formula=y~x+I(x^2), se=F, size=3, color='black') +
  xlab('Evenness Difference') + ylab('|ANPP Difference|') +
  theme(legend.position='none') +
  annotate('text', x=-0.3, y=2.4, label=expression(paste('partial ', R^2,'=0.076',sep='')), size=8, hjust='left') +
  coord_cartesian(ylim=c(0,2.5))
rankPlot <- ggplot(data=subset(correSEMdataTrt, treatment_year<11), aes(x=rank_diff, y=abs(anpp_pdiff))) +
  geom_smooth(method='lm', se=F, color='grey', aes(group=interaction(site_project_comm, treatment))) +
  geom_smooth(method='lm', formula=y~x+I(x^2), se=F, size=3, color='black') +
  xlab('Rank Difference') + ylab('|ANPP Difference|') +
  theme(legend.position='none') +
  annotate('text', x=0.01, y=2.4, label=expression(paste('partial ', R^2,'=0.201',sep='')), size=8, hjust='left') +
  coord_cartesian(ylim=c(0,2.5))
gainsPlot <- ggplot(data=subset(correSEMdataTrt, treatment_year<11), aes(x=gains, y=abs(anpp_pdiff))) +
  geom_smooth(method='lm', se=F, color='grey', aes(group=interaction(site_project_comm, treatment))) +
  # geom_smooth(method='lm', formula=y~x+I(x^2), se=F, size=3, color='black') +
  xlab('Species Gains') + ylab('|ANPP Difference|') +
  theme(legend.position='none') +
  annotate('text', x=0.05, y=2.4, label='n.s.', size=8, hjust='left') +
  coord_cartesian(ylim=c(0,2.5))
lossesPlot <- ggplot(data=subset(correSEMdataTrt, treatment_year<11), aes(x=losses, y=abs(anpp_pdiff))) +
  geom_smooth(method='lm', se=F, color='grey', aes(group=interaction(site_project_comm, treatment))) +
  geom_smooth(method='lm', formula=y~x+I(x^2), se=F, size=3, color='black') +
  xlab('Species Losses') + ylab('|ANPP Difference|') +
  theme(legend.position='none') +
  annotate('text', x=0.13, y=2.4, label=expression(paste('partial ', R^2,'=0.564',sep='')), size=8, hjust='left') +
  coord_cartesian(ylim=c(0,2.5))
richnessPlot <- ggplot(data=subset(correSEMdataTrt, treatment_year<11), aes(x=richness_diff, y=abs(anpp_pdiff))) +
  geom_smooth(method='lm', se=F, color='grey', aes(group=interaction(site_project_comm, treatment))) +
  geom_smooth(method='lm', formula=y~x+I(x^2), se=F, size=3, color='black') +
  xlab('Richness Difference') + ylab('|ANPP Difference|') +
  theme(legend.position='none') +
  annotate('text', x=-0.85, y=2.4, label=expression(paste('partial ', R^2,'=0.157',sep='')), size=8, hjust='left') +
  coord_cartesian(ylim=c(0,2.5))

pushViewport(viewport(layout=grid.layout(2,2)))
print(evennessPlot, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
print(rankPlot, vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
print(gainsPlot, vp=viewport(layout.pos.row = 2, layout.pos.col = 1))
print(lossesPlot, vp=viewport(layout.pos.row = 2, layout.pos.col = 2))
#export at 2000x1000


###richness only
modelRich <- lmer(abs(anpp_pdiff) ~ poly(richness_diff,2) + (1|site_project_comm/treatment),
                 data=subset(correSEMdataTrt, treatment_year<11&!is.na(richness_diff)))
anova(modelRich, type='I')
r.squaredGLMM(modelRich) #marginal R2=0.1980784, conditional R2=0.1980784
AIC(modelRich)

ggplot(data=subset(correSEMdataTrt, treatment_year<11), aes(x=richness_diff, y=abs(anpp_pdiff))) +
  # geom_point(color='grey') +
  geom_smooth(method='lm', se=F, color='grey',aes(group=interaction(site_project_comm, treatment))) +
  geom_smooth(method='lm', formula=y~x+I(x^2), se=F, size=3, color='black') +
  xlab('Richness Difference') + ylab('|ANPP Difference|') +
  theme(legend.position='none') +
  annotate('text', x=-0.85, y=2.2, label=expression(paste(R^2,'=0.198',sep='')), size=8, hjust='left')
#export at 1000x1000




###exploratory correlations and histograms (all variables compare treatment to control plots)
dataVis <- correSEMdataTrt%>%
  select(anpp_pdiff, composition_diff, richness_diff, evenness_diff, rank_diff, species_diff) #make visualization dataframe
chart.Correlation(dataVis, histogram=T, pch=19)

dataVis <- correSEMdataTrt%>%
  select(anpp_pdiff, composition_change, richness_change, evenness_change, rank_change, gains, losses) #make visualization dataframe
chart.Correlation(dataVis, histogram=T, pch=19)





##################
### By treatment types
##################

###n treatments only
Ntrt <- correSEMdataTrt%>%
  filter(trt_type=='N'&n>0)

###sequential sum of squares from regression (N trt first, then community metrics)
NmodelAll <- lmer(anpp_pdiff_transform ~ n + 
                 evenness_diff + rank_difference + species_difference + richness_difference + 
                 (1|site_code/project_name),
               data=Ntrt)
summary(NmodelAll)
anova(NmodelAll, type='I') #richness
r.squaredGLMM(NmodelAll)
#R2m (marginal)=0.1683844, R2c (conditional)=0.1842924


###p treatments only
Ptrt <- correSEMdataTrt%>%
  filter(trt_type=='P'&p>0)

###sequential sum of squares from regression (P trt first, then community metrics)
PmodelAll <- lmer(anpp_pdiff_transform ~ p + 
                    evenness_diff + rank_difference + species_difference + richness_difference + 
                    (1|site_code/project_name),
                  data=Ptrt)
summary(PmodelAll)
anova(PmodelAll, type='I')
r.squaredGLMM(PmodelAll)
#R2m (marginal)=0.06608303, R2c (conditional)=0.06608303


###irr treatments only
irrTrt <- correSEMdataTrt%>%
  filter(trt_type=='irr')

###sequential sum of squares from regression (P trt first, then community metrics)
irrModelAll <- lmer(anpp_pdiff_transform ~ precip + 
                    evenness_diff + rank_difference + species_difference + richness_difference + 
                    (1|site_code/project_name),
                  data=irrTrt)
summary(irrModelAll)
anova(irrModelAll, type='I')
r.squaredGLMM(irrModelAll) #irrigation, rank
#R2m (marginal)=0.3853351, R2c (conditional)=0.3853351


###drought treatments only
droTrt <- correSEMdataTrt%>%
  filter(trt_type=='drought')

###sequential sum of squares from regression (P trt first, then community metrics)
droModelAll <- lmer(anpp_pdiff_transform ~ precip + 
                      evenness_diff + rank_difference + species_difference + richness_difference + 
                      (1|site_code/project_name),
                    data=droTrt)
summary(droModelAll)
anova(droModelAll, type='I')
r.squaredGLMM(droModelAll)
#R2m (marginal)=0.2083582, R2c (conditional)=0.2083582


###temperature treatments only
tempTrt <- correSEMdataTrt%>%
  filter(trt_type=='temp')

###sequential sum of squares from regression (P trt first, then community metrics)
tempModelAll <- lmer(anpp_pdiff_transform ~ temp + 
                      evenness_diff + rank_difference + species_difference + richness_difference + 
                      (1|site_code/project_name),
                    data=tempTrt)
summary(tempModelAll)
anova(tempModelAll, type='I') #evenness, species diff (marginally)
r.squaredGLMM(tempModelAll)
#R2m (marginal)=0.3581464, R2c (conditional)=0.5777553
