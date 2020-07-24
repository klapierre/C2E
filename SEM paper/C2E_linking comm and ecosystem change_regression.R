library(tidyverse)
library(grid)
library(PerformanceAnalytics)
library(lme4)
library(lmerTest)
library(MuMIn)

#kim's desktop
setwd('C:\\Users\\komatsuk\\Dropbox (Smithsonian)\\working groups\\CoRRE\\converge_diverge\\datasets\\LongForm')
#kim's laptop
setwd('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\CoRRE\\converge_diverge\\datasets\\LongForm')

#source data management code  -- if on desktop, change source code
source('C:\\Users\\lapie\\Desktop\\R files laptop\\C2E\\SEM paper\\C2E_SEM_data processing.R')

###all treatments
###using absolute value of anpp percent difference
#temporal trends in ANPP
modelTime <- lmer(abs(anpp_pdiff) ~ treatment_year + (1|site_code/project_name),
                  data=subset(correSEMdataTrt, treatment_year<11))
summary(modelTime)
anova(modelTime, type='I') #composition difference
r.squaredGLMM(modelTime)
#R2m (marginal)=0.007671168, R2c (conditional)=0.1539671; ***marginal R2 estimate is based on only fixed factors, while conditional R2 estimate is based on full model (fixed and random effects)
ggplot(data=subset(correSEMdataTrt, treatment_year<11), aes(x=treatment_year, y=abs(anpp_pdiff))) +
  geom_point(color='grey') +
  geom_smooth(method='lm', formula=y~x+I(x^2), se=F, color='grey',aes(group=site_project_comm)) +
  geom_smooth(method='lm', formula=y~x+I(x^2), se=F, size=3, color='black') +
  xlab('Treatment Year') + ylab('|ANPP Difference|') +
  theme(legend.position='none')


#overall compositional difference
modelComp <- lmer(abs(anpp_pdiff) ~ composition_diff + (1|site_code/project_name),
                 data=subset(correSEMdataTrt, treatment_year<11))
summary(modelComp)
anova(modelComp, type='I') #composition difference
r.squaredGLMM(modelComp)
#R2m (marginal)=0.1826489, R2c (conditional)=0.318324; ***marginal R2 estimate is based on only fixed factors, while conditional R2 estimate is based on full model (fixed and random effects)

ggplot(data=subset(correSEMdataTrt, treatment_year<11), aes(x=composition_diff, y=abs(anpp_pdiff))) +
  geom_point(color='grey') +
  geom_smooth(method='lm', se=F, color='grey',aes(group=site_project_comm)) +
  geom_smooth(method='lm', se=F, size=3, color='black') +
  xlab('Composition Difference') + ylab('|ANPP Difference|') +
  theme(legend.position='none')


#specific community difference metrics
modelAll <- lmer(abs(anpp_pdiff) ~ evenness_diff + rank_difference + species_difference + richness_difference + (1|site_code/project_name),
                  data=subset(correSEMdataTrt, treatment_year<11))
summary(modelAll)
anova(modelAll, type='I') #rank, richness, and spp diff
r.squaredGLMM(modelAll)
#R2m (marginal)=0.2096342, R2c (conditional)=0.3422267; ***marginal R2 estimate is based on only fixed factors, while conditional R2 estimate is based on full model (fixed and random effects)

evennessPlot <- ggplot(data=subset(correSEMdataTrt, treatment_year<11), aes(x=evenness_diff, y=abs(anpp_pdiff))) +
  geom_point(color='grey') +
  geom_smooth(method='lm', se=F, color='grey',aes(group=site_project_comm)) +
  geom_smooth(method='lm', se=F, size=3, color='black') +
  xlab('Evenness Difference') + ylab('|ANPP Difference|') +
  theme(legend.position='none')
rankPlot <- ggplot(data=subset(correSEMdataTrt, treatment_year<11), aes(x=rank_difference, y=abs(anpp_pdiff))) +
  geom_point(color='grey') +
  geom_smooth(method='lm', se=F, color='grey',aes(group=site_project_comm)) +
  geom_smooth(method='lm', se=F, size=3, color='black') +
  xlab('Rank Difference') + ylab('|ANPP Difference|') +
  theme(legend.position='none')
sppdiffPlot <- ggplot(data=subset(correSEMdataTrt, treatment_year<11), aes(x=species_difference, y=abs(anpp_pdiff))) +
  geom_point(color='grey') +
  geom_smooth(method='lm', se=F, color='grey',aes(group=site_project_comm)) +
  geom_smooth(method='lm', se=F, size=3, color='black') +
  xlab('Species Difference') + ylab('|ANPP Difference|') +
  theme(legend.position='none')
richnessPlot <- ggplot(data=subset(correSEMdataTrt, treatment_year<11), aes(x=richness_difference, y=abs(anpp_pdiff))) +
  geom_point(color='grey') +
  geom_smooth(method='lm', se=F, color='grey',aes(group=site_project_comm)) +
  geom_smooth(method='lm', se=F, size=3, color='black') +
  xlab('Richness Difference') + ylab('|ANPP Difference|') +
  theme(legend.position='none')

pushViewport(viewport(layout=grid.layout(2,2)))
print(evennessPlot, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
print(rankPlot, vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
print(sppdiffPlot, vp=viewport(layout.pos.row = 2, layout.pos.col = 1))
print(richnessPlot, vp=viewport(layout.pos.row = 2, layout.pos.col = 2))












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
