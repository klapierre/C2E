library(tidyverse)
library(grid)
library(PerformanceAnalytics)
library(lme4)
library(lmerTest)
library(MuMIn)

#kim's desktop
setwd('C:\\Users\\komatsuk\\Dropbox (Smithsonian)\\working groups\\CoRRE\\converge_diverge\\datasets\\LongForm')

#source data management code  -- if on desktop, change source code
source('C:\\Users\\komatsuk\\Desktop\\R files\\C2E\\SEM paper\\C2E_SEM_data processing.R')

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
#R2m (marginal)=0.06608303, R2c (conditional)=0.06608303; ***marginal R2 estimate is based on only fixed factors, while conditional R2 estimate is based on full mode (fixed and random effects)


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
anova #evenness, species diff (marginally)
r.squaredGLMM(tempModelAll)
#R2m (marginal)=0.3581464, R2c (conditional)=0.5777553
