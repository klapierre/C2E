library(tidyverse)
library(lme4)
library(car)

dat<-read.csv("~/Dropbox/C2E/Products/CommunityChange/March2018 WG/CORRE_RACS_Subset_Perm.csv")

pplots<-dat%>%
  filter(site_project_comm == "KNZ_pplots_0"& treatment %in% c("N1P0", "N2P0"))

#repeated measures anova
summary(rich<-lmer(richness_change~as.factor(treatment_year)*treatment + (1|plot_id), data = pplots))
Anova(rich)

rank<-lmer(rank_change~as.factor(treatment_year)*treatment + (1|plot_id), data = pplots)
anova(rank)

