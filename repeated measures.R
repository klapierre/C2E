library(tidyverse)
library(lme4)
library(car)

dat<-read.csv("C:\\Users\\wilco\\Dropbox\\C2E\\Products\\CommunityChange\\March2018 WG\\CORRE_RACS_Subset_Perm.csv")

pplots<-dat%>%
  filter(site_project_comm == "CDR_e001_D"& treatment %in% c(6, 9)) %>%
  mutate(abs_evenness_change = abs(evenness_change))

#repeated measures anova
summary(rich<-lmer(richness_change~as.factor(treatment_year)*treatment + (1|plot_id), data = pplots))
Anova(rich)

rank<-lmer(rank_change~as.factor(treatment_year)*treatment + (1|plot_id), data = pplots)
Anova(rank)

evenness<-lmer(abs_evenness_change~as.factor(treatment_year)*treatment + (1|plot_id), data = pplots)
Anova(evenness)
ggplot(pplots, aes(x=treatment_year, y=abs_evenness_change, col=factor(treatment) )) +
  geom_point()+
  geom_smooth()
