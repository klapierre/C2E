library(tidyverse)
library(lme4)
library(car)

dat<-read.csv("~/Dropbox/C2E/Products/CommunityChange/March2018 WG/CORRE_RACS_Subset_Perm.csv")

pplots<-dat%>%
  filter(site_project_comm == "KNZ_pplots_0"& treatment %in% c("N1P0", "N2P0"))

#repeated measures anova
summary(rich<-lmer(cbind(richness_change+rank_change)~as.factor(treatment_year)*treatment + (1|plot_id), data = pplots))
Anova(rich)

rank<-lmer(rank_change~as.factor(treatment_year)*treatment + (1|plot_id), data = pplots)
anova(rank)


##looping through to run RM-anova for each c-t comparision
spc<-unique(dat$site_project_comm)
rm.anova_output<-data.frame()

for (i in 1:length(spc)){
  
  subset<-dat%>%
    filter(site_project_comm==spc[i])
  
  control<-subset%>%
    filter(plot_mani==0)
  
  treat_list<-unique(subset(subset, plot_mani>0)$treatment)
  
  for(j in 1:length(treat_list)) {
    spc1<-spc[i]
    trt <- treat_list[j]
    
    subset_trt<-subset%>%
      filter(treatment==treat_list[j])
    
    #dataset of two treatments    
    subset_t12<-rbind(control, subset_trt)
    
    out<-lmer(cbind(richness_change+evenness_change+rank_change+gains+losses)~as.factor(treatment_year)*treatment + (1|plot_id), data = subset_t12)
    
    output<-data.frame(site_project_comm = spc1,
                       treatment = trt,
                       p.value = out$p.value)
    
    ttest_output <- rbind(ttest_output, output)
  }}
