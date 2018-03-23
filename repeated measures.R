library(tidyverse)
library(lme4)
library(car)

dat<-read.csv("C:\\Users\\wilco\\Dropbox\\C2E\\Products\\CommunityChange\\March2018 WG\\CORRE_RACS_Subset_Perm.csv")
dat.b<-
pplots<-dat%>%
  filter(site_project_comm == "CDR_e001_D"& treatment %in% c(6, 9)) %>%
  mutate(abs_evenness_change = abs(evenness_change))

#repeated measures anova
summary(rich<-lmer(cbind(richness_change+rank_change)~as.factor(treatment_year)*treatment + (1|plot_id), data = pplots))
Anova(rich)

rank<-lmer(rank_change~as.factor(treatment_year)*treatment + (1|plot_id), data = pplots)
Anova(rank)


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
    #repeated measures for richness
    rm_s<-lmer(richness_change~as.factor(treatment_year)*treatment + (1|plot_id), data = subset_t12)
        a_s2<-Anova(rm_s)
        S_out<-data.frame(site_project_comm = spc1,
                       treatment = trt,
                       S_year = a_s2[1,3],
                       S_trt = a_s2[2,3],
                       S_trt.yr = a_s2[3,3])
      #repeated measuers for evenness
        rm_e<-lmer(evenness_change~as.factor(treatment_year)*treatment + (1|plot_id), data = subset_t12)
        a_e2<-Anova(rm_e)
        E_out<-data.frame(site_project_comm = spc1,
                          treatment = trt,
                          E_year = a_e2[1,3],
                          E_trt = a_e2[2,3],
                          E_trt.yr = a_e2[3,3])
        #repeated measuers for Rank
        rm_r<-lmer(rank_change~as.factor(treatment_year)*treatment + (1|plot_id), data = subset_t12)
        a_r2<-Anova(rm_r)
        R_out<-data.frame(site_project_comm = spc1,
                          treatment = trt,
                          R_year = a_r2[1,3],
                          R_trt = a_r2[2,3],
                          R_trt.yr = a_r2[3,3])
        #repeated measuers for gain
        rm_g<-lmer(gains~as.factor(treatment_year)*treatment + (1|plot_id), data = subset_t12)
        a_g2<-Anova(rm_g)
        G_out<-data.frame(site_project_comm = spc1,
                          treatment = trt,
                          G_year = a_g2[1,3],
                          G_trt = a_g2[2,3],
                          G_trt.yr = a_g2[3,3])
        #repeated measuers for evenness
        rm_l<-lmer(losses~as.factor(treatment_year)*treatment + (1|plot_id), data = subset_t12)
        a_l2<-Anova(rm_l)
        L_out<-data.frame(site_project_comm = spc1,
                          treatment = trt,
                          L_year = a_l2[1,3],
                          L_trt = a_l2[2,3],
                          L_trt.yr = a_l2[3,3])
    
        rmresults<-S_out%>%
          left_join(E_out)%>%
          left_join(R_out)%>%
          left_join(G_out)%>%
          left_join(L_out)
        
        rm.anova_output<-rbind(rm.anova_output, rmresults)
        
  }}

#pull out data on trt only

trt_effects<-rm.anova_output%>%
  mutate(richness = ifelse(S_trt<0.05, 1, 0),
         evenness = ifelse(E_trt<0.05, 1, 0),
         rank = ifelse(R_trt<0.05, 1, 0),
         gain = ifelse(G_trt<0.05, 1, 0),
         loss = ifelse(L_trt<0.05, 1, 0))%>%
  select(site_project_comm, treatment, richness, evenness, rank, gain, loss )%>%
  mutate(site_treatment = paste(site_project_comm, treatment, sep = "_"))%>%
  gather(metric, value, richness:loss)

summed<-trt_effects%>%
  group_by(metric)%>%
  summarize(sum=sum(value))

ggplot(trt_effects, aes(y = site_treatment, x = metric))+
  geom_tile(aes(fill = as.factor(value)))+
  scale_fill_brewer(type = "qual", 
                    labels = c("C and T not different", "C and T different", "NA"), 
                    name = NULL)+
  xlab("Metric")+
  ylab("Site and Treatment")
ggsave(filename ="C:\\Users\\wilco\\Dropbox\\C2E\\Products\\CommunityChange\\March2018 WG\\rmanoval_figure_trt.pdf",
       height = 14, 
       width = 7, 
       units = "in")

trt_effects_trt_yr<-rm.anova_output%>%
  mutate(richness = ifelse(S_trt.yr<0.05, 1, 0),
         evenness = ifelse(E_trt.yr<0.05, 1, 0),
         rank = ifelse(R_trt.yr<0.05, 1, 0),
         gain = ifelse(G_trt.yr<0.05, 1, 0),
         loss = ifelse(L_trt.yr<0.05, 1, 0))%>%
  select(site_project_comm, treatment, richness, evenness, rank, gain, loss )%>%
  mutate(site_treatment = paste(site_project_comm, treatment, sep = "_"))%>%
  gather(metric, value, richness:loss)

summed<-trt_effects_trt_yr%>%
  group_by(metric)%>%
  summarize(sum=sum(value))

ggplot(trt_effects_trt_yr, aes(y = site_treatment, x = metric))+
  geom_tile(aes(fill = as.factor(value)))+
  scale_fill_brewer(type = "qual", 
                    labels = c("C and T not different", "C and T different", "NA"), 
                    name = NULL)+
  xlab("Metric")+
  ylab("Site and Treatment")
ggsave(filename ="C:\\Users\\wilco\\Dropbox\\C2E\\Products\\CommunityChange\\March2018 WG\\rmanoval_figure_trt_yr.pdf",
height = 14, 
width = 7, 
units = "in")
