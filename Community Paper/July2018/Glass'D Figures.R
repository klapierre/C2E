library(tidyverse)
library(ggthemes)
library(grid)
library(gridExtra)

#meghan's working directory
setwd("C:\\Users\\megha\\Dropbox\\")
setwd("~/Dropbox/")

theme_set(theme_bw(12))

### Read in data 

sig<-read.csv("C2E/Products/CommunityChange/Summer2018_Results/gam_comparison_table.csv")%>%
  select(site_proj_comm, treatment, response_var, final_treatment_year, sig_diff_cntrl_trt)


change_metrics <- read.csv("C2E\\Products\\CommunityChange\\March2018 WG\\MetricsTrts_July2018.csv") %>%
  mutate(abs_richness_change = abs(richness_change),
         abs_evenness_change = abs(evenness_change))

subset<-read.csv("C2E\\Products\\CommunityChange\\March2018 WG\\experiment_trt_subset.csv")

### Control data
change_control <- change_metrics %>%
  filter(plot_mani==0) %>%
  dplyr::select(treatment_year, treatment_year2, abs_richness_change, abs_evenness_change, 
                rank_change, gains, losses, site_project_comm, treatment, plot_mani) %>%
  rename(abs_richness_change_ctrl = abs_richness_change,
         abs_evenness_change_ctrl = abs_evenness_change,
         rank_change_ctrl = rank_change,
         gains_ctrl = gains,
         losses_ctrl = losses
  ) %>%
  group_by(site_project_comm, treatment_year2) %>%
  summarise_at(vars(abs_richness_change_ctrl:losses_ctrl), funs(mean, sd), na.rm=T)

change_glass_d <- change_metrics %>%
  filter(plot_mani != 0) %>%
  group_by(site_project_comm, treatment, treatment_year2, plot_mani) %>%
  summarise(abs_richness_change = mean(abs_richness_change,na.rm=T),
            abs_evenness_change = mean(abs_evenness_change, na.rm=T),
            rank_change = mean(rank_change, na.rm=T),
            gains = mean(gains, na.rm=T),
            losses = mean(losses, na.rm=T)) %>%
  left_join(change_control, by=c("site_project_comm","treatment_year2")) %>%
  mutate(abs_richness_glass = (abs_richness_change-abs_richness_change_ctrl_mean)/abs_richness_change_ctrl_sd,
         abs_evenness_glass = (abs_evenness_change-abs_evenness_change_ctrl_mean)/abs_evenness_change_ctrl_sd,
         rank_glass = (rank_change-rank_change_ctrl_mean)/rank_change_ctrl_sd,
         gains_glass = (gains-gains_ctrl_mean)/gains_ctrl_sd,
         losses_glass = (losses-losses_ctrl_mean)/losses_ctrl_sd
  )%>%
  select(site_project_comm, treatment, treatment_year2, plot_mani, abs_richness_glass, abs_evenness_glass, rank_glass, gains_glass, losses_glass)

#change_glass_d is the thing that we want

## replace Inf with NAs in change_glass_d and select only treatments I want
change_glass_d <- change_glass_d %>%
  mutate(abs_richness_glass=replace(abs_richness_glass, abs_richness_glass=="Inf"|abs_richness_glass=="NaN", NA)) %>%
  mutate(abs_evenness_glass=replace(abs_evenness_glass, abs_evenness_glass=="Inf"|abs_evenness_glass=="NaN", NA)) %>%
  mutate(rank_glass=replace(rank_glass, rank_glass=="Inf"|rank_glass=="NaN", NA)) %>%
  mutate(gains_glass=replace(gains_glass, gains_glass=="Inf"|gains_glass=="NaN", NA)) %>%
  mutate(losses_glass=replace(losses_glass, losses_glass=="Inf"|losses_glass=="NaN", NA))

# read in treatment variables for subsetting later
info.trt=read.csv("C2E\\Products\\CommunityChange\\March2018 WG\\treatment interactions_July2018.csv") 

### calculate mean change through time and combine with predictor variables
GlassD<-change_glass_d%>%
  select(-plot_mani)%>%
  right_join(subset)%>%
  rename(richness_change_abs=abs_richness_glass,
         evenness_change_abs=abs_evenness_glass,
         rank_change=rank_glass,
         gains=gains_glass,
         losses=losses_glass)%>%
  gather(response_var, glassd, richness_change_abs:losses)%>%
  left_join(sig)%>%
  left_join(info.trt)%>%
  select(-site_proj_comm, -site_code, -project_name, -community_type, -trt_type, -press)

###doing with all years
allyears_all <- GlassD %>%
  group_by(site_project_comm, treatment, response_var, trt_type2, use) %>%
  summarise(mglassd=mean(glassd, na.rm=T))

allyears_sigonly <- GlassD %>%
  filter(sig_diff_cntrl_trt=="yes")%>%
  group_by(site_project_comm, treatment, response_var, trt_type2, use) %>%
  summarise(mglassd=mean(glassd, na.rm=T))

###graphing this ALL YEARS A - all data, B - sig only

##A
glassD_trta<-allyears_all%>%
  group_by(response_var, trt_type2)%>%
  summarize(mean=mean(mglassd),
            sd=sd(mglassd),
            num=length(mglassd))%>%
  mutate(se=sd/sqrt(num))%>%
  filter(trt_type2=="N"|trt_type2=="Mult. Nuts."|trt_type2=="Irrigation"|trt_type2=="CO2"|trt_type2=="P"|trt_type2=="Temperature"|trt_type2=="Irr + Temp")

glassD_alla<-allyears_all%>%
  ungroup()%>%
  group_by(response_var)%>%
  summarize(mean=mean(mglassd),
            sd=sd(mglassd),
            num=length(mglassd))%>%
  mutate(se=sd/sqrt(num))%>%
  mutate(trt_type2="All Trts")

glassD_alldata<-glassD_alla%>%
  bind_rows(glassD_trta)

a<-ggplot(data=glassD_alldata, aes(x=trt_type2, y=mean, fill=trt_type2))+
  geom_bar(position=position_dodge(), stat="identity")+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),position= position_dodge(0.9), width=0.2)+
  ylab("Glass's D")+
  xlab("")+
  scale_x_discrete(limits=c("All Trts","CO2","Irrigation","Temperature","Irr + Temp","N","P","Mult. Nuts."), labels=c("All Trts", "CO2","Irrigation", "Temp","Irr + Temp","Nitrogen","Phosphorus","Mult Nuts"))+
  scale_fill_manual(values=c("black","green3",'purple','blue','darkorange','orange','gold3','red'))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.5))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")+
  geom_vline(xintercept = 1.5, size = 1)+
  facet_wrap(~response_var, ncol=1)+
  ggtitle("All Data")

##B
glassD_trtb<-allyears_sigonly%>%
  group_by(response_var, trt_type2)%>%
  summarize(mean=mean(mglassd),
            sd=sd(mglassd),
            num=length(mglassd))%>%
  mutate(se=sd/sqrt(num))%>%
  filter(trt_type2=="N"|trt_type2=="Mult. Nuts."|trt_type2=="Irrigation"|trt_type2=="CO2"|trt_type2=="P"|trt_type2=="Temperature"|trt_type2=="Irr + Temp")

glassD_allb<-allyears_sigonly%>%
  ungroup()%>%
  group_by(response_var)%>%
  summarize(mean=mean(mglassd),
            sd=sd(mglassd),
            num=length(mglassd))%>%
  mutate(se=sd/sqrt(num))%>%
  mutate(trt_type2="All Trts")

glassD_alldatb<-glassD_allb%>%
  bind_rows(glassD_trtb)

b<-ggplot(data=glassD_alldatb, aes(x=trt_type2, y=mean, fill=trt_type2))+
  geom_bar(position=position_dodge(), stat="identity")+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),position= position_dodge(0.9), width=0.2)+
  ylab("Glass's D")+
  xlab("")+
  scale_x_discrete(limits=c("All Trts","CO2","Irrigation","Temperature","Irr + Temp","N","P","Mult. Nuts."), labels=c("All Trts", "CO2","Irrigation", "Temp","Irr + Temp","Nitrogen","Phosphorus","Mult Nuts"))+
  scale_fill_manual(values=c("black","green3",'purple','blue','darkorange','orange','gold3','red'))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.5))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")+
  geom_vline(xintercept = 1.5, size = 1)+
  facet_wrap(~response_var, ncol=1)+
  ggtitle("Sig Only")

grid.arrange(a, b, ncol=2)


###doing with 5th year only
five_all <- GlassD %>%
  filter(treatment_year2==5)%>%
  group_by(site_project_comm, treatment, response_var, trt_type2, use) %>%
  summarise(mglassd=mean(glassd, na.rm=T))

five_sigonly <- GlassD %>%
  filter(sig_diff_cntrl_trt=="yes")%>%
  filter(treatment_year2==5)%>%
  group_by(site_project_comm, treatment, response_var, trt_type2, use) %>%
  summarise(mglassd=mean(glassd, na.rm=T))

###graphing this 5th year only A - all data, B - sig only

##A
glassD_trta<-five_all%>%
  group_by(response_var, trt_type2)%>%
  summarize(mean=mean(mglassd),
            sd=sd(mglassd),
            num=length(mglassd))%>%
  mutate(se=sd/sqrt(num))%>%
  filter(trt_type2=="N"|trt_type2=="Mult. Nuts."|trt_type2=="Irrigation"|trt_type2=="CO2"|trt_type2=="P"|trt_type2=="Temperature"|trt_type2=="Irr + Temp")

glassD_alla<-five_all%>%
  ungroup()%>%
  group_by(response_var)%>%
  summarize(mean=mean(mglassd),
            sd=sd(mglassd),
            num=length(mglassd))%>%
  mutate(se=sd/sqrt(num))%>%
  mutate(trt_type2="All Trts")

glassD_alldata<-glassD_alla%>%
  bind_rows(glassD_trta)

a<-ggplot(data=glassD_alldata, aes(x=trt_type2, y=mean, fill=trt_type2))+
  geom_bar(position=position_dodge(), stat="identity")+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),position= position_dodge(0.9), width=0.2)+
  ylab("Glass's D")+
  xlab("")+
  scale_x_discrete(limits=c("All Trts","CO2","Irrigation","Temperature","Irr + Temp","N","P","Mult. Nuts."), labels=c("All Trts", "CO2","Irrigation", "Temp","Irr + Temp","Nitrogen","Phosphorus","Mult Nuts"))+
  scale_fill_manual(values=c("black","green3",'purple','blue','darkorange','orange','gold3','red'))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.5))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")+
  geom_vline(xintercept = 1.5, size = 1)+
  facet_wrap(~response_var, ncol=1)+
  ggtitle("All Data")

##B
glassD_trtb<-five_sigonly%>%
  group_by(response_var, trt_type2)%>%
  summarize(mean=mean(mglassd),
            sd=sd(mglassd),
            num=length(mglassd))%>%
  mutate(se=sd/sqrt(num))%>%
  filter(trt_type2=="N"|trt_type2=="Mult. Nuts."|trt_type2=="Irrigation"|trt_type2=="CO2"|trt_type2=="P"|trt_type2=="Temperature"|trt_type2=="Irr + Temp")

glassD_allb<-five_sigonly%>%
  ungroup()%>%
  group_by(response_var)%>%
  summarize(mean=mean(mglassd),
            sd=sd(mglassd),
            num=length(mglassd))%>%
  mutate(se=sd/sqrt(num))%>%
  mutate(trt_type2="All Trts")

glassD_alldatb<-glassD_allb%>%
  bind_rows(glassD_trtb)

b<-ggplot(data=glassD_alldatb, aes(x=trt_type2, y=mean, fill=trt_type2))+
  geom_bar(position=position_dodge(), stat="identity")+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),position= position_dodge(0.9), width=0.2)+
  ylab("Glass's D")+
  xlab("")+
  scale_x_discrete(limits=c("All Trts","CO2","Irrigation","Temperature","Irr + Temp","N","P","Mult. Nuts."), labels=c("All Trts", "CO2","Irrigation", "Temp","Irr + Temp","Nitrogen","Phosphorus","Mult Nuts"))+
  scale_fill_manual(values=c("black","green3",'purple','blue','darkorange','orange','gold3','red'))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.5))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")+
  geom_vline(xintercept = 1.5, size = 1)+
  facet_wrap(~response_var, ncol=1)+
  ggtitle("Sig Only")

grid.arrange(a, b, ncol=2)

###doing with last year only
last_all <- GlassD %>%
  group_by(site_project_comm, treatment)%>%
  mutate(maxyear=max(treatment_year2))%>%
  filter(treatment_year2==maxyear)%>%
  group_by(site_project_comm, treatment, response_var, trt_type2, use) %>%
  summarise(mglassd=mean(glassd, na.rm=T))

last_sigonly <- GlassD %>%
  filter(sig_diff_cntrl_trt=="yes")%>%
  group_by(site_project_comm, treatment)%>%
  mutate(maxyear=max(treatment_year2))%>%
  filter(treatment_year2==maxyear)%>%
  group_by(site_project_comm, treatment, response_var, trt_type2, use) %>%
  summarise(mglassd=mean(glassd, na.rm=T))

###graphing this last year only A - all data, B - sig only

##A
glassD_trta<-last_all%>%
  group_by(response_var, trt_type2)%>%
  summarize(mean=mean(mglassd),
            sd=sd(mglassd),
            num=length(mglassd))%>%
  mutate(se=sd/sqrt(num))%>%
  filter(trt_type2=="N"|trt_type2=="Mult. Nuts."|trt_type2=="Irrigation"|trt_type2=="CO2"|trt_type2=="P"|trt_type2=="Temperature"|trt_type2=="Irr + Temp")

glassD_alla<-last_all%>%
  ungroup()%>%
  group_by(response_var)%>%
  summarize(mean=mean(mglassd),
            sd=sd(mglassd),
            num=length(mglassd))%>%
  mutate(se=sd/sqrt(num))%>%
  mutate(trt_type2="All Trts")

glassD_alldata<-glassD_alla%>%
  bind_rows(glassD_trta)

a<-ggplot(data=glassD_alldata, aes(x=trt_type2, y=mean, fill=trt_type2))+
  geom_bar(position=position_dodge(), stat="identity")+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),position= position_dodge(0.9), width=0.2)+
  ylab("Glass's D")+
  xlab("")+
  scale_x_discrete(limits=c("All Trts","CO2","Irrigation","Temperature","Irr + Temp","N","P","Mult. Nuts."), labels=c("All Trts", "CO2","Irrigation", "Temp","Irr + Temp","Nitrogen","Phosphorus","Mult Nuts"))+
  scale_fill_manual(values=c("black","green3",'purple','blue','darkorange','orange','gold3','red'))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.5))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")+
  geom_vline(xintercept = 1.5, size = 1)+
  facet_wrap(~response_var, ncol=1)+
  ggtitle("All Data")

##B
glassD_trtb<-last_sigonly%>%
  group_by(response_var, trt_type2)%>%
  summarize(mean=mean(mglassd),
            sd=sd(mglassd),
            num=length(mglassd))%>%
  mutate(se=sd/sqrt(num))%>%
  filter(trt_type2=="N"|trt_type2=="Mult. Nuts."|trt_type2=="Irrigation"|trt_type2=="CO2"|trt_type2=="P"|trt_type2=="Temperature"|trt_type2=="Irr + Temp")

glassD_allb<-last_sigonly%>%
  ungroup()%>%
  group_by(response_var)%>%
  summarize(mean=mean(mglassd),
            sd=sd(mglassd),
            num=length(mglassd))%>%
  mutate(se=sd/sqrt(num))%>%
  mutate(trt_type2="All Trts")

glassD_alldatb<-glassD_allb%>%
  bind_rows(glassD_trtb)

b<-ggplot(data=glassD_alldatb, aes(x=trt_type2, y=mean, fill=trt_type2))+
  geom_bar(position=position_dodge(), stat="identity")+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),position= position_dodge(0.9), width=0.2)+
  ylab("Glass's D")+
  xlab("")+
  scale_x_discrete(limits=c("All Trts","CO2","Irrigation","Temperature","Irr + Temp","N","P","Mult. Nuts."), labels=c("All Trts", "CO2","Irrigation", "Temp","Irr + Temp","Nitrogen","Phosphorus","Mult Nuts"))+
  scale_fill_manual(values=c("black","green3",'purple','blue','darkorange','orange','gold3','red'))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.5))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")+
  geom_vline(xintercept = 1.5, size = 1)+
  facet_wrap(~response_var, ncol=1)+
  ggtitle("Sig Only")

grid.arrange(a, b, ncol=2)


##All stats. I will do this once we decide which we like best
t.test(change_glass_d_mean$abs_richness_glass, mu=0)#over all sig.
t.test(change_glass_d_mean$abs_evenness_glass, mu=0)#overall sig.
t.test(change_glass_d_mean$rank_glass, mu=0)#not sig
t.test(change_glass_d_mean$gains_glass, mu=0)#not sig
t.test(change_glass_d_mean$losses_glass, mu=0)#sig.

#water
irr<-subset(change_glass_d_mean, trt_type2=="Irrigation")
t.test(irr$abs_richness_glass, mu=0)#not sig
t.test(irr$abs_evenness_glass, mu=0)#not sig
t.test(irr$rank_glass, mu=0)#not sig
t.test(irr$gains_glass, mu=0)#not sig
t.test(irr$losses_glass, mu=0)#not sig

#N
N<-subset(change_glass_d_mean, trt_type2=="N")
t.test(N$abs_richness_glass, mu=0)#sig
t.test(N$abs_evenness_glass, mu=0)#sig
t.test(N$rank_glass, mu=0)# sig
t.test(N$gains_glass, mu=0)# sig
t.test(N$losses_glass, mu=0)# sig

#P
P<-subset(change_glass_d_mean, trt_type2=="P")
t.test(P$abs_richness_glass, mu=0)#not sig
t.test(P$abs_evenness_glass, mu=0)#not sig
t.test(P$rank_glass, mu=0)# not sig
t.test(P$gains_glass, mu=0)# not sig
t.test(P$losses_glass, mu=0)#not sig

#CO2
CO2<-subset(change_glass_d_mean, trt_type2=="CO2")
t.test(CO2$abs_richness_glass, mu=0)#not sig
t.test(CO2$abs_evenness_glass, mu=0)#not sig
t.test(CO2$rank_glass, mu=0)# not sig
t.test(CO2$gains_glass, mu=0)# not sig
t.test(CO2$losses_glass, mu=0)#not sig

#mult nuts
multnuts<-subset(change_glass_d_mean, trt_type2=="Mult. Nuts.")
t.test(multnuts$abs_richness_glass, mu=0)# sig
t.test(multnuts$abs_evenness_glass, mu=0)# sig
t.test(multnuts$rank_glass, mu=0)# not sig
t.test(multnuts$gains_glass, mu=0)# not sig
t.test(multnuts$losses_glass, mu=0)# sig

#temperature
temp<-subset(change_glass_d_mean, trt_type2=="Temperature")
t.test(temp$abs_richness_glass, mu=0)# not sig
t.test(temp$abs_evenness_glass, mu=0)# not sig
t.test(temp$rank_glass, mu=0)# not sig
t.test(temp$gains_glass, mu=0)# not sig
t.test(temp$losses_glass, mu=0)# not sig

#water and temperature
irrtemp<-subset(change_glass_d_mean, trt_type2=="Irr + Temp")
t.test(irrtemp$abs_richness_glass, mu=0)# not sig
t.test(irrtemp$abs_evenness_glass, mu=0)# not sig
t.test(irrtemp$rank_glass, mu=0)# not sig
t.test(irrtemp$gains_glass, mu=0)# not sig
t.test(irrtemp$losses_glass, mu=0)# not sig
