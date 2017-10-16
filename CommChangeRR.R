library(tidyverse)
library(gridExtra)
library(grid)
library(gtable)
library(codyn)
library(vegan)
library(Kendall)


dat<-read.csv("C:\\Users\\megha\\Dropbox\\converge_diverge\\datasets\\LongForm\\CORRE_RAC_Metrics_Oct2017_compareyears.csv")%>%
  select(-X)

trt<-read.csv("C:\\Users\\megha\\Dropbox\\converge_diverge\\datasets\\LongForm\\ExperimentInformation_May2017.csv")%>%
  select(site_code, project_name, community_type, treatment,plot_mani)%>%
  unique()%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))

dat2<-merge(dat, trt, by=c("site_project_comm","treatment"))

control<-dat2%>%
  filter(plot_mani==0)%>%
  mutate(C_S=S, C_even=even, C_gain=gain, C_loss=loss, C_MRSc=MRSc)%>%
  select(site_project_comm, treatment_year, C_S, C_even, C_gain, C_loss, C_MRSc)

# logRR<-merge(dat2, control, by=c("site_project_comm","treatment_year"))%>%
#   filter(plot_mani!=0)%>%
#   mutate(lrS=log(S/C_S), lrE=log(even/C_even), lrG=log(gain/C_gain), lrL=log(loss/C_loss), lrR=log(MRSc/C_MRSc))



RR<-merge(dat2, control, by=c("site_project_comm","calendar_year"))%>%
  filter(plot_mani!=0)%>%
  mutate(rS=((S-C_S)/C_S), rE=((even-C_even)/C_even), rG=((gain-C_gain)/C_gain), rL=((loss-C_loss)/C_loss), rR=((MRSc-C_MRSc)/C_MRSc))%>%
  select(site_project_comm, treatment_year, treatment,rS, rE, rG, rL, rR)%>%
  mutate(id=paste(site_project_comm, treatment, sep="::"))

S<-ggplot(data=RR, aes(x=treatment_year, y=rS))+
  geom_line(aes(group=id))+
  geom_smooth(method="lm", se=F, color="red", size=3)+
  ggtitle("Richness")+
  geom_hline(yintercept=0, size=1)
even<-ggplot(data=RR, aes(x=treatment_year, y=rE))+
  geom_line(aes(group=id))+
  geom_smooth(method="lm", se=F, color="red", size=3)+
  ggtitle("Evenness")+
  geom_hline(yintercept=0, size=1)
gain<-ggplot(data=RR, aes(x=treatment_year, y=rG))+
  geom_line(aes(group=id))+
  geom_smooth(method="lm", se=F, color="red", size=3)+
  ggtitle("Gains")+
  geom_hline(yintercept=0, size=1)
loss<-ggplot(data=RR, aes(x=treatment_year, y=rL))+
  geom_line(aes(group=id))+
  geom_smooth(method="lm", se=F, color="red", size=3)+
  ggtitle("Losses")+
  geom_hline(yintercept=0, size=1)
mrsc<-ggplot(data=RR, aes(x=treatment_year, y=rR))+
  geom_line(aes(group=id))+
  geom_smooth(method="lm", se=F, color="red", size=3)+
  ggtitle("Reordering")+
  geom_hline(yintercept=0, size=1)
  
grid.arrange(S, even, gain, loss, mrsc)