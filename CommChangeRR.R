library(tidyverse)
library(gridExtra)
library(grid)
library(gtable)
library(codyn)
library(vegan)
library(Kendall)

theme_set(theme_bw(12))

#read in the data
dat1<-read.csv("C:\\Users\\megha\\Dropbox\\converge_diverge\\datasets\\LongForm\\CORRE_RAC_Metrics_Oct2017_compareyears.csv")%>%
  select(-X)

dat1<-read.csv("~/Dropbox/converge_diverge/datasets/LongForm/CORRE_RAC_Metrics_Oct2017_compareyears.csv")%>%
  select(-X)

trt<-read.csv("C:\\Users\\megha\\Dropbox\\converge_diverge\\datasets\\LongForm\\ExperimentInformation_May2017.csv")%>%
  select(site_code, project_name, community_type, treatment,plot_mani)%>%
  unique()%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))

dat2<-merge(dat, trt, by=c("site_project_comm","treatment"))

##overall what is the relationsip between the metrics
#graphing this
dat1<-dat1%>%
  na.omit

panel.pearson <- function(x, y, ...) {
  horizontal <- (par("usr")[1] + par("usr")[2]) / 2; 
  vertical <- (par("usr")[3] + par("usr")[4]) / 2; 
  text(horizontal, vertical, format(cor(x,y), digits=3, cex=10)) 
}

pairs(dat1[,c(5:11)], upper.panel = panel.pearson)


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

###doing the horizonal line comparison
dat<-read.csv("~/Dropbox/converge_diverge/datasets/LongForm/CORRE_ContTreat_Compare_OCT2017.csv")%>%
  select(-X)%>%
  mutate(id=paste(site_project_comm, treatment, sep="::"),
         plotmani=as.factor(plot_mani))

trt<-read.csv("~/Dropbox/converge_diverge/datasets/LongForm/ExperimentInformation_May2017.csv")%>%
  select(-calendar_year, -treatment_year, -X, -plot_mani)%>%
  unique()%>%
  mutate(nut=ifelse(n!=0|p!=0|k!=0, 1, 0),
         co2=ifelse(CO2!=0, 1, 0), 
         wat=ifelse(precip!=0, 1, 0),
         other=ifelse(temp!=0|mow_clip!=0|burn!=0|herb_removal!=0|management!=0|plant_mani!=0,1,0),
         trt=nut+co2+wat+other,
         gtreat=ifelse(trt==1&nut==1, "nutrient",ifelse(trt==1&co2==1, "carbon", ifelse(trt==1&wat==1,"water", ifelse(trt==1&other==1,"non resource", ifelse(trt>1,"multiple", "control"))))))%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))%>%
  select(site_project_comm, treatment, gtreat)

dat3<-merge(dat, trt, by=c("site_project_comm","treatment"))

%>%
  filter(pulse!=1, successional!=1)

S<-
  ggplot(data=dat, aes(x=treatment_year, y=PCSdiff, color=plotmani))+
  geom_line(aes(group=id))+
  geom_smooth(method="lm", se=F, size=2, aes(group=plotmani, color=plotmani))+
  ggtitle("Richness")+
  geom_hline(yintercept=0, size=1)

even<-
  ggplot(data=dat, aes(x=treatment_year, y=PCEvendiff, color=plotmani))+
  geom_line(aes(group=id))+
  geom_smooth(method="lm", se=F,size=2, aes(group=plotmani, color=plotmani))+
  ggtitle("Evenness")+
  geom_hline(yintercept=0, size=1)

loss<-
  ggplot(data=dat, aes(x=treatment_year, y=spdiffc, color=plotmani))+
  geom_line(aes(group=id))+
  geom_smooth(method="lm", se=F, size=2, aes(group=plotmani, color=plotmani))+
  ggtitle("Sp Difference")+
  geom_hline(yintercept=0, size=1)

mrsc<-
  ggplot(data=dat, aes(x=treatment_year, y=MRSc_diff, color=plotmani))+
  geom_line(aes(group=id))+
  geom_smooth(method="lm", se=F, size=2, aes(group=plotmani, color=plotmani))+
  ggtitle("Reordering")+
  geom_hline(yintercept=0, size=1)

grid.arrange(S, even, loss, mrsc)

#figure for kevin
S<-
  ggplot(data=dat, aes(x=treatment_year, y=PCSdiff))+
  geom_line(aes(group=id))+
  geom_smooth(method="lm", se=F, size=2, color="red")+
  ggtitle("Richness")+
  geom_hline(yintercept=0, size=1)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ylab("% Change in Richness")+
  xlab("Treatment Year")

even<-
  ggplot(data=dat, aes(x=treatment_year, y=PCEvendiff))+
  geom_line(aes(group=id))+
  geom_smooth(method="lm", se=F,size=2, color="red")+
  ggtitle("Evenness")+
  geom_hline(yintercept=0, size=1)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ylab("% Change in Evenness")+
  xlab("Treatment Year")

loss<-
  ggplot(data=dat, aes(x=treatment_year, y=spdiffc))+
  geom_line(aes(group=id))+
  geom_smooth(method="lm", se=F, size=2, color="red")+
  ggtitle("Species Difference")+
  geom_hline(yintercept=0, size=1)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ylab("Proportion of Species Different")+
  xlab("Treatment Year")

mrsc<-
  ggplot(data=dat, aes(x=treatment_year, y=MRSc_diff))+
  geom_line(aes(group=id))+
  geom_smooth(method="lm", se=F, size=2, color="red")+
  ggtitle("Rank Differences")+
  geom_hline(yintercept=0, size=1)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ 
  ylab("Mean Rank Shift")+
  xlab("Treatment Year")

grid.arrange(S, even, loss, mrsc)

