library(tidyverse)
library(gridExtra)
library(grid)
library(gtable)
library(codyn)
library(vegan)
library(Kendall)
library(relaimpo)
library(MASS)
library(lme4)

theme_set(theme_bw(12))
#path work
C:\\Users\\megha\\Dropbox\\converge_diverge\\datasets\\LongForm\\
#path home
~/Dropbox/converge_diverge/datasets/LongForm/

#read in the data
dat<-read.csv("~/Dropbox/converge_diverge/datasets/LongForm/CORRE_RAC_Metrics_Feb2018_allReplicates.csv")%>%
  select(-X)

dat_mult<-read.csv("~/Dropbox/converge_diverge/datasets/LongForm/CORRE_Mult_Metrics_Feb2018.csv")%>%
  select(-X)

plotid<-read.csv("~/Dropbox/converge_diverge/datasets/LongForm/SpeciesRelativeAbundance_Oct2017.csv")%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))

plotid2<-plotid%>%
  select(site_project_comm, plot_id, treatment)%>%
  unique()

trt<-read.csv("~/Dropbox/converge_diverge/datasets/LongForm/ExperimentInformation_Nov2017.csv")%>%
  select(site_code, project_name, community_type, treatment,plot_mani)%>%
  unique()%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))

dat1<-dat%>%
  left_join(plotid)%>%
  left_join(trt)

##overall what is the relationsip between the metrics
#graphing this
dat2<-dat1%>%
  na.omit

panel.pearson <- function(x, y, ...) {
  horizontal <- (par("usr")[1] + par("usr")[2]) / 2; 
  vertical <- (par("usr")[3] + par("usr")[4]) / 2; 
  text(horizontal, vertical, format(cor(x,y), digits=3, cex=10)) 
}

pairs(dat2[,c(3:7)], upper.panel = panel.pearson)

ave<-dat2%>%
  group_by(calendar_year_pair, site_project_comm, treatment, plot_mani)%>%
  summarize(rich=mean(richness_change),
            even=mean(evenness_change),
            rank=mean(rank_change),
            gain=mean(gains),
            loss=mean(losses))%>%
  mutate(trt=ifelse(plot_mani==0, "control","treatment"))

merged_data<-ave%>%
  left_join(dat_mult)


###doing multiple regression
stepAIC(lm(composition_change~rich+even+rank+gain+loss, data=merged_data))
summary(m3<-lm(composition_change~rich+even+rank+gain+loss, data=merged_data))
calc.relimp(m3, type="lmg", rela=F)

#mixed model - not sure of the code here.
m2 <- lmer(composition_change ~ rich+even+rank+gain+loss+
             (composition_change | site_project_comm),
           data = merged_data)
summary(m2)
m2
###datasets 10 yrs of longer
datalen<-dat_mult%>%
  ungroup()%>%
  select(calendar_year_pair, site_project_comm)%>%
  unique()%>%
  group_by(site_project_comm)%>%
  summarize(len=length(calendar_year_pair))%>%
  mutate(keep=ifelse(len>7, 1, 0))

long<-merged_data%>%
  left_join(datalen)%>%
  filter(keep==1)

control<-long%>%
  filter(plot_mani==0)%>%
  ungroup()%>%
  mutate(crich=rich, ceven=even, crank=rank, cgain=gain, closs=loss)%>%
  select(site_project_comm, calendar_year_pair, crich, ceven, crank, cgain, closs)

trt<-long%>%
  filter(plot_mani>0)

cont_trt<-trt%>%
  left_join(control)%>%
  mutate(rS=((rich-crich)/crich), rE=((even-ceven)/ceven), rG=((gain-cgain)/cgain), rL=((loss-closs)/closs), rR=((rank-crank)/crank),  trt=ifelse(plot_mani==0, "control","treatment"))%>%
  separate(calendar_year_pair, into=c("year", "year1"), sep="-")%>%
  filter(!is.infinite(rS), !is.infinite(rE), !is.infinite(rR), !is.infinite(rG), !is.infinite(rL),
         !is.nan(rS), !is.nan(rE), !is.nan(rR), !is.nan(rG), !is.nan(rL))

###graphing data add line for controls and line for treatments]

  ggplot(data=cont_trt, aes(x=year, y=rS))+
    geom_point(size=3)+
    geom_smooth(method="lm", se=F, color="red", size=3)+
    facet_wrap(~site_project_comm, scales="free")

  ggplot(data=cont_trt, aes(x=year, y=rE))+
    geom_point(size=3)+
    geom_smooth(method="lm", se=F, color="red", size=3)+
    facet_wrap(~site_project_comm, scales="free")
  
  ggplot(data=cont_trt, aes(x=year, y=rR))+
    geom_point(size=3)+
    geom_smooth(method="lm", se=F, color="red", size=3)+
    facet_wrap(~site_project_comm, scales="free")
  
  ggplot(data=cont_trt, aes(x=year, y=rL))+
    geom_point(size=3)+
    geom_smooth(method="lm", se=F, color="red", size=3)+
    facet_wrap(~site_project_comm, scales="free")
  
  ggplot(data=cont_trt, aes(x=year, y=rG))+
    geom_point(size=3)+
    geom_smooth(method="lm", se=F, color="red", size=3)+
    facet_wrap(~site_project_comm, scales="free")
  

even<-ggplot(data=dat1, aes(x=treatment_year, y=evenness_change))+
  geom_line(aes(group=id, color=trt), size=0.1)+
  geom_smooth(method="lm", se=F, size=2, aes(color=trt))+
  ggtitle("Evenness Change")
gain<-ggplot(data=dat1, aes(x=treatment_year, y=gains))+
  geom_line(aes(group=id, color=trt), size=0.1)+
  geom_smooth(method="lm", se=F, size=2, aes(color=trt))+
  ggtitle("Gains")
loss<-ggplot(data=dat1, aes(x=treatment_year, y=losses))+
  geom_line(aes(group=id, color=trt), size=0.1)+
  geom_smooth(method="lm", se=F, size=2, aes(color=trt))+
  ggtitle("Losses")
mrsc<-ggplot(data=dat1, aes(x=treatment_year, y=rank_change))+
  geom_line(aes(group=id, color=trt), size=0.1)+
  geom_smooth(method="lm", se=F,  size=3, aes(color=trt))+
  ggtitle("Reordering")

grid.arrange(S, even, gain, loss, mrsc)

###do some kind of ratio to controls
control<-dat1%>%
  filter(plot_mani==0)%>%
  mutate(C_S=S_diff, C_even=even_diff, C_gain=gain, C_loss=loss, C_MRSc=MRSc)%>%
  select(site_project_comm, treatment_year, C_S, C_even, C_gain, C_loss, C_MRSc)

logRR<-merge(dat1, control, by=c("site_project_comm","treatment_year"))%>%
  filter(plot_mani!=0)%>%
  mutate(lrS=log((S_diff+.01)/(C_S+.01)), lrE=log((even_diff+.01)/(C_even+.01)), lrG=log((gain+.01)/(C_gain+.01)), lrL=log((loss+.01)/(C_loss+.01)), lrR=log((MRSc+.01)/(C_MRSc+.01)))%>%
  na.omit%>%
  mutate(id=paste(site_project_comm, treatment, sep="::"))

S<-ggplot(data=logRR, aes(x=treatment_year, y=lrS))+
  geom_line(aes(group=id))+
  geom_smooth(method="lm", se=F, color="red", size=3)+
  ggtitle("Richness")+
  geom_hline(yintercept=0, size=1)
even<-ggplot(data=logRR, aes(x=treatment_year, y=lrE))+
  geom_line(aes(group=id))+
  geom_smooth(method="lm", se=F, color="red", size=3)+
  ggtitle("Evenness")+
  geom_hline(yintercept=0, size=1)
gain<-ggplot(data=logRR, aes(x=treatment_year, y=lrG))+
  geom_line(aes(group=id))+
  q
  ggtitle("Gains")+
  geom_hline(yintercept=0, size=1)
loss<-ggplot(data=logRR, aes(x=treatment_year, y=lrL))+
  geom_line(aes(group=id))+
  geom_smooth(method="lm", se=F, color="red", size=3)+
  ggtitle("Losses")+
  geom_hline(yintercept=0, size=1)
mrsc<-ggplot(data=logRR, aes(x=treatment_year, y=lrR))+
  geom_line(aes(group=id))+
  geom_smooth(method="lm", se=F, color="red", size=3)+
  ggtitle("Reordering")+
  geom_hline(yintercept=0, size=1)

grid.arrange(S, even, gain, loss, mrsc)

# 
# RR<-merge(dat2, control, by=c("site_project_comm","treatment_year"))%>%
#   filter(plot_mani!=0)%>%
#   mutate(rS=((S-C_S)/C_S), rE=((even-C_even)/C_even), rG=((gain-C_gain)/C_gain), rL=((loss-C_loss)/C_loss), rR=((MRSc-C_MRSc)/C_MRSc))%>%
#   select(site_project_comm, treatment_year, treatment,rS, rE, rG, rL, rR)%>%
#   mutate(id=paste(site_project_comm, treatment, sep="::"))
# 
# S<-ggplot(data=RR, aes(x=treatment_year, y=rS))+
#   geom_line(aes(group=id))+
#   geom_smooth(method="lm", se=F, color="red", size=3)+
#   ggtitle("Richness")+
#   geom_hline(yintercept=0, size=1)
# even<-ggplot(data=RR, aes(x=treatment_year, y=rE))+
#   geom_line(aes(group=id))+
#   geom_smooth(method="lm", se=F, color="red", size=3)+
#   ggtitle("Evenness")+
#   geom_hline(yintercept=0, size=1)
# gain<-ggplot(data=RR, aes(x=treatment_year, y=rG))+
#   geom_line(aes(group=id))+
#   geom_smooth(method="lm", se=F, color="red", size=3)+
#   ggtitle("Gains")+
#   geom_hline(yintercept=0, size=1)
# loss<-ggplot(data=RR, aes(x=treatment_year, y=rL))+
#   geom_line(aes(group=id))+
#   geom_smooth(method="lm", se=F, color="red", size=3)+
#   ggtitle("Losses")+
#   geom_hline(yintercept=0, size=1)
# mrsc<-ggplot(data=RR, aes(x=treatment_year, y=rR))+
#   geom_line(aes(group=id))+
#   geom_smooth(method="lm", se=F, color="red", size=3)+
#   ggtitle("Reordering")+
#   geom_hline(yintercept=0, size=1)

###doing the horizonal line comparison
dath<-read.csv("C:\\Users\\megha\\Dropbox\\converge_diverge\\datasets\\LongForm\\CORRE_ContTreat_Compare_OCT2017.csv")%>%
  select(-X)%>%
  mutate(id=paste(site_project_comm, treatment, sep="::"),
         plotmani=as.factor(plot_mani))


S<-  ggplot(data=dath, aes(x=treatment_year, y=PCSdiff))+
  geom_line(aes(group=id))+
  geom_smooth(method="lm", se=F, size=2, color="red")+
  ggtitle("Richness")+
  geom_hline(yintercept=0, size=1)
even<-  ggplot(data=dath, aes(x=treatment_year, y=PCEvendiff))+
  geom_line(aes(group=id))+
  geom_smooth(method="lm", se=F,size=2, color="red")+
  ggtitle("Evenness")+
  geom_hline(yintercept=0, size=1)
loss<-  ggplot(data=dath, aes(x=treatment_year, y=spdiffc))+
  geom_line(aes(group=id))+
  geom_smooth(method="lm", se=F, size=2,color="red")+
  ggtitle("Sp Difference")+
  geom_hline(yintercept=0, size=1)
mrsc<-  ggplot(data=dath, aes(x=treatment_year, y=MRSc_diff))+
  geom_line(aes(group=id))+
  geom_smooth(method="lm", se=F, size=2, color="red")+
  ggtitle("Reordering")+
  geom_hline(yintercept=0, size=1)

grid.arrange(S, even, loss, mrsc)

#trying to color by treatmetns
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


