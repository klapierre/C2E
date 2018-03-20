library(tidyverse)
library(gridExtra)
library(grid)
library(gtable)
library(codyn)
library(vegan)
library(Kendall)


#read in data and subset only experiemnts that are 10 or more years

dat<-read.csv("C:\\Users\\megha\\Dropbox\\converge_diverge\\datasets\\LongForm\\CORRE_RAC_Metrics_Oct2017_compareyears.csv")%>%
  select(-X)

dat<-read.csv("~/Dropbox/converge_diverge/datasets/LongForm/CORRE_RAC_Metrics_Oct2017_allyears_2.csv")%>%
  select(-X)

trt<-read.csv("~/Dropbox/converge_diverge/datasets/LongForm/ExperimentInformation_May2017.csv")%>%
  select(site_code, project_name, community_type, treatment,plot_mani)%>%
  unique()%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))


longset<-dat%>%
  select(site_project_comm, treatment_year)%>%
  group_by(site_project_comm)%>%
  summarize(len=max(treatment_year))%>%
  filter(len>9)%>%
  select(-len)

datsub<-merge(dat, longset, by="site_project_comm")%>%
  mutate(id=paste(site_project_comm, treatment, sep="::"))%>%
  na.omit #3 averages for evenness have NAs and need to omit

##LOOP THIS

list<-unique(datsub$id)

vp_output<-data.frame()

for (i in 1:length(list)){
  
  #subset to get a single treatment for a site_proj_comm
  subset<-datsub%>%
    filter(id==list[i])
  
  #do vp and pull out what we want
  vp<-varpart(subset$mean_change, ~even_diff, ~gain, ~loss, ~MRSc, data=subset)
  
  standdev<-sd(subset$mean_change)
  
  adjR.temp <- data.frame(id=unique(subset$id),
                          treatment=unique(subset$treatment),
                          metric=c("even","gain","loss","MRSc"),
                          adj.r2=vp[["part"]][["indfract"]][c(1:4),"Adj.R.square"],
                          resid=vp[["part"]][["indfract"]][16,"Adj.R.square"],
                          sd_meanchange=standdev)
  
  vp_output<-rbind(vp_output, adjR.temp)
}

vp_output2<-vp_output%>%
  separate(id, into=c("site_project_comm", "trt"), sep="::")%>%
  mutate(treat=as.factor(treatment))%>%
  na.omit

vp_output3<-merge(vp_output2, trt, by=c("site_project_comm","treatment"))

theme_set(theme_bw(12))
ggplot(data=vp_output2, aes(x=metric, y=adj.r2, group=treatment))+
        geom_bar(stat="identity", position=position_dodge(),aes(fill=treat))+
        theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank())+
        xlab("Metric")+
        ylab("Adj R2")+
        #geom_text(label=resid)+
        #geom_text(label=sttdev)+
        facet_wrap(~site_project_comm)+
  theme(legend.position="none")

#didn't work for KNZ_BGP; ASGA_clonal, DCGS_gap, ORNL_Face
probelmatic<-datsub%>%
  filter(site_project_comm=="ASGA_clonal_0"|site_project_comm=="KNZ_BGP_0"|site_project_comm=="dcgs_gap_0"|site_project_comm=="ORNL_FACE_0")

test<-probelmatic%>%
  filter(id=="ASGA_clonal_0::non-clonal_UN")

vp<-varpart(test$mean_change, ~even, ~gain, ~loss, ~MRSc, data=test)





# looking at this
 plot(vp)
# plot(test_c$calendar_year,test_c$mean_change)
# points(test_np$calendar_year, test_np$mean_change, col="red", pch=19)
