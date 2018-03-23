### NMDS for select experiments

### Set up workspace

library(tidyverse)
library(vegan)
library(ggthemes)
corredat1<-read.csv("~/Dropbox/converge_diverge/datasets/LongForm/SpeciesRelativeAbundance_Oct2017.csv")%>%
  select(-X)%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))%>%
  filter(site_project_comm!="GVN_FACE_0", site_project_comm!="AZI_NitPhos_0", site_project_comm!="JRN_study278_0", site_project_comm!="KNZ_GFP_4F", site_project_comm!="Saskatchewan_CCD_0")

##several studies only have two measurments of a plot. I am dropping those plots
azi<-read.csv("~/Dropbox/converge_diverge/datasets/LongForm/SpeciesRelativeAbundance_Oct2017.csv")%>%
  select(-X)%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))%>%
  filter(site_code=="AZI")%>%
  filter(plot_id!=11&plot_id!=15&plot_id!=35&plot_id!=37)

jrn<-read.csv("~/Dropbox/converge_diverge/datasets/LongForm/SpeciesRelativeAbundance_Oct2017.csv")%>%
  select(-X)%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))%>%
  filter(site_project_comm=="JRN_study278_0")%>%
  filter(plot_id!=211&plot_id!=210)

knz<-read.csv("~/Dropbox/converge_diverge/datasets/LongForm/SpeciesRelativeAbundance_Oct2017.csv")%>%
  select(-X)%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))%>%
  filter(site_project_comm=="KNZ_GFP_4F")%>%
  filter(plot_id!="7_1_1"&plot_id!="7_2_1")

sak<-read.csv("~/Dropbox/converge_diverge/datasets/LongForm/SpeciesRelativeAbundance_Oct2017.csv")%>%
  select(-X)%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))%>%
  filter(site_project_comm=="Saskatchewan_CCD_0")%>%
  filter(plot_id!=2)

corredat_raw <-rbind(corredat1, azi, jrn, knz, sak)

treatment_info<-read.csv("~/Dropbox/converge_diverge/datasets/LongForm/ExperimentInformation_Nov2017.csv")%>%
  dplyr::select(site_code, project_name, community_type, treatment,plot_mani)%>%
  unique()%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))

corredat <- corredat_raw %>%
  filter(site_project_comm != "GVN_FACE_0") %>%
  left_join(treatment_info, by=c( "site_code","project_name","community_type", "treatment","site_project_comm"))%>%
  filter(site_project_comm == "SERC_TMECE_SP"|site_project_comm=="KNZ_GFP_K1B"|site_project_comm=="Alberta_CCD_0")




###running through all NMDS
#this might not work this was earlier code
site_project_comm_vec <- unique(corredat$site_project_comm)

 nmds=data.frame()

 for (i in 1:length(site_project_comm_vec)){
   subset<-corredat%>%
     filter(site_project_comm==site_project_comm_vec[i])%>%
    spread(genus_species, relcov, fill=0)

  info<-subset[,1:10]

  mds<-metaMDS(subset[,11:ncol(subset)],autotransform=FALSE, shrink=FALSE)

  scores <- data.frame(scores(mds, display="sites"))
  scores2<- cbind(info, scores)

  nmds<-rbind(nmds, scores2)
  }

nmds_mean<-nmds%>%
  group_by(site_project_comm, calendar_year, treatment)%>%
  summarize(n=length(NMDS1),
            m1=mean(NMDS1),
            m2=mean(NMDS2),
            sd1=sd(NMDS1),
            sd2=sd(NMDS2),
            se1=sd1/sqrt(n),
            se2=sd2/sqrt(n))
theme_set(theme_bw(12))


for (i in 1:length(site_project_comm_vec)){

  tograph<-subset(nmds_mean, site_project_comm==site_project_comm_vec[i])

  print(ggplot(data=tograph, aes(x=m1, y=m2, color=treatment))+
    geom_point(size=5)+
    geom_errorbar(aes(ymin=m2-se2, ymax=m2+se2), color="black")+
    geom_errorbarh(aes(xmin=m1-se1, xmax=m1+se1),color="black")+
    theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank())+
    xlab("NMDS1")+
    ylab("NMDS2")+
    ggtitle(site_project_comm_vec[i])+
    facet_wrap(~calendar_year))

  }
