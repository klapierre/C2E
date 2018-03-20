### NMDS for select experiments


### Set up workspace
setwd("C:\\Users\\wilco\\Dropbox\\")
library(tidyverse)
library(vegan)
library(ggthemes)
corredat1<-read.csv("C2E\\Products\\CommunityChange\\March2018 WG\\data\\SpeciesRelativeAbundance_Oct2017.csv")%>%
  select(-X)%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))%>%
  filter(site_project_comm!="GVN_FACE_0", site_project_comm!="AZI_NitPhos_0", site_project_comm!="JRN_study278_0", site_project_comm!="KNZ_GFP_4F", site_project_comm!="Saskatchewan_CCD_0")

##several studies only have two measurments of a plot. I am dropping those plots
azi<-read.csv("C2E\\Products\\CommunityChange\\March2018 WG\\data\\SpeciesRelativeAbundance_Oct2017.csv")%>%
  select(-X)%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))%>%
  filter(site_code=="AZI")%>%
  filter(plot_id!=11&plot_id!=15&plot_id!=35&plot_id!=37)

jrn<-read.csv("C2E\\Products\\CommunityChange\\March2018 WG\\data\\SpeciesRelativeAbundance_Oct2017.csv")%>%
  select(-X)%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))%>%
  filter(site_project_comm=="JRN_study278_0")%>%
  filter(plot_id!=211&plot_id!=210)

knz<-read.csv("C2E\\Products\\CommunityChange\\March2018 WG\\data\\SpeciesRelativeAbundance_Oct2017.csv")%>%
  select(-X)%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))%>%
  filter(site_project_comm=="KNZ_GFP_4F")%>%
  filter(plot_id!="7_1_1"&plot_id!="7_2_1")

sak<-read.csv("C2E\\Products\\CommunityChange\\March2018 WG\\data\\SpeciesRelativeAbundance_Oct2017.csv")%>%
  select(-X)%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))%>%
  filter(site_project_comm=="Saskatchewan_CCD_0")%>%
  filter(plot_id!=2)

corredat_raw <-rbind(corredat1, azi, jrn, knz, sak)

treatment_info<-read.csv("C2E\\Products\\CommunityChange\\March2018 WG\\data\\ExperimentInformation_Nov2017.csv")%>%
  dplyr::select(site_code, project_name, community_type, treatment,plot_mani)%>%
  unique()%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))

corredat <- corredat_raw %>%
  filter(site_project_comm != "GVN_FACE_0") %>%
  left_join(treatment_info, by=c( "site_code","project_name","community_type", "treatment","site_project_comm"))

site_project_comm_vec <- unique(corredat$site_project_comm)
