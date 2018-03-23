library(tidyverse)

dat <- read.csv("~/Dropbox/C2E/Products/CommunityChange/March2018 WG/CORRE_RACS_Subset_Bayes.csv")

plotid<-read.csv("~/Dropbox/converge_diverge/datasets/LongForm/SpeciesRelativeAbundance_Oct2017.csv")%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))

plotid2<-plotid%>%
  select(site_project_comm, plot_id, treatment)%>%
  unique()

trt<-read.csv("~/Dropbox/converge_diverge/datasets/LongForm/ExperimentInformation_Nov2017.csv")%>%
  select(site_code, project_name, community_type, treatment,plot_mani)%>%
  unique()%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))


max <- dat%>%
  group_by(site_project_comm, treatment, plot_id)%>%
  filter(treatment_year<11)%>%
  mutate(S = max(richness_change),
         E = max(evenness_change),
         R = max(rank_change),
         G = max(gains),
         L = max(losses),
         Smax = ifelse(S == richness_change, treatment_year, 0),
         Emax = ifelse(E == evenness_change, treatment_year, 0),
         Rmax = ifelse(R == rank_change, treatment_year, 0),
         Gmax = ifelse(G == gains, treatment_year, 0),
         Lmax = ifelse(L == losses, treatment_year, 0))%>%
  gather(max_metric, max_value, Smax:Lmax)%>%
  select(site_project_comm, treatment, plot_id, max_metric, max_value)%>%
  mutate(merge = ifelse(max_metric == "Smax","Richness",ifelse(max_metric=="Emax","Evenness",ifelse(max_metric=="Rmax","Rank",ifelse(max_metric=="Gmax","Gains", "Losses")))))%>%
  filter(max_value!=0)

change_value<-dat%>%
  group_by(site_project_comm, treatment, plot_id)%>%
  filter(treatment_year<11)%>%
  mutate(S = max(richness_change),
         E = max(evenness_change),
         R = max(rank_change),
         G = max(gains),
         L = max(losses))%>%
  select(site_project_comm, treatment, plot_id, S, E, R, G, L)%>%
  unique()%>%
  gather(change_metric, change_value, S:L)%>%
  mutate(merge = ifelse(change_metric == "S","Richness",ifelse(change_metric=="E","Evenness",ifelse(change_metric=="R","Rank",ifelse(change_metric=="G","Gains", "Losses")))))


max_toplot<-max%>%
  left_join(plotid2)%>%
  left_join(trt)%>%
  mutate(trt = ifelse(plot_mani == 0, "C", "T"))%>%
  left_join(change_value)

ggplot(data = max_toplot, aes(x = max_metric, y = max_value, color = change_value))+
  geom_point()+
  geom_jitter()+
  
  facet_wrap(~trt)
  