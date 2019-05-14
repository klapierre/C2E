library(tidyverse)
setwd("C:\\Users\\wilco\\Dropbox\\C2E\\Products\\CommunityChange\\March2018 WG\\")

#meghan's computer
setwd("C:\\Users\\megha\\Dropbox\\C2E\\Products\\CommunityChange\\March2018 WG")
setwd("~/Dropbox/C2E/Products/CommunityChange/March2018 WG")

dat <- read.csv("MetricsTrts_July2018.csv")%>%
  select(-X)

sig_com<-read.csv('../Summer2018_Results/gam_com_sig_change.csv')%>%
  mutate(site_project_comm = site_proj_comm)

plotid<-read.csv("SpeciesRelativeAbundance_Oct2017.csv")%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))

plotid2<-plotid%>%
  select(site_project_comm, plot_id, treatment)%>%
  unique()

trt<-read.csv("ExperimentInformation_Nov2017.csv")%>%
  group_by(site_code, project_name, community_type) %>%
  mutate(exp_length = max(treatment_year)) %>%
  ungroup() %>%
  select(site_code, project_name, community_type, treatment,plot_mani, exp_length)%>%
  unique()%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))

sigdat<-dat%>%
  left_join(sig_com)%>%
  filter(keep == "yes"|plot_mani ==0)

max <- sigdat%>%
  group_by(site_project_comm, treatment, plot_id)%>%
#  filter(exp_length>7)%>%
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

# change_value<-dat%>%
#   group_by(site_project_comm, treatment, plot_id)%>%
# #  filter(exp_length>7)%>%
#   mutate(S = max(abs(richness_change)),
#          E = max(abs(evenness_change)),
#          R = max(rank_change),
#          G = max(gains),
#          L = max(losses))%>%
#   select(site_project_comm, treatment, plot_id, S, E, R, G, L)%>%
#   unique()%>%
#   gather(change_metric, change_value, S:L)%>%
#   mutate(merge = ifelse(change_metric == "S","Richness",ifelse(change_metric=="E","Evenness",ifelse(change_metric=="R","Rank",ifelse(change_metric=="G","Gains", "Losses")))))


max_toplot<-max%>%
  left_join(plotid2)%>%
  left_join(trt)%>%
  mutate(trt = ifelse(plot_mani == 0, "Control", "Treatment"))%>%
#  left_join(change_value) %>%
  group_by(site_project_comm, treatment, max_metric, trt, exp_length) %>%
  summarise_at(vars(max_value), funs(mean))


longterm_plot <- ggplot(data = subset(max_toplot,exp_length>7), aes(x = max_metric, y = max_value))+
  geom_jitter()+
  geom_boxplot(alpha=.1) +
  facet_wrap(~trt) +
  theme_bw() +
  xlab("Change Metric") +
  ylab("Treatment Year") +
  theme(axis.text=element_text(size=12, color="black"), strip.text.x=element_text(size=12)) +
  scale_x_discrete(limits=c("Smax","Emax","Rmax","Gmax","Lmax"),labels=c("Rich","Even","Rank","Gain","Loss"))

pdf(paste0("..\\Summer2018_Results\\Figures\\maxvalue change boxplots_exp gt 7_", Sys.Date(), ".pdf"), height=4, width=7, useDingbats=F)
print(longterm_plot)
dev.off()

allexp_plot <- ggplot(data = max_toplot, aes(x = max_metric, y = max_value))+
  geom_jitter()+
  geom_boxplot(alpha=0.1) +
  facet_wrap(~trt) +
  #theme_few() +
  xlab("Change Metric") +
  ylab("Treatment Year") +
  coord_flip() +
  theme_bw() +
  theme(axis.text=element_text(size=12, color="black"), strip.text.x=element_text(size=12)) +
  scale_x_discrete(limits=c("Smax","Emax","Rmax","Gmax","Lmax"),labels=c("Rich","Even","Rank","Gain","Loss"))

pdf(paste0("..\\Summer2018_Results\\Figures\\maxvalue change boxplots_all exp_", Sys.Date(), ".pdf"), height=4, width=7, useDingbats=F)
print(allexp_plot)
dev.off()

write.csv(max_toplot, 
          paste0("../Summer2018_Results/year of max SERGL change_only significant comm change included_",Sys.Date(),".csv"),
                 row.names=F)

