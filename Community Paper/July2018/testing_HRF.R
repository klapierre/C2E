#meghan's working directory
setwd("C:\\Users\\megha\\Dropbox\\")

library(tidyverse)
library(gridExtra)

change_metrics <- read.csv("C2E\\Products\\CommunityChange\\March2018 WG\\MetricsTrts_July2018.csv") %>%
  mutate(richness_change_abs = abs(richness_change),
         evenness_change_abs = abs(evenness_change))%>%
  select(-X, -richness_change, -evenness_change, -composition_change, -site_code, -project_name, -community_type, -trt_type, -use, -composition_change, -dispersion_change)

subset_studies<-read.csv("C:\\Users\\megha\\Dropbox\\C2E\\Products\\CommunityChange\\March2018 WG\\experiment_trt_subset.csv")

#not doing this.
# sig_com<-read.csv('C2E\\Products\\CommunityChange\\Summer2018_Results\\gam_com_sig_change.csv')%>%
#   mutate(site_project_comm = site_proj_comm)

metrics_sig<-read.csv("C2E\\Products\\CommunityChange\\Summer2018_Results\\gam_metrics_sig_change.csv")%>%
  mutate(site_project_comm=site_proj_comm)%>%
  select(-site_proj_comm)%>%
  filter(response_var!="richness_change_abs")%>%
  right_join(subset_studies)

### Control data
change_control <- change_metrics %>%
  filter(plot_mani==0) %>%
  select(-treatment, -plot_mani)%>%
  rename(richness_change_abs_ctrl = richness_change_abs,
         evenness_change_abs_ctrl = evenness_change_abs,
         rank_change_ctrl = rank_change,
         gains_ctrl = gains,
         losses_ctrl = losses) %>%
  group_by(site_project_comm, treatment_year2) %>%
  summarise_at(vars(c(rank_change_ctrl:losses_ctrl, richness_change_abs_ctrl, evenness_change_abs_ctrl)), funs(mean, sd), na.rm=T)

change_glass_d <- change_metrics %>%
  filter(plot_mani != 0) %>%
  group_by(site_project_comm, treatment, treatment_year2, plot_mani) %>%
  summarise(richness_change_abs = mean(richness_change_abs,na.rm=T),
            evenness_change_abs = mean(evenness_change_abs, na.rm=T),
            rank_change = mean(rank_change, na.rm=T),
            gains = mean(gains, na.rm=T),
            losses = mean(losses, na.rm=T)) %>%
  left_join(change_control, by=c("site_project_comm","treatment_year2")) %>%
  mutate(abs_richness_glass = (richness_change_abs-richness_change_abs_ctrl_mean)/richness_change_abs_ctrl_sd,
         abs_evenness_glass = (evenness_change_abs-evenness_change_abs_ctrl_mean)/evenness_change_abs_ctrl_sd,
         rank_glass = (rank_change-rank_change_ctrl_mean)/rank_change_ctrl_sd,
         gains_glass = (gains-gains_ctrl_mean)/gains_ctrl_sd,
         losses_glass = (losses-losses_ctrl_mean)/losses_ctrl_sd
  ) %>%
  dplyr::select(site_project_comm:plot_mani, abs_richness_glass:losses_glass) %>%
  ungroup()

#change_glass_d is the thing that we want

## replace Inf with 0 in change_glass_d. I am replaceing with 0 because I am ranking it all and this way it wont' rank.
change_glass_d <- change_glass_d %>%
  mutate(abs_richness_glass=replace(abs_richness_glass, abs_richness_glass=="Inf", NA)) %>%
  mutate(rank_glass=replace(rank_glass, rank_glass=="Inf", 0)) %>%
  mutate(gains_glass=replace(gains_glass, gains_glass=="Inf", 0)) %>%
  mutate(losses_glass=replace(losses_glass, losses_glass=="Inf", 0))%>%
  mutate(abs_evenness_glass=replace(abs_evenness_glass, abs_evenness_glass=="NaN", 0))%>%
  mutate(rank_glass=replace(rank_glass, rank_glass=="NaN", 0)) %>%
  mutate(gains_glass=replace(gains_glass, gains_glass=="NaN", 0))%>%
  mutate(losses_glass=replace(losses_glass, losses_glass=="NaN", 0))

##getting max glass D for each metric
max <- change_glass_d%>%
  group_by(site_project_comm, treatment)%>%
  mutate(S = max(abs_richness_glass),
         E = max(abs_evenness_glass),
         R = max(rank_glass),
         G = max(gains_glass),
         L = max(losses_glass),
         Smax = ifelse(S == abs_richness_glass, treatment_year2, 0),
         Emax = ifelse(E == abs_evenness_glass, treatment_year2, 0),
         Rmax = ifelse(R == rank_glass, treatment_year2, 0),
         Gmax = ifelse(G == gains_glass, treatment_year2, 0),
         Lmax = ifelse(L == losses_glass, treatment_year2, 0))%>%
  gather(max_metric, max_value, Smax:Lmax)%>%
  select(site_project_comm, treatment, max_metric, max_value)%>%
  filter(max_value!=0)

# ##ranking everything and then seeing when changed - not doing it this way
# rank<-max%>%
#   group_by(site_project_comm, treatment)%>%
#   filter(max_metric!="Smax")%>%
#   mutate(rank=rank(max_value, ties.method="average"))%>%
#   mutate(response_var = ifelse(max_metric=="Emax","evenness_change_abs",ifelse(max_metric=="Rmax","rank_change",ifelse(max_metric=="Gmax","gains", "losses"))))%>%
#   right_join(metrics_sig)
#   
# rank_mean<-rank%>%
#   group_by(response_var)%>%
#   summarize(mrank=mean(rank), sdrank=sd(rank), n=length(response_var))%>%
#   mutate(se = sdrank/sqrt(n))


  

##dropping what didn't change and then ranking everything - I think this is the way to do it.
rank_sig<-max%>%
  mutate(response_var = ifelse(max_metric=="Emax","evenness_change_abs",ifelse(max_metric=="Rmax","rank_change",ifelse(max_metric=="Gmax","gains", "losses"))))%>%
  right_join(metrics_sig)%>%
  group_by(site_project_comm, treatment)%>%
  filter(max_metric!="Smax")%>%
  mutate(rank=rank(max_value, ties.method="average"))

rank_sig_mean<-rank_sig%>%
  group_by(response_var)%>%
  summarize(mrank=mean(rank), sdrank=sd(rank), n=length(response_var))%>%
  mutate(se = sdrank/sqrt(n))

summary(aov(rank~response_var, data=rank_sig))

theme_set(theme_bw(12))
# ggplot(data=rank_mean, aes(x=response_var, y = mrank))+
#   geom_bar(stat="identity", position=position_dodge(0.9))+
#   geom_errorbar(aes(ymin=mrank-se, ymax=mrank+se), position=position_dodge(0.9), width=0.2)+
#   scale_x_discrete(limits=c("evenness_change_abs", "rank_change", "gains", "losses"), labels=c("Evenness","Rank","Gains","Losses"))+
#   ylab("Average Rank")+
#   xlab("Communtiy Change Metric")

ggplot(data=rank_sig_mean, aes(x=response_var, y = mrank))+
  geom_bar(stat="identity", position=position_dodge(0.9))+
  geom_errorbar(aes(ymin=mrank-se, ymax=mrank+se), position=position_dodge(0.9), width=0.2)+
  scale_x_discrete(limits=c("evenness_change_abs", "rank_change", "gains", "losses"), labels=c("Evenness","Rank","Gains","Losses"))+
  ylab("Average Rank")+
  xlab("Communtiy Change Metric")




##THS IS MAKING THE OLD PLOT WITH GLASS'S D
##need to only select those that saw sig change.
ggplot(data = subset(max), aes(x = max_metric, y = max_value))+
  geom_jitter()+
  geom_boxplot(alpha=.1) +
  xlab("Change Metric") +
  ylab("Treatment Year") +
  theme(axis.text=element_text(size=12, color="black"), strip.text.x=element_text(size=12)) +
  scale_x_discrete(limits=c("Smax","Emax","Rmax","Gmax","Lmax"),labels=c("Rich","Even","Rank","Gain","Loss"))+
  coord_flip()



