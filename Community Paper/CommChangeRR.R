################################################################################
##  CommChangesRR.R: This script does too much. I am focusing on the multiple regression models.
##
##  Author: Meghan Avolio (meghan.avolio@gmail.com)
##  Date: March 19, 2018
##  Update: July 30, 2018
################################################################################


library(tidyverse)
library(gridExtra)
library(grid)
library(gtable)
library(codyn)
library(vegan)

library(lme4)
library(relaimpo)

theme_set(theme_bw(12))

#read in the data
dat<-read.csv("~/Dropbox/C2E/Products/CommunityChange/March2018 WG/CORRE_RAC_Metrics_July2018_trtyr.csv")%>%
  select(-X)

dat_mult<-read.csv("~/Dropbox/C2E/Products/CommunityChange/March2018 WG/CORRE_Mult_Metrics_July2018.csv")%>%
select(-X)

trts_interactions<-read.csv("~/Dropbox/C2E/Products/CommunityChange/March2018 WG/treatment interactions_July2018.csv")%>%
  select(-site_project_comm)

sig_exp_bayes <- read.csv("~/Dropbox/C2E/Products/CommunityChange/March2018 WG/treatments_sig mult diff.csv")%>%
  select(site_project_comm, treatment)

unique(sig_exp_bayes$site_project_comm)

sig_exp_perm <- read.csv("~/Dropbox/C2E/Products/CommunityChange/March2018 WG/permanova out.csv")%>%
  filter(tot_pval>2)%>%
  select(site_project_comm, treatment)

unique(sig_exp_perm$site_project_comm)


plotid<-read.csv("~/Dropbox/converge_diverge/datasets/LongForm/SpeciesRelativeAbundance_Oct2017.csv")%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))

plotid2<-plotid%>%
  select(site_project_comm, plot_id, treatment)%>%
  unique()

trt<-read.csv("~/Dropbox/converge_diverge/datasets/LongForm/ExperimentInformation_Nov2017.csv")%>%
  select(site_code, project_name, community_type, treatment,plot_mani)%>%
  unique()%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))

# #filtering down to experiments and treatments that have sig change
# dat2<-dat%>%
#   left_join(plotid2)%>%
#   left_join(trt)%>%
#   right_join(sig_exp_bayes)
# 
# sig_exp2<-sig_exp_bayes%>%
#   select(-treatment)%>%
#   unique()
# 
# controls<-dat%>%
#   left_join(plotid2)%>%
#   left_join(trt)%>%
#   filter(plot_mani==0)%>%
#   right_join(sig_exp2)
# 
# sig_exp_dat<-rbind(controls, dat2)
# 
# write.csv(sig_exp_dat, "~/Dropbox/C2E/Products/CommunityChange/March2018 WG/CORRE_RACS_Subset_Perm.csv")


dat1<-dat%>%
  left_join(plotid2)%>%
  left_join(trt)%>%
  left_join(trts_interactions)%>%
  left_join(dat_mult)

pairs(dat1[4:8])

##overall what is the relationship
summary(mall<-lm(composition_change~richness_change+evenness_change+rank_change+gains+losses, data=dat1))
calc.relimp(mall, type="lmg", rela=F)

controls<-dat1%>%
  filter(plot_mani==0)%>%
  left_join(dat_mult)

##what is the relationship for the controls - rank change is the most important
summary(mc<-lm(composition_change~richness_change+evenness_change+rank_change+gains+losses, data=controls))
calc.relimp(mc, type="lmg", rela=F)

trts<-dat1%>%
  filter(plot_mani>0)%>%
  left_join(dat_mult)

##what is the relationship for the treatments - rank change is the most important
summary(mt<-lm(composition_change~richness_change+evenness_change+rank_change+gains+losses, data=trts))
calc.relimp(mt, type="lmg", rela=F)


##what is the relationship for the N trts - rank change is the most important
summary(mN<-lm(composition_change~richness_change+evenness_change+rank_change+gains+losses, data=subset(trts, trt_type == "N")))
calc.relimp(mN, type="lmg", rela=F)

##what is the relationship for the irg trts - rank change is the most important
summary(mi<-lm(composition_change~richness_change+evenness_change+rank_change+gains+losses, data=subset(trts, trt_type == "irr")))
calc.relimp(mi, type="lmg", rela=F)

##what is the relationship for the other trts - rank change is the most important
summary(mo<-lm(composition_change~richness_change+evenness_change+rank_change+gains+losses, data=subset(trts, trt_type == "other")))
calc.relimp(mo, type="lmg", rela=F)



####get average of all replicates in a treatment.
ave<-sig_exp_dat_perm%>%
  mutate_at(vars(richness_change, evenness_change, rank_change, gains, losses), abs)%>%
  group_by(treatment_year, treatment_year2, site_project_comm, treatment, plot_mani)%>%
  summarize_at(vars(richness_change, evenness_change, rank_change, gains, losses), funs(mean), na.rm = T)%>%
  mutate(trt=ifelse(plot_mani==0, "control","treatment"))


control<-ave%>%
  filter(plot_mani==0)%>%
  ungroup()%>%
  mutate(C_S=richness_change, C_even=evenness_change, C_gain=gains, C_loss=losses, C_rank=rank_change)%>%
  select(site_project_comm, treatment_year, treatment_year2, C_S, C_even, C_gain, C_loss, C_rank)


###do some log reponse ratio (natural log(T/C)) and add a constant for when 0 is present.
# how much to add for when it is zero
#rich
  rich <- ave$richness_change
  rich <- rich[rich !=0]
  min(rich) # 0.01333333
  rich_const<-0.001
#even
  even <- ave$evenness_change
  even <- even[even !=0]
  even <- na.omit(even)
  min(even) # 5.17201e-05
  even_const<-1e-06
#rank
  rank <- ave$rank_change
  rank <- rank[rank !=0]
  min(rank) # 0.01111111
  rank_const<-0.001
#gain
  gain <- ave$gains
  gain <- gain[gain !=0]
  min(gain) # 0.01086957
  gain_const<-0.001
#loss  
  loss <- ave$losses
  loss <- loss[loss !=0]
  min(loss) # 0.009259259
  loss_const<-0.0001
  
logRR<-merge(ave, control, by=c("site_project_comm","treatment_year", "treatment_year2"))%>%
  ungroup()%>%
  filter(plot_mani!=0)%>%
  mutate(richness_change1=ifelse(richness_change==0, 0.001, richness_change),
        evenness_change1=ifelse(evenness_change==0, 1e-06, evenness_change),
        rank_change1=ifelse(rank_change==0, 0.001, rank_change),
        gains1=ifelse(gains == 0, 0.001, gains),
        losses1=ifelse(losses==0,0.0001, losses),
        C_S1=ifelse(C_S==0, 0.001, C_S),
        C_even1=ifelse(C_even==0, 1e-06, C_even),
        C_rank1=ifelse(C_rank==0, 0.001, C_rank),
        C_gain1=ifelse(C_gain == 0, 0.001, C_gain),
        C_loss1=ifelse(C_loss==0,0.0001, C_loss))%>%
  mutate(lrS=log((richness_change1/C_S1)),
         lrE=log((evenness_change1/C_even1)),
         lrR = log((rank_change1/C_rank1)),
         lrL= log((losses1/C_loss1)),
         lrG= log((gains1/C_gain1)))%>%
 select(site_project_comm, treatment_year, treatment_year2, treatment,lrS, lrE, lrR, lrL, lrG)


###we do not like this because there are so many outliers when there is a zero in the dataset it really pulls it the data adding a lot of noise.
write.csv(logRR, "~/Dropbox/C2E/Products/CommunityChange/March2018 WG/CORRE_RAC_LogRR_March2018_trtyr.csv")

#####doing difference
diff<-merge(ave, control, by=c("site_project_comm","treatment_year", "treatment_year2"))%>%
  ungroup()%>%
  filter(plot_mani!=0)%>%
  mutate(D_Rich = richness_change - C_S,
         D_Even = evenness_change - C_even,
         D_Rank = rank_change - C_rank,
         D_Loss= losses - C_loss,
         D_Gain= gains - C_gain)%>%
  select(site_project_comm, treatment_year, treatment_year2, treatment,D_Rich, D_Even, D_Rank, D_Loss, D_Gain)%>%
  mutate(id = paste(site_project_comm, treatment, sep = "::"))

###making graphs of all treatments

S<-ggplot(data=diff, aes(x=treatment_year, y=D_Rich))+
  geom_line(aes(group=id))+
  geom_smooth(method="loess", se=F, color="red", size=1)+
  ggtitle("Richness")
even<-ggplot(data=diff, aes(x=treatment_year, y=D_Even))+
  geom_line(aes(group=id))+
  geom_smooth(method="loess", se=F, color="red", size=1)+
  ggtitle("Evenness")
gain<-ggplot(data=diff, aes(x=treatment_year, y=D_Gain))+
  geom_line(aes(group=id))+
  geom_smooth(method="loess", se=F, color="red", size=1)+
  ggtitle("Gains")
loss<-ggplot(data=diff, aes(x=treatment_year, y=D_Loss))+
  geom_line(aes(group=id))+
  geom_smooth(method="loess", se=F, color="red", size=1)+
  ggtitle("Losses")
rank<-ggplot(data=diff, aes(x=treatment_year, y=D_Rank))+
  geom_line(aes(group=id))+
  geom_smooth(method="loess", se=F, color="red", size=1)+
  ggtitle("Rank")

grid.arrange(S, even, gain, loss, rank)

#trying to color by treatments
trt2<-read.csv("~/Dropbox/converge_diverge/datasets/LongForm/ExperimentInformation_Nov2017.csv")%>%
  select(-calendar_year, -treatment_year, -X)%>%
  unique()

##code to add to the end of these steps to figure out how many experimetn-treatment compbinations are in each manipulation type.

# %>%
#   select(site_project_comm, treatment, trt_type)%>%
#   unique()%>%
#   group_by(trt_type)%>%
#   summarize(num=length(trt_type))

singleResource <- diff%>%
  left_join(trt2)%>%
  #filter pretrt data
  filter(treatment_year!=0)%>%
  #filter just single resource manipulations
  filter(resource_mani==1, plot_mani==1)%>%
  #set CEH Megarich nutrient values to 0 (added to all megaliths, not a treatment)
  mutate(n2=ifelse(site_code=='CEH', 0, n), p2=ifelse(site_code=='CEH', 0, p), k2=ifelse(site_code=='CEH', 0, k))%>%
  #drop lime added, as only one trt does this
  filter(other_trt!='lime added')%>%
  #create categorical treatment type column
  mutate(trt_type=ifelse(n2>0, 'N', ifelse(p2>0, 'P', ifelse(k2>0, 'K', ifelse(precip<0, 'drought', ifelse(precip>0, 'irr', ifelse(CO2>0, 'CO2', 'precip_vari')))))))%>%
  #keep just relevent column names for this analysis
  select(site_project_comm, treatment_year, treatment_year2, treatment, id, trt_type, D_Rich, D_Even, D_Rank, D_Loss, D_Gain)%>%
  gather(metric, value, D_Rich:D_Gain)

#analysis 2: single non-resource
singleNonresource <- diff%>%
  left_join(trt2)%>%
  #filter pretrt data
  filter(treatment_year!=0)%>%
  #filter just single resource manipulations
  filter(resource_mani==0, plot_mani==1)%>%
  #drop tilled, stone, and fungicide as only one trt each do these; drop 'current pattern' of rainfall, as this is basically a control
  filter(other_trt!='tilled', other_trt!='shallow soil', other_trt!='fungicide added', other_trt!='current pattern')%>%
  #create categorical treatment type column
  mutate(trt_type=ifelse(burn==1, 'burn', ifelse(mow_clip==1, 'mow_clip', ifelse(herb_removal==1, 'herb_rem', ifelse(temp>0, 'temp', ifelse(plant_trt==1, 'plant_mani', 'other'))))))%>%
  #keep just relevent column names for this analysis
  select(site_project_comm, treatment_year, treatment_year2, treatment, id, trt_type, D_Rich, D_Even, D_Rank, D_Loss, D_Gain)%>%
  gather(metric, value, D_Rich:D_Gain)


#analysis 3: 2-way interactions
twoWay <- diff%>%
  left_join(trt2)%>%
  #filter pretrt data
  filter(treatment_year!=0)%>%
  #filter just single resource manipulations
  filter(plot_mani==2)%>%
  #drop tilled, stone, and fungicide as only one trt each do these; drop 'current pattern' of rainfall, as this is basically a control
  filter(other_trt!='tilled', other_trt!='shallow soil', other_trt!='fungicide added', other_trt!='current pattern')%>%
  #create categorical treatment type column
  mutate(trt_type=ifelse(resource_mani==1&burn==1, 'R*burn', ifelse(resource_mani==1&mow_clip==1, 'R*mow_clip', ifelse(resource_mani==1&herb_removal==1, 'R*herb_rem', ifelse(resource_mani==1&temp>0, 'R*temp', ifelse(resource_mani==1&plant_trt==1, 'R*plant_mani', ifelse(resource_mani==1&other_trt!=0, 'R*other', ifelse(n>0&p>0, 'R*R', ifelse(n>0&CO2>0, 'R*R', ifelse(n>0&precip!=0, 'R*R', ifelse(p>0&k>0, 'R*R', ifelse(CO2>0&precip!=0, 'R*R', 'N*N'))))))))))))%>%
  #drop R*herb_removal (single rep)
  filter(trt_type!='R*herb_rem')%>%
  #keep just relevent column names for this analysis
  select(site_project_comm, treatment_year, treatment_year2, treatment, id, trt_type, D_Rich, D_Even, D_Rank, D_Loss, D_Gain)%>%
  gather(metric, value, D_Rich:D_Gain)

#analysis 4: 3+ way interactions
threeWay <- diff%>%
  left_join(trt2)%>%
  #filter pretrt data
  filter(treatment_year!=0)%>%
  #filter just single resource manipulations
  filter(plot_mani>2, plot_mani<6)%>%
  #drop tilled, stone, and fungicide as only one trt each do these; drop 'current pattern' of rainfall, as this is basically a control
  filter(other_trt!='tilled', other_trt!='shallow soil', other_trt!='fungicide added', other_trt!='current pattern')%>%
  #create categorical treatment type column
  mutate(trt_type=ifelse(burn==0&mow_clip==0&herb_removal==0&temp==0&plant_trt==0, 'all_resource', ifelse(n==0&p==0&k==0&CO2==0&precip==0, 'all_nonresource', 'both')))%>%
  #drop single all-nonresource treatment (NIN herbdiv 5NF)
  filter(trt_type!='all_nonresource')%>%
  #keep just relevent column names for this analysis
  select(site_project_comm, treatment_year, treatment_year2, treatment, id, trt_type, D_Rich, D_Even, D_Rank, D_Loss, D_Gain)%>%
  gather(metric, value, D_Rich:D_Gain)


###figures for the treatments
##single resources
ggplot(data=singleResource, aes(x=treatment_year, y=value)) +
  geom_line(aes(group = id)) +
  geom_smooth(method='loess', se=T, aes(color=trt_type)) +
  xlab('Treatment Year') + ylab('Diff (Trt-Cont)') +
  geom_hline(yintercept = 0)+
  facet_wrap(~metric, scales = "free")+
  ggtitle("Single Resource Raw")
##single nonresources
ggplot(data=singleNonresource, aes(x=treatment_year, y=value)) +
  geom_line(aes(group = id)) +
  geom_smooth(method='loess', se=T, aes(color=trt_type)) +
  xlab('Treatment Year') + ylab('Diff (Trt-Cont)') +
  geom_hline(yintercept = 0)+
  facet_wrap(~metric, scales = "free")+
  ggtitle("Single Non Resource Raw")
##twoWay
ggplot(data=twoWay, aes(x=treatment_year, y=value)) +
  geom_line(aes(group = id)) +
  geom_smooth(method='loess', se=T, aes(color=trt_type)) +
  xlab('Treatment Year') + ylab('Diff (Trt-Cont)') +
  geom_hline(yintercept = 0)+
  facet_wrap(~metric, scales = "free")+
  ggtitle("Two Manip Raw")
##three way
ggplot(data=threeWay, aes(x=treatment_year, y=value)) +
  geom_line(aes(group = id)) +
  geom_smooth(method='loess', se=T, aes(color=trt_type)) +
  xlab('Treatment Year') + ylab('Diff (Trt-Cont)') +
  geom_hline(yintercept = 0)+
  facet_wrap(~metric, scales = "free")+
  ggtitle("Two + Manip Raw")


#####doing this for cumulative
change_cumsum <- sig_exp_dat_perm %>%
  arrange(site_project_comm, treatment, plot_id, treatment_year)%>%
  group_by(site_project_comm, treatment, plot_id) %>%
  mutate(richness_change_abs = abs(richness_change)) %>%
  mutate(evenness_change_abs = abs(evenness_change)) %>%
  mutate_at(vars(richness_change, richness_change_abs, evenness_change,evenness_change_abs, rank_change, gains, losses), funs(cumsum) ) %>%
  mutate(control = ifelse(plot_mani==0,"control","treatment"))

ave_cum<-change_cumsum%>%
  mutate_at(vars(richness_change, evenness_change, rank_change, gains, losses), abs)%>%
  group_by(treatment_year, treatment_year2, site_project_comm, treatment, plot_mani)%>%
  summarize_at(vars(richness_change, evenness_change, rank_change, gains, losses), funs(mean), na.rm = T)%>%
  mutate(trt=ifelse(plot_mani==0, "control","treatment"))


control_cum<-ave_cum%>%
  filter(plot_mani==0)%>%
  ungroup()%>%
  mutate(C_S=richness_change, C_even=evenness_change, C_gain=gains, C_loss=losses, C_rank=rank_change)%>%
  select(site_project_comm, treatment_year, treatment_year2, C_S, C_even, C_gain, C_loss, C_rank)

diff_cum<-merge(ave_cum, control_cum, by=c("site_project_comm","treatment_year", "treatment_year2"))%>%
  ungroup()%>%
  filter(plot_mani!=0)%>%
  mutate(D_Rich = richness_change - C_S,
         D_Even = evenness_change - C_even,
         D_Rank = rank_change - C_rank,
         D_Loss= losses - C_loss,
         D_Gain= gains - C_gain)%>%
  select(site_project_comm, treatment_year, treatment_year2, treatment,D_Rich, D_Even, D_Rank, D_Loss, D_Gain)%>%
  mutate(id = paste(site_project_comm, treatment, sep = "::"))

singleResource_cumulative <- diff_cum%>%
  left_join(trt2)%>%
  #filter pretrt data
  filter(treatment_year!=0)%>%
  #filter just single resource manipulations
  filter(resource_mani==1, plot_mani==1)%>%
  #set CEH Megarich nutrient values to 0 (added to all megaliths, not a treatment)
  mutate(n2=ifelse(site_code=='CEH', 0, n), p2=ifelse(site_code=='CEH', 0, p), k2=ifelse(site_code=='CEH', 0, k))%>%
  #drop lime added, as only one trt does this
  filter(other_trt!='lime added')%>%
  #create categorical treatment type column
  mutate(trt_type=ifelse(n2>0, 'N', ifelse(p2>0, 'P', ifelse(k2>0, 'K', ifelse(precip<0, 'drought', ifelse(precip>0, 'irr', ifelse(CO2>0, 'CO2', 'precip_vari')))))))%>%
  #keep just relevent column names for this analysis
  select(site_project_comm, treatment_year, treatment_year2, treatment, id, trt_type, D_Rich, D_Even, D_Rank, D_Loss, D_Gain)%>%
  gather(metric, value, D_Rich:D_Gain)

#analysis 2: single non-resource
singleNonresource_cumulative <- diff_cum%>%
  left_join(trt2)%>%
  #filter pretrt data
  filter(treatment_year!=0)%>%
  #filter just single resource manipulations
  filter(resource_mani==0, plot_mani==1)%>%
  #drop tilled, stone, and fungicide as only one trt each do these; drop 'current pattern' of rainfall, as this is basically a control
  filter(other_trt!='tilled', other_trt!='shallow soil', other_trt!='fungicide added', other_trt!='current pattern')%>%
  #create categorical treatment type column
  mutate(trt_type=ifelse(burn==1, 'burn', ifelse(mow_clip==1, 'mow_clip', ifelse(herb_removal==1, 'herb_rem', ifelse(temp>0, 'temp', ifelse(plant_trt==1, 'plant_mani', 'other'))))))%>%
  #keep just relevent column names for this analysis
  select(site_project_comm, treatment_year, treatment_year2, treatment, id, trt_type, D_Rich, D_Even, D_Rank, D_Loss, D_Gain)%>%
  gather(metric, value, D_Rich:D_Gain)


#analysis 3: 2-way interactions
twoWay_cumulative <- diff_cum%>%
  left_join(trt2)%>%
  #filter pretrt data
  filter(treatment_year!=0)%>%
  #filter just single resource manipulations
  filter(plot_mani==2)%>%
  #drop tilled, stone, and fungicide as only one trt each do these; drop 'current pattern' of rainfall, as this is basically a control
  filter(other_trt!='tilled', other_trt!='shallow soil', other_trt!='fungicide added', other_trt!='current pattern')%>%
  #create categorical treatment type column
  mutate(trt_type=ifelse(resource_mani==1&burn==1, 'R*burn', ifelse(resource_mani==1&mow_clip==1, 'R*mow_clip', ifelse(resource_mani==1&herb_removal==1, 'R*herb_rem', ifelse(resource_mani==1&temp>0, 'R*temp', ifelse(resource_mani==1&plant_trt==1, 'R*plant_mani', ifelse(resource_mani==1&other_trt!=0, 'R*other', ifelse(n>0&p>0, 'R*R', ifelse(n>0&CO2>0, 'R*R', ifelse(n>0&precip!=0, 'R*R', ifelse(p>0&k>0, 'R*R', ifelse(CO2>0&precip!=0, 'R*R', 'N*N'))))))))))))%>%
  #drop R*herb_removal (single rep)
  filter(trt_type!='R*herb_rem')%>%
  #keep just relevent column names for this analysis
  select(site_project_comm, treatment_year, treatment_year2, treatment, id, trt_type, D_Rich, D_Even, D_Rank, D_Loss, D_Gain)%>%
  gather(metric, value, D_Rich:D_Gain)

#analysis 4: 3+ way interactions
threeWay_cumulative <- diff_cum%>%
  left_join(trt2)%>%
  #filter pretrt data
  filter(treatment_year!=0)%>%
  #filter just single resource manipulations
  filter(plot_mani>2, plot_mani<6)%>%
  #drop tilled, stone, and fungicide as only one trt each do these; drop 'current pattern' of rainfall, as this is basically a control
  filter(other_trt!='tilled', other_trt!='shallow soil', other_trt!='fungicide added', other_trt!='current pattern')%>%
  #create categorical treatment type column
  mutate(trt_type=ifelse(burn==0&mow_clip==0&herb_removal==0&temp==0&plant_trt==0, 'all_resource', ifelse(n==0&p==0&k==0&CO2==0&precip==0, 'all_nonresource', 'both')))%>%
  #drop single all-nonresource treatment (NIN herbdiv 5NF)
  filter(trt_type!='all_nonresource')%>%
  #keep just relevent column names for this analysis
  select(site_project_comm, treatment_year, treatment_year2, treatment, id, trt_type, D_Rich, D_Even, D_Rank, D_Loss, D_Gain)%>%
  gather(metric, value, D_Rich:D_Gain)


###figures for the treatments
##single resources
ggplot(data=singleResource_cumulative, aes(x=treatment_year, y=value)) +
  geom_line(aes(group = id)) +
  geom_smooth(method='loess', se=T, aes(color=trt_type)) +
  xlab('Treatment Year') + ylab('Diff (Trt-Cont)') +
  geom_hline(yintercept = 0)+
  facet_wrap(~metric, scales = "free")+
  ggtitle("Single Resouce Cumulative")
##single nonresources
ggplot(data=singleNonresource_cumulative, aes(x=treatment_year, y=value)) +
  geom_line(aes(group = id)) +
  geom_smooth(method='loess', se=T, aes(color=trt_type)) +
  xlab('Treatment Year') + ylab('Diff (Trt-Cont)') +
  geom_hline(yintercept = 0)+
  facet_wrap(~metric, scales = "free")+
  ggtitle("Single Non Resource Cumulative")
##twoWay
ggplot(data=twoWay_cumulative, aes(x=treatment_year, y=value)) +
  geom_line(aes(group = id)) +
  geom_smooth(method='loess', se=T, aes(color=trt_type)) +
  xlab('Treatment Year') + ylab('Diff (Trt-Cont)') +
  geom_hline(yintercept = 0)+
  facet_wrap(~metric, scales = "free")+
  ggtitle("Two Manip Cumulative")
##three way
ggplot(data=threeWay_cumulative, aes(x=treatment_year, y=value)) +
  geom_line(aes(group = id)) +
  geom_smooth(method='loess', se=T, aes(color=trt_type)) +
  xlab('Treatment Year') + ylab('Diff (Trt-Cont)') +
  geom_hline(yintercept = 0)+
  facet_wrap(~metric, scales = "free")+
  ggtitle("Two + Manip Cumulative")


###investigateing crd2a - it is real

sub<-sig_exp_dat_perm%>%
  filter(site_project_comm == "CDR_e002_A")

ggplot(data=sub, aes(x = as.factor(treatment_year), y = rank_change))+
  geom_boxplot()+
  facet_grid(~treatment)

###historgrams
hist<-dat%>%
  gather(metric, value, richness_change:losses)

ggplot(data=subset(hist, metric == "evenness_change"), aes(x = value))+
  geom_histogram()+
  stat_bin(binwidth = 0.0001)+
  xlim(-0.01, 0.01)

ggplot(data=hist, aes(x = value))+
  geom_histogram()+
  facet_wrap(~metric, scales = "free")

hist_evar<-delta_rac%>%
  gather(metric, value, richness_change:losses)

ggplot(data=hist_evar, aes(x = value))+
  geom_histogram()+
  facet_wrap(~metric, scales = "free")
