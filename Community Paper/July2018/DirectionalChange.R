### Calculating directional change for all exps
###
### Authors:  Sally Koerner (sekoerne@uncg.edu) & Meghan Avolio (meghan.avolio@jhu.edu)
### Last updated: Nov 1, 2018

library(tidyverse)
library(gridExtra)
library(grid)
library(gtable)
library(codyn)
library(vegan)
#home
setwd("~/Dropbox/")

#Files from home
corredat<-read.csv("converge_diverge/datasets/LongForm/SpeciesRelativeAbundance_Oct2017.csv")

#gvn face - only 2 years of data so will only have one point for the dataset, therefore we are removing this dataset from these analyses.
corredat1<-corredat%>%
  select(-X)%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))%>%
  filter(site_project_comm!="GVN_FACE_0", site_project_comm!="AZI_NitPhos_0", site_project_comm!="JRN_study278_0", site_project_comm!="KNZ_GFP_4F", site_project_comm!="Saskatchewan_CCD_0", project_name!="e001", project_name!="e002")

##several studies only have two measurments of a plot. I am dropping those plots
azi<-corredat%>%
  select(-X)%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))%>%
  filter(site_code=="AZI")%>%
  filter(plot_id!=11&plot_id!=15&plot_id!=35&plot_id!=37)

jrn<-corredat%>%
  select(-X)%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))%>%
  filter(site_project_comm=="JRN_study278_0")%>%
  filter(plot_id!=211&plot_id!=210)

knz<-corredat%>%
  select(-X)%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))%>%
  filter(site_project_comm=="KNZ_GFP_4F")%>%
  filter(plot_id!="7_1_1"&plot_id!="7_2_1")

sak<-corredat%>%
  select(-X)%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))%>%
  filter(site_project_comm=="Saskatchewan_CCD_0")%>%
  filter(plot_id!=2)

###remove extra treatments from CDR e001 and e002
cdr <- corredat%>%
  select(-X)%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))%>%
  filter(project_name=="e001"|project_name=="e002")%>%
  filter(treatment==1|treatment==6|treatment==8|treatment==9|treatment=='1_f_u_n'|treatment=='6_f_u_n'|treatment=='8_f_u_n'|treatment=='9_f_u_n')


###final dataset to use
corredat_raw<-rbind(corredat1, azi, jrn, knz, sak, cdr)

treatment_info<-read.csv("converge_diverge/datasets/LongForm/ExperimentInformation_Nov2017.csv")%>%
  dplyr::select(site_code, project_name, community_type, treatment,plot_mani)%>%
  unique()%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))

corredat <- corredat_raw %>%
  # left_join(treatment_info, by=c( "site_code","project_name","community_type", "treatment","site_project_comm"))%>%
              mutate(site_project_comm_trt=paste(site_project_comm, treatment, sep="::"))

numyears<- corredat%>%
  select(site_project_comm, calendar_year)%>%
  unique()%>%
  group_by(site_project_comm)%>%
  summarize(num=length(calendar_year))%>%
  filter(num>3)

corredat_sub<-corredat%>%
  right_join(numyears)
  

### need to subset to only long experiments (8 or more years)

#### now finalized dataset is ready and we will use corredat

#####look at directional change and slope - Dont use this because it takes forever AND we need the p-value anyways
# spct<-unique(corredat$site_project_comm_trt)
# rate_change<-data.frame()
# 
# for (i in 1:length(spct)){
#   
#   subset<-corredat%>%
#     filter(site_project_comm_trt==spct[i])
#   
#   out<-rate_change(subset, time.var = 'calendar_year', species.var = "genus_species", abundance.var = 'relcov', replicate.var = 'plot_id')
#   out$site_project_comm<-spct[i]
#   
#   rate_change<-rbind(rate_change, out)
# }
# 
# rate_change_4yrs<-rate_change%>%
#   mutate(spct=site_project_comm)%>%
#   separate (spct, into=c("site_project_comm", "treatment"), sep="::")%>%
#   right_join(numyears)
# write.csv(rate_change_4yrs, "C2E/Products/CommunityChange/Summer2018_Results/rate_change_interval_4orMoreYrs.csv", row.names=F)


### to get all data time lags
spct<-unique(corredat$site_project_comm_trt)

rate_change_interval<-data.frame()

for (i in 1:length(spct)){
  
  subset<-corredat%>%
    filter(site_project_comm_trt==spct[i])
  
  out<-rate_change_interval(subset, time.var = 'calendar_year', species.var = "genus_species", abundance.var = 'relcov', replicate.var = 'plot_id')
  out$site_project_comm<-spct[i]
  
  rate_change_interval<-rbind(rate_change_interval, out)
}
###only need this if we rerun with corredat_sub
rate_change_interval_4yrs<-rate_change_interval%>%
  mutate(spct=site_project_comm)%>%
  separate (spct, into=c("site_project_comm", "treatment"), sep="::")%>%
  right_join(numyears)%>%
  mutate(site_project_comm_trt=paste(site_project_comm, treatment, sep="::"))

##to get slopes and p-values for each treatment through time (i.e., putting in all plot to plot combinations into a single linear regression across replicates)
##get slopes for each treatment including controls
spct<-unique(rate_change_interval_4yrs$site_project_comm_trt)

lm.slopes<-data.frame()
for (i in 1:length(spct)){
  subset<-rate_change_interval_4yrs%>%
    filter(site_project_comm_trt==spct[i])
  test.lm<-lm(distance~interval, data=subset)
  output.lm<-data.frame(site_project_comm=unique(subset$site_project_comm), 
                        treatment=unique(subset$treatment), 
                        est=summary(test.lm)$coef["interval", c("Estimate")], 
                        st.er=summary(test.lm)$coef["interval", c("Std. Error")], 
                        p.val=summary(test.lm)$coef["interval","Pr(>|t|)"])
  lm.slopes<-rbind(lm.slopes, output.lm)
}

lm.slopes2<-lm.slopes%>%
  mutate(sig=ifelse(p.val<0.0501, 1, 0))%>%
  left_join(treatment_info)



##test for sig diff between trt-control slopes
##there are so few differences that not going to pay attention to this.

# spc2<-unique(anpp_precip$site_project_comm)
# test.lm<-data.frame()
# for (i in 1:length(spc2)){
#   subset<-anpp_precip%>%
#     filter(site_project_comm==spc2[i])
#   control<-subset%>%
#     filter(plot_mani==0)
#   treat<-subset%>%
#   filter(plot_mani!=0)
# trt_list<-unique(treat$treatment)
# for (i in 1:length(trt_list)){
#   subset2<-treat%>%
#     filter(treatment==trt_list[i])
#   trt<-trt_list[i]
#   ct<-rbind(subset2, control)
#   ct.lm<-lm(anpp~precip_mm*trt, data=ct)
#   output.lm<-data.frame(site_project_comm=unique(subset$site_project_comm), 
#                         treatment=trt, 
#                         est=summary(ct.lm)$coef["precip_mm:trtT", c("Estimate")],
#                         val=summary(ct.lm)$coef["precip_mm:trtT","Pr(>|t|)"])
#   test.lm<-rbind(test.lm, output.lm)
# }
# }
#write.csv(rate_change_interval, "converge_diverge/datasets/rate_change_interval.csv", 
         #row.names=F)

#Merge with exeriment info
#Using rate_change_interval - recalculate slopes and get slope and p-value (0.5)
#How often are we seeing directional change? Does it differe treatment vs control? (plotmani=0)
#Does it different by trt type