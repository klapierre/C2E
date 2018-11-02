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
#Meghan Work
setwd("C:\\Users\\megha\\Dropbox")

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


# ### to get all data time lags
# spct<-unique(corredat_sub$site_project_comm_trt)
# 
# rate_change_interval<-data.frame()
# 
# for (i in 1:length(spct)){
#   
#   subset<-corredat%>%
#     filter(site_project_comm_trt==spct[i])
#   
#   out<-rate_change_interval(subset, time.var = 'calendar_year', species.var = "genus_species", abundance.var = 'relcov', replicate.var = 'plot_id')
#   out$site_project_comm_trt<-spct[i]
#   
#   rate_change_interval<-rbind(rate_change_interval, out)
# }
# 
# write.csv(rate_change_interval, "C2E/Products/CommunityChange/Summer2018_Results/rate_change_interval.csv", row.names=F)

rate_change_interval <- read.csv("C2E/Products/CommunityChange/Summer2018_Results/rate_change_interval.csv")

##to get slopes and p-values for each treatment through time (i.e., putting in all plot to plot combinations into a single linear regression across replicates)
##get slopes for each treatment including controls

spct<-unique(rate_change_interval$site_project_comm_trt)

lm.slopes<-data.frame()
for (i in 1:length(spct)){
  subset<-rate_change_interval%>%
    filter(site_project_comm_trt==spct[i])
  test.lm<-lm(distance~interval, data=subset)
  output.lm<-data.frame(site_project_comm_trt=unique(subset$site_project_comm), 
                        est=summary(test.lm)$coef["interval", c("Estimate")], 
                        st.er=summary(test.lm)$coef["interval", c("Std. Error")],
                        p.val=summary(test.lm)$coef["interval","Pr(>|t|)"])
  lm.slopes<-rbind(lm.slopes, output.lm)
}

lm.slopes2<-lm.slopes%>%
  separate(site_project_comm_trt, into=c("site_project_comm","treatment"), sep = "::")%>%
  mutate(sig=ifelse(p.val<0.0501, 1, 0))%>%
  left_join(treatment_info)%>%
  mutate(trt=ifelse(plot_mani==0, "C", "T"))

prop_sig<-lm.slopes2%>%
  group_by(trt)%>%
  summarise(sum=sum(sig), n=length(sig))%>%
  mutate(prop_sig=sum/n)
# 74% of the time the controls are experiencing directional change and 73% of the time the treatments are experience directional change.


##test for sig diff between trt-control slopes
##there are so few differences that not going to pay attention to this.
rci2<-rate_change_interval%>%
  separate(site_project_comm_trt, into=c("site_project_comm","treatment"), sep = "::")%>%
  left_join(treatment_info)%>%
  mutate(trt=ifelse(plot_mani==0, "C", "T"))

spc<-unique(rci2$site_project_comm)
test.lm<-data.frame()
for (i in 1:length(spc)){
  subset<-rci2%>%
    filter(site_project_comm==spc[i])
  control<-subset%>%
    filter(plot_mani==0)
  treat<-subset%>%
  filter(plot_mani!=0)
trt_list<-unique(treat$treatment)
for (i in 1:length(trt_list)){
  subset2<-treat%>%
    filter(treatment==trt_list[i])
  trt<-trt_list[i]
  ct<-rbind(subset2, control)
  ct.lm<-lm(distance~interval*trt, data=ct)
  output.lm<-data.frame(site_project_comm=unique(subset$site_project_comm),
                        treatment=trt,
                        est=summary(ct.lm)$coef["interval:trtT", c("Estimate")],
                        pval=summary(ct.lm)$coef["interval:trtT","Pr(>|t|)"])
  test.lm<-rbind(test.lm, output.lm)
}
}

test.lm2<-test.lm%>%
  mutate(sig=ifelse(pval<0.0501, 1, 0))

#29% of the time the contorls and treatments have different slopes.
prop_diff<-sum(test.lm2$sig)/289


##import treatments
trts_interactions<-read.csv("C2E/Products/CommunityChange/March2018 WG/treatment interactions_July2018.csv")%>%
  select(-site_project_comm)%>%
  mutate(site_project_comm = paste(site_code, project_name, community_type, sep = "_"))

##of the trt-control differences, what treatments?
test_trt<-test.lm2%>%
  left_join(trts_interactions)%>%
  filter(use==1)

prop_sig_trt<-test_trt%>%
  group_by(trt_type)%>%
  summarise(sum=sum(sig), n=length(sig))%>%
  mutate(prop_sig=sum/n)

#test differences chi-sq.
chisq<-test_trt%>%
  mutate(sign=ifelse(sig==0, "num_nonsig", "num_sig"))%>%
  group_by(trt_type, sign)%>%
  summarise(sum=length(sig))%>%
  spread(sign, sum, fill=0)%>%
  filter(trt_type!="drought"&trt_type!="N+irr")#dropping drought bc only too few replicates.

#chi-sq
prop.test(x=as.matrix(chisq[,c('num_sig', 'num_nonsig')]), alternative='two.sided')
#p = 0.070 there is no difference between treatmetns in whether directional change occurs.


#Merge with exeriment info
#Using rate_change_interval - recalculate slopes and get slope and p-value (0.5)
#How often are we seeing directional change? OFTEN Does it differe treatment vs control? (plotmani=0) NO
#Does it different by trt type YES. Higher for N + P or N+PK over N, CO2 and P.