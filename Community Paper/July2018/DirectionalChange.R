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
setwd("C:\\Users\\mavolio2\\Dropbox")

#to subset only the treatments I want
subset_studies<-read.csv("C2E/Products/CommunityChange/March2018 WG/experiment_trt_subset_May2019.csv")

#Files from home
corredat<-read.csv("converge_diverge/datasets/LongForm/SpeciesRelativeAbundance_March2019.csv")%>%
  select(-X)

#gvn face - only 2 years of data so will only have one point for the dataset, therefore we are removing this dataset from these analyses.
corredat1<-corredat%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))%>%
  filter(site_project_comm!="GVN_FACE_0", site_project_comm!="AZI_NitPhos_0", site_project_comm!="JRN_study278_0", site_project_comm!="KNZ_GFP_4F", site_project_comm!="Saskatchewan_CCD_0", project_name!="e001", project_name!="e002")

##several studies only have two measurments of a plot. I am dropping those plots
azi<-corredat%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))%>%
  filter(site_code=="AZI")%>%
  filter(plot_id!=11&plot_id!=15&plot_id!=35&plot_id!=37)

jrn<-corredat%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))%>%
  filter(site_project_comm=="JRN_study278_0")%>%
  filter(plot_id!=211&plot_id!=210)

knz<-corredat%>%
    mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))%>%
  filter(site_project_comm=="KNZ_GFP_4F")%>%
  filter(plot_id!="7_1_1"&plot_id!="7_2_1")

sak<-corredat%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))%>%
  filter(site_project_comm=="Saskatchewan_CCD_0")%>%
  filter(plot_id!=2)

###remove extra treatments from CDR e001 and e002
cdr <- corredat%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))%>%
  filter(project_name=="e001"|project_name=="e002")%>%
  filter(treatment==1|treatment==6|treatment==8|treatment==9|treatment=='1_f_u_n'|treatment=='6_f_u_n'|treatment=='8_f_u_n'|treatment=='9_f_u_n')

##remove one of 2 pre-treatment years in edge for CHY, SGS, and HAYS
edge<-corredat%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))%>%
  filter(site_code=="CHY"|site_code=="SGS"|site_code=="HYS"&project_name=="EDGE")%>%
  filter(calendar_year!=2012)
  


###final dataset to use
corredat_raw<-rbind(corredat1, azi, jrn, knz, sak, cdr, edge)

treatment_info<-read.csv("converge_diverge/datasets/LongForm/ExperimentInformation_March2019.csv")%>%
  select(site_code, project_name, community_type, treatment,plot_mani, trt_type)%>%
  unique()%>%
  #filter(plot_mani!=0)%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))%>%
  mutate(use=ifelse(trt_type=="N"|trt_type=="P"|trt_type=="CO2"|trt_type=="irr"|trt_type=="temp"|trt_type=="N*P"|trt_type=="mult_nutrient"|trt_type=='precip_vari', 1, 0))%>%
  mutate(trt_type2=ifelse(trt_type=="N"|trt_type=="control","N", 
                          ifelse(trt_type=="P", "P", 
                                 ifelse(trt_type=="CO2", "CO2",
                                        ifelse(trt_type=="irr", "Irrigation",
                                               ifelse(trt_type=="temp", "Temperature", 
                                                      ifelse(trt_type=="N*P"|trt_type=="mult_nutrient", "Mult. Nuts.", 
                                                             ifelse(trt_type=="drought", "drought", 
                                                                    ifelse(trt_type=="CO2*temp", "CO2*temp", 
                                                                           ifelse(trt_type=="drought*temp", "drought*temp", 
                                                                                  ifelse(trt_type=="irr*temp", "irr*temp",
                                                                                         ifelse(trt_type=="irr*CO2*temp"|trt_type=="N*CO2*temp"|trt_type=="N*irr*temp"|trt_type=="N*irr*CO2*temp", "mult_res*temp", 
                                                                                                ifelse(trt_type=="irr*herb_removal"|trt_type=="irr*plant_mani"|trt_type=="irr*plant_mani*herb_removal", "irr*NR", 
                                                                                                       ifelse(trt_type=="herb_removal"|trt_type=="till"|trt_type=="mow_clip"|trt_type=="burn"|trt_type=="plant_mani"|trt_type=="stone"|trt_type=="graze"|trt_type=="burn*graze"|trt_type=="fungicide"|trt_type=="plant_mani*herb_removal"|trt_type=="burn*mow_clip", "NR", 
                                                                                                              ifelse(trt_type=="precip_vari", "Precip. Vari.",  
                                                                                                                     ifelse(trt_type=="N*plant_mani"|trt_type=="N*burn"|trt_type=="N*mow_clip"|trt_type=="N*till"|trt_type=="N*stone"|trt_type=="N*burn*graze"|trt_type=="N*burn*mow_clip", "N*NR", 
                                                                                                                            ifelse(trt_type=="N*temp", "N*temp", 
                                                                                                                                   ifelse(trt_type=="N*CO2", "N*CO2",
                                                                                                                                          ifelse(trt_type=="irr*CO2", "irr*CO2",
                                                                                                                                                 ifelse(trt_type=="N*irr", "N*irr",
                                                                                                                                                        ifelse(trt_type=="mult_nutrient*herb_removal"|trt_type=="mult_nutrient*fungicide"|trt_type=="N*P*burn*graze"|trt_type=="N*P*burn"|trt_type=="*P*mow_clip"|trt_type=="N*P*burn*mow_clip"|trt_type=="N*P*mow_clip", "mult_nutrients*NR",
                                                                                                                                                               ifelse(trt_type=="P*mow_clip"|trt_type=="P*burn"|trt_type=="P*burn*graze"|trt_type=="P*burn*mow_clip", "P*NR", 
                                                                                                                                                                      ifelse(trt_type=="precip_vari*temp", "precip_vari*temp", 
                                                                                                                                                                             ifelse(trt_type=="N*irr*CO2", "mult_res", 999))))))))))))))))))))))))

corredat <- corredat_raw %>%
  # left_join(treatment_info, by=c( "site_code","project_name","community_type", "treatment","site_project_comm"))%>%
              mutate(site_project_comm_trt=paste(site_project_comm, treatment, sep="::"))

info.trt2<-read.csv("C2E\\Products\\CommunityChange\\March2018 WG\\ForAnalysis_allAnalysisAllDatasets_04082019.csv") %>%
  mutate(site_project_comm = paste(site_code, project_name, community_type, sep="_"))%>%
  select(site_project_comm, treatment, trt_type)%>%
  unique()%>%
  mutate(trt_type3=ifelse(trt_type=="drought"|trt_type=="irr"|trt_type=="N"|trt_type=="precip_vari"|trt_type=="P"|trt_type=="CO2"|trt_type=="other_resource", "Res.", ifelse(trt_type=="mow_clip"|trt_type=="temp"|trt_type=="plant_mani"|trt_type=="other_nonresource"|trt_type=="herb_removal"|trt_type=="NxN", "Non-Res.", ifelse(trt_type=="RxR"|trt_type=="RxRxR", "Mult. Res.", ifelse(trt_type=="threeway"|trt_type=="RxN", "Res.+Non-Res.", "oops")))))%>%
  select(-trt_type)

numyears<- corredat%>%
  select(site_project_comm, treatment_year)%>%
  unique()%>%
  filter(treatment_year!=0)%>%
  group_by(site_project_comm)%>%
  summarize(num=length(treatment_year))%>%
  filter(num>4)%>%
  mutate(use=1)

corredat_sub<-corredat_raw%>%
  right_join(numyears)%>%
  mutate(site_project_comm_trt=paste(site_project_comm, treatment, sep="::"))
  

### need to subset to only long experiments (8 or more years)

#### now finalized dataset is ready and we will use corredat

#####look at directional change and slope - Dont use this because it takes forever AND we need the p-value anyways
# spct<-unique(corredat_sub$site_project_comm_trt)
# rate_change<-data.frame()
# 
# for (i in 1:length(spct)){
# 
#   subset<-corredat_sub%>%
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
spct<-unique(corredat_sub$site_project_comm_trt)

rate_change_interval<-data.frame()

for (i in 1:length(spct)){

  subset<-corredat_sub%>%
    filter(site_project_comm_trt==spct[i])

  out<-rate_change_interval(subset, time.var = 'calendar_year', species.var = "genus_species", abundance.var = 'relcov', replicate.var = 'plot_id')
  out$site_project_comm_trt<-spct[i]

  rate_change_interval<-rbind(rate_change_interval, out)
}

# write.csv(rate_change_interval, "C2E/Products/CommunityChange/Summer2018_Results/rate_change_interval_may2019.csv", row.names=F)

rate_change_interval <- read.csv("C2E/Products/CommunityChange/Summer2018_Results/rate_change_interval_may2019.csv")

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

cslopes<-lm.slopes2%>%
  filter(trt=="C")%>%
  rename(c_est=est, 
         c_pval=p.val,
         c_sig=sig)%>%
  select(site_project_comm, c_est, c_pval, c_sig)

diffslopes<-lm.slopes2%>%
  filter(trt=="T")%>%
  right_join(subset_studies)%>%
  left_join(cslopes)%>%
  mutate(diff=est-c_est)

#write.csv(diffslopes, "C2E/Products/CommunityChange/Summer2018_Results/diff_directinal_slopes.csv", row.names=F)

conly<-cslopes<-lm.slopes2%>%
  filter(trt=="C")
prop_sig<-lm.slopes2%>%
  filter(trt=="T")%>%
  right_join(subset_studies)%>%
  bind_rows(conly)%>%
  group_by(trt)%>%
  summarise(sum=sum(sig), n=length(sig))%>%
  mutate(prop_sig=sum/n)
# 77% of the time the controls are experiencing directional change and 84% of the time the treatments are experience directional change - this is not the best way to do this.


##test for sig diff between trt-control slopes
##there are so few differences that not going to pay attention to this.
rci2<-rate_change_interval%>%
  separate(site_project_comm_trt, into=c("site_project_comm","treatment"), sep = "::")%>%
  left_join(treatment_info)%>%
  mutate(trt=ifelse(plot_mani==0, "C", "T"))


#plot this - need to sqrt the x-axis.
ggplot(data=subset(rci2, site_project_comm=='SERC_CXN_0'&treatment=="t1"|treatment=="t4"), aes(x = sqrt(interval), y = distance, color=trt))+
  #geom_point()+
  geom_smooth(method = "loess")+
  geom_jitter()

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
for (j in 1:length(trt_list)){
  subset2<-treat%>%
    filter(treatment==trt_list[j])
  trt<-trt_list[j]
  ct<-rbind(subset2, control)
  ct.lm<-lm(distance~sqrt(interval)*trt, data=ct)
  output.lm<-data.frame(site_project_comm=unique(subset$site_project_comm),
                        treatment=trt,
                        dest=summary(ct.lm)$coef["sqrt(interval):trtT", c("Estimate")],
                        dpval=summary(ct.lm)$coef["sqrt(interval):trtT","Pr(>|t|)"])
  test.lm<-rbind(test.lm, output.lm)
}
}


# ###loop making pictures
# spc<-unique(rci2$site_project_comm)
# test.lm<-data.frame()
# for (i in 1:length(spc)){
#   subset<-rci2%>%
#     filter(site_project_comm==spc[i])
#   control<-subset%>%
#     filter(plot_mani==0)
#   treat<-subset%>%
#     filter(plot_mani!=0)
#   trt_list<-unique(treat$treatment)
#   for (j in 1:length(trt_list)){
#     subset2<-treat%>%
#       filter(treatment==trt_list[j])
#     trt<-trt_list[j]
#     ct<-rbind(subset2, control)
# 
#     print(ggplot(data=ct, aes(x = sqrt(interval), y = distance, color=trt))+
#       geom_point()+
#       geom_smooth(method = "loess", se=F)+
#       geom_smooth(method = "lm", se=F, linetype=2)+
#       ggtitle(spc[i]))
#   }}
    
    

test.lm2<-test.lm%>%
  mutate(dsig=ifelse(dpval<0.0501, 1, 0))%>%
  right_join(subset_studies)

#32% of the time the contorls and treatments have different slopes.
prop_diff<-sum(test.lm2$dsig)/219


diffslopes2<-test.lm2%>%
  left_join(diffslopes)%>%
  select(site_project_comm, treatment, dsig, est, sig, c_est, c_sig, diff)

###binning the reuslts to 3 categoreis, no diff, C>T, C<T"
tally_all<-diffslopes2%>%
  mutate(resp=ifelse(dsig==0, "c", ifelse(dsig==1&c_est>est, "b", ifelse(dsig==1&c_est<est, "a", 999))))%>%
  group_by(resp)%>%
  summarize(num=length(resp))%>%
  mutate(trt_type2="All GCDs",
         pct = num/219)

diffslope_trt<-diffslopes2%>%
  left_join(treatment_info)%>%
  filter(use==1)%>%
  select(site_project_comm, treatment, dsig, est, sig, c_est, c_sig, diff, trt_type2)%>%
  mutate(resp=ifelse(dsig==0, "c", ifelse(dsig==1&c_est>est, "b", ifelse(dsig==1&c_est<est, "a", 999))))%>%
  group_by(trt_type2, resp)%>%
  summarize(num=length(resp))%>%
  mutate(pct = ifelse(trt_type2=="CO2", num/7, ifelse(trt_type2=="Irrigation",num/12, ifelse(trt_type2=="Precip. Vari.", num/8, ifelse(trt_type2=="Temperature", num/7, ifelse(trt_type2=="N" , num/31, ifelse(trt_type2=="P", num/9, ifelse(trt_type2=="Mult. Nuts.", num/51,999))))))))

diffslope_trt3<-diffslopes2%>%
  left_join(info.trt2)%>%
  select(site_project_comm, treatment, dsig, est, sig, c_est, c_sig, diff, trt_type3)%>%
  mutate(resp=ifelse(dsig==0, "c", ifelse(dsig==1&c_est>est, "b", ifelse(dsig==1&c_est<est, "a", 999))))%>%
  group_by(trt_type3, resp)%>%
  summarize(num=length(resp))%>%
  mutate(pct = ifelse(trt_type3=="Res.", num/70, ifelse(trt_type3=="Mult. Res.",num/60, ifelse(trt_type3=="Res.+Non-Res.", num/58, ifelse(trt_type3=="Non-Res.", num/31, 999)))))%>%
  rename(trt_type2=trt_type3)
  
tograph<-diffslope_trt%>%
  bind_rows(tally_all, diffslope_trt3)

ggplot(tograph, aes(x = trt_type2, y = pct, fill = resp)) +
  geom_col(width = 0.7) +
  coord_flip() +
  theme_minimal()+
  scale_fill_manual(name = "Response", limits=c("c", "b", "a"),labels = c("C=T", "C>T", "C<T"), values=c("light gray", "skyblue",'gold')) +
  scale_x_discrete(limits=c( 'Mult. Nuts.','P', 'N',"Temperature",'Precip. Vari.', 'Irrigation','CO2', 'Res.+Non-Res.', "Mult. Res.", "Res.", 'Non-Res.',  'All GCDs'), labels=c("Mult. Nuts.", "Phosphorus","Nitrogen","Temperature" , "Precip. Vari.","Irrigation","CO2", "Res.+Non-Res.","Multiple Res.","Single Res.","Non-Res." , "All GCDs"))+
  labs(x = "Treatment", y = "Proportion of communities") +
  theme(legend.position = "top")+
  geom_vline(xintercept = 7.5, linetype="dashed")+
  geom_vline(xintercept = 11.5, linetype="dashed")
  # geom_text(x=1, y = 0.05, label="61%", size=4,check_overlap = TRUE)+
  # geom_text(x=2, y = 0.05, label="11%", size=4,check_overlap = TRUE)+
  # geom_text(x=3, y = 0.05, label="19%", size=4,check_overlap = TRUE)+
  # geom_text(x=4, y = 0.05, label="29%", size=4,check_overlap = TRUE)+
  # geom_text(x=5, y = 0.05, label="38%", size=4,check_overlap = TRUE)+
  # geom_text(x=6, y = 0.05, label="25%", size=4,check_overlap = TRUE)+
  # geom_text(x=7, y = 0.05, label="14%", size=4,check_overlap = TRUE)+
  # geom_text(x=8, y = 0.05, label="32%", size=4, check_overlap = T)


##this bins to 5 categories
# tally_all<-diffslopes2%>%
#   mutate(resp=ifelse(dsig==0, "e", ifelse(c_sig==0&sig==1, "b", ifelse(c_sig==1&sig==0, "d", ifelse(c_sig==1&sig==1&c_est>est, "c", ifelse(c_sig==1&sig==1&c_est<est, "a", 999))))))%>%
#   group_by(resp)%>%
#   summarize(num=length(resp))%>%
#   mutate(trt_type2="All GCDs",
#          pct = num/219)
# 
# diffslope_trt<-diffslopes2%>%
#   left_join(treatment_info)%>%
#   filter(use==1)%>%
#   select(site_project_comm, treatment, dsig, est, sig, c_est, c_sig, diff, trt_type2)%>%
#   mutate(resp=ifelse(dsig==0, "e", ifelse(c_sig==0&sig==1, "b", ifelse(c_sig==1&sig==0, "d", ifelse(c_sig==1&sig==1&c_est>est, "c", ifelse(c_sig==1&sig==1&c_est<est, "a", 999))))))%>%
#   group_by(trt_type2, resp)%>%
#   summarize(num=length(resp))%>%
#   mutate(pct = ifelse(trt_type2=="CO2", num/7, ifelse(trt_type2=="Irrigation",num/12, ifelse(trt_type2=="Precip. Vari.", num/8, ifelse(trt_type2=="Temperature", num/7, ifelse(trt_type2=="N" , num/31, ifelse(trt_type2=="P", num/9, ifelse(trt_type2=="Mult. Nuts.", num/51,999))))))))%>%
#   bind_rows(tally_all)
# 
# ggplot(diffslope_trt, aes(x = trt_type2, y = pct, fill = resp)) +
#   geom_col(width = 0.7) +
#   coord_flip() +
#   theme_minimal()+
#   scale_fill_manual(name = "Response", limits=c("e", "d", "c", "b", "a"),labels = c("C=T", "CdTnd","C>T", "CndTd", "C<T"), values=c("light gray", "skyblue","steelblue", "khaki1",'gold')) +
#   scale_x_discrete(limits=c( 'Mult. Nuts.','P', 'N',"Temperature",'Precip. Vari.', 'Irrigation','CO2',  'All GCDs'))+
#   labs(x = "Treatment", y = "Proportion of communities") +
#   theme(legend.position = "top")+
#   geom_vline(xintercept = 7.5, linetype="dashed")+
#   geom_text(x=1, y = 0.05, label="61%", size=4,check_overlap = TRUE)+
#   geom_text(x=2, y = 0.05, label="11%", size=4,check_overlap = TRUE)+
#   geom_text(x=3, y = 0.05, label="19%", size=4,check_overlap = TRUE)+
#   geom_text(x=4, y = 0.05, label="29%", size=4,check_overlap = TRUE)+
#   geom_text(x=5, y = 0.05, label="38%", size=4,check_overlap = TRUE)+
#   geom_text(x=6, y = 0.05, label="25%", size=4,check_overlap = TRUE)+
#   geom_text(x=7, y = 0.05, label="14%", size=4,check_overlap = TRUE)+
#   geom_text(x=8, y = 0.05, label="32%", size=4, check_overlap = T)

##of the trt-control differences, what treatments?
test_trt<-test.lm2%>%
  left_join(treatment_info)%>%
  filter(use==1)

test_trt3<-test.lm2%>%
  left_join(info.trt2)


alldir<-data.frame(trt_type2=c('All Trts.', 'All Trts.'), n=c(219,219), value=c('prop_sig', 'pnotsig'), sig=c(0.324, 0.676))

prop_sig_trt<-test_trt%>%
  group_by(trt_type2)%>%
  summarise(sum=sum(dsig), n=length(dsig))%>%
  mutate(prop_sig=sum/n)%>%
    filter(trt_type2!="Irr + Temp")%>%
  mutate(pnotsig=1-prop_sig)%>%
  select(-sum)%>%
  gather(value, sig, prop_sig:pnotsig)

#test differences chi-sq.
chisq<-test_trt%>%
  mutate(sign=ifelse(dsig==0, "num_nonsig", "num_sig"))%>%
  group_by(trt_type2, sign)%>%
  summarise(sum=length(dsig))%>%
  spread(sign, sum, fill=0)%>%
  filter(trt_type2!="Irr + Temp")

#chi-sq
prop.test(x=as.matrix(chisq[,c('num_sig', 'num_nonsig')]), alternative='two.sided')
#p = 0.0002

prop_sig_trt3<-test_trt3%>%
  group_by(trt_type3)%>%
  summarise(sum=sum(dsig), n=length(dsig))%>%
  mutate(prop_sig=sum/n)%>%
  mutate(pnotsig=1-prop_sig)%>%
  select(-sum)%>%
  gather(value, sig, prop_sig:pnotsig)

chisq2<-test_trt3%>%
  mutate(sign=ifelse(dsig==0, "num_nonsig", "num_sig"))%>%
  group_by(trt_type3, sign)%>%
  summarise(sum=length(dsig))%>%
  spread(sign, sum, fill=0)

prop.test(x=as.matrix(chisq2[,c('num_sig', 'num_nonsig')]), alternative='two.sided')
#p <0.001

# tograph<-prop_sig_trt%>%
#   bind_rows(alldir)
# 
# ggplot(tograph, aes(x = trt_type2, y = sig, fill = value)) +
#   geom_col(width = 0.7) +
#   coord_flip() +
#   theme_minimal()+
#   scale_fill_brewer(name = "", labels = c("Not significant", "Significant")) +
#   scale_x_discrete(limits=c("Temperature",'Precip. Vari.', 'P', 'N','Mult. Nuts.', 'Irrigation','CO2',  'All Trts.'))+
#   labs(x = "Treatment", y = "Proportion of communities") +
#   theme(legend.position = "top")+
#   geom_hline(yintercept = 0.5)+
#   geom_vline(xintercept = 7.5, linetype="dashed")#+
#   # geom_text(x=1, y = 0.05, label="n = 7", size=4,check_overlap = TRUE)+
#   # geom_text(x=2, y = 0.06, label="n = 14", size=4,check_overlap = TRUE)+
#   # geom_text(x=3, y = 0.06, label="n = 52", size=4,check_overlap = TRUE)+
#   # geom_text(x=4, y = 0.06, label="n = 34", size=4,check_overlap = TRUE)+
#   # geom_text(x=5, y = 0.05, label="n = 9", size=4,check_overlap = TRUE)+
#   # geom_text(x=6, y = 0.05, label="n = 7", size=4,check_overlap = TRUE)+
#   # geom_text(x=7, y = 0.07, label="n = 218", size=4,check_overlap = TRUE)
# 
# 
# #Merge with exeriment info
# #Using rate_change_interval - recalculate slopes and get slope and p-value (0.5)
# #How often are we seeing directional change? OFTEN Does it differe treatment vs control? (plotmani=0) NO
# #Does it different by trt type YES. Higher for N + P or N+PK over N, CO2 and P.