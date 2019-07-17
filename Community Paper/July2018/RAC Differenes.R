################################################################################
##  RAC Changes.R: This script creates the RAC differences and multivariate differences metrics for the codyn dataset 
##
##  Author: Meghan Avolio (meghan.avolio@gmail.com)
##  Date: Oct 30, 2018
################################################################################

library(tidyverse)
library(gridExtra)
library(grid)
library(gtable)
library(codyn)
library(vegan)
#home
setwd("~/Dropbox/")
#work
setwd("C:\\Users\\megha\\Dropbox\\")

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

treatment_info<-read.csv("C2E\\Products\\CommunityChange\\March2018 WG\\data\\ExperimentInformation_Nov2017.csv")%>%
  dplyr::select(site_code, project_name, community_type, treatment,plot_mani)%>%
  unique()%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))

corredat <- corredat_raw %>%
 left_join(treatment_info, by=c( "site_code","project_name","community_type", "treatment","site_project_comm"))

#####CALCULATING RAC differences
spc<-unique(corredat$site_project_comm)
diff_rac<-data.frame()

for (i in 1:length(spc)){
  subset<-corredat%>%
    filter(site_project_comm==spc[i])
  
  ref_trt <- unique(subset(subset, plot_mani==0)$treatment)
  
  out<-RAC_difference(subset, time.var = 'calendar_year', species.var = "genus_species", abundance.var = 'relcov', replicate.var = 'plot_id', treatment.var = "treatment", reference.treatment = ref_trt)
  
  out$site_project_comm<-spc[i]
  
  diff_rac<-rbind(diff_rac, out)
}

write.csv(diff_rac, "C2E\\Products\\CommunityChange\\March2018 WG\\CORRE_RAC_Diff_Metrics_Oct2018.csv", row.names = F)


#####CALCULATING RAC differences CONTROLS ONLY
corredat_control<-corredat%>%
  filter(plot_mani==0)

diff_rac_c<-data.frame()

for (i in 1:length(spc)){
  subset<-corredat_control%>%
    filter(site_project_comm==spc[i])
  
  out<-RAC_difference(subset, time.var = 'calendar_year', species.var = "genus_species", abundance.var = 'relcov', replicate.var = 'plot_id', treatment.var = "treatment")
  
  out$site_project_comm<-spc[i]
  
  diff_rac_c<-rbind(diff_rac_c, out)
}

write.csv(diff_rac_c, "C2E\\Products\\CommunityChange\\March2018 WG\\CORRE_RAC_Diff_control_Metrics_Oct2018.csv", row.names = F)

#####CALCULATING multivariate differences
spc<-unique(corredat$site_project_comm)
diff_mult<-data.frame()

for (i in 1:length(spc)){
  subset<-corredat%>%
    filter(site_project_comm==spc[i])
  
  ref_trt <- unique(subset(subset, plot_mani==0)$treatment)
  
  out<-multivariate_difference(subset, time.var = 'calendar_year', species.var = "genus_species", abundance.var = 'relcov', replicate.var = 'plot_id', treatment.var = "treatment", reference.treatment = ref_trt)
  
  out$site_project_comm<-spc[i]
  
  diff_mult<-rbind(diff_mult, out)
}

control<-treatment_info%>%
  filter(plot_mani==0)%>%
  mutate(control=treatment)%>%
  select(site_project_comm, control)

diff_mult2 <- diff_mult%>%
  left_join(control)%>%
  mutate(trt_greater_disp=as.character(trt_greater_disp),
         control=as.character(control))%>%
  mutate(greater_disp=ifelse(trt_greater_disp == control, "C","T"))

write.csv(diff_mult2, "C2E\\Products\\CommunityChange\\March2018 WG\\CORRE_Mult_diff_Metrics_Oct2018.csv", row.names = F)