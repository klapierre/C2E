################################################################################
##  RAC Changes.R: This script creates the RAC differences and multivariate differences metrics for the codyn dataset 
##
##  Author: Meghan Avolio (meghan.avolio@jhu.edu)
##  Date: Oct 13, 2021
################################################################################

library(tidyverse)
library(gridExtra)
library(grid)
library(gtable)
library(devtools)
library(codyn)
library(vegan)
#home
setwd("~/Dropbox/")
#work
setwd("C:\\Users\\mavolio2\\Dropbox\\")

#Read in sp. rel. abund.
corredat<-read.csv("converge_diverge/datasets/LongForm/SpeciesRelativeAbundance_Nov2019.csv")

#gvn face - only 2 years of data so will only have one point for the dataset, therefore we are removing this dataset from these analyses.
corredat1<-corredat%>%
  select(-X)%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))%>%
  filter(site_project_comm!="GVN_FACE_0", site_project_comm!="AZI_NitPhos_0", site_project_comm!="JRN_study278_0", site_project_comm!="KNZ_GFP_4F", site_project_comm!="Saskatchewan_CCD_0", project_name!="e001", project_name!="e002",site_project_comm!="CHY_EDGE_0", site_project_comm!="HYS_EDGE_0", site_project_comm!="SGS_EDGE_0")

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

##remove one of 2 pre-treatment years in edge for CHY, SGS, and HAYS
edge<-corredat%>%
  select(-X)%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))%>%
  filter(site_code=="CHY"|site_code=="SGS"|site_code=="HYS"&project_name=="EDGE")%>%
  filter(calendar_year!=2012)


###final dataset to use
corredat_raw<-rbind(corredat1, azi, jrn, knz, sak, cdr, edge)

treatment_info<-read.csv("converge_diverge/datasets/LongForm/ExperimentInformation_March2019.csv")%>%
  dplyr::select(site_code, project_name, community_type, treatment,plot_mani)%>%
  unique()%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))

corredat <- corredat_raw %>%
 left_join(treatment_info, by=c( "site_code","project_name","community_type", "treatment","site_project_comm"))


########
noraresp<-corredat%>%
  filter(plot_mani==0)%>%
  group_by(site_project_comm, genus_species)%>%
  summarize(mrelcov=mean(relcov))%>%
  filter(mrelcov<0.01)%>%
  select(-mrelcov)

ggplot(data=subset(noraresp, mrelcov<0.1), aes(x=mrelcov))+
  geom_histogram()

corredat_norare <- corredat %>%
  right_join(noraresp)

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

write.csv(diff_rac, "C2E\\Products\\Testing Hypots\\CORRE_RAC_Diff_Metrics_Oct2021.csv", row.names = F)


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

write.csv(diff_rac_c, "C2E\\Products\\Testing Hypots\\CORRE_RAC_Diff_control_Metrics_Oct2021.csv", row.names = F)

####NO RARES
#####CALCULATING RAC differences
spc<-unique(corredat_norare$site_project_comm)
diff_rac_norare<-data.frame()

for (i in 1:length(spc)){
  subset<-corredat_norare%>%
    filter(site_project_comm==spc[i])
  
  ref_trt <- unique(subset(subset, plot_mani==0)$treatment)
  
  out<-RAC_difference(subset, time.var = 'calendar_year', species.var = "genus_species", abundance.var = 'relcov', replicate.var = 'plot_id', treatment.var = "treatment", reference.treatment = ref_trt)
  
  out$site_project_comm<-spc[i]
  
  diff_rac_norare<-rbind(diff_rac_norare, out)
}

write.csv(diff_rac_norare, "C2E\\Products\\Testing Hypots\\CORRE_RAC_Diff_Metrics_norares_Oct2021.csv", row.names = F)

##no RARE
#####CALCULATING RAC differences CONTROLS ONLY
corredat_control_nr<-corredat_norare%>%
  filter(plot_mani==0)%>%
  filter(site_project_comm!="Finse_WarmNut_0")

spc<-unique(corredat_norare$site_project_comm)

diff_rac_c_nr<-data.frame()

for (i in 1:length(spc)){
  subset<-corredat_control_nr%>%
    filter(site_project_comm==spc[i])
  
  out<-RAC_difference(subset, time.var = 'calendar_year', species.var = "genus_species", abundance.var = 'relcov', replicate.var = 'plot_id', treatment.var = "treatment")
  
  out$site_project_comm<-spc[i]
  
  diff_rac_c_nr<-rbind(diff_rac_c_nr, out)
}

write.csv(diff_rac_c_nr, "C2E\\Products\\Testing Hypots\\CORRE_RAC_Diff_control_Metrics__norares_Oct2021.csv", row.names = F)

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

write.csv(diff_mult2, "C2E\\Products\\Testing Hypots\\CORRE_Mult_diff_Metrics_Oct2021.csv", row.names = F)




#####CALCULATING multivariate differences no rares

spc<-unique(corredat_norare$site_project_comm)
diff_mult_norare<-data.frame()

for (i in 1:length(spc)){
  subset<-corredat_norare%>%
    filter(site_project_comm==spc[i])
  
  ref_trt <- unique(subset(subset, plot_mani==0)$treatment)
  
  out<-multivariate_difference(subset, time.var = 'calendar_year', species.var = "genus_species", abundance.var = 'relcov', replicate.var = 'plot_id', treatment.var = "treatment", reference.treatment = ref_trt)
  
  out$site_project_comm<-spc[i]
  
  diff_mult_norare<-rbind(diff_mult_norare, out)
}

control<-treatment_info%>%
  filter(plot_mani==0)%>%
  mutate(control=treatment)%>%
  select(site_project_comm, control)

diff_mult_norare2 <- diff_mult_norare%>%
  left_join(control)%>%
  mutate(trt_greater_disp=as.character(trt_greater_disp),
         control=as.character(control))%>%
  mutate(greater_disp=ifelse(trt_greater_disp == control, "C","T"))

write.csv(diff_mult_norare2, "C2E\\Products\\CommunityChange\\March2018 WG\\CORRE_Mult_diff_Metrics_norare_Jun2019.csv", row.names = F)

####species differences
spc<-unique(corredat$site_project_comm)
diff_abund<-data.frame()

for (i in 1:length(spc)){
  subset<-corredat%>%
    filter(site_project_comm==spc[i])
  
  ref_trt <- unique(subset(subset, plot_mani==0)$treatment)
  
  out<-abundance_difference(subset, time.var = 'calendar_year', species.var = "genus_species", abundance.var = 'relcov', replicate.var = 'plot_id', treatment.var = "treatment", reference.treatment = ref_trt)
  
  out$site_project_comm<-spc[i]
  
  diff_abund<-rbind(diff_abund, out)
}

write.csv(diff_abund, "C2E\\Products\\CommunityChange\\March2018 WG\\CORRE_Abund_Diff_Nov2021.csv", row.names = F)
