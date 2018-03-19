library(tidyverse)
library(gridExtra)
library(reldist)
library(grid)
library(gtable)
library(codyn)
library(vegan)
library(Kendall)


#Files from home

corredat<-read.csv("~/Dropbox/converge_diverge/datasets/LongForm/SpeciesRelativeAbundance_Oct2017.csv")%>%
  select(-X)%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))%>%
  filter(site_project_comm!="GVN_FACE_0")

corredat1<-read.csv("~/Dropbox/converge_diverge/datasets/LongForm/SpeciesRelativeAbundance_Oct2017.csv")%>%
  select(-X)%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))%>%
  filter(site_project_comm!="GVN_FACE_0", site_project_comm!="AZI_NitPhos_0", site_project_comm!="JRN_study278_0", site_project_comm!="KNZ_GFP_4F", site_project_comm!="Saskatchewan_CCD_0")

##several studies only have two measurments of a plot. I am dropping those plots
azi<-read.csv("~/Dropbox/converge_diverge/datasets/LongForm/SpeciesRelativeAbundance_Oct2017.csv")%>%
  select(-X)%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))%>%
  filter(site_code=="AZI")%>%
  filter(plot_id!=11&plot_id!=15&plot_id!=35&plot_id!=37)

jrn<-read.csv("~/Dropbox/converge_diverge/datasets/LongForm/SpeciesRelativeAbundance_Oct2017.csv")%>%
  select(-X)%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))%>%
  filter(site_project_comm=="JRN_study278_0")%>%
  filter(plot_id!=211&plot_id!=210)

knz<-read.csv("~/Dropbox/converge_diverge/datasets/LongForm/SpeciesRelativeAbundance_Oct2017.csv")%>%
  select(-X)%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))%>%
  filter(site_project_comm=="KNZ_GFP_4F")%>%
  filter(plot_id!="7_1_1"&plot_id!="7_2_1")

sak<-read.csv("~/Dropbox/converge_diverge/datasets/LongForm/SpeciesRelativeAbundance_Oct2017.csv")%>%
  select(-X)%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))%>%
  filter(site_project_comm=="Saskatchewan_CCD_0")%>%
  filter(plot_id!=2)

corredat<-rbind(corredat1, azi, jrn, knz, sak)


####files from work
corredat<-read.csv("C:\\Users\\megha\\Dropbox\\converge_diverge\\datasets\\LongForm\\SpeciesRelativeAbundance_Oct2017.csv")%>%
  select(-X)%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))%>%
  filter(site_project_comm!="GVN_FACE_0")

corredat1<-read.csv("C:\\Users\\megha\\Dropbox\\converge_diverge\\datasets\\LongForm\\SpeciesRelativeAbundance_Oct2017.csv")%>%
  select(-X)%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))%>%
  filter(site_project_comm!="GVN_FACE_0", site_project_comm!="AZI_NitPhos_0", site_project_comm!="JRN_study278_0", site_project_comm!="KNZ_GFP_4F", site_project_comm!="Saskatchewan_CCD_0")

##several studies only have two measurments of a plot. I am dropping those plots
azi<-read.csv("C:\\Users\\megha\\Dropbox\\converge_diverge\\datasets\\LongForm\\SpeciesRelativeAbundance_Oct2017.csv")%>%
  select(-X)%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))%>%
  filter(site_code=="AZI")%>%
  filter(plot_id!=11&plot_id!=15&plot_id!=35&plot_id!=37)

jrn<-read.csv("C:\\Users\\megha\\Dropbox\\converge_diverge\\datasets\\LongForm\\SpeciesRelativeAbundance_Oct2017.csv")%>%
  select(-X)%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))%>%
  filter(site_project_comm=="JRN_study278_0")%>%
  filter(plot_id!=211&plot_id!=210)

knz<-read.csv("C:\\Users\\megha\\Dropbox\\converge_diverge\\datasets\\LongForm\\SpeciesRelativeAbundance_Oct2017.csv")%>%
  select(-X)%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))%>%
  filter(site_project_comm=="KNZ_GFP_4F")%>%
  filter(plot_id!="7_1_1"&plot_id!="7_2_1")

sak<-read.csv("C:\\Users\\megha\\Dropbox\\converge_diverge\\datasets\\LongForm\\SpeciesRelativeAbundance_Oct2017.csv")%>%
  select(-X)%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))%>%
  filter(site_project_comm=="Saskatchewan_CCD_0")%>%
  filter(plot_id!=2)

corredat<-rbind(corredat1, azi, jrn, knz, sak)%>%
  filter(!is.na(genus_species))##get rid of problem with RIO

#problems
#gvn face - only 2 years of data so will only have one point for the dataset.


plotinfo<-corredat%>%
  select(site_project_comm, calendar_year, plot_id, treatment, treatment_year)%>%
  unique()


#######Doing difference

#problems
#Sakatchewan, says Error in mapply(FUN = f, ..., SIMPLIFY = FALSE)
#zero-length inputs cannot be mixed with those of non-zero length 

#####CALCULATING RAC differences with block

blocked<-corredat%>%
  filter(block!=0)%>%
  filter(site_project_comm!="ARC_MNT_0"&site_project_comm!="BAY_LIND_0"&site_project_comm!="dcgs_gap_0"&site_project_comm!="JRN_study278_0"&site_project_comm!="KLU_KGFert_0"&site_project_comm!="KLU_BFFert_0"&site_project_comm!=""&site_project_comm!="LATNJA_CLIP_Heath"&site_project_comm!="LATNJA_CLIP_Meadow"&site_project_comm!="NWT_bowman_DryBowman"&site_project_comm!="NWT_bowman_WetBowman"&site_project_comm!="NWT_snow_0"&site_project_comm!="TRA_Lovegrass_0")

spc<-unique(blocked$site_project_comm)
diff_rac_block<-data.frame()

for (i in 1:length(spc)){
  subset<-blocked%>%
    filter(site_project_comm==spc[i])
  
  out<-RAC_difference(subset, time.var = 'calendar_year', species.var = "genus_species", abundance.var = 'relcov', replicate.var = 'plot_id', block.var = 'block', treatment.var = 'treatment')
  out$site_project_comm<-spc[i]
  
  diff_rac_block<-rbind(diff_rac_block, out)
}

#####CALCULATING RAC differences without blocks pooling up to treatment

trt_control<-corredat%>%
  filter(block==0)

spc<-unique(trt_control$site_project_comm)
diff_rac_ct<-data.frame()

for (i in 1:length(spc)){
  subset<-trt_control%>%
    filter(site_project_comm==spc[i])
  
  out<-RAC_difference(subset, time.var = 'calendar_year', species.var = "genus_species", abundance.var = 'relcov', replicate.var = 'plot_id', treatment.var = 'treatment', pool = "YES")
  out$site_project_comm<-spc[i]
  
  diff_rac_ct<-rbind(diff_rac_ct, out)
}

#####CALCULATING abundance differences with block
spc<-unique(blocked$site_project_comm)
diff_abund_block<-data.frame()

for (i in 1:length(spc)){
  subset<-blocked%>%
    filter(site_project_comm==spc[i])
  
  out<-abundance_difference(subset, time.var = 'calendar_year', species.var = "genus_species", abundance.var = 'relcov', replicate.var = 'plot_id', block.var = 'block', treatment.var = 'treatment')
  out$site_project_comm<-spc[i]
  
  diff_abund_block<-rbind(diff_abund_block, out)
}
#####CALCULATING abundance differences without blocks pooling up to treatment

trt_control<-corredat%>%
  filter(block==0)

spc<-unique(trt_control$site_project_comm)
diff_abund_ct<-data.frame()

for (i in 1:length(spc)){
  subset<-trt_control%>%
    filter(site_project_comm==spc[i])
  
  out<-abundance_difference(subset, time.var = 'calendar_year', species.var = "genus_species", abundance.var = 'relcov', replicate.var = 'plot_id', treatment.var = 'treatment', pool = "YES")
  out$site_project_comm<-spc[i]
  
  diff_abund_ct<-rbind(diff_abund_ct, out)
}
#####CALCULATING RAC differences without blocks pooling up to treatment for all datasets
spc<-unique(corredat$site_project_comm)
diff_rac_all<-data.frame()

for (i in 1:length(spc)){
  subset<-corredat%>%
    filter(site_project_comm==spc[i])
  
  out<-RAC_difference(subset, time.var = 'calendar_year', species.var = "genus_species", abundance.var = 'relcov', replicate.var = 'plot_id', treatment.var = 'treatment', pool = "YES")
  out$site_project_comm<-spc[i]
  
  diff_rac_all<-rbind(diff_rac_all, out)
}


##calculating multivariate differences

spc<-unique(corredat$site_project_comm)
diff_mult<-data.frame()

for (i in 1:length(spc)){
  subset<-corredat%>%
    filter(site_project_comm==spc[i])
  
  out<-multivariate_difference(subset, time.var = 'calendar_year', species.var = "genus_species", abundance.var = 'relcov', replicate.var = 'plot_id', treatment='treatment')
  out$site_project_comm<-spc[i]
  
  diff_mult<-rbind(diff_mult, out)
}

#####CALCULATING curve differences with block
spc<-unique(blocked$site_project_comm)
diff_curve_block<-data.frame()

for (i in 1:length(spc)){
  subset<-blocked%>%
    filter(site_project_comm==spc[i])
  
  out<-curve_difference(subset, time.var = 'calendar_year', species.var = "genus_species", abundance.var = 'relcov', replicate.var = 'plot_id', block.var = 'block', treatment.var = 'treatment')
  out$site_project_comm<-spc[i]
  
  diff_curve_block<-rbind(diff_curve_block, out)
}

#####CALCULATING curve differences without blocks pooling up to treatment
spc<-unique(corredat$site_project_comm)
diff_curve_ct<-data.frame()

for (i in 1:length(spc)){
  subset<-corredat%>%
    filter(site_project_comm==spc[i])
  
  out<-curve_difference(subset, time.var = 'calendar_year', species.var = "genus_species", abundance.var = 'relcov', replicate.var = 'plot_id', treatment.var = 'treatment', pool = "YES")
  out$site_project_comm<-spc[i]
  
  diff_curve_ct<-rbind(diff_curve_ct, out)
}


# merge1<-merge(div_diff, reordering_ct, by=c("site_project_comm","calendar_year","treatment"))
# all_Cont_Treat_Compare<-merge(merge1, corre_braycurtis_control_treat,by=c("site_project_comm","calendar_year","treatment"))
# 
# write.csv(all_Cont_Treat_Compare, "~/Dropbox/converge_diverge/datasets/LongForm/CORRE_ContTreat_Compare_OCT2017.csv")
# 


write.csv(all_Cont_Treat_Compare, "C:\\Users\\megha\\Dropbox\\converge_diverge\\datasets\\LongForm\\CORRE_ContTreat_Compare_Nov2017.csv")
