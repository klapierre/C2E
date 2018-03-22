library(devtools)
install_github("mavolio/codyn", ref = "RACs_cleaner")
library(tidyverse)
library(gridExtra)
library(grid)
library(gtable)
library(codyn)
library(vegan)


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

#problems
#gvn face - only 2 years of data so will only have one point for the dataset.


plotinfo<-corredat%>%
  select(site_project_comm, calendar_year, plot_id, treatment, treatment_year)%>%
  unique()


#####CALCULATING DIVERSITY METRICs

spc<-unique(corredat$site_project_comm)
div_eq<-data.frame()

for (i in 1:length(spc)){
  subset<-corredat%>%
    filter(site_project_comm==spc[i])
  
  out<-community_structure(subset, time.var = 'calendar_year', abundance.var = 'relcov', replicate.var = 'plot_id')
  out$site_project_comm<-spc[i]
  
  div_eq<-rbind(div_eq, out)
}

#####CALCULATING RAC changes
spc<-unique(corredat$site_project_comm)
delta_rac<-data.frame()

for (i in 1:length(spc)){
  subset<-corredat%>%
    filter(site_project_comm==spc[i])
  
  out<-RAC_change(subset, time.var = 'calendar_year', species.var = "genus_species", abundance.var = 'relcov', replicate.var = 'plot_id')
  out$site_project_comm<-spc[i]
  
  delta_rac<-rbind(delta_rac, out)
}

write.csv(delta_rac, "C:\\Users\\megha\\Dropbox\\converge_diverge\\datasets\\LongForm\\CORRE_RAC_Metrics_Feb2018_allReplicates.csv")


write.csv(delta_rac, "~/Dropbox/C2E/Products/CommunityChange/March2018 WG/CORRE_RAC_Metrics_March2018.csv")

####based on treatment_year
spc<-unique(corredat$site_project_comm)
delta_rac<-data.frame()

for (i in 1:length(spc)){
  subset<-corredat%>%
    filter(site_project_comm==spc[i])
  
  out<-RAC_change(subset, time.var = 'treatment_year', species.var = "genus_species", abundance.var = 'relcov', replicate.var = 'plot_id')
  out$site_project_comm<-spc[i]
  
  delta_rac<-rbind(delta_rac, out)
}
write.csv(delta_rac, "~/Dropbox/C2E/Products/CommunityChange/March2018 WG/CORRE_RAC_Metrics_March2018_trtyr.csv")

###calculating multivariate changes
spc<-unique(corredat$site_project_comm)
delta_mult<-data.frame()

for (i in 1:length(spc)){
  subset<-corredat%>%
    filter(site_project_comm==spc[i])
  
  out<-centroid_change(subset, time.var = 'treatment_year', species.var = "genus_species", abundance.var = 'relcov', replicate.var = 'plot_id', treatment.var = "treatment")
  out$site_project_comm<-spc[i]
  
  delta_mult<-rbind(delta_mult, out)
}

write.csv(delta_mult, "~/Dropbox/C2E/Products/CommunityChange/March2018 WG/CORRE_Mult_Metrics_March2018.csv")

# ##getting the average for each treatment in a year
# 
# corre_diversity<-merge(plotinfo, diversity, by=c("site_project_comm","calendar_year","plot_id"))%>%
#   group_by(site_project_comm, calendar_year, treatment_year, treatment)%>%
#   summarize(S_diff=mean(S_diff), even_diff=mean(E_diff, na.rm=T))
# 
# corre_gainloss<-merge(plotinfo, gain_loss, by=c("site_project_comm","calendar_year","plot_id"))%>%
#   group_by(site_project_comm, calendar_year, treatment_year, treatment)%>%
#   summarize(gain=mean(appearance), loss=mean(disappearance))
# 
# corre_reordering<-merge(plotinfo, reordering, by=c("site_project_comm","calendar_year","plot_id"))%>%
#   group_by(site_project_comm, calendar_year, treatment_year, treatment)%>%
#   summarize(MRSc=mean(MRSc))
# 
# ####MERGING TO A SINGE DATASET and exporting
# 
# merge1<-merge(corre_diversity, corre_gainloss, by=c("site_project_comm","calendar_year","treatment_year","treatment"), all=T)
# merge2<-merge(merge1, corre_reordering, by=c("site_project_comm","calendar_year","treatment_year","treatment"), all=T)
# all_metrics<-merge(merge2, corre_braycurtis, by=c("site_project_comm","calendar_year","treatment"), all=T)
# 
# write.csv(all_metrics, "C:\\Users\\megha\\Dropbox\\converge_diverge\\datasets\\LongForm\\CORRE_RAC_Metrics_Oct2017_allyears_2.csv")
# 
# 
# write.csv(all_metrics, "~/Dropbox/converge_diverge/datasets/LongForm/CORRE_RAC_Metrics_Oct2017_allyears_2.csv")


# ###Getting b-C distnace of each plot to itself comparing t1 to t2.
# 
# corredat$expplot<-paste(corredat$site_project_comm, corredat$plot_id, sep="::")
# 
# exp_plot_list<-unique(corredat$expplot)
# 
# 
# #makes an empty dataframe
# bray_curtis_dissim=data.frame(site_project_comm_plot=c(), calendar_year=c(), bc_dissim=c()) 
# 
# ##calculating bray-curtis mean change and disperison differecnes
# for(i in 1:length(exp_plot_list)) {
#   
#   #subsets out each dataset
#   subset=corredat%>%
#     filter(expplot==exp_plot_list[i])%>%
#     select(site_project_comm, treatment, calendar_year, genus_species, relcov, plot_id)
#   
#   #get years
#   experiment_years<-sort(unique(subset$calendar_year))
#   
#   #transpose data
#   species=subset%>%
#     spread(genus_species, relcov, fill=0)
#   
#   #calculate bray-curtis dissimilarities
#   bc=as.matrix(vegdist(species[,5:ncol(species)], method="bray"))
#   
#   ###experiment_year is year x+1
#   bc_dis=data.frame(site_project_comm_plot=exp_plot_list[i],
#                     calendar_year=experiment_years[2:length(experiment_years)],
#                     bc_dissim=diag(bc[2:nrow(bc),1:(ncol(bc)-1)]))
#   
#   #pasting dispersions into the dataframe made for this analysis
#   bray_curtis_dissim=rbind(bc_dis, bray_curtis_dissim)  
# }
# 
# corre_braycurtis<-bray_curtis_dissim%>%
#   separate(site_project_comm_plot, into=c("site_project_comm","plot_id"), sep="::")
# 
# ###merging to a single dataset and adding treatment information
# merge1<-merge(gain_loss, diversity, by=c("site_project_comm","calendar_year","plot_id"))
# merge2<-merge(merge1, reordering,by=c("site_project_comm","calendar_year","plot_id")) 
# merge3<-merge(merge2, corre_braycurtis, by=c("site_project_comm","calendar_year","plot_id"))
# corre_all<-merge(plotinfo, merge3, by=c("site_project_comm","calendar_year","plot_id"))
# 
# write.csv(corre_all, "C:\\Users\\megha\\Dropbox\\converge_diverge\\datasets\\LongForm\\CORRE_RAC_Metrics_Oct2017_allReplicates.csv")
# 
# write.csv(corre_all, "~/Dropbox/converge_diverge/datasets/LongForm/CORRE_RAC_Metrics_Oct2017_allReplicates_2.csv")


write.csv(corre_all, "C:\\Users\\megha\\Dropbox\\converge_diverge\\datasets\\LongForm\\CORRE_RAC_Metrics_NOV2017_allReplicates.csv")

