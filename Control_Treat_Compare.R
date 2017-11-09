library(tidyverse)
library(gridExtra)
library(grid)
library(gtable)
library(codyn)
library(vegan)
library(Kendall)

#function to calculate richness
#' @x the vector of abundances of each species
S<-function(x){
  x1<-x[x!=0]
  length(x1)
}

#function to calculate EQ evenness from Smith and Wilson 1996
#' @x the vector of abundances of each species
#' if all abundances are equal it returns a 1
E_q<-function(x){
  x1<-x[x!=0]
  if (length(x1)==1) {
    return(NA)
  }
  if (abs(max(x1) - min(x1)) < .Machine$double.eps^0.5) {##bad idea to test for zero, so this is basically doing the same thing testing for a very small number
    return(1)
  }
  r<-rank(x1, ties.method = "average")
  r_scale<-r/max(r)
  x_log<-log(x1)
  fit<-lm(r_scale~x_log)
  b<-fit$coefficients[[2]]
  2/pi*atan(b)
}

#read in the data FIX THE PATH LATER

corredat<-read.csv("~/Dropbox/converge_diverge/datasets/LongForm/SpeciesRelativeAbundance_Oct2017.csv")%>%
  select(-X)%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))

plotinfo<-read.csv("~/Dropbox/converge_diverge/datasets/LongForm/ExperimentInformation_May2017.csv")%>%
  select(site_code, project_name, community_type, calendar_year, treatment, plot_mani)%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))



corredat<-read.csv("C:\\Users\\megha\\Dropbox\\converge_diverge\\datasets\\LongForm\\SpeciesRelativeAbundance_Oct2017.csv")%>%
  select(-X)%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))%>%
  filter(site_code!="RIO")

plotinfo<-read.csv("C:\\Users\\megha\\Dropbox\\converge_diverge\\datasets\\LongForm\\ExperimentInformation_May2017.csv")%>%
  select(site_code, project_name, community_type, calendar_year, treatment, plot_mani)%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))

#problems
#Sakatchewan, says Error in mapply(FUN = f, ..., SIMPLIFY = FALSE)
#zero-length inputs cannot be mixed with those of non-zero length 
# RIO has NA in species_genus

##not sure why i had to do this for the diversity measures.
# ##fill in the zeros
# explist<-unique(corredat$site_project_comm)
# 
# 
# corredat_sppool<-data.frame()
# 
# for (i in 1:length(explist)){
#   ##get zero abundances to be filled in for all species.
#   ##this works the first time only
#   subset<-corredat%>%
#     filter(site_project_comm==explist[i])%>%
#     spread(genus_species, relcov, fill=0)
#   
#   ##make long and get averages of each species by treatment
#   long<-subset%>%
#     gather(genus_species, relcov, 10:ncol(subset))%>%
#     group_by(site_project_comm, calendar_year, treatment, treatment_year, genus_species)%>%
#     summarize(relcov=mean(relcov))%>%
#     ungroup
#   
#   corredat_sppool<-rbind(corredat_sppool, long)
# }



###calculate species differences and reordering and evenness and richness

#label the control versus treatment plots

corredat_treat_control<-merge(plotinfo, corredat,by=c("site_code","project_name","community_type","calendar_year",'treatment',"site_project_comm"))%>%
  mutate(expyear=paste(site_project_comm, calendar_year, sep="::"))
  

reordering_ct=data.frame(site_project_comm=c(), treatment=c(), calendar_year=c(), Sd=c(), Ed=c(), Rd=c(), spd=c())

explist<-unique(corredat_treat_control$site_project_comm)

for (i in 1:length(explist)){
  ##get zero abundances to be filled in for all species.
  ##this works the first time only
  subset<-corredat_treat_control%>%
    filter(site_project_comm==explist[i])%>%
    spread(genus_species, relcov, fill=0)
  
  spc<-explist[i]
  
  ##make long and get averages of each species by treatment
  long<-subset%>%
    gather(genus_species, relcov, 12:ncol(subset))%>%
    group_by(site_project_comm, calendar_year, treatment, plot_mani, genus_species)%>%
    summarize(relcov=mean(relcov))
  
  ##add ranks dropping zeros
  rank_pres<-long%>%
    filter(relcov!=0)%>%
    tbl_df()%>%
    group_by(site_project_comm, calendar_year, treatment, plot_mani)%>%
    mutate(rank=rank(-relcov, ties.method = "average"))%>%
    tbl_df()
  
  ###make zero abundant species have the rank S+1 (the size of the species pool plus 1)
  ##pull out zeros
  zeros<-long%>%
    filter(relcov==0)
  ##get species richness for each year
  rich<-group_by(long, site_project_comm, calendar_year, treatment, plot_mani)%>%
    summarize(S=S(relcov))
  ##merge together make zero abundances rank S+1
  zero_rank<-merge(zeros, rich, by=c("site_project_comm","calendar_year", "treatment", "plot_mani"))%>%
    mutate(rank=S+1)%>%
    select(-S)%>%
    tbl_df()
  ##combine all
  rank<-rbind(rank_pres, zero_rank)
  
  timestep<-sort(unique(rank$calendar_year)) 
  
  for(i in 1:(length(timestep))){
    
    time<-rank%>%
      filter(calendar_year==timestep[i])
    
    time_id<-timestep[i]
    
    #fitler out control plots
    control<-time%>%
      filter(plot_mani==0)
    
    treat_list<-unique(subset(time, plot_mani!=0)$treatment)
    
    for (i in 1:length(treat_list)){
      treat<-time%>%
        filter(treatment==treat_list[i])
      
      treat_id<-treat_list[i]
      
      subset_ct<-merge(control, treat, by=c("site_project_comm", "calendar_year","genus_species"), all=T)%>%
        filter(relcov.x!=0|relcov.y!=0)
      
      MRSc_diff<-mean(abs(subset_ct$rank.x-subset_ct$rank.y))/nrow(subset_ct)
      
      spdiff<-subset_ct%>%
        filter(relcov.x==0|relcov.y==0)
      
      spdiffc<-nrow(spdiff)/nrow(subset_ct)
      
      ##eveness richness
      s_c <- S(subset_ct$relcov.x)
      e_c <- E_q(subset_ct$relcov.x)
      s_t <- S(subset_ct$relcov.y)
      e_t <- E_q(subset_ct$relcov.y)
      
      sdiff<-abs(s_c-s_t)/nrow(subset_ct)
      ediff<-abs(e_c-e_t)/nrow(subset_ct)
      
      metrics<-data.frame(site_project_comm=spc, treatment=treat_id, calendar_year=time_id, Sd=sdiff, Ed=ediff, Rd=MRSc_diff, spd=spdiffc)#spc_id
      ##calculate differences for these year comparison and rbind to what I want.
      
      reordering_ct=rbind(metrics, reordering_ct)  
    }
  }
}

###mean change and dispersion

#####Calculating Bray-Curtis both comparing the mean community change between treatment and control plots in a time step

###first, get bray curtis dissimilarity values for each all years within each experiment between all combinations of plots
###second, get distance of each plot to its year centroid 
###third: mean_change is the distance the centroids of consequtive years
####fourth: dispersion_diff is the average dispersion of plots within a treatment to treatment centriod then compared between consequtive years

###list of all treats within an experiment

exp_year<-unique(corredat_treat_control$expyear)

#makes an empty dataframe
bray_curtis=data.frame() 
##calculating bray-curtis mean change and disperison differecnes
for(i in 1:length(exp_year)) {
  
  #subsets out each dataset
  subset<-corredat_treat_control%>%
    filter(expyear==exp_year[i])%>%
    select(site_project_comm, treatment, calendar_year, genus_species, relcov, plot_id, plot_mani)

  #need this to keep track of plot mani
  labels=subset%>%
    select(plot_mani, treatment)%>%
    unique()
  
  #transpose data
  species=subset%>%
    spread(genus_species, relcov, fill=0)
  
  #calculate bray-curtis dissimilarities
  bc=vegdist(species[,6:ncol(species)], method="bray")
  
  #calculate distances of each plot to treatment centroid (i.e., dispersion)
  disp=betadisper(bc, species$treatment, type="centroid")
  
  #getting distances between centroids over years; these centroids are in BC space, so that's why this uses euclidean distances
  cent_dist=as.data.frame(as.matrix(vegdist(disp$centroids, method="euclidean")))
  
  #extracting only the distances we need and adding labels for the comparisons;
  cent_C_T=data.frame(site_project_comm_year=exp_year[i],
                      treatment=row.names(cent_dist),
                      mean_change=t(cent_dist[names(cent_dist)==labels$treatment[labels$plot_mani==0],]))
  
  #renaming column
  colnames(cent_C_T)[3]<-"mean_change"
  
  #collecting and labeling distances to centroid from betadisper to get a measure of dispersion and then take the mean for a treatment
  disp2=data.frame(site_project_comm_year=exp_year[i],
                   treatment=species$treatment,
                   plot_mani=species$plot_mani,
                   plot_id=species$plot_id,
                   dist=disp$distances)%>%
    tbl_df%>%
    group_by(site_project_comm_year, treatment, plot_mani)%>%
    summarize(dispersion=mean(dist))
  
  control<-disp2$dispersion[disp2$plot_mani==0]
  
  ##subtract control from treatments
  disp_treat=disp2%>%
    mutate(disp_diff=dispersion-control)%>%
    select(-dispersion)
  
  #merge together change in mean and dispersion data
  distances<-merge(cent_C_T, disp_treat, by=c("site_project_comm_year","treatment"))
  
  #pasting dispersions into the dataframe made for this analysis
  bray_curtis=rbind(bray_curtis, distances)  
}

corre_braycurtis_control_treat<-bray_curtis%>%
  separate(site_project_comm_year, into=c("site_project_comm","calendar_year"), sep="::")%>%
  filter(plot_mani!=0)

all_Cont_Treat_Compare<-merge(reordering_ct, corre_braycurtis_control_treat,by=c("site_project_comm","calendar_year","treatment"))

write.csv(all_Cont_Treat_Compare, "C:\\Users\\megha\\Dropbox\\converge_diverge\\datasets\\LongForm\\CORRE_ContTreat_Compare_Nov2017.csv")
