library(tidyverse)
library(gridExtra)
library(reldist)
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

corredat<-read.csv("~/Documents/SpeciesRelativeAbundance_May2017.csv")%>%
  select(-X)%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))%>%
  filter(site_project_comm!="IMGERS_Yu_0"&site_project_comm!="Saskatchewan_CCD_0"&site_project_comm!="GVN_FACE_0")

corredat<-read.csv("C:\\Users\\megha\\Dropbox\\converge_diverge\\datasets\\LongForm\\SpeciesRelativeAbundance_May2017.csv")%>%
  select(-X)%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))%>%
  filter(site_project_comm!="IMGERS_Yu_0"&site_project_comm!="Saskatchewan_CCD_0"&site_project_comm!="GVN_FACE_0")

plotinfo<-read.csv("C:\\Users\\megha\\Dropbox\\converge_diverge\\datasets\\LongForm\\ExperimentInformation_May2017.csv")%>%
  select(site_code, project_name, community_type, calendar_year, treatment, plot_mani)

#problems
#IMGERS_Yu, plot 303 is there twice, it should be plot 304 where treatment=N5
#Sakatchewan, says Error in mapply(FUN = f, ..., SIMPLIFY = FALSE)
#zero-length inputs cannot be mixed with those of non-zero length 
#gvn face - only 2 years of data so will only have one point for the dataset.

##fill in the zeros
explist<-unique(corredat$site_project_comm)


corredat_sppool<-data.frame()

for (i in 1:length(explist)){
  ##get zero abundances to be filled in for all species.
  ##this works the first time only
  subset<-corredat%>%
    filter(site_project_comm==explist[i])%>%
    spread(genus_species, relcov, fill=0)
  
  ##make long and get averages of each species by treatment
  long<-subset%>%
    gather(genus_species, relcov, 10:ncol(subset))%>%
    group_by(site_project_comm, calendar_year, treatment, genus_species)%>%
    summarize(relcov=mean(relcov))
  
  corredat_sppool<-rbind(corredat_sppool, long)
}


###richness and evenness
diversity <- group_by(sp_pooled, site_project_comm, calendar_year,treatment_year, treatment) %>% 
  summarize(S=S(relcov),
            Even=E_q(relcov))%>%
  tbl_df()


###calculate species differences and reordering

#label the control versus treatment plots

corredat_treat_control<-merge(plotinfo, corredat,by=c("site_code","project_name","community_type","calendar_year",'treatment'))
  

reordering_ct=data.frame(site_project_comm=c(), treatment=c(), calendar_year=c(), MRSc_diff=c(), spdiffc=c())

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
    gather(genus_species, relcov, 11:ncol(subset))%>%
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
      
      metrics<-data.frame(site_project_comm=spc, treatment=treat_id, calendar_year=time_id, MRSc_diff=MRSc_diff, spdiffc=spdiffc)#spc_id
      ##calculate differences for these year comparison and rbind to what I want.
      
      reordering_ct=rbind(metrics, reordering_ct)  
    }
  }
}

