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
  filter(site_project_comm!="IMGERS_Yu_0"&site_project_comm!="Saskatchewan_CCD_0")


#problems
#IMGERS_Yu, plot 303 is there twice, it should be plot 304 where treatment=N5
#Sakatchewan, says Error in mapply(FUN = f, ..., SIMPLIFY = FALSE) : 
#zero-length inputs cannot be mixed with those of non-zero length 

plotinfo<-corredat%>%
  select(site_project_comm, calendar_year, plot_id, treatment, treatment_year)%>%
  unique()


#####CALCULATING DIVERSITY METRICS WITHIN A TIME STEP FOR EACH REPLICATE 

##need to get this working with NAs for mean calculations
diversity <- group_by(corredat, site_project_comm, calendar_year, plot_id) %>% 
  summarize(S=S(relcov),
            E_q=E_q(relcov))%>%
            tbl_df()

#####CALCULATING DIVERSITY METRICS ACROSS CONSECUTIVE TIME STEPS FOR EACH REPLICATE
explist<-unique(corredat$site_project_comm)

gain_loss<-data.frame()

for (i in 1:length(explist)){
  subset<-corredat%>%
    filter(site_project_comm==explist[i])

loss<-turnover(df=subset, time.var="calendar_year", species.var="genus_species", abundance.var="relcov", replicate.var="plot_id", metric="disappearance")
gain<-turnover(df=subset, time.var="calendar_year", species.var="genus_species", abundance.var="relcov", replicate.var="plot_id", metric="appearance")

gain$site_project_comm<-explist[i]
loss$site_project_comm<-explist[i]

gl<-merge(gain, loss, by=c("site_project_comm","plot_id","calendar_year"))

gain_loss<-rbind(gain_loss, gl)
}

####New appraoch to Rank Shifts
###ranks - taking into account that all speices are not always present.
###Give all species with zerio abundace the S+1 rank for that year.
###includes species that are not present in year X but appear in year X+1 or are present in year X and disappear in year X+1

##add ranks 
ranks<-corredat%>%
  filter(relcov!=0)%>%
  tbl_df()%>%
  group_by(site_project_comm, calendar_year, plot_id)%>%
  mutate(rank=rank(-relcov, ties.method = "average"))%>%
  tbl_df()

#adding zeros
addzero<-data.frame()

for (i in 1:length(explist)){
  subset<-corredat%>%
    filter(site_project_comm==explist[i])

  wide<-subset%>%
    spread(key=genus_species, value=relcov, fill=0)
  
  long<-wide%>%
    gather(key=genus_species, value=relcov, 10:ncol(wide))
  
  addzero<-rbind(addzero, long)
}

###make zero abundant species have the rank S+1 (the size of the species pool plus 1)
##pull out zeros
corre_zeros<-addzero%>%
  filter(relcov==0)
##get species richness for each year
corre_S<-group_by(corredat, site_project_comm, calendar_year, plot_id)%>%
  summarize(S=S(relcov))
##merge together make zero abundances rank S+1
corre_zero_rank<-merge(corre_zeros, corre_S, by=c("site_project_comm","calendar_year","plot_id"))%>%
  mutate(rank=S+1)%>%
  select(-S)%>%
  tbl_df()
##combine all
corre_rank<-rbind(ranks, corre_zero_rank)%>%
  mutate(exp_plot=paste(site_project_comm, plot_id, sep="::"))

##calculate reordering between time steps by mean ranks shifts corrected for the size of the speceis pool

reordering=data.frame(id=c(), calendar_year=c(), MRSc=c())#expeiment year is year of timestep2

exp_plot_list<-unique(corre_rank$exp_plot)

for (i in 1:length(exp_plot_list)){
  subset<-corre_rank%>%
    filter(exp_plot==exp_plot_list[i])
  id<-exp_plot_list[i]
  
  splist<-subset%>%
    select(genus_species)%>%
    unique()
  sppool<-length(splist$genus_species)
  
  #now get all timestep within an experiment
  timestep<-sort(unique(subset$calendar_year))    
  
  for(i in 1:(length(timestep)-1)) {#minus 1 will keep me in year bounds NOT WORKING
    subset_t1<-subset%>%
      filter(calendar_year==timestep[i])
    
    subset_t2<-subset%>%
      filter(calendar_year==timestep[i+1])
    
    subset_t12<-merge(subset_t1, subset_t2, by=c("genus_species","site_project_comm","site_code","project_name","community_type"), all=T)%>%
      filter(relcov.x!=0|relcov.y!=0)
    
    MRSc<-mean(abs(subset_t12$rank.x-subset_t12$rank.y))/nrow(subset_t12)
    
    metrics<-data.frame(id=id, experiment_year=timestep[i+1], MRSc=MRSc)#spc_id
    ##calculate differences for these year comparison and rbind to what I want.
    
    reordering=rbind(metrics, reordering)  
  }
}

corre_reordering<-reordering%>%
  separate(id, c("site_project_comm","plot_id"), sep="::")

#####Calculating Bray-Curtis both comparing the mean community change between consequtive time steps and the change in dispersion between two time steps.
##Doing this for all years of an experiment at one time point because want to ensure all points are in the same space.
###first, get bray curtis dissimilarity values for each all years within each experiment between all combinations of plots
###second, get distance of each plot to its year centroid 
###third: mean_change is the distance the centroids of consequtive years
####fourth: dispersion_diff is the average dispersion of plots within a treatment to treatment centriod then compared between consequtive years

#makes an empty dataframe
bray_curtis=data.frame(site_project_comm=c(), calendar_year=c(), bc_mean_change=c(), bc_dispersion_diff=c()) 
##calculating bray-curtis mean change and disperison differecnes
for(i in 1:length(explist)) {
  
  #subsets out each dataset
  subset=corredat%>%
    filter(site_project_comm==explist[i])%>%
    select(site_project_comm, calendar_year, genus_species, relcov, plot_id)
  
  #get years
  experiment_years<-sort(unique(subset$calendar_year))
  
  #transpose data
  species=subset%>%
    spread(genus_species, relcov, fill=0)
  
  #calculate bray-curtis dissimilarities
  bc=vegdist(species[,4:ncol(species)], method="bray")
  
  #calculate distances of each plot to year centroid (i.e., dispersion)
  disp=betadisper(bc, species$calendar_year, type="centroid")
  
  #getting distances between centroids over years; these centroids are in BC space, so that's why this uses euclidean distances
  cent_dist=as.matrix(vegdist(disp$centroids, method="euclidean"))
  
  ##extracting only the comparisions we want year x to year x=1.
  ###(experiment_year is year x+1
  cent_dist_yrs=data.frame(site_project_comm=explist[i],
                           calendar_year=experiment_years[2:length(experiment_years)],
                           mean_change=diag(cent_dist[2:nrow(cent_dist),1:(ncol(cent_dist)-1)]))
  
  #collecting and labeling distances to centroid from betadisper to get a measure of dispersion and then take the mean for a year
  disp2=data.frame(site_project_comm=site_project_comm_u[i],
                   experiment_year=species$experiment_year,
                   plot_id=species$plot_id,
                   dist=disp$distances)%>%
    tbl_df%>%
    group_by(site_project_comm, experiment_year)%>%
    summarize(dispersion=mean(dist))
  
  ##subtract consequtive years subtracts year x+1 - x. So if it is positive there was greater dispersion in year x+1 and if negative less dispersion in year x+1
  disp_yrs=data.frame(site_project_comm=site_project_comm_u[i],
                      experiment_year=experiment_years[2:length(experiment_years)],
                      dispersion_diff=diff(disp2$dispersion))
  
  #merge together change in mean and dispersion data
  distances<-merge(cent_dist_yrs, disp_yrs, by=c("site_project_comm","experiment_year"))
  
  #pasting dispersions into the dataframe made for this analysis
  bray_curtis=rbind(distances, bray_curtis)  
}

codyndat_braycurtis<-bray_curtis


####MERGING TO A SINGE DATASET
merge1<-merge(codyndat_diversity, codyndat_gains_loss, by=c("site_project_comm","experiment_year"))
merge2<-merge(merge1, codyndat_reorder, by=c("site_project_comm","experiment_year"))
merge3<-merge(merge2, codyndat_braycurtis, by=c("site_project_comm","experiment_year"))
merge4<-merge(merge3, codyndat_dstar, by=c("site_project_comm","experiment_year"))
codyndat_allmetrics<-merge(merge4, codyndat_info, by="site_project_comm")

