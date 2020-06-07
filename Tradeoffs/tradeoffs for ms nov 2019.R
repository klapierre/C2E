#to do: 

# TRY TO FIGURE OUT WHO THESE WINNERS ARE:
# subset to only the species that are common in at least 5? site_proj_comm x treatment combinations, do correlation matrix to see if they are winners everywhere or only sometimes. but species names differ across experiments. 
# subset to only winners, see if any occur across multiple sites? 
# subset to common species only (>1% abund) and look at those winners (to avoid the problem that some winners are super rare and thus irrelevant)
# or just look at average responses (across all treatments) for each species. but species names differ across experiments.

# pick out a few winners haphazardly and discuss them. 
# systematically can we find big winners (and/or losers) across all studies/treatment combinations?


setwd("/Volumes/GoogleDrive/My Drive/RESEARCH/C2E tradeoffs")
setwd("C:\\Users\\wilco\\Dropbox\\shared working groups\\C2E\\Tradeoffs- Adam\\") ## KW HP laptop

###
### Set up workspace
###

## Load ibraries
library(tidyverse)
library(reshape2)
library(ggpubr)
library(lme4)
library(grid)

## Graphing parameters
theme_set(theme_bw()) 
theme_eg=theme_update(panel.grid.minor=element_blank(), panel.grid.major=element_blank(), strip.background=element_rect(color="white", fill="white"))

### 
### Read in and clean data
###

abun_long <- read.csv("SpeciesRawAbundance_March2019.csv") %>%
  dplyr::select(-X) %>% # remove excess columns
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_")) # create site_proj_comm column

exp_info <- read.csv("ExperimentInformationTradeoff.csv")
projlist=unique(abun_long$site_project_comm)

###
### Loop through experiments and add zeros for species found in each experiment, but not in a particular plot
###
abun_long_with_zeros <- {}

for(PROJECT in 1:length(projlist)){
  abun_long_temp <- filter(abun_long, site_project_comm == projlist[PROJECT]) %>% # Subsets to PROJECT
    mutate(genus_species=factor(genus_species)) %>% # stops R from thinking that genus_species levels include all species within CoRRE
    spread(key=genus_species, value=abundance) %>% # similar to melt function
    replace(is.na(.), 0) %>% # Replace NAs with zeros
    gather(key=genus_species, value=abundance, -site_code:-site_project_comm) # convert to long form
  abun_long_with_zeros <- rbind(abun_long_with_zeros, abun_long_temp) # combine with master data frame
  rm(abun_long_temp)
}

abun_long_with_zeros <- abun_long_with_zeros %>%
  rename(species=genus_species)

## Write file 
write.csv(abun_long_with_zeros, paste0("SpeciesRawAbundance_WithZeros_", Sys.Date(),".csv"), row.names=F)

###
### Merging with adam's hand coded file that finds single factor treatments of multifactor studies and gives them a nice label
###

exp_select <- filter(exp_info, pulse==0 & !expts$TradeoffTrt %in% c("Combined treatments", "X", "herb_removal"))#removes pulse experiments, weird treatments (plant manipulations, precipitation variability, etc), multi factor treatments of multifactor studies (eg N+P addition), and herbivore removal experiments (not a GCD))

abun_long_merged <- abun_long_with_zeros %>%
  full_join(exp_select, by=c("site_code", "project_name", "calendar_year", "treatment_year", "treatment", "community_type")) %>%
  filter(treatment_year > 0) # Drop treatment year zero
  
useadamlong.1m=melt(abun_long_merged, id=c("site_project_comm", "calendar_year", "species"), measure="TradeoffTrt")

###
### Which sites/experiments have which treatments?
###

trtsinsite=dcast(useadamlong.1m, site_project_comm~value, fun=length)

#subsetting to sites with at least 3 plot types (control + 2 trts)
trtsinsite$usesite=rowSums(trtsinsite[,c(3:length(trtsinsite)),]>0)
usethesesites=data.frame(site_project_comm=trtsinsite$site_project_comm[trtsinsite$usesite>2])
useadamlong.1s=merge(useadamlong.1, usethesesites, by="site_project_comm")
#write.csv(useadamlong.1s, "files for adam/only sites with at least 3 treatments may2020.csv", row.names=F)

#average across replicate plots within a year
useadamlong.1s.m=melt(useadamlong.1s, id=c("site_project_comm", "calendar_year", "treatment_year", "plot_id", "species", "TradeoffTrt"), measure="abundance")
useadam.byyear=dcast(useadamlong.1s.m, site_project_comm + calendar_year + treatment_year + species~TradeoffTrt, fun=mean)

#average over time for a species 
useadam.byyear.m=melt(useadam.byyear, id=c("site_project_comm", "calendar_year", "treatment_year", "species"))
useadam.bysp=dcast(useadam.byyear.m, site_project_comm + species ~ variable, fun=mean, na.rm=T)


#----------GENERATING NULL MODEL


#permuting the community to get null distribution. thank you Kevin!

### Read in data and choose relavent columns -- PART OF THIS IS UNNECCESARY, WOULD BE NICE TO CLEAN THIS UP A BIT
full_df <- useadamlong.1s[,c("site_project_comm", "treatment_year", "TradeoffTrt", "plot_id", "species", "abundance")] %>%
  mutate(site_proj_comm_yr_sp=as.factor(paste(full_df$site_project_comm, full_df$treatment_year, full_df$species, sep="::")))

### Create data frame with number of plots for each site_proj_comm
group_length <- full_df %>%
  dplyr::select(site_project_comm, plot_id) %>%
  unique() %>%
  group_by(site_project_comm) %>%
  mutate(num_of_plots = length(levels(factor(plot_id)))) %>%
  dplyr::select(-plot_id) %>%
  unique()

### Create treatment data frame
treatment_info <- useadamlong.1s %>%
  dplyr::select(site_project_comm, treatment, plot_id, nutrients:other_manipulation, trt_details:TradeoffTrt) %>%
  unique() %>%
  arrange(site_project_comm)

real_plot_id_and_treatment <- treatment_info %>%
  dplyr::select(site_project_comm, treatment, plot_id) %>%
  rename(actual_treatment=treatment)

treatment_info <- treatment_info %>% select(-plot_id) ## This data frame will be combined with randomized vector, then sorted, and recombined with plot id to randomize treatments for plot id's

### Loop that creates n randomized dataframes
averaged_df_master <- {}

for(perm in 1:999){

  ### Create a random vector that randomly selects numbers from 1 to the number of plots in each site_project_comm
  rand_vector_master <- {}
  for(i in 1:nrow(group_length)){
    rand_vector_temp <- data.frame(site_project_comm = as.character(group_length$site_project_comm[i]), 
                                   rand_plot_id = sample(1:group_length$num_of_plots[i], group_length$num_of_plots[i]))
    rand_vector_master <- rbind(rand_vector_master, rand_vector_temp)
    rm(rand_vector_temp)
  }

  ### Append random vector to treatment information, sort to randomize, add back to original plot ids
  randomized_trt_info <- treatment_info %>%
    bind_cols(
      rand_vector_master %>% arrange(site_project_comm)) %>%
    arrange(site_project_comm, rand_plot_id) %>%
    bind_cols(
      real_plot_id_and_treatment
    ) %>%
    rename(randomized_treatment = treatment)

  ### Combine randomized treatment information with species abundance (based on actual plot id)
  randomized_trt_with_abundances <- full_df %>%
    dplyr::select(-TradeoffTrt, -site_proj_comm_yr_sp) %>%
    full_join(
      randomized_trt_info %>% dplyr::select(site_project_comm:TradeoffTrt, plot_id),
      by=c("site_project_comm", "plot_id")
    ) %>%
    arrange(site_project_comm, plot_id, treatment_year) %>%
    mutate(permutation = perm)
  
  averaged_df <- randomized_trt_with_abundances %>%
    group_by(permutation, site_project_comm, TradeoffTrt, species) %>%
    summarise(abundance = mean(abundance))

  averaged_df_master <- rbind(averaged_df_master, averaged_df)
  rm(averaged_df, randomized_trt_with_abundances, randomized_trt_info, rand_vector_master)
}

write.csv(averaged_df_master, paste0("999 permutations",Sys.Date(),".csv"), row.names=F)

#why do we have prmutations for these sites when we don't have real data??
#KLU_BFFert_0, NIN_herbdiv_0, SR_Water_A, SR_Water_P, TRA_Lovegrass_0


#--------MERGE REAL AND PERMUTED DATA 

#permuted=read.csv("999 permutations july 2019.csv")
#permuted=as.data.frame(averaged_df_master)

permuted.m=melt(permuted, id=c("permutation", "site_project_comm", "TradeoffTrt", "species"))
permuted.c=dcast(permuted.m, permutation+site_project_comm+species~TradeoffTrt)
permuted.c$herb_removal=NULL

real.c=useadam.bysp
real.c$permutation="real"

fulldat=rbind(real.c, permuted.c)


#----------WHICH TRADEOFFS DO WE HAVE THE DATA TO TEST?


trtsinsite.s=trtsinsite[trtsinsite$usesite>2,]
trtsinsite.s$usesite=NULL
trtsinsite.s$control=NULL

trtsinsite.s$BMCTxCO2=ifelse(trtsinsite.s$burnmowcliptill>0 & trtsinsite.s$CO2>0, 1, 0)
trtsinsite.s$BMCTxdrought=ifelse(trtsinsite.s$burnmowcliptill>0 & trtsinsite.s$drought>0, 1, 0)
trtsinsite.s$BMCTxirr=ifelse(trtsinsite.s$burnmowcliptill>0 & trtsinsite.s$irr>0, 1, 0)
trtsinsite.s$BMCTxmult_nutrient=ifelse(trtsinsite.s$burnmowcliptill>0 & trtsinsite.s$mult_nutrient>0, 1, 0)
trtsinsite.s$BMCTxN=ifelse(trtsinsite.s$burnmowcliptill>0 & trtsinsite.s$N>0, 1, 0)
trtsinsite.s$BMCTxP=ifelse(trtsinsite.s$burnmowcliptill>0 & trtsinsite.s$P>0, 1, 0)
trtsinsite.s$BMCTxtemp=ifelse(trtsinsite.s$burnmowcliptill>0 & trtsinsite.s$temp>0, 1, 0)
trtsinsite.s$CO2xdrought=ifelse(trtsinsite.s$CO2>0 & trtsinsite.s$drought>0, 1, 0)
trtsinsite.s$CO2xirr=ifelse(trtsinsite.s$CO2>0 & trtsinsite.s$irr>0, 1, 0)
trtsinsite.s$CO2xmult_nutrient=ifelse(trtsinsite.s$CO2>0 & trtsinsite.s$mult_nutrient>0, 1, 0)
trtsinsite.s$CO2xN=ifelse(trtsinsite.s$CO2>0 & trtsinsite.s$N>0, 1, 0)
trtsinsite.s$CO2xP=ifelse(trtsinsite.s$CO2>0 & trtsinsite.s$P>0, 1, 0)
trtsinsite.s$CO2xtemp=ifelse(trtsinsite.s$CO2>0 & trtsinsite.s$temp>0, 1, 0)
trtsinsite.s$droughtxirr=ifelse(trtsinsite.s$drought>0 & trtsinsite.s$irr>0, 1, 0)
trtsinsite.s$droughtxmult_nutrient=ifelse(trtsinsite.s$drought>0 & trtsinsite.s$mult_nutrient>0, 1, 0)
trtsinsite.s$droughtxN=ifelse(trtsinsite.s$drought>0 & trtsinsite.s$N>0, 1, 0)
trtsinsite.s$droughtxP=ifelse(trtsinsite.s$drought>0 & trtsinsite.s$P>0, 1, 0)
trtsinsite.s$droughtxtemp=ifelse(trtsinsite.s$drought>0 & trtsinsite.s$temp>0, 1, 0)
trtsinsite.s$irrxmult_nutrient=ifelse(trtsinsite.s$irr>0 & trtsinsite.s$mult_nutrient>0, 1, 0)
trtsinsite.s$irrxN=ifelse(trtsinsite.s$irr>0 & trtsinsite.s$N>0, 1, 0)
trtsinsite.s$irrxP=ifelse(trtsinsite.s$irr>0 & trtsinsite.s$P>0, 1, 0)
trtsinsite.s$irrxtemp=ifelse(trtsinsite.s$irr>0 & trtsinsite.s$temp>0, 1, 0)
trtsinsite.s$mult_nutrientxN=ifelse(trtsinsite.s$mult_nutrient>0 & trtsinsite.s$N>0, 1, 0)
trtsinsite.s$mult_nutrientxP=ifelse(trtsinsite.s$mult_nutrient>0 & trtsinsite.s$P>0, 1, 0)
trtsinsite.s$mult_nutrientxtemp=ifelse(trtsinsite.s$mult_nutrient>0 & trtsinsite.s$temp>0, 1, 0)
trtsinsite.s$NxP=ifelse(trtsinsite.s$N>0 & trtsinsite.s$P>0, 1, 0)
trtsinsite.s$Nxtemp=ifelse(trtsinsite.s$N>0 & trtsinsite.s$temp>0, 1, 0)
trtsinsite.s$Pxtemp=ifelse(trtsinsite.s$P>0 & trtsinsite.s$temp>0, 1, 0)
#write.csv(trtsinsite.s[,1:9], "which sites have which treatments.csv", row.names=F)

#how many studies for each tradeoff? formatting for use later in loops
temp=subset(trtsinsite.s, select=-c(site_project_comm, burnmowcliptill, CO2, drought, irr, mult_nutrient, N, P, temp))
temp2=colSums(temp)
temp2=data.frame(studies=temp2)
temp2$tradeoffs=row.names(temp2)
trtcombos.avail=temp2[temp2$studies>0,]
temp=matrix(unlist(strsplit(as.character(unique(trtcombos.avail$tradeoffs)), "x")), ncol=2, byrow=T)
trtcombos.avail$var1=temp[,1]
trtcombos.avail$var1pretty=as.factor(trtcombos.avail$var1)
levels(trtcombos.avail$var1pretty)=c("Disturbance", "CO2", "Drought", "Irrigation", "Multiple Nutrients", "Nitrogen")
trtcombos.avail$var2=temp[,2]
trtcombos.avail$var2pretty=as.factor(trtcombos.avail$var2)
levels(trtcombos.avail$var2pretty)=c("Drought", "Irrigation", "Multiple Nutrients", "Nitrogen", "Phosphorus", "Temperature")
trtcombos.avail.m=melt(trtcombos.avail, id="tradeoffs", measure=c("var1", "var2"))
names(trtcombos.avail.m)=c("tradeoff", "variable", "trt")
ntradeoffs.all=unique(trtcombos.avail.m$tradeoff)

#identifying which tradeoffs can be tested in which site(s)
trtcombos.m=melt(trtsinsite.s, id="site_project_comm", measure=c(row.names(trtcombos.avail)))
trtcombos.m=trtcombos.m[trtcombos.m$value>0,]

#preparing for tradeoffs that we can test at only one site:
trtcombos.1site=melt(trtcombos.avail[trtcombos.avail$studies==1,], id="tradeoffs", measure=c("var1", "var2"))
names(trtcombos.1site)[3]="trt"
ntradeoffs.1site=unique(trtcombos.1site$tradeoffs)

#preparing for tradeoffs that we can test at more than one site:
trtcombos.manysite=melt(trtcombos.avail[trtcombos.avail$studies>1,], id="tradeoffs", measure=c("var1", "var2"))
names(trtcombos.manysite)[3]="trt"
ntradeoffs.manysite=unique(trtcombos.manysite$tradeoffs)


#----------FINDING SPECIES RELATIVE ABUNDANCE IN CONTROL PLOTS


realdat.tot=dcast(fulldat[fulldat$permutation=="real",], site_project_comm~., value.var="control", fun=sum)
names(realdat.tot)[2]="real.control.tot"
sprelabund=merge(fulldat[fulldat$permutation=="real",], realdat.tot, by="site_project_comm")
sprelabund=sprelabund[,c("site_project_comm", "species", "control", "real.control.tot")]
sprelabund$relabund.control=sprelabund$control/sprelabund$real.control.tot
write.csv(sprelabund, "files for adam/species relative abundances in control plots.csv", row.names=F)
length(sprelabund$relabund.control[sprelabund$relabund.control<=0.01]) #number of rare experiment-species
length(sprelabund$relabund.control[sprelabund$relabund.control>0.01]) #number of common experiment-species

#----------CALCULATE EFFECT SIZES: USING E ONLY THROUGHTOUT THIS WHOLE FILE (SEE tradeoffs for ms, pre-nov 2019.R for LRR, Eo, %stim)

#E (standardizes by mean abundance across all trts)

E.f=fulldat[,c("site_project_comm", "species", "permutation")]
E.f$BMCT=(fulldat$burnmowcliptill-fulldat$control)/((fulldat$burnmowcliptill+fulldat$control))
E.f$CO2=(fulldat$CO2-fulldat$control)/((fulldat$CO2+fulldat$control))
E.f$drought=(fulldat$drought-fulldat$control)/((fulldat$drought+fulldat$control))
E.f$irr=(fulldat$irr-fulldat$control)/((fulldat$irr+fulldat$control))
E.f$mult_nutrient=(fulldat$mult_nutrient-fulldat$control)/((fulldat$mult_nutrient+fulldat$control))
E.f$N=(fulldat$N-fulldat$control)/((fulldat$N+fulldat$control))
E.f$P=(fulldat$P-fulldat$control)/((fulldat$P+fulldat$control))
E.f$temp=(fulldat$temp-fulldat$control)/((fulldat$temp+fulldat$control))
E.f.m=melt(E.f, id=c("site_project_comm", "species", "permutation"))

#summary statistics on effect sizes for rare vs common species:
E.f.m2=E.f.m[E.f.m$permutation=="real" & !is.na(E.f.m$value),]
E.f.m2=merge(E.f.m2, sprelabund, by=c("site_project_comm", "species"))
E.f.m2$resp=as.factor(ifelse(E.f.m2$value>0, "pos", ifelse(E.f.m2$value<0, "neg", "null")))
E.f.m2$abund=as.factor(ifelse(E.f.m2$relabund.control>0.01, "abundant", "rare"))
allsphist=qplot(E.f.m2$value, ylab="Number of species", xlab="", main="a) All species", xlim=c(-1.1, 1.1))
dcast(E.f.m2, .~resp, value.var="value", fun=length)
1830/(1830+61+1900) #neg 0.48
1900/(1830+61+1900) #pos 0.50
pos.neg.by.trt=dcast(E.f.m2, variable~resp, value.var="value", fun=length)
pos.neg.by.trt$totalsp=pos.neg.by.trt$neg + pos.neg.by.trt$null + pos.neg.by.trt$pos
pos.neg.by.trt$prop.pos.allsp=pos.neg.by.trt$pos/pos.neg.by.trt$totalsp
pos.neg.by.trt$prop.neg.allsp=pos.neg.by.trt$neg/pos.neg.by.trt$totalsp

E.f.m.raresp=E.f.m2[E.f.m2$relabund.control<=0.01,]
raresphist=qplot(E.f.m.raresp$value, ylab="Number of species", xlab="", main="b) Rare species", xlim=c(-1.1, 1.1))
dcast(E.f.m.raresp, .~resp, value.var="value", fun=length)
1129/(1129+59+1434) #neg 0.43
1434/(1129+59+1434) #pos 0.55
raresp.pos.neg.by.trt=dcast(E.f.m.raresp, variable~resp, value.var="value", fun=length)
raresp.pos.neg.by.trt$raresp=raresp.pos.neg.by.trt$neg + raresp.pos.neg.by.trt$null + raresp.pos.neg.by.trt$pos
raresp.pos.neg.by.trt$prop.pos.raresp=raresp.pos.neg.by.trt$pos/raresp.pos.neg.by.trt$raresp
raresp.pos.neg.by.trt$prop.neg.raresp=raresp.pos.neg.by.trt$neg/raresp.pos.neg.by.trt$raresp

E.f.m.commonsp=E.f.m2[E.f.m2$relabund.control>0.01,]
commonsphist=qplot(E.f.m.commonsp$value, ylab="Number of species", xlab="Effect size (E)", main="c) Abundant species", xlim=c(-1.1, 1.1))
dcast(E.f.m.commonsp, .~resp, value.var="value", fun=length)
701/(701+2+466) #neg 0.60
466/(701+2+466) #pos 0.40
commonsp.pos.neg.by.trt=dcast(E.f.m.commonsp, variable~resp, value.var="value", fun=length)
commonsp.pos.neg.by.trt$commonsp=commonsp.pos.neg.by.trt$neg + commonsp.pos.neg.by.trt$null + commonsp.pos.neg.by.trt$pos
commonsp.pos.neg.by.trt$prop.pos.commonsp=commonsp.pos.neg.by.trt$pos/commonsp.pos.neg.by.trt$commonsp
commonsp.pos.neg.by.trt$prop.neg.commonsp=commonsp.pos.neg.by.trt$neg/commonsp.pos.neg.by.trt$commonsp

pos.neg.by.trt=merge(pos.neg.by.trt, raresp.pos.neg.by.trt, by="variable")
pos.neg.by.trt=merge(pos.neg.by.trt, commonsp.pos.neg.by.trt, by="variable")
pos.neg.by.trt=pos.neg.by.trt[,c("variable", "totalsp", "prop.pos.allsp", "prop.neg.allsp", "raresp", "prop.pos.raresp", "prop.neg.raresp", "commonsp", "prop.pos.commonsp", "prop.neg.commonsp")]
write.csv(pos.neg.by.trt, "supplementary figs/Table S2, rare and common species responses by treatment.csv", row.names=F)

ggarrange(allsphist, raresphist, commonsphist, ncol=1)
ggsave("supplementary figs/Fig S2, histograms of E.pdf", height=7, width=4)


# COMPARING E TO LRR AND % STIM FOR REAL DATA ONLY

#LRR; fills NA whenever abundance in control or treatment plots is zero
LRR.f=fulldat[,c("site_project_comm", "species", "permutation")]
LRR.f$BMCT=ifelse(fulldat$control>0 & fulldat$burnmowcliptill>0, log(fulldat$burnmowcliptill/fulldat$control), NA)
LRR.f$CO2=ifelse(fulldat$control>0 & fulldat$CO2>0, log(fulldat$CO2/fulldat$control), NA) 
LRR.f$drought=ifelse(fulldat$control>0 & fulldat$drought>0, log(fulldat$drought/fulldat$control), NA)
LRR.f$irr=ifelse(fulldat$control>0 & fulldat$irr>0, log(fulldat$irr/fulldat$control), NA)
LRR.f$mult_nutrient=ifelse(fulldat$control>0 & fulldat$mult_nutrient>0, log(fulldat$mult_nutrient/fulldat$control), NA)
LRR.f$N=ifelse(fulldat$control>0 & fulldat$N>0, log(fulldat$N/fulldat$control), NA)
LRR.f$P=ifelse(fulldat$control>0 & fulldat$P>0, log(fulldat$P/fulldat$control), NA)
LRR.f$temp=ifelse(fulldat$control>0 & fulldat$temp>0, log(fulldat$temp/fulldat$control), NA)
LRR.f.m=melt(LRR.f, id=c("site_project_comm", "species", "permutation"))
LRR.f.m2=LRR.f.m[LRR.f.m$permutation=="real" & !is.na(LRR.f.m$value),]

#percent stimulation; fills NA whenever abundance in control plots is zero
PS.f=fulldat[,c("site_project_comm", "species", "permutation")]
PS.f$BMCT=ifelse(fulldat$control>0, (fulldat$burnmowcliptill-fulldat$control)/(fulldat$control), NA)
PS.f$CO2=ifelse(fulldat$control>0, (fulldat$CO2-fulldat$control)/(fulldat$control), NA)
PS.f$drought=ifelse(fulldat$control>0, (fulldat$drought-fulldat$control)/(fulldat$control), NA)
PS.f$irr=ifelse(fulldat$control>0, (fulldat$irr-fulldat$control)/(fulldat$control), NA)
PS.f$mult_nutrient=ifelse(fulldat$control>0, (fulldat$mult_nutrient-fulldat$control)/(fulldat$control), NA)
PS.f$N=ifelse(fulldat$control>0, (fulldat$N-fulldat$control)/(fulldat$control), NA)
PS.f$P=ifelse(fulldat$control>0, (fulldat$P-fulldat$control)/(fulldat$control), NA)
PS.f$temp=ifelse(fulldat$control>0, (fulldat$temp-fulldat$control)/(fulldat$control), NA)
PS.f.m=melt(PS.f, id=c("site_project_comm", "species", "permutation"))
PS.f.m2=PS.f.m[PS.f.m$permutation=="real" & !is.na(PS.f.m$value),]

compare.metrics=merge(E.f.m2, LRR.f.m2, by=c("site_project_comm", "species", "variable", "permutation"))
compare.metrics=merge(compare.metrics, PS.f.m2, by=c("site_project_comm", "species", "variable", "permutation"))

EvsLRR=qplot(value.x, value.y, data=compare.metrics, xlab="E", ylab="LRR", alpha=I(0.1))
cor.test(compare.metrics$value.x, compare.metrics$value.y, method="spearman")

EvsPS=qplot(value.x, value, data=compare.metrics, xlab="E", ylab="Percent stimulation", alpha=I(0.1), ylim=c(-10, 200))
cor.test(compare.metrics$value.x, compare.metrics$value, method="spearman")

ggarrange(EvsLRR, EvsPS, ncol=1)
ggsave("supplementary figs/Fig S1, comparison of E with LRR and percent stim.pdf", width=3, height=5)


#---------COMPARING SIMULATED DATA TO REAL DATA: all species


#looping through all tradeoffs we want
perm.list=unique(fulldat$permutation)
output.lm.E.allsp=numeric(0)
datallreal=numeric(0)

for (h in 1:length(ntradeoffs.all)) {
	trdf=as.data.frame(trtcombos.m[trtcombos.m$variable==as.character(ntradeoffs.all[h]), c("site_project_comm")])
	names(trdf)="site_project_comm"
	dath=merge(E.f.m, trdf, by="site_project_comm")
	names(dath)[4:5]=c("trt", "eff")
	trdf.vars.m=trtcombos.avail.m[trtcombos.avail.m$tradeoff==as.character(ntradeoffs.all[h]),]
	dath=merge(dath, trdf.vars.m, by="trt")
	want=dcast(dath, site_project_comm+species+tradeoff+permutation~variable, value.var="eff")
	wantreal=want[want$permutation=="real",]
	datallreal=rbind(datallreal, wantreal)
	nsites=unique(want$site_project_comm)
	temp.i=numeric(0)
	
	for(i in 1:length(perm.list)) {
		dati=want[want$permutation==as.character(perm.list[i]),]
		temp.j=numeric(0)
		
		for(j in 1:length(nsites)) {
			datj=dati[dati$site_project_comm==as.character(nsites[j]),]
			datj=datj[complete.cases(datj[,c("var1", "var2")]),]
			datj$QI=ifelse(datj$var1>0 & datj$var2>0, 1, 0)
			datj$QII=ifelse(datj$var1<0 & datj$var2>0, 1, 0)
			datj$QIII=ifelse(datj$var1<0 & datj$var2<0, 1, 0)
			datj$QIV=ifelse(datj$var1>0 & datj$var2<0, 1, 0)
			nspecies=length(datj$species)
			mod=lm(var2~var1, data=datj)
			temp.j2=data.frame(site_project_comm=as.character(nsites[j]), intercept=mod$coef[1], slope=mod$coef[2], R2=summary(mod)$r.squared, prop.QI=sum(datj$QI)/nspecies, prop.QII=sum(datj$QII)/nspecies, prop.QIII=sum(datj$QIII)/nspecies, prop.QIV=sum(datj$QIV)/nspecies)
			temp.j=rbind(temp.j, temp.j2)
			#for each site (in that permutation), get the relationship between var1 and var2 and export stuff, labeling with the site name
			}

		temp.i2=data.frame(temp.j, permutation=as.character(perm.list[i]))
		temp.i=rbind(temp.i, temp.i2)
		#collect all that for each permutation and label with the permutation
	}
	
	#loop through sites to rank observed slope against all permuted slopes
	temp.h=numeric(0)
	for(g in 1:length(nsites)) {
	dat=temp.i[temp.i$site_project_comm==as.character(nsites[g]),]
	dat$prop.correspondence=dat$prop.QI + dat$prop.QIII
	dat$prop.tradeoffs=dat$prop.QII + dat$prop.QIV
	dat$slope.ranks=rank(dat$slope)
	dat$QI.ranks=rank(dat$prop.QI)
	dat$QII.ranks=rank(dat$prop.QII)
	dat$QIII.ranks=rank(dat$prop.QIII)
	dat$QIV.ranks=rank(dat$prop.QIV)	
	dat$correspondence.ranks=rank(dat$prop.correspondence)
	dat$tradeoffs.ranks=rank(dat$prop.tradeoffs)
	temp.h=rbind(temp.h, dat)
	
}

	temp.h2=data.frame(temp.h, tradeoff=as.character(ntradeoffs.all[h]))
	output.lm.E.allsp=rbind(output.lm.E.allsp, temp.h2)
	#collect all that for each tradeoff and label with the tradeoff
}

output.m.E.allsp=melt(output.lm.E.allsp, id=c("site_project_comm", "permutation", "tradeoff"))

#getting nice labels for plotting later:
tradeoff.labels=data.frame(tradeoff=unique(output.m.E.allsp$tradeoff))
tradeoff.labels$tradeoff2=tradeoff.labels$tradeoff
levels(tradeoff.labels$tradeoff2)=c("Disturbance x Drought", "Disturbance x Irrigation", "Disturbance x Nitrogen", "Disturbance x Phosphorus", "Disturbance x Temperature", "CO2 x Irrigation", "CO2 x Nitrogen", "CO2 x Temperature", "Drought x Irrigation", "Drought x Nitrogen", "Drought x Temperature", "Irrigation x Multiple Nutrients", "Irrigation x Nitrogen", "Irrigation x Phosphorus", "Irrigation x Temperature", "Multiple Nutrients x Temperature", "Nitrogen x Phosphorus", "Nitrogen x Temperature")

#slopes
slopesforadam=output.m.E.allsp[output.m.E.allsp$variable=="slope",]
slopesforadam=merge(slopesforadam, tradeoff.labels, by="tradeoff")
write.csv(slopesforadam, "files for adam/slopes of real and permuted communities for each tradeoff.csv", row.names=F)

output.realslopes=dcast(output.m.E.allsp[output.m.E.allsp$permutation=="real" & output.m.E.allsp$variable=="slope",], site_project_comm+tradeoff~variable)
names(output.realslopes)[3]="realslope"
output.simslopes=dcast(output.m.E.allsp[output.m.E.allsp$variable=="slope" & !output.m.E.allsp$permutation=="real",], site_project_comm+tradeoff~variable, mean)
names(output.simslopes)[3]="simulatedslope"
output.slopes=merge(output.realslopes, output.simslopes, by=c("site_project_comm", "tradeoff"))
output.slopes$dif=output.slopes$realslope-output.slopes$simulatedslope

#getting p-values on our slopes
nperms=length(perm.list)
output.realsloperanks=dcast(output.m.E.allsp[output.m.E.allsp$permutation=="real" & output.m.E.allsp$variable=="slope.ranks",], site_project_comm+tradeoff~variable)
output.realsloperanks=merge(output.slopes, output.realsloperanks, by=c("site_project_comm", "tradeoff"))
output.realsloperanks$p=ifelse(output.realsloperanks$slope.ranks<=nperms/2, output.realsloperanks$slope.ranks/nperms, (nperms-output.realsloperanks$slope.ranks)/nperms)*2
write.csv(output.realsloperanks, "files for adam/slope differences and significances for each tradeoff.csv", row.names=F)

#plotting histogram of slopes
realslope.hist=qplot(realslope, data=output.realsloperanks, xlab="", ylab="Count", xlim=c(-1.1, 1.1), main="a) Observed slopes")
simslope.hist=qplot(simulatedslope, data=output.realsloperanks, xlab="", ylab="Count", xlim=c(-1.1, 1.1), main="b) Simulated slopes")
difslope.hist=qplot(dif, data=output.realsloperanks, xlab="Slope", ylab="Count", xlim=c(-1.1, 1.1), main="c) Correspondence index\n(Observed slope - simulated slope)")
ggarrange(realslope.hist, simslope.hist, difslope.hist, ncol=1)
ggsave("supplementary figs/Fig S4, slope distributions.pdf", width=4, height=8)

#proportion of species in 4 quadrants or with more useful bins
output.realprop=dcast(output.m.E.allsp[output.m.E.allsp$permutation=="real" & output.m.E.allsp$variable %in% c("prop.QI", "prop.QII", "prop.QIII", "prop.QIV", "prop.correspondence", "prop.tradeoffs"),], site_project_comm+tradeoff+variable~.)
names(output.realprop)[4]="obsprop"
output.simprop=dcast(output.m.E.allsp[output.m.E.allsp$variable %in% c("prop.QI", "prop.QII", "prop.QIII", "prop.QIV", "prop.correspondence", "prop.tradeoffs") & !output.m.E.allsp$permutation=="real",], site_project_comm+tradeoff+variable~., mean)
names(output.simprop)[4]="simulatedprop"
temp=dcast(output.simprop, tradeoff~variable, value.var="simulatedprop", fun=mean); temp #why do quadrants I and III have more species than II and IV???
output.proportions=merge(output.realprop, output.simprop, by=c("site_project_comm", "tradeoff", "variable"))
names(output.proportions)[3]="quadrant"
output.proportions$dif.in.prop=output.proportions$obsprop-output.proportions$simulatedprop
write.csv(output.proportions, "files for adam/proportions of species in dif quadrants, all species.csv", row.names=F)

#getting p-values on our proportions
nperms=length(perm.list)
output.realpropranks=dcast(output.m.E.allsp[output.m.E.allsp$permutation=="real" & output.m.E.allsp$variable %in% c("QI.ranks", "QII.ranks", "QIII.ranks", "QIV.ranks", "correspondence.ranks", "tradeoffs.ranks"),], site_project_comm+tradeoff+variable~.)
names(output.realpropranks)[4]="quad.rank"
output.realpropranks$quadrant=output.realpropranks$variable
levels(output.realpropranks$quadrant)=c("intercept", "slope", "R2", "prop.QI", "prop.QII", "prop.QIII", "prop.QIV", "prop.correspondence", "prop.tradeoffs", "slope.ranks", "QI", "QII", "QIII", "QIV", "correspondence", "tradeoffs")
temp=output.proportions
levels(temp$quadrant)=c("intercept", "slope", "R2", "QI", "QII", "QIII", "QIV", "correspondence", "tradeoffs", "slope.ranks", "QI.ranks", "QII.ranks", "QIII.ranks", "QIV.ranks", "correspondence.ranks", "tradeoffs.ranks")
output.realpropranks=merge(temp, output.realpropranks, by=c("site_project_comm", "tradeoff", "quadrant"))
output.realpropranks$p=ifelse(output.realpropranks$quad.rank<=nperms/2, output.realpropranks$quad.rank/nperms, (nperms-output.realpropranks$quad.rank)/nperms)*2

#testing to see if proportions in quadrants correlates with slope difference (to see if our 2 metrics/approaches agree)
ranks.m=melt(output.realpropranks, id=c("site_project_comm", "tradeoff", "quadrant"), measure="dif.in.prop")
ranks.c=dcast(ranks.m, site_project_comm + tradeoff~quadrant)
for.hypoth.of.causality=merge(output.realsloperanks[,c("site_project_comm", "tradeoff", "dif")], ranks.c, by=c("site_project_comm", "tradeoff"))
cor.test(for.hypoth.of.causality$dif, for.hypoth.of.causality$correspondence, method="spearman")
qplot(correspondence, dif, data=for.hypoth.of.causality, xlab="proportion of species showing correspondence", ylab="difference between simulated and observed slope") + facet_wrap(~tradeoff)
ggsave("other figs/correlation of slope difference with correspondence.pdf", height=6, width=8)
cor.test(for.hypoth.of.causality$dif, for.hypoth.of.causality$tradeoffs, method="spearman")
qplot(tradeoffs, dif, data=for.hypoth.of.causality, xlab="proportion of species showing tradeoffs", ylab="difference between simulated and observed slope") + facet_wrap(~tradeoff)
ggsave("other figs/correlation of slope difference with tradeoffs.pdf", height=6, width=8)
qplot(tradeoffs, dif, data=for.hypoth.of.causality, xlab="proportion of species showing tradeoffs", ylab="difference between simulated and observed slope")


#--APPROACH 1 QUADRANTS------PLOTTING: SPECIES RESPONSES TO TREATMENT PAIRS:

#first cleaning up long list of species responses to pairs of treatments for each tradeoff:
datallreal=datallreal[complete.cases(datallreal[,c("var1", "var2")]),]
datallreal=merge(datallreal, tradeoff.labels, by="tradeoff")

#adding mean abund in control plots:
sprelabund2=sprelabund[,c("site_project_comm", "species", "relabund.control"),]
datallreal=merge(datallreal, sprelabund2, by=c("site_project_comm", "species"))

write.csv(datallreal, "files for adam/species responses to pairs of treatments.csv", row.names=F)

#all sites together on one plot:
qplot(var1, var2, data=datallreal, size=relabund.control, alpha=I(0.2), xlab="Effect of treatment 1 (E)", ylab="Effect of treatment 2 (E)") + facet_wrap(~tradeoff2, ncol=3) + geom_vline(xintercept=0) + geom_hline(yintercept=0) + geom_smooth(method="lm", color="black") + theme(legend.position="none")
ggsave("supplementary figs/fig S3, species responses to pairs of treatments, by treatment combo.pdf", height=12, width=6.5)

qplot(var1, var2, size=relabund.control, alpha=I(0.2), data=datallreal) + geom_vline(xintercept=0) + geom_hline(yintercept=0) + geom_smooth(method="lm", color="black") + xlab("Effect of treatment 1 (E)") + ylab("Effect of treatment 2 (E)") + theme(legend.position="none")
ggsave("ms figs/fig 2, species responses to pairs of treatments, all sites together.pdf", height=6, width=6)

#some summary statistics on that plot:
datallreal$quad=as.factor(ifelse(datallreal$var1<0 & datallreal$var2<0, "loser", ifelse(datallreal$var1>0 & datallreal$var2>0, "winner", "mixed")))
summary(datallreal$quad) #number of species in each quadrant
summary(datallreal$quad[datallreal$relabund.control<=0.01]) #number of rare species in each quadrant
summary(datallreal$quad[datallreal$relabund.control>0.01]) #number of common species in each quadrant


#--APPROACH 1 QUADRANTS-------PLOTTING: BOXPLOT OF TRADEOFFS VS WINNERS/LOSERS/MIXED:

#plotting significance of the proportion winners/losers/mixed, sequencing by proportion mixed (=tradeoffs):
output.realWLMranks=output.realpropranks[output.realpropranks$quadrant %in% c("QI", "QIII", "tradeoffs"),]
output.realWLMranks=merge(output.realWLMranks, tradeoff.labels, by="tradeoff")
output.realWLMranks$colorcode=ifelse(output.realWLMranks$p<0.05 & output.realWLMranks$quadrant=="QI", "red", ifelse(output.realWLMranks$p<0.05 & output.realWLMranks$quadrant=="QIII", "red", ifelse(output.realWLMranks$p<0.05 & output.realWLMranks$quadrant=="tradeoffs", "blue", "NA")))
output.realWLMranks$quadrant=droplevels(output.realWLMranks$quadrant)
levels(output.realWLMranks$quadrant)=c("Winners", "Losers", "Mixed")

#adding means of proportions for each tradeoff for plotting:
props.m=melt(output.realWLMranks, id=c("site_project_comm", "tradeoff2", "quadrant"), measure="dif.in.prop")
meanprops=dcast(props.m, tradeoff2+quadrant~., fun=mean)
names(meanprops)[3]="average"
write.csv(meanprops, "files for adam/average proportion of winners, losers, mixed for quilt.csv", row.names=F)
meanprops=meanprops[meanprops$quadrant=="Mixed",]
meanprops$quadrant=NULL
output.realWLMranks=merge(output.realWLMranks, meanprops, by="tradeoff2")

temp=output.realWLMranks
levels(temp$quadrant)=c("a) Winners", "b) Losers", "c) Mixed")
ggplot(data=temp, aes(x=dif.in.prop, y=reorder(tradeoff2, average))) + facet_wrap(~quadrant, nrow=1) + geom_point(shape=I(21), fill=output.realWLMranks$colorcode) + xlab("Observed - expected proportion of species") + ylab("") + geom_vline(xintercept=0) + theme(strip.text = element_text(hjust = 0))
ggsave("ms figs/fig 3, boxplot of species that are winners, losers, mixed.pdf", height=4, width=5)

#compiling info for correlation analyses later:
WLMranks.allsp=output.realWLMranks
WLMranks.allsp$average=NULL
WLMranks.allsp$spgroup="All species"

#--APPROACH 2 SLOPES-------PLOTTING: BOXPLOT OF TRADEOFFS VS WINNERS/LOSERS FOR ALL TRADEOFFS: 

output.realsloperanks=merge(output.realsloperanks, tradeoff.labels, by="tradeoff")
output.realsloperanks$colorcode=ifelse(output.realsloperanks$p<0.05 & output.realsloperanks$dif>0, "red", ifelse(output.realsloperanks$p<0.05 & output.realsloperanks$dif<0, "blue", "NA"))

#adding means of differences for each tradeoff for plotting:
slopes.m=melt(output.realsloperanks, id=c("site_project_comm", "tradeoff2"), measure="dif")
meanslopes=dcast(slopes.m, tradeoff2~., fun=mean)
names(meanslopes)[2]="average"
output.realsloperanks=merge(output.realsloperanks, meanslopes, by="tradeoff2")

#compiling all that for one big figure later:
big.boxplot.E.allsp=output.realsloperanks
big.boxplot.E.allsp$spgroup="All species"
big.boxplot.E.means.allsp=meanslopes
big.boxplot.E.means.allsp$spgroup="All species"

ggplot(data=output.realsloperanks, aes(x=dif, y=reorder(tradeoff2, average))) + geom_point(shape=I(21), fill=output.realsloperanks$colorcode) + xlab("Observed - expected slope") + ylab("") + geom_vline(xintercept=0)
ggsave("other figs/boxplot of all species slopes E, without means.pdf", height=4, width=5)

ggplot(data=output.realsloperanks) + geom_point(shape=I(21), aes(dif, y=reorder(tradeoff2, average)), fill=output.realsloperanks$colorcode) + xlab("Observed - expected slope") + ylab("") + theme(legend.position="none") + geom_vline(xintercept=0) + geom_point(data=meanslopes, aes(average, tradeoff2, size=I(2)))
ggsave("other figs/boxplot of all species slopes E, with means.pdf", height=4, width=5)


#--APPROACH 2 SLOPES-------PLOTTING: REGRESSIONS AND HISTOGRAMS FOR EACH TRADEOFF, EACH SITE: 

#wow: http://www.sthda.com/english/articles/32-r-graphics-essentials/126-combine-multiple-ggplots-in-one-graph/

#looping through to create a single figure containing the regression + histogram for each site, each tradeoff

for (h in 1:length(ntradeoffs.all)) {
	trdf=as.data.frame(trtcombos.m[trtcombos.m$variable==as.character(ntradeoffs.all[h]), c("site_project_comm")])
	names(trdf)="site_project_comm"
	dath=merge(E.f.m, trdf, by="site_project_comm")
	names(dath)[4:5]=c("trt", "eff")
	trdf.vars.m=trtcombos.avail.m[trtcombos.avail.m$tradeoff==as.character(ntradeoffs.all[h]),]
	dath=merge(dath, trdf.vars.m, by="trt")
	want=dcast(dath, site_project_comm+species+tradeoff+permutation~variable, value.var="eff")
	nsites=unique(want$site_project_comm)
	
	for(g in 1:length(nsites)) {
	dat=want[want$site_project_comm==as.character(nsites[g]),]

		xlabname=paste("Species response to", trtcombos.avail$var1pretty[trtcombos.avail$tradeoffs==as.character(ntradeoffs.all[h])], "(E)", sep=" ")
		ylabname=paste("Species response to", trtcombos.avail$var2pretty[trtcombos.avail$tradeoffs==as.character(ntradeoffs.all[h])], "(E)", sep=" ")
		reg=ggplot(dat[!dat$permutation=="real",], aes(var1, var2)) + geom_line(stat="smooth", method="lm", color="gray", alpha=0.1, aes(group=permutation)) + geom_point(data=dat[dat$permutation=="real",], alpha=0.9) + geom_line(data=dat[dat$permutation=="real",], stat="smooth", method="lm") + xlab(xlabname) + ylab(ylabname)

		dat2=output.lm.E.allsp[output.lm.E.allsp$tradeoff==as.character(ntradeoffs.all[h]) & output.lm.E.allsp$site_project_comm==as.character(nsites[g]),]
		xlabname2=paste("Slope of species responses to\n", trtcombos.avail$var1pretty[trtcombos.avail$tradeoffs==as.character(ntradeoffs.all[h])], "and", trtcombos.avail$var2pretty[trtcombos.avail$tradeoffs==as.character(ntradeoffs.all[h])], sep=" ")
		his=ggplot(data=dat2[!dat2$permutation=="real",], aes(slope, fill=I("gray"))) + geom_histogram() + geom_vline(data=dat2[dat2$permutation=="real",], aes(xintercept=slope)) + xlab(xlabname2)

		figure=ggarrange(reg, his, ncol=1)
		sitelabel=paste("Site: ", as.character(nsites[g]), sep="")
		annotate_figure(figure, top=text_grob(sitelabel))
		ggsave(paste("supplementary figs/reg and hist/", as.character(ntradeoffs.all[h]), ",", as.character(nsites[g]), ".pdf", sep=""), width=4, height=6)
	}
}

#----Nxtemp for example figure

trdf=as.data.frame(trtcombos.m[trtcombos.m$variable=="Nxtemp", c("site_project_comm")])
names(trdf)="site_project_comm"
dath=merge(E.f.m, trdf, by="site_project_comm")
names(dath)[4:5]=c("trt", "eff")
trdf.vars.m=trtcombos.avail.m[trtcombos.avail.m$tradeoff=="Nxtemp",]
dath=merge(dath, trdf.vars.m, by="trt")
want=dcast(dath, site_project_comm+species+tradeoff+permutation~variable, value.var="eff")
nsites=unique(want$site_project_comm)
nsites.pretty=c("JSP", "NWT", "SEV")
textsize=c(0.01, 0.025, 0.03)
figure=list()

for(g in 1:length(nsites)) {
	dat=want[want$site_project_comm==as.character(nsites[g]),]
	sprelabund2.want=sprelabund2[sprelabund2$site_project_comm==as.character(nsites[g]),]
	dat=merge(dat, sprelabund2.want, by=c("site_project_comm", "species"))	
	
	regraw=ggplot(dat[dat$permutation=="real",], aes(var1, var2)) + geom_hline(yintercept=0, color="gray") + geom_vline(xintercept=0, color="gray") + geom_line(stat="smooth", method="lm", color="gray", alpha=0.1, aes(group=permutation)) + geom_point(aes(size=relabund.control)) + geom_text(aes(label=species, hjust=-0.1, vjust=0, size=textsize[g])) + geom_line(data=dat[dat$permutation=="real",], stat="smooth", method="lm") + xlab("Nitrogen effect (E)") + ylab("Temperature effect (E)") + ggtitle(as.character(nsites.pretty[g])) + theme(legend.position = "none")
	regsim=ggplot(dat[!dat$permutation=="real",], aes(var1, var2)) + geom_line(stat="smooth", method="lm", color="gray", alpha=0.1, aes(group=permutation)) + geom_line(data=dat[dat$permutation=="real",], stat="smooth", method="lm") + xlab("Nitrogen effect (E)") + ylab("Temperature effect (E)")

	dat2=output.lm.E.allsp[output.lm.E.allsp$tradeoff=="Nxtemp" & output.lm.E.allsp$site_project_comm==as.character(nsites[g]),]
	his=ggplot(data=dat2[!dat2$permutation=="real",], aes(slope, fill=I("gray"))) + geom_histogram() + geom_vline(data=dat2[dat2$permutation=="real",], size=1, aes(xintercept=slope)) + xlab("Slope") + ylab("Number of communities")

	figure[[g]]=ggarrange(regraw, regsim, his, ncol=1)
	ggarrange(regraw, regsim, his, ncol=1)
	ggsave(paste("ms figs/fig 4, all species E ", as.character(nsites[g]), ".pdf", sep=""), width=3, height=7)
	}
ggarrange(figure[[1]], figure[[2]], figure[[3]], nrow=1)	
ggsave("ms figs/fig 4, all species E all panels with labels.pdf", width=10, height=10)

#same thing without species labels
trdf=as.data.frame(trtcombos.m[trtcombos.m$variable=="Nxtemp", c("site_project_comm")])
names(trdf)="site_project_comm"
dath=merge(E.f.m, trdf, by="site_project_comm")
names(dath)[4:5]=c("trt", "eff")
trdf.vars.m=trtcombos.avail.m[trtcombos.avail.m$tradeoff=="Nxtemp",]
dath=merge(dath, trdf.vars.m, by="trt")
want=dcast(dath, site_project_comm+species+tradeoff+permutation~variable, value.var="eff")
nsites=unique(want$site_project_comm)
nsites.pretty=c("JSP", "NWT", "SEV")
figure=list()

for(g in 1:length(nsites)) {
	dat=want[want$site_project_comm==as.character(nsites[g]),]
	sprelabund2.want=sprelabund2[sprelabund2$site_project_comm==as.character(nsites[g]),]
	dat=merge(dat, sprelabund2.want, by=c("site_project_comm", "species"))	
	
	regraw=ggplot(dat[dat$permutation=="real",], aes(var1, var2)) + geom_hline(yintercept=0, color="gray") + geom_vline(xintercept=0, color="gray") + geom_line(stat="smooth", method="lm", color="gray", alpha=0.1, aes(group=permutation)) + geom_point(aes(size=relabund.control)) + geom_line(data=dat[dat$permutation=="real",], stat="smooth", method="lm") + xlab("Nitrogen effect (E)") + ylab("Temperature effect (E)") + ggtitle(as.character(nsites.pretty[g])) + theme(legend.position = "none")
	regsim=ggplot(dat[!dat$permutation=="real",], aes(var1, var2)) + geom_line(stat="smooth", method="lm", color="gray", alpha=0.1, aes(group=permutation)) + geom_line(data=dat[dat$permutation=="real",], stat="smooth", method="lm") + xlab("Nitrogen effect (E)") + ylab("Temperature effect (E)")

	dat2=output.lm.E.allsp[output.lm.E.allsp$tradeoff=="Nxtemp" & output.lm.E.allsp$site_project_comm==as.character(nsites[g]),]
	his=ggplot(data=dat2[!dat2$permutation=="real",], aes(slope, fill=I("gray"))) + geom_histogram() + geom_vline(data=dat2[dat2$permutation=="real",], size=1, aes(xintercept=slope)) + xlab("Slope") + ylab("Number of communities")

	figure[[g]]=ggarrange(regraw, regsim, his, ncol=1)
	}
ggarrange(figure[[1]], figure[[2]], figure[[3]], nrow=1)	
ggsave("ms figs/fig 4, all species E all panels no labels.pdf", width=10, height=10)

#---------SUBSETTING TO ONLY COMMON SPECIES (MORE THAN 1% MEAN ABUND ACROSS CONTROL PLOTS)


commonsp=sprelabund[sprelabund$relabund.control>0.01,c("site_project_comm", "species")]

length(fulldat$species[fulldat$permutation=="real"])
length(commonsp$species)
#cuts dataset from 1967 species/site_project_comm combinations to 471

E.c.m=merge(E.f.m, commonsp, by=c("site_project_comm", "species"))


#---------COMPARING SIMULATED DATA TO REAL DATA: only common species


#looping through all tradeoffs we want
perm.list=unique(fulldat$permutation)
output.lm.E.commonsp=numeric(0)

for (h in 1:length(ntradeoffs.all)) {
	trdf=as.data.frame(trtcombos.m[trtcombos.m$variable==as.character(ntradeoffs.all[h]), c("site_project_comm")])
	names(trdf)="site_project_comm"
	dath=merge(E.c.m, trdf, by="site_project_comm")
	names(dath)[4:5]=c("trt", "eff")
	trdf.vars.m=trtcombos.avail.m[trtcombos.avail.m$tradeoff==as.character(ntradeoffs.all[h]),]
	dath=merge(dath, trdf.vars.m, by="trt")
	want=dcast(dath, site_project_comm+species+tradeoff+permutation~variable, value.var="eff")
	nsites=unique(want$site_project_comm)
	temp.i=numeric(0)
	
	for(i in 1:length(perm.list)) {
		dati=want[want$permutation==as.character(perm.list[i]),]
		temp.j=numeric(0)
		
		for(j in 1:length(nsites)) {
			datj=dati[dati$site_project_comm==as.character(nsites[j]),]
			datj=datj[complete.cases(datj[,c("var1", "var2")]),]
			datj$QI=ifelse(datj$var1>0 & datj$var2>0, 1, 0)
			datj$QII=ifelse(datj$var1<0 & datj$var2>0, 1, 0)
			datj$QIII=ifelse(datj$var1<0 & datj$var2<0, 1, 0)
			datj$QIV=ifelse(datj$var1>0 & datj$var2<0, 1, 0)
			nspecies=length(datj$species)
			mod=lm(var2~var1, data=datj)
			temp.j2=data.frame(site_project_comm=as.character(nsites[j]), intercept=mod$coef[1], slope=mod$coef[2], R2=summary(mod)$r.squared, prop.QI=sum(datj$QI)/nspecies, prop.QII=sum(datj$QII)/nspecies, prop.QIII=sum(datj$QIII)/nspecies, prop.QIV=sum(datj$QIV)/nspecies)
			temp.j=rbind(temp.j, temp.j2)
			#for each site (in that permutation), get the relationship between var1 and var2 and export stuff, labeling with the site name
			}

		temp.i2=data.frame(temp.j, permutation=as.character(perm.list[i]))
		temp.i=rbind(temp.i, temp.i2)
		#collect all that for each permutation and label with the permutation
	}
	
	#loop through sites to rank observed slope against all permuted slopes
	temp.h=numeric(0)
	for(g in 1:length(nsites)) {
	dat=temp.i[temp.i$site_project_comm==as.character(nsites[g]),]
	dat$prop.correspondence=dat$prop.QI + dat$prop.QIII
	dat$prop.tradeoffs=dat$prop.QII + dat$prop.QIV
	dat$slope.ranks=rank(dat$slope)
	dat$QI.ranks=rank(dat$prop.QI)
	dat$QII.ranks=rank(dat$prop.QII)
	dat$QIII.ranks=rank(dat$prop.QIII)
	dat$QIV.ranks=rank(dat$prop.QIV)	
	dat$correspondence.ranks=rank(dat$prop.correspondence)
	dat$tradeoffs.ranks=rank(dat$prop.tradeoffs)
	temp.h=rbind(temp.h, dat)
}

	temp.h2=data.frame(temp.h, tradeoff=as.character(ntradeoffs.all[h]))
	output.lm.E.commonsp=rbind(output.lm.E.commonsp, temp.h2)
	#collect all that for each tradeoff and label with the tradeoff
}

output.m.E.commonsp=melt(output.lm.E.commonsp, id=c("site_project_comm", "permutation", "tradeoff"))

#slopes
output.realslopes=dcast(output.m.E.commonsp[output.m.E.commonsp$permutation=="real" & output.m.E.commonsp$variable=="slope",], site_project_comm+tradeoff~variable)
names(output.realslopes)[3]="realslope"
output.simslopes=dcast(output.m.E.commonsp[output.m.E.commonsp$variable=="slope" & !output.m.E.commonsp$permutation=="real",], site_project_comm+tradeoff~variable, mean)
names(output.simslopes)[3]="simulatedslope"
output.slopes=merge(output.realslopes, output.simslopes, by=c("site_project_comm", "tradeoff"))
output.slopes$dif=output.slopes$realslope-output.slopes$simulatedslope

#getting p-values on our slopes
nperms=length(perm.list)
output.realsloperanks=dcast(output.m.E.commonsp[output.m.E.commonsp$permutation=="real" & output.m.E.commonsp$variable=="slope.ranks",], site_project_comm+tradeoff~variable)
output.realsloperanks=merge(output.slopes, output.realsloperanks, by=c("site_project_comm", "tradeoff"))
output.realsloperanks$p=ifelse(output.realsloperanks$slope.ranks<=nperms/2, output.realsloperanks$slope.ranks/nperms, (nperms-output.realsloperanks$slope.ranks)/nperms)*2

#proportion of species in 4 quadrants or with more useful bins
output.realprop=dcast(output.m.E.commonsp[output.m.E.commonsp$permutation=="real" & output.m.E.commonsp$variable %in% c("prop.QI", "prop.QII", "prop.QIII", "prop.QIV", "prop.correspondence", "prop.tradeoffs"),], site_project_comm+tradeoff+variable~.)
names(output.realprop)[4]="obsprop"
output.simprop=dcast(output.m.E.commonsp[output.m.E.commonsp$variable %in% c("prop.QI", "prop.QII", "prop.QIII", "prop.QIV", "prop.correspondence", "prop.tradeoffs") & !output.m.E.commonsp$permutation=="real",], site_project_comm+tradeoff+variable~., mean)
names(output.simprop)[4]="simulatedprop"
output.proportions=merge(output.realprop, output.simprop, by=c("site_project_comm", "tradeoff", "variable"))
names(output.proportions)[3]="quadrant"
output.proportions$dif.in.prop=output.proportions$obsprop-output.proportions$simulatedprop
write.csv(output.proportions, "files for adam/proportions of species in dif quadrants, common species.csv", row.names=F)

#getting p-values on our proportions
nperms=length(perm.list)
output.realpropranks=dcast(output.m.E.commonsp[output.m.E.commonsp$permutation=="real" & output.m.E.commonsp$variable %in% c("QI.ranks", "QII.ranks", "QIII.ranks", "QIV.ranks", "correspondence.ranks", "tradeoffs.ranks"),], site_project_comm+tradeoff+variable~.)
names(output.realpropranks)[4]="quad.rank"
output.realpropranks$quadrant=output.realpropranks$variable
levels(output.realpropranks$quadrant)=c("intercept", "slope", "R2", "prop.QI", "prop.QII", "prop.QIII", "prop.QIV", "prop.correspondence", "prop.tradeoffs", "slope.ranks", "QI", "QII", "QIII", "QIV", "correspondence", "tradeoffs")
temp=output.proportions
levels(temp$quadrant)=c("intercept", "slope", "R2", "QI", "QII", "QIII", "QIV", "correspondence", "tradeoffs", "slope.ranks", "QI.ranks", "QII.ranks", "QIII.ranks", "QIV.ranks", "correspondence.ranks", "tradeoffs.ranks")
output.realpropranks=merge(temp, output.realpropranks, by=c("site_project_comm", "tradeoff", "quadrant"))
output.realpropranks$p=ifelse(output.realpropranks$quad.rank<=nperms/2, output.realpropranks$quad.rank/nperms, (nperms-output.realpropranks$quad.rank)/nperms)*2


#--APPROACH 1 QUADRANTS-------PLOTTING: BOXPLOT OF TRADEOFFS VS WINNERS/LOSERS/MIXED:


#getting the proportion winners/losers/mixed, sequencing by proportion mixed (=tradeoffs):
output.realWLMranks=output.realpropranks[output.realpropranks$quadrant %in% c("QI", "QIII", "tradeoffs"),]
output.realWLMranks=merge(output.realWLMranks, tradeoff.labels, by="tradeoff")
output.realWLMranks$colorcode=ifelse(output.realWLMranks$p<0.05 & output.realWLMranks$quadrant=="QI", "red", ifelse(output.realWLMranks$p<0.05 & output.realWLMranks$quadrant=="QIII", "red", ifelse(output.realWLMranks$p<0.05 & output.realWLMranks$quadrant=="tradeoffs", "blue", "NA")))
output.realWLMranks$quadrant=droplevels(output.realWLMranks$quadrant)
levels(output.realWLMranks$quadrant)=c("Winners", "Losers", "Mixed")
output.realWLMranks=output.realWLMranks[!output.realWLMranks$site_project_comm=="SERC_CXN_0",] #removing SERC because only 3 species there

#compiling info for correlation analyses later:
WLMranks.commonsp=output.realWLMranks
WLMranks.commonsp$spgroup="Common species"


#--APPROACH 2 SLOPES-------PLOTTING: BOXPLOT OF TRADEOFFS VS WINNERS/LOSERS FOR ALL TRADEOFFS: 

output.realsloperanks=merge(output.realsloperanks, tradeoff.labels, by="tradeoff")
output.realsloperanks$colorcode=ifelse(output.realsloperanks$p<0.05 & output.realsloperanks$dif>0, "red", ifelse(output.realsloperanks$p<0.05 & output.realsloperanks$dif<0, "blue", "NA"))
output.realsloperanks=output.realsloperanks[!output.realsloperanks$site_project_comm=="SERC_CXN_0",] #removing SERC because only 3 species there

#adding means of differences for each tradeoff for plotting:
slopes.m=melt(output.realsloperanks, id=c("site_project_comm", "tradeoff2"), measure="dif")
meanslopes=dcast(slopes.m, tradeoff2~., fun=mean)
names(meanslopes)[2]="average"
output.realsloperanks=merge(output.realsloperanks, meanslopes, by="tradeoff2")

#compiling all that for one big figure later:
big.boxplot.E.commonsp=output.realsloperanks
big.boxplot.E.commonsp$spgroup="Common species"
big.boxplot.E.means.commonsp=meanslopes
big.boxplot.E.means.commonsp$spgroup="Common species"

ggplot(data=output.realsloperanks, aes(x=dif, y=reorder(tradeoff2, average))) + geom_point(shape=I(21), fill=output.realsloperanks$colorcode) + xlab("Observed - expected slope") + ylab("") + geom_vline(xintercept=0)
ggsave("other figs/boxplot of common species slopes E, without means.pdf", height=4, width=5)

ggplot(data=output.realsloperanks) + geom_point(shape=I(21), aes(dif, y=reorder(tradeoff2, average)), fill=output.realsloperanks$colorcode) + xlab("Observed - expected slope") + ylab("") + theme(legend.position="none") + geom_vline(xintercept=0) + geom_point(data=meanslopes, aes(average, tradeoff2, size=I(2)))
ggsave("other figs/boxplot of common species slopes E, with means.pdf", height=4, width=5)
	
	
#---------SUBSETTING TO ONLY RARE SPECIES (MEAN LESS THAN OR EQUAL TO 1% ABUND ACROSS CONTROL PLOTS)


raresp=sprelabund[sprelabund$relabund.control<=0.01,c("site_project_comm", "species")]

length(fulldat$species[fulldat$permutation=="real"])
length(raresp$species)
#cuts dataset from 1967 species/site_project_comm combinations to 1496

E.r.m=merge(E.f.m, raresp, by=c("site_project_comm", "species"))


#---------COMPARING SIMULATED DATA TO REAL DATA: only rare species


#looping through all tradeoffs we want
perm.list=unique(fulldat$permutation)
output.lm.E.raresp=numeric(0)

for (h in 1:length(ntradeoffs.all)) {
	trdf=as.data.frame(trtcombos.m[trtcombos.m$variable==as.character(ntradeoffs.all[h]), c("site_project_comm")])
	names(trdf)="site_project_comm"
	dath=merge(E.r.m, trdf, by="site_project_comm")
	names(dath)[4:5]=c("trt", "eff")
	trdf.vars.m=trtcombos.avail.m[trtcombos.avail.m$tradeoff==as.character(ntradeoffs.all[h]),]
	dath=merge(dath, trdf.vars.m, by="trt")
	want=dcast(dath, site_project_comm+species+tradeoff+permutation~variable, value.var="eff")
	nsites=unique(want$site_project_comm)
	temp.i=numeric(0)
	
	for(i in 1:length(perm.list)) {
		dati=want[want$permutation==as.character(perm.list[i]),]
		temp.j=numeric(0)
		
		for(j in 1:length(nsites)) {
			datj=dati[dati$site_project_comm==as.character(nsites[j]),]
			datj=datj[complete.cases(datj[,c("var1", "var2")]),]
			datj$QI=ifelse(datj$var1>0 & datj$var2>0, 1, 0)
			datj$QII=ifelse(datj$var1<0 & datj$var2>0, 1, 0)
			datj$QIII=ifelse(datj$var1<0 & datj$var2<0, 1, 0)
			datj$QIV=ifelse(datj$var1>0 & datj$var2<0, 1, 0)
			nspecies=length(datj$species)
			mod=lm(var2~var1, data=datj)
			temp.j2=data.frame(site_project_comm=as.character(nsites[j]), intercept=mod$coef[1], slope=mod$coef[2], R2=summary(mod)$r.squared, prop.QI=sum(datj$QI)/nspecies, prop.QII=sum(datj$QII)/nspecies, prop.QIII=sum(datj$QIII)/nspecies, prop.QIV=sum(datj$QIV)/nspecies)
			temp.j=rbind(temp.j, temp.j2)
			#for each site (in that permutation), get the relationship between var1 and var2 and export stuff, labeling with the site name
			}

		temp.i2=data.frame(temp.j, permutation=as.character(perm.list[i]))
		temp.i=rbind(temp.i, temp.i2)
		#collect all that for each permutation and label with the permutation
	}
	
	#loop through sites to rank observed slope against all permuted slopes
	temp.h=numeric(0)
	for(g in 1:length(nsites)) {
	dat=temp.i[temp.i$site_project_comm==as.character(nsites[g]),]
	dat$prop.correspondence=dat$prop.QI + dat$prop.QIII
	dat$prop.tradeoffs=dat$prop.QII + dat$prop.QIV
	dat$slope.ranks=rank(dat$slope)
	dat$QI.ranks=rank(dat$prop.QI)
	dat$QII.ranks=rank(dat$prop.QII)
	dat$QIII.ranks=rank(dat$prop.QIII)
	dat$QIV.ranks=rank(dat$prop.QIV)	
	dat$correspondence.ranks=rank(dat$prop.correspondence)
	dat$tradeoffs.ranks=rank(dat$prop.tradeoffs)
	temp.h=rbind(temp.h, dat)
	
}

	temp.h2=data.frame(temp.h, tradeoff=as.character(ntradeoffs.all[h]))
	output.lm.E.raresp=rbind(output.lm.E.raresp, temp.h2)
	#collect all that for each tradeoff and label with the tradeoff
}

#comparing mean of simulated communities to real communities
output.m.E.raresp=melt(output.lm.E.raresp, id=c("site_project_comm", "permutation", "tradeoff"))

#slopes
output.realslopes=dcast(output.m.E.raresp[output.m.E.raresp$permutation=="real" & output.m.E.raresp$variable=="slope",], site_project_comm+tradeoff~variable)
names(output.realslopes)[3]="realslope"
output.simslopes=dcast(output.m.E.raresp[output.m.E.raresp$variable=="slope" & !output.m.E.raresp$permutation=="real",], site_project_comm+tradeoff~variable, mean)
names(output.simslopes)[3]="simulatedslope"
output.slopes=merge(output.realslopes, output.simslopes, by=c("site_project_comm", "tradeoff"))
output.slopes$dif=output.slopes$realslope-output.slopes$simulatedslope

#getting p-values on our slopes
nperms=length(perm.list)
output.realsloperanks=dcast(output.m.E.raresp[output.m.E.raresp$permutation=="real" & output.m.E.raresp$variable=="slope.ranks",], site_project_comm+tradeoff~variable)
output.realsloperanks=merge(output.slopes, output.realsloperanks, by=c("site_project_comm", "tradeoff"))
output.realsloperanks$p=ifelse(output.realsloperanks$slope.ranks<=nperms/2, output.realsloperanks$slope.ranks/nperms, (nperms-output.realsloperanks$slope.ranks)/nperms)*2

#proportion of species in 4 quadrants or with more useful bins
output.realprop=dcast(output.m.E.raresp[output.m.E.raresp$permutation=="real" & output.m.E.raresp$variable %in% c("prop.QI", "prop.QII", "prop.QIII", "prop.QIV", "prop.correspondence", "prop.tradeoffs"),], site_project_comm+tradeoff+variable~.)
names(output.realprop)[4]="obsprop"
output.simprop=dcast(output.m.E.raresp[output.m.E.raresp$variable %in% c("prop.QI", "prop.QII", "prop.QIII", "prop.QIV", "prop.correspondence", "prop.tradeoffs") & !output.m.E.raresp$permutation=="real",], site_project_comm+tradeoff+variable~., mean)
names(output.simprop)[4]="simulatedprop"
output.proportions=merge(output.realprop, output.simprop, by=c("site_project_comm", "tradeoff", "variable"))
names(output.proportions)[3]="quadrant"
output.proportions$dif.in.prop=output.proportions$obsprop-output.proportions$simulatedprop
write.csv(output.proportions, "files for adam/proportions of species in dif quadrants, rare species.csv", row.names=F)

#getting p-values on our proportions
nperms=length(perm.list)
output.realpropranks=dcast(output.m.E.raresp[output.m.E.raresp$permutation=="real" & output.m.E.raresp$variable %in% c("QI.ranks", "QII.ranks", "QIII.ranks", "QIV.ranks", "correspondence.ranks", "tradeoffs.ranks"),], site_project_comm+tradeoff+variable~.)
names(output.realpropranks)[4]="quad.rank"
output.realpropranks$quadrant=output.realpropranks$variable
levels(output.realpropranks$quadrant)=c("intercept", "slope", "R2", "prop.QI", "prop.QII", "prop.QIII", "prop.QIV", "prop.correspondence", "prop.tradeoffs", "slope.ranks", "QI", "QII", "QIII", "QIV", "correspondence", "tradeoffs")
temp=output.proportions
levels(temp$quadrant)=c("intercept", "slope", "R2", "QI", "QII", "QIII", "QIV", "correspondence", "tradeoffs", "slope.ranks", "QI.ranks", "QII.ranks", "QIII.ranks", "QIV.ranks", "correspondence.ranks", "tradeoffs.ranks")
output.realpropranks=merge(temp, output.realpropranks, by=c("site_project_comm", "tradeoff", "quadrant"))
output.realpropranks$p=ifelse(output.realpropranks$quad.rank<=nperms/2, output.realpropranks$quad.rank/nperms, (nperms-output.realpropranks$quad.rank)/nperms)*2


#--APPROACH 1 QUADRANTS-------PLOTTING: BOXPLOT OF TRADEOFFS VS WINNERS/LOSERS/MIXED:


#getting the proportion winners/losers/mixed, sequencing by proportion mixed (=tradeoffs):
output.realWLMranks=output.realpropranks[output.realpropranks$quadrant %in% c("QI", "QIII", "tradeoffs"),]
output.realWLMranks=merge(output.realWLMranks, tradeoff.labels, by="tradeoff")
output.realWLMranks$colorcode=ifelse(output.realWLMranks$p<0.05 & output.realWLMranks$quadrant=="QI", "red", ifelse(output.realWLMranks$p<0.05 & output.realWLMranks$quadrant=="QIII", "red", ifelse(output.realWLMranks$p<0.05 & output.realWLMranks$quadrant=="tradeoffs", "blue", "NA")))
output.realWLMranks$quadrant=droplevels(output.realWLMranks$quadrant)
levels(output.realWLMranks$quadrant)=c("Winners", "Losers", "Mixed")

#compiling info for correlation analyses later:
WLMranks.raresp=output.realWLMranks
WLMranks.raresp$spgroup="Rare species"


#--APPROACH 2 SLOPES-------PLOTTING: BOXPLOT OF TRADEOFFS VS WINNERS/LOSERS FOR ALL TRADEOFFS: 

output.realsloperanks=merge(output.realsloperanks, tradeoff.labels, by="tradeoff")
output.realsloperanks$colorcode=ifelse(output.realsloperanks$p<0.05 & output.realsloperanks$dif>0, "red", ifelse(output.realsloperanks$p<0.05 & output.realsloperanks$dif<0, "blue", "NA"))

#adding means of differences for each tradeoff for plotting:
slopes.m=melt(output.realsloperanks, id=c("site_project_comm", "tradeoff2"), measure="dif")
meanslopes=dcast(slopes.m, tradeoff2~., fun=mean)
names(meanslopes)[2]="average"
output.realsloperanks=merge(output.realsloperanks, meanslopes, by="tradeoff2")

#compiling all that for one big figure later:
big.boxplot.E.raresp=output.realsloperanks
big.boxplot.E.raresp$spgroup="Rare species"
big.boxplot.E.means.raresp=meanslopes
big.boxplot.E.means.raresp$spgroup="Rare species"

ggplot(data=output.realsloperanks, aes(x=dif, y=reorder(tradeoff2, average))) + geom_point(shape=I(21), fill=output.realsloperanks$colorcode) + xlab("Observed - expected slope") + ylab("") + geom_vline(xintercept=0)
ggsave("other figs/boxplot of rare species slopes E, without means.pdf", height=4, width=5)

ggplot(data=output.realsloperanks) + geom_point(shape=I(21), aes(dif, y=reorder(tradeoff2, average)), fill=output.realsloperanks$colorcode) + xlab("Observed - expected slope") + ylab("") + theme(legend.position="none") + geom_vline(xintercept=0) + geom_point(data=meanslopes, aes(average, tradeoff2, size=I(2)))
ggsave("other figs/boxplot of rare species slopes E, with means.pdf", height=4, width=5)


#-------------------OVERALL COMPARISONS OF RARE, COMMON, ALL SPECIES:


#--APPROACH 1: PROPORTIONS OF SPECIES IN QUADRANTS:


WLMranks=rbind(WLMranks.allsp, WLMranks.commonsp, WLMranks.raresp)

temp=WLMranks[WLMranks$spgroup=="All species",]
temp=merge(temp, meanprops, by="tradeoff2")
levels(temp$quadrant)=c("a) Winners", "b) Losers", "c) Mixed")
fig3a=ggplot(data=temp, aes(x=dif.in.prop, y=reorder(tradeoff2, average))) + facet_wrap(~quadrant, nrow=1) + geom_point(shape=I(21), fill=temp$colorcode, size=I(3)) + xlab("Observed - expected proportion of species") + ylab("") + geom_vline(xintercept=0) + theme(strip.text = element_text(hjust = 0))
ggplot(data=temp, aes(y=dif.in.prop, x=reorder(tradeoff2, average))) + facet_wrap(~quadrant, ncol=1) + geom_point(shape=I(21), fill=temp$colorcode, size=I(3)) + ylab("Observed - expected proportion of species") + xlab("") + geom_hline(yintercept=0) + theme(strip.text = element_text(hjust = 0), axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
ggsave("ms figs/fig 3 transposed.pdf", width=4, height=10)

WLMranks=WLMranks[!(WLMranks$site_project_comm=="SERC_CXN_0" & WLMranks$spgroup=="Common species"),] #dropping SERC from common sp (already absent from rare sp)
wlmranks=melt(WLMranks, id=c("tradeoff", "tradeoff2", "site_project_comm", "quadrant", "spgroup"), measure=c("obsprop", "simulatedprop", "dif.in.prop"))
wlmranks[wlmranks$site_project_comm=="SERC_CXN_0",] #confirming that SERC removed for rare and common sp

# is the proportion of losers higher amongst rare than abundant sp? proportion of winners higher amongst common sp?

meanwlmranks=dcast(wlmranks, quadrant+spgroup+variable~., fun=mean)
names(meanwlmranks)[4]="average"
meanwlmranks.sd=dcast(wlmranks, quadrant+spgroup+variable~., fun=sd)
names(meanwlmranks.sd)[4]="stdev"
meanwlmranks.n=dcast(wlmranks, quadrant+spgroup+variable~., fun=length)
names(meanwlmranks.n)[4]="n"
meanwlmranks=merge(meanwlmranks, meanwlmranks.sd, by=c("quadrant", "spgroup", "variable"))
meanwlmranks=merge(meanwlmranks, meanwlmranks.n, by=c("quadrant", "spgroup", "variable"))
meanwlmranks$se=meanwlmranks$stdev/sqrt(meanwlmranks$n)
meanwlmranks$spgroup=as.factor(meanwlmranks$spgroup)
levels(meanwlmranks$spgroup)=c("All", "Abundant", "Rare")
levels(meanwlmranks$variable)=c("Observed", "Expected", "Observed-Expected")

qplot(spgroup, average, data=meanwlmranks, ylab="Mean proportion +/- SE", xlab="Species group") + facet_grid(variable~quadrant, scales="free") + geom_errorbar(aes(ymin=average-se, ymax=average+se, width=0.1))
ggsave("other figs/average proportions of winners, losers, mixed.pdf", width=6, height=5)

levels(meanwlmranks$spgroup)=c("a) All species", "b) Abundant species", "c) Rare species")
qplot(quadrant, average, data=meanwlmranks[meanwlmranks$variable=="Observed-Expected",], ylab="Observed-expected proportion", xlab="Species responses to pairs of treatments") + facet_wrap(~spgroup, ncol=3) + geom_errorbar(aes(ymin=average-se, ymax=average+se, width=0.1)) + theme(strip.text = element_text(hjust = 0)) + geom_hline(yintercept=0)
ggsave("other figs/potential Fig 4, mean quadrant abund for all, abundant, rare sp.pdf", width=6, height=3)

# comparing raresp to commonsp in terms of obs-exp proportions winners and losers. Of the common species, is there a higher proportion losers than would be expected? Similarly, among rare species, is there a higher proportion losers? And is this proportion higher in rare or common sp?

dif.in.proportions=dcast(wlmranks[wlmranks$variable=="dif.in.prop",], tradeoff+site_project_comm+quadrant~spgroup)
names(dif.in.proportions)[4:6]=c("allsp", "commonsp", "raresp")
dif.in.proportions$rare.exceed.common=as.factor(ifelse(dif.in.proportions$raresp>dif.in.proportions$commonsp, "more.rare", ifelse(dif.in.proportions$raresp<dif.in.proportions$commonsp, "more.common", "equal")))
dcast(dif.in.proportions, quadrant~rare.exceed.common, value.var="rare.exceed.common", fun=length) #NAs are coming from SERC, which was removed
levels(dif.in.proportions$quadrant)=c("d) Winners", "e) Losers", "f) Mixed")
qplot(allsp, commonsp, data=dif.in.proportions) + facet_wrap(~quadrant) + geom_abline(slope=1, intercept=0)
qplot(allsp, raresp, data=dif.in.proportions) + facet_wrap(~quadrant) + geom_abline(slope=1, intercept=0)
fig3b.1=qplot(commonsp, raresp, data=dif.in.proportions, xlab="Observed-expected proportion amongst abundant species", ylab="Observed-expected proportion\namongst rare species") + facet_wrap(~quadrant) + geom_abline(slope=1, intercept=0) + theme(strip.text = element_text(hjust = 0))
fig3b=ggarrange(grid.rect(gp=gpar(col=NA)), fig3b.1, nrow=1, widths=c(1, 6))
ggarrange(fig3a, fig3b, nrow=2, heights=c(2, 1))
ggsave("ms figs/fig 3 with two panels.pdf", height=8, width=9)


#--APPROACH 2: SLOPES----------boxplot of slopes:


big.boxplot.E=rbind(big.boxplot.E.allsp, big.boxplot.E.commonsp, big.boxplot.E.raresp)
big.boxplot.E[big.boxplot.E$site_project_comm=="SERC_CXN_0",]
big.boxplot.E.means=rbind(big.boxplot.E.means.allsp, big.boxplot.E.means.commonsp, big.boxplot.E.means.raresp)

assign.seq=big.boxplot.E.means[big.boxplot.E.means$spgroup=="All species",c("tradeoff2", "average")]
assign.seq$tradeoff.ordered=factor(assign.seq$tradeoff2, levels=unique(assign.seq$tradeoff2[order(assign.seq$average)]), ordered=T)
assign.seq$average=NULL

big.boxplot.E=merge(big.boxplot.E, assign.seq, by="tradeoff2")
big.boxplot.E$spgroup=as.factor(big.boxplot.E$spgroup)
levels(big.boxplot.E$spgroup)=c("a) All species", "b) Abundant species", "c) Rare species")

big.boxplot.E.means=merge(big.boxplot.E.means, assign.seq, by="tradeoff2")
big.boxplot.E.means$spgroup=as.factor(big.boxplot.E.means$spgroup)
levels(big.boxplot.E.means$spgroup)=c("a) All species", "b) Abundant species", "c) Rare species")

ggplot(data=big.boxplot.E) + geom_point(aes(dif, y=tradeoff.ordered, shape=I(21)), fill=big.boxplot.E$colorcode) + facet_wrap(~spgroup, scales="free_x") + xlab("Correspondence index") + ylab("") + theme(legend.position="none", strip.text = element_text(hjust = 0)) + geom_vline(xintercept=0) + geom_point(data=big.boxplot.E.means, aes(average, tradeoff.ordered, size=I(2)))

ggsave("ms figs/fig 6, slopes boxplot with means, E.pdf", height=4, width=7)

# correlation coef of dif comparing allsp to raresp and to commonsp

slopediff=melt(big.boxplot.E, id=c("tradeoff", "tradeoff2", "site_project_comm", "spgroup"), measure=c("realslope", "simulatedslope", "dif"))

simulated.slopes=dcast(slopediff[slopediff$variable=="simulatedslope",], tradeoff+site_project_comm~spgroup)
names(simulated.slopes)[3:5]=c("allsp", "commonsp", "raresp")
qplot(allsp, commonsp, data=simulated.slopes) + geom_abline(slope=1, intercept=0)
qplot(allsp, raresp, data=simulated.slopes) + geom_abline(slope=1, intercept=0)
cor(simulated.slopes[,3:5], use="complete.obs")

dif.in.slopes=dcast(slopediff[slopediff$variable=="dif",], tradeoff+site_project_comm~spgroup)
names(dif.in.slopes)[3:5]=c("allsp", "commonsp", "raresp")
qplot(allsp, commonsp, data=dif.in.slopes) + geom_abline(slope=1, intercept=0)
qplot(allsp, raresp, data=dif.in.slopes) + geom_abline(slope=1, intercept=0)
cor(dif.in.slopes[3:5], use="complete.obs")




