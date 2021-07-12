################################################################################
##  C2E_data formatting.R: Compiling and formatting data for analysis.
##
##  Author: Kimberly Komatsu
##  Date created: July 12, 2021
################################################################################

library(codyn)
library(grid)
library(PerformanceAnalytics)
library(tidyverse)

#kim's laptop
setwd('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\CoRRE\\CoRRE_database\\Data\\CompiledData')


##### import community data #####
trt <- read.csv('ExperimentInfo.csv')%>%
  select(-X)%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep='::'))

relCover <- read.csv('RelativeCover.csv')%>%
  select(-X)%>%
  filter(treatment_year!=0)%>%
  filter(site_code!='ANR')%>% #drop ANR fert1, repeat entries for spp and treatmnet_year wrong
  group_by(site_code, project_name, community_type, treatment_year, calendar_year, treatment, plot_id, genus_species)%>%
  summarise(relcov=mean(relcov))%>% #need to fix duplicate entires in several experiments (AZI EELplot; Bt DrougthNet; CDR BioCON; CDR e001; NIN HerbDiv; Sil NASH; SIU TON)
  ungroup()%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep='::'))%>%
  left_join(trt)%>%
  mutate(plot_mani=ifelse(site_project_comm %in% c('Bt::DroughtNet::0', 'Glen::Fert::0') & treatment %in% c('control', 'C'), 0,
                          ifelse(site_project_comm %in% c('Bt::NPKDNet::0', 'Glen::Fert::0') & treatment %in% c('drought*fert', 'NP'), 2,
                                 ifelse(site_project_comm %in% c('Bt::DroughtNet::0', 'Bt::NPKDNet::0', 'Glen::Fert::0', 'NGBER::gb::0') & treatment %in% c('drought', 'fert', 'N', 'P', 'AMBIENT'), 1, plot_mani))))%>%
  filter(site_project_comm!='Bt::NPKDNet::0', site_project_comm!='HAYS::Precip::0', site_project_comm!='KUFS::E2::0', site_project_comm!='SORBAS::CLIMARID::0')%>% #experiment needs controls repeated in data; trt names messed up
  filter(site_project_comm!='DCGS::gap::0')%>% #pulse experiment
  mutate(treatment_year=ifelse(site_project_comm=='DL::NSFC::0'&calendar_year==2013, 9,
                               ifelse(site_project_comm=='DL::NSFC::0'&calendar_year==2014, 10,
                                      ifelse(site_project_comm=='DL::NSFC::0'&calendar_year==2015, 11,
                                             ifelse(site_project_comm=='DL::NSFC::0'&calendar_year==2016, 12,
                                                    treatment_year)))))


##### calculating community differences #####
#make a new dataframe with just the label
exp_year=relCover%>%
  select(site_project_comm)%>%
  unique()

#makes an empty dataframe
correCommDiff=data.frame(row.names=1) 

for(i in 1:length(exp_year$site_project_comm)) {
  
  #creates a dataset for each unique year, trt, exp combo
  subset <- relCover[relCover$site_project_comm==as.character(exp_year$site_project_comm[i]),]%>%
    select(site_project_comm, treatment_year, treatment, plot_mani, genus_species, relcov, plot_id)%>%
    mutate(treatment2=ifelse(plot_mani==0, 'TRUECONTROL', as.character(treatment)))
  
  #need this to keep track of treatment
  treatments <- subset%>%
    select(plot_id, treatment)%>%
    unique()
  
  #calculating composition difference and abs(dispersion difference)
  multivariate <- multivariate_difference(subset, time.var = 'treatment_year', species.var = "genus_species", abundance.var = 'relcov', replicate.var = 'plot_id', treatment.var='treatment2', reference.treatment='TRUECONTROL')%>%
    rename(treatment=treatment22)%>%
    select(-treatment2, -trt_greater_disp)
  
  #calculating univariate community metrics for each plot
  univariate <- RAC_difference(subset, time.var = 'treatment_year', species.var = "genus_species", abundance.var = 'relcov', replicate.var = 'plot_id', treatment.var='treatment2', reference.treatment='TRUECONTROL')%>%
    rename(treatment=treatment22)%>%
    select(-treatment2, -plot_id2)%>%
    group_by(treatment_year, treatment)%>%
    summarise(richness_diff=mean(richness_diff), evenness_diff=mean(evenness_diff), rank_diff=mean(rank_diff), species_diff=mean(species_diff))%>%
    ungroup()
  
  #merge multivariate and univariate metrics
  all <- univariate%>%
    full_join(multivariate)%>%
    mutate(site_project_comm=exp_year$site_project_comm[i])
  
  #pasting dispersions into the dataframe made for this analysis
  correCommDiff=rbind(all, correCommDiff)  
}



##### calculating community changes #####
#makes an empty dataframe
correCommChange=data.frame(row.names=1) 

exp_year2 <- subset(exp_year, site_project_comm!='GVN::FACE::0')

for(i in 1:length(exp_year2$site_project_comm)) {
  
  #creates a dataset for each unique year, trt, exp combo
  subset <- relCover[relCover$site_project_comm==as.character(exp_year2$site_project_comm[i]),]%>%
    select(site_project_comm, treatment_year, treatment, plot_mani, genus_species, relcov, plot_id)%>%
    mutate(treatment2=ifelse(plot_mani==0, 'TRUECONTROL', as.character(treatment)))
  
  #need this to keep track of treatment
  treatments <- subset%>%
    select(plot_id, treatment)%>%
    unique()
  
  referenceTime=min(subset$treatment_year)
  
  #calculating composition difference and abs(dispersion difference)
  multivariate <- multivariate_change(subset, time.var = 'treatment_year', species.var = "genus_species", abundance.var = 'relcov', replicate.var = 'plot_id', treatment.var='treatment2', reference.time=referenceTime)%>%
    select(-treatment_year)%>%
    rename(treatment_year=treatment_year2,
           treatment=treatment2)
  
  #calculating univariate community metrics for each plot
  univariate <- RAC_change(subset, time.var = 'treatment_year', species.var = "genus_species", abundance.var = 'relcov', replicate.var = 'plot_id', reference.time=referenceTime)%>%
    select(-treatment_year)%>%
    rename(treatment_year=treatment_year2)%>%
    left_join(treatments)%>%
    group_by(treatment_year, treatment)%>%
    summarise(richness_change=mean(richness_change), evenness_change=mean(evenness_change), rank_change=mean(rank_change), gains=mean(gains), losses=mean(losses))%>%
    ungroup()
  
  #merge multivariate and univariate metrics
  all <- univariate%>%
    full_join(multivariate)%>%
    mutate(site_project_comm=exp_year2$site_project_comm[i])
  
  #pasting dispersions into the dataframe made for this analysis
  correCommChange=rbind(all, correCommChange)  
}



##### import anpp data #####
correANPP <- read.csv('ANPP2020.csv')%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep='::'))%>%
  select(-X, -calendar_year)%>%
  filter(!is.na(anpp))%>%
  left_join(trt)%>%
  mutate(treatment_year_2=ifelse(site_project_comm=='SEV::Nfert::0'&calendar_year==2004, 10, ifelse(site_project_comm=='SEV::Nfert::0'&calendar_year==2005, 11, ifelse(site_project_comm=='SEV::Nfert::0'&calendar_year==2006, 12, ifelse(site_project_comm=='SEV::Nfert::0'&calendar_year==2007, 13, ifelse(site_project_comm=='SEV::Nfert::0'&calendar_year==2008, 14, ifelse(site_project_comm=='SEV::Nfert::0'&calendar_year==2009, 15, ifelse(site_project_comm=='SEV::Nfert::0'&calendar_year==2010, 16, ifelse(site_project_comm=='SEV::Nfert::0'&calendar_year==2011, 17, ifelse(site_project_comm=='SEV::Nfert::0'&calendar_year==2012, 18, treatment_year))))))))))%>%
  select(-treatment_year)%>%mutate(treatment_year=treatment_year_2)%>%select(-treatment_year_2)%>%
  #remove subset of CDR trts and KNZ BGP mowed treatment (they herbicided)
  mutate(drop=ifelse(site_code=='CDR'&treatment=='2', 1, ifelse(site_code=='CDR'&treatment=='3', 1, ifelse(site_code=='CDR'&treatment=='4', 1, ifelse(site_code=='CDR'&treatment=='5', 1, ifelse(site_code=='CDR'&treatment=='7', 1, ifelse(site_code=='CDR'&treatment=='8', 1, ifelse(site_code=='CDR'&treatment=='2_f_u_n', 1, ifelse(site_code=='CDR'&treatment=='3_f_u_n', 1, ifelse(site_code=='CDR'&treatment=='4_f_u_n', 1, ifelse(site_code=='CDR'&treatment=='5_f_u_n', 1, ifelse(site_code=='CDR'&treatment=='7_f_u_n', 1, ifelse(site_code=='CDR'&treatment=='8_f_u_n', 1, ifelse(project_name=='BGP'&treatment=='u_m_n', 1, ifelse(project_name=='BGP'&treatment=='u_m_p', 1, ifelse(project_name=='BGP'&treatment=='u_m_b', 1, ifelse(project_name=='BGP'&treatment=='u_m_c', 1, ifelse(project_name=='BGP'&treatment=='b_m_n', 1, ifelse(project_name=='BGP'&treatment=='b_m_p', 1, ifelse(project_name=='BGP'&treatment=='b_m_b', 1, ifelse(project_name=='BGP'&treatment=='b_m_c', 1, 0)))))))))))))))))))))%>%
  #remove NANT wet because it only has ANPP in one year and has a much higher rate of N added (67.2 gm-2)
  filter(site_code!='NANT')%>%
  filter(drop==0)%>%
  ####NOTE: check outliers with data providers and omit this step when data is correct
  mutate(outlier=ifelse(site_code=='ORNL'&anpp>1080, 1, ifelse(project_name=='snow'&anpp>1010, 1, ifelse(project_name=='T7'&anpp>1860, 1, 0))))%>%
  filter(outlier==0)%>%
  filter(anpp>0)%>%
  select(-drop, -outlier)


# #checking site level data for outliers in anpp
# ggplot(correANPP, aes(anpp)) + geom_histogram() + facet_wrap(~project_name, scales='free')


##### calculating percent anpp difference #####
#anpp ctl data
correANPPctl <- correANPP%>%
  filter(plot_mani==0)%>%
  group_by(site_code, project_name, community_type)%>%
  summarise(anpp_ctl=mean(anpp))%>%
  ungroup()
#anpp difference
correANPPdiff <- correANPP%>%
  filter(plot_mani!=0)%>%
  group_by(site_code, project_name, community_type, treatment_year, treatment)%>%
  summarise(anpp=mean(anpp))%>%
  ungroup()%>%
  left_join(correANPPctl)%>%
  #calculate anpp diff as percent diff from ctl in each year
  mutate(anpp_pdiff=(anpp-anpp_ctl)/anpp_ctl)%>%
  select(-anpp, -anpp_ctl)%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep='::'))

#merge
correAllData <- correCommDiff%>%
  full_join(correCommChange)%>%
  filter(!is.na(richness_diff))%>%
  mutate(evenness_diff=ifelse(is.na(evenness_diff), 0, evenness_diff),
         evenness_change=ifelse(is.na(evenness_change), 0, evenness_change))%>%
  full_join(correANPPdiff)%>%
  filter(!is.na(anpp_pdiff))%>%
  left_join(trt)%>%
  #transform to improve normality
  mutate(anpp_pdiff_transform=log10(anpp_pdiff+1-min(anpp_pdiff)), composition_diff_transform=(composition_diff)^(1/3))


##### exploratory correlations and histograms (all variables compare treatment to control plots) #####
#difference plots
dataVis <- correAllData%>%
  select(anpp_pdiff, composition_diff, richness_diff, evenness_diff, rank_diff, species_diff)
chart.Correlation(dataVis, histogram=T, pch=19)

#change plots
dataVis <- correAllData%>%
  select(anpp_pdiff, composition_change, richness_change, evenness_change, rank_change, gains, losses) #make visualization dataframe
chart.Correlation(dataVis, histogram=T, pch=19)

# qqnorm(correSEMdata$anpp_pdiff_transform)
# qqline(correSEMdata$anpp_pdiff_transform, col="red")
# 
# qqnorm(correSEMdata$composition_diff_transform)
# qqline(correSEMdata$composition_diff_transform, col="red")



# #detour: are sites with lower gamma diversity more prone to high Ed?
# siteInfo <- read.csv('SiteExperimentDetails_Dec2016.csv')%>%
#   left_join(correSEMdata)
# with(siteInfo, plot(Ed~rrich)) #yes! but we can transform.



#keep all data
correAllDataTrt <- correAllData%>%
  left_join(trt)%>%
  #drop successional and pulse treatments
  filter(successional!=1, pulse!=1)%>%
  #split precip into irrigation and drought
  mutate(irr=ifelse(precip>0, precip, 0), drought=ifelse(precip<0, precip, 0))%>%
  #make binary treatments just in case we want them later
  mutate(n_trt=ifelse(n>0, 1, 0), p_trt=ifelse(p>0, 1, 0), k_trt=ifelse(k>0, 1, 0), CO2_trt=ifelse(CO2>0, 1, 0), irr_trt=ifelse(precip>0, 1, 0), drought_trt=ifelse(precip<0, 1, 0), temp_trt=n_trt+p_trt+k_trt+CO2_trt+irr_trt+drought_trt, other_trt=ifelse((plot_mani-temp_trt)>0, 1, 0))%>%
  # na.omit()%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep='_'), site_project_comm_trt=paste(site_code, project_name, community_type, treatment, sep='_'))%>%
  mutate(trt_type_2=ifelse(trt_type %in% c('CO2', 'irr', 'mult_nutrient', 'N', 'N*CO2', 'P', 'N*P', 'N*irr'), 'added', ifelse(trt_type %in% c('drought'), 'removed', ifelse(trt_type %in% c('fungicide', 'mow_clip', 'precip_vari', 'temp'), 'other', 'multiple'))))

# write.csv(correAllDataTrt, 'CoRRE_comm and anpp diff_07122021.csv')