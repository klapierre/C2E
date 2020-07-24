library(tidyverse)
library(grid)
library(PerformanceAnalytics)

#kim's desktop
setwd('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\CoRRE\\converge_diverge\\datasets\\LongForm')


theme_set(theme_bw())
theme_update(axis.title.x=element_text(size=40, vjust=-0.35, margin=margin(t=15)), axis.text.x=element_text(size=34, color='black'),
             axis.title.y=element_text(size=40, angle=90, vjust=0.5, margin=margin(r=15)), axis.text.y=element_text(size=34, color='black'),
             plot.title = element_text(size=40, vjust=2),
             panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
             legend.title=element_blank(), legend.text=element_text(size=20))


####data input
###add treatment info to get binary treatments of various manipulation types
trt <- read.csv('ExperimentInformation_March2019.csv')%>%
  select(-X)%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep='::'))%>%
  select(-site_code, -project_name, -community_type, -public)

###community data
correCommDiff <- read.csv('corre_community_difference_July2019.csv')%>%
  select(-treatment)%>%
  rename(treatment=treatment2)%>%
  rename(treatment_year=time)

###anpp data
correANPP <- read.csv('ANPP_Oct2017.csv')%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep='::'))%>%
  select(-X, -calendar_year)%>%
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

# test <- correANPP%>%filter(is.na(trt_type))%>%select(site_project_comm, treatment, treatment_year)%>%unique()


# #checking site level data for outliers in anpp
# ggplot(correANPP, aes(anpp)) + geom_histogram() + facet_wrap(~project_name, scales='free')


#calculating anpp difference
#anpp ctl data
correANPPctl <- correANPP%>%
  filter(plot_mani==0)%>%
  group_by(site_code, project_name, community_type)%>%
  summarise(anpp_ctl=mean(anpp))%>%
  ungroup()
#anpp change
correANPPdiff <- correANPP%>%
  filter(plot_mani!=0)%>%
  group_by(site_code, project_name, community_type, treatment_year, treatment)%>%
  summarise(anpp=mean(anpp))%>%
  ungroup()%>%
  left_join(correANPPctl)%>%
  #calculate anpp change as percent change from ctl in each year
  mutate(anpp_pdiff=(anpp-anpp_ctl)/anpp_ctl)%>%
  select(-anpp, -anpp_ctl)%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep='::'))

#merge
correSEMdata <- correANPPdiff%>%
  left_join(correCommDiff)%>%
  left_join(trt)%>%
  na.omit()%>%
  #transform to improve normality
  mutate(anpp_pdiff_transform=log10(anpp_pdiff+1-min(anpp_pdiff)), composition_diff_transform=(composition_diff)^(1/3))
#all of these compare each year to year 1!


###exploratory correlations and histograms (all variables compare treatment to control plots)
dataVis <- correSEMdata%>%
  select(anpp_pdiff, composition_diff, richness_difference, evenness_diff, rank_difference, species_difference) #make visualization dataframe
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
correSEMdataTrt <- correSEMdata%>%
  left_join(trt)%>%
  #drop successional and pulse treatments
  filter(successional!=1, pulse!=1)%>%
  #split precip into irrigation and drought
  mutate(irr=ifelse(precip>0, precip, 0), drought=ifelse(precip<0, precip, 0))%>%
  #make binary treatments just in case we want them later
  mutate(n_trt=ifelse(n>0, 1, 0), p_trt=ifelse(p>0, 1, 0), k_trt=ifelse(k>0, 1, 0), CO2_trt=ifelse(CO2>0, 1, 0), irr_trt=ifelse(precip>0, 1, 0), drought_trt=ifelse(precip<0, 1, 0), temp_trt=n_trt+p_trt+k_trt+CO2_trt+irr_trt+drought_trt, other_trt=ifelse((plot_mani-temp_trt)>0, 1, 0))%>%
  na.omit()%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep='_'), site_project_comm_trt=paste(site_code, project_name, community_type, treatment, sep='_'))%>%
  mutate(trt_type_2=ifelse(trt_type %in% c('CO2', 'irr', 'mult_nutrient', 'N', 'N*CO2', 'P', 'N*P', 'N*irr'), 'added', ifelse(trt_type %in% c('drought'), 'removed', ifelse(trt_type %in% c('fungicide', 'mow_clip', 'precip_vari', 'temp'), 'other', 'multiple'))))

rm(list=setdiff(ls(), "correSEMdataTrt"))
# write.csv(correSEMdataTrt, 'CoRRE_comm and anpp diff_07160219.csv')