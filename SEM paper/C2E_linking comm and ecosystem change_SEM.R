library(tidyverse)
library(ggplot2)
library(grid)
library(PerformanceAnalytics)
library(lavaan)
library(lavaan.survey)
library(semPlot)

#kim's desktop
setwd('C:\\Users\\la pierrek\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\datasets\\LongForm')

#kim's laptop
setwd('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\datasets\\LongForm')


theme_set(theme_bw())
theme_update(axis.title.x=element_text(size=40, vjust=-0.35, margin=margin(t=15)), axis.text.x=element_text(size=34, color='black'),
             axis.title.y=element_text(size=40, angle=90, vjust=0.5, margin=margin(r=15)), axis.text.y=element_text(size=34, color='black'),
             plot.title = element_text(size=40, vjust=2),
             panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
             legend.title=element_blank(), legend.text=element_text(size=20))


#######MODEL STRUCTURE: use differences between treatment and control plots in each year (horizontal arrows); include treatments as factors at the bottom, with multi-factor treatments having more than one greater than 0 (and make them categorical response variables)-------------------
####data input
###add treatment info to get binary treatments of various manipulation types
trt <- read.csv('ExperimentInformation_anpp_Dec2017.csv')%>%
  select(-X)%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep='_'))%>%
  select(-site_code, -project_name, -community_type, -public)

###community data
correCommChange <- read.csv('CORRE_ContTreat_Compare_Nov2017.csv')%>%
  select(-X)

###anpp data
correANPP <- read.csv('ANPP_Dec2017.csv')%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep='_'))%>%
  select(-X)%>%
  mutate(treatment_year_2=ifelse(site_project_comm=='SEV_Nfert_0'&calendar_year==2004, 10, ifelse(site_project_comm=='SEV_Nfert_0'&calendar_year==2005, 11, ifelse(site_project_comm=='SEV_Nfert_0'&calendar_year==2006, 12, ifelse(site_project_comm=='SEV_Nfert_0'&calendar_year==2007, 13, ifelse(site_project_comm=='SEV_Nfert_0'&calendar_year==2008, 14, ifelse(site_project_comm=='SEV_Nfert_0'&calendar_year==2009, 15, ifelse(site_project_comm=='SEV_Nfert_0'&calendar_year==2010, 16, ifelse(site_project_comm=='SEV_Nfert_0'&calendar_year==2011, 17, ifelse(site_project_comm=='SEV_Nfert_0'&calendar_year==2012, 18, treatment_year))))))))))%>%
  select(-treatment_year)%>%mutate(treatment_year=treatment_year_2)%>%select(-treatment_year_2)%>%
  #remove subset of CDR trts and KNZ BGP mowed treatment (they herbicided)
  mutate(drop=ifelse(site_code=='CDR'&treatment=='2', 1, ifelse(site_code=='CDR'&treatment=='3', 1, ifelse(site_code=='CDR'&treatment=='4', 1, ifelse(site_code=='CDR'&treatment=='5', 1, ifelse(site_code=='CDR'&treatment=='7', 1, ifelse(site_code=='CDR'&treatment=='2_f_u_n', 1, ifelse(site_code=='CDR'&treatment=='3_f_u_n', 1, ifelse(site_code=='CDR'&treatment=='4_f_u_n', 1, ifelse(site_code=='CDR'&treatment=='5_f_u_n', 1, ifelse(site_code=='CDR'&treatment=='7_f_u_n', 1, ifelse(project_name=='BGP'&treatment=='u_m_n', 1, ifelse(project_name=='BGP'&treatment=='u_m_p', 1, ifelse(project_name=='BGP'&treatment=='u_m_b', 1, ifelse(project_name=='BGP'&treatment=='u_m_c', 1, ifelse(project_name=='BGP'&treatment=='b_m_n', 1, ifelse(project_name=='BGP'&treatment=='b_m_p', 1, ifelse(project_name=='BGP'&treatment=='b_m_b', 1, ifelse(project_name=='BGP'&treatment=='b_m_c', 1, 0)))))))))))))))))))%>%
  #remove NANT wet because it only has ANPP in one year and has a much higher rate of N added (67.2 gm-2)
  filter(site_code!='NANT')%>%
  filter(drop==0)%>%
  ####NOTE: check outliers with data providers and omit this step when data is correct
  mutate(outlier=ifelse(site_code=='ORNL'&anpp>1080, 1, ifelse(project_name=='snow'&anpp>1010, 1, ifelse(project_name=='T7'&anpp>1860, 1, 0))))%>%
  filter(outlier==0)%>%
  filter(anpp>0)%>%
  select(-drop, -outlier)%>%
  left_join(trt)


# #checking site level data for outliers in anpp
# ggplot(correANPP, aes(anpp)) + geom_histogram() + facet_wrap(~project_name, scales='free')


#calculating anpp difference
#anpp ctl data
correANPPctl <- correANPP%>%
  filter(plot_mani==0)%>%
  select(site_code, project_name, community_type, calendar_year, treatment_year, anpp)%>%
  group_by(site_code, project_name, community_type, calendar_year, treatment_year)%>%
  summarise(anpp_ctl=mean(anpp))%>%
  ungroup()
#anpp change
correANPPchange <- correANPP%>%
  filter(plot_mani!=0)%>%
  group_by(site_code, project_name, community_type, calendar_year, treatment_year, treatment)%>%
  summarise(anpp=mean(anpp))%>%
  ungroup()%>%
  left_join(correANPPctl)%>%
  #calculate anpp change as percent change from ctl in each year
  mutate(anpp_PC=(anpp-anpp_ctl)/anpp_ctl)%>%
  select(-anpp, -anpp_ctl)%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep='_'))

#merge
correSEMdata <- correANPPchange%>%
  left_join(correCommChange)%>%
  na.omit()%>%
  #transform to improve normality
  mutate(mean_change_transform=sqrt(mean_change), Ed_transform=log10(as.numeric(Ed)), Sd_transform=log10(as.numeric(Sd)+1), anpp_PC_transform=log10(anpp_PC+1-min(anpp_PC)))
#all of these compare treatment to control!


###exploratory correlations and histograms (all variables compare treatment to control plots)
dataVis <- correSEMdata%>%
  select(anpp_PC_transform, mean_change_transform, Sd_transform, Ed_transform, Rd, spd) #make visualization dataframe
chart.Correlation(dataVis, histogram=T, pch=19)

# #detour: are sites with lower gamma diversity more prone to high Ed?
# siteInfo <- read.csv('SiteExperimentDetails_Dec2016.csv')%>%
#   left_join(correSEMdata)
# with(siteInfo, plot(Ed~rrich)) #yes! but we can transform.



#keep all data
correSEMdataTrt <- correSEMdata%>%
  left_join(trt)%>%
  #drop successional and pulse treatments
  filter(successional!=1, pulse!=1)%>%
  mutate(n_trt=ifelse(n>0, 1, 0), p_trt=ifelse(p>0, 1, 0), k_trt=ifelse(k>0, 1, 0), CO2_trt=ifelse(CO2>0, 1, 0), irr_trt=ifelse(precip>0, 1, 0), drought_trt=ifelse(precip<0, 1, 0), temp_trt=n_trt+p_trt+k_trt+CO2_trt+irr_trt+drought_trt, other_trt=ifelse((plot_mani-temp_trt)>0, 1, 0))%>%
  na.omit()

#get treatment types
trtInteractions <- read.csv('treatment interactions_11152017.csv')%>%
  select(site_code, project_name, community_type, treatment, trt_type)
trtInteractionsANPP <- correSEMdataTrt%>%
  left_join(trtInteractions)

##anpp change through time by trt type
#all yrs
ggplot(data=trtInteractionsANPP, aes(x=treatment_year, y=anpp_PC_transform)) +
  geom_smooth(method='lm', formula=y~x+I(x^2), se=F, aes(color=trt_type)) +
  geom_smooth(method='lm', formula=y~x+I(x^2), size=3, color='black') +
  xlab('Treatment Year') + ylab('log ANPP (%) Difference')
ggplot(data=subset(trtInteractionsANPP, treatment_year<11), aes(x=treatment_year, y=anpp_PC_transform)) +
  geom_smooth(method='lm', formula=y~x+I(x^2), se=F, aes(color=trt_type)) +
  geom_smooth(method='lm', formula=y~x+I(x^2), size=3, color='black') +
  xlab('Treatment Year') + ylab('log ANPP (%) Difference')
#export at 1000 x 800



###anpp change correlations
#with mean_change
anpp_mean_y1 <- ggplot(data=subset(correSEMdataTrt, treatment_year==1), aes(x=mean_change_transform, y=anpp_PC_transform)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('ANPP Difference') + xlab('Community Difference') + ylim(0,0.7) + xlim(0,1)
anpp_mean_y2 <- ggplot(data=subset(correSEMdataTrt, treatment_year==2), aes(x=mean_change_transform, y=anpp_PC_transform)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('ANPP Difference') + xlab('Community Difference') + ylim(0,0.7) + xlim(0,1)
anpp_mean_y3 <- ggplot(data=subset(correSEMdataTrt, treatment_year==3), aes(x=mean_change_transform, y=anpp_PC_transform)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('ANPP Difference') + xlab('Community Difference') + ylim(0,0.7) + xlim(0,1)
anpp_mean_y4 <- ggplot(data=subset(correSEMdataTrt, treatment_year==4), aes(x=mean_change_transform, y=anpp_PC_transform)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('ANPP Difference') + xlab('Community Difference') + ylim(0,0.7) + xlim(0,1)
anpp_mean_y5 <- ggplot(data=subset(correSEMdataTrt, treatment_year==5), aes(x=mean_change_transform, y=anpp_PC_transform)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('ANPP Difference') + xlab('Community Difference') + ylim(0,0.7) + xlim(0,1)
anpp_mean_y6 <- ggplot(data=subset(correSEMdataTrt, treatment_year==6), aes(x=mean_change_transform, y=anpp_PC_transform)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('ANPP Difference') + xlab('Community Difference') + ylim(0,0.7) + xlim(0,1)
anpp_mean_y7 <- ggplot(data=subset(correSEMdataTrt, treatment_year==7), aes(x=mean_change_transform, y=anpp_PC_transform)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('ANPP Difference') + xlab('Community Difference') + ylim(0,0.7) + xlim(0,1)
anpp_mean_y8 <- ggplot(data=subset(correSEMdataTrt, treatment_year==8), aes(x=mean_change_transform, y=anpp_PC_transform)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('ANPP Difference') + xlab('Community Difference') + ylim(0,0.7) + xlim(0,1)
anpp_mean_y9 <- ggplot(data=subset(correSEMdataTrt, treatment_year==9), aes(x=mean_change_transform, y=anpp_PC_transform)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('ANPP Difference') + xlab('Community Difference') + ylim(0,0.7) + xlim(0,1)
anpp_mean_y10 <- ggplot(data=subset(correSEMdataTrt, treatment_year==10), aes(x=mean_change_transform, y=anpp_PC_transform)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('ANPP Difference') + xlab('Community Difference') + ylim(0,0.7) + xlim(0,1)
#with N
anpp_n_y1 <- ggplot(data=subset(correSEMdataTrt, treatment_year==1), aes(x=n, y=anpp_PC_transform)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('ANPP Difference') + xlab('N') + ylim(0,0.7) + xlim(0,60)
anpp_n_y2 <- ggplot(data=subset(correSEMdataTrt, treatment_year==2), aes(x=n, y=anpp_PC_transform)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('ANPP Difference') + xlab('N') + ylim(0,0.7) + xlim(0,60)
anpp_n_y3 <- ggplot(data=subset(correSEMdataTrt, treatment_year==3), aes(x=n, y=anpp_PC_transform)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('ANPP Difference') + xlab('N') + ylim(0,0.7) + xlim(0,60)
anpp_n_y4 <- ggplot(data=subset(correSEMdataTrt, treatment_year==4), aes(x=n, y=anpp_PC_transform)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('ANPP Difference') + xlab('N') + ylim(0,0.7) + xlim(0,60)
anpp_n_y5 <- ggplot(data=subset(correSEMdataTrt, treatment_year==5), aes(x=n, y=anpp_PC_transform)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('ANPP Difference') + xlab('N') + ylim(0,0.7) + xlim(0,60)
anpp_n_y6 <- ggplot(data=subset(correSEMdataTrt, treatment_year==6), aes(x=n, y=anpp_PC_transform)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('ANPP Difference') + xlab('N') + ylim(0,0.7) + xlim(0,60)
anpp_n_y7 <- ggplot(data=subset(correSEMdataTrt, treatment_year==7), aes(x=n, y=anpp_PC_transform)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('ANPP Difference') + xlab('N') + ylim(0,0.7) + xlim(0,60)
anpp_n_y8 <- ggplot(data=subset(correSEMdataTrt, treatment_year==8), aes(x=n, y=anpp_PC_transform)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('ANPP Difference') + xlab('N') + ylim(0,0.7) + xlim(0,60)
anpp_n_y9 <- ggplot(data=subset(correSEMdataTrt, treatment_year==9), aes(x=n, y=anpp_PC_transform)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('ANPP Difference') + xlab('N') + ylim(0,0.7) + xlim(0,60)
anpp_n_y10 <- ggplot(data=subset(correSEMdataTrt, treatment_year==10), aes(x=n, y=anpp_PC_transform)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('ANPP Difference') + xlab('N') + ylim(0,0.7) + xlim(0,60)

pushViewport(viewport(layout=grid.layout(10,2)))
print(anpp_mean_y1, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
print(anpp_mean_y2, vp=viewport(layout.pos.row = 2, layout.pos.col = 1))
print(anpp_mean_y3, vp=viewport(layout.pos.row = 3, layout.pos.col = 1))
print(anpp_mean_y4, vp=viewport(layout.pos.row = 4, layout.pos.col = 1))
print(anpp_mean_y5, vp=viewport(layout.pos.row = 5, layout.pos.col = 1))
print(anpp_mean_y6, vp=viewport(layout.pos.row = 6, layout.pos.col = 1))
print(anpp_mean_y7, vp=viewport(layout.pos.row = 7, layout.pos.col = 1))
print(anpp_mean_y8, vp=viewport(layout.pos.row = 8, layout.pos.col = 1))
print(anpp_mean_y9, vp=viewport(layout.pos.row = 9, layout.pos.col = 1))
print(anpp_mean_y10, vp=viewport(layout.pos.row = 10, layout.pos.col = 1))
print(anpp_n_y1, vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
print(anpp_n_y2, vp=viewport(layout.pos.row = 2, layout.pos.col = 2))
print(anpp_n_y3, vp=viewport(layout.pos.row = 3, layout.pos.col = 2))
print(anpp_n_y4, vp=viewport(layout.pos.row = 4, layout.pos.col = 2))
print(anpp_n_y5, vp=viewport(layout.pos.row = 5, layout.pos.col = 2))
print(anpp_n_y6, vp=viewport(layout.pos.row = 6, layout.pos.col = 2))
print(anpp_n_y7, vp=viewport(layout.pos.row = 7, layout.pos.col = 2))
print(anpp_n_y8, vp=viewport(layout.pos.row = 8, layout.pos.col = 2))
print(anpp_n_y9, vp=viewport(layout.pos.row = 9, layout.pos.col = 2))
print(anpp_n_y10, vp=viewport(layout.pos.row = 10, layout.pos.col = 2))
#export at 2000 x 4000


###Community Change correlations
#with mean_change
mean_Rd_y1 <- ggplot(data=subset(correSEMdataTrt, treatment_year==1), aes(x=Rd, y=mean_change_transform)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
ylab('Community Difference') + xlab('Rd') + ylim(0,1) + xlim(0,0.41)
mean_Rd_y2 <- ggplot(data=subset(correSEMdataTrt, treatment_year==2), aes(x=Rd, y=mean_change_transform)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('Community Difference') + xlab('Rd') + ylim(0,1) + xlim(0,0.41)
mean_Rd_y3 <- ggplot(data=subset(correSEMdataTrt, treatment_year==3), aes(x=Rd, y=mean_change_transform)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('Community Difference') + xlab('Rd') + ylim(0,1) + xlim(0,0.41)
mean_Rd_y4 <- ggplot(data=subset(correSEMdataTrt, treatment_year==4), aes(x=Rd, y=mean_change_transform)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('Community Difference') + xlab('Rd') + ylim(0,1) + xlim(0,0.41)
mean_Rd_y5 <- ggplot(data=subset(correSEMdataTrt, treatment_year==5), aes(x=Rd, y=mean_change_transform)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('Community Difference') + xlab('Rd') + ylim(0,1) + xlim(0,0.41)
mean_Rd_y6 <- ggplot(data=subset(correSEMdataTrt, treatment_year==6), aes(x=Rd, y=mean_change_transform)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('Community Difference') + xlab('Rd') + ylim(0,1) + xlim(0,0.41)
mean_Rd_y7 <- ggplot(data=subset(correSEMdataTrt, treatment_year==7), aes(x=Rd, y=mean_change_transform)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('Community Difference') + xlab('Rd') + ylim(0,1) + xlim(0,0.41)
mean_Rd_y8 <- ggplot(data=subset(correSEMdataTrt, treatment_year==8), aes(x=Rd, y=mean_change_transform)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('Community Difference') + xlab('Rd') + ylim(0,1) + xlim(0,0.41)
mean_Rd_y9 <- ggplot(data=subset(correSEMdataTrt, treatment_year==9), aes(x=Rd, y=mean_change_transform)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('Community Difference') + xlab('Rd') + ylim(0,1) + xlim(0,0.41)
mean_Rd_y10 <- ggplot(data=subset(correSEMdataTrt, treatment_year==10), aes(x=Rd, y=mean_change_transform)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('Community Difference') + xlab('Rd') + ylim(0,1) + xlim(0,0.41)
#with N
mean_n_y1 <- ggplot(data=subset(correSEMdataTrt, treatment_year==1), aes(x=n, y=mean_change_transform)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('Community Difference') + xlab('N') + ylim(0,1) + xlim(0,60)
mean_n_y2 <- ggplot(data=subset(correSEMdataTrt, treatment_year==2), aes(x=n, y=mean_change_transform)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('Community Difference') + xlab('N') + ylim(0,1) + xlim(0,60)
mean_n_y3 <- ggplot(data=subset(correSEMdataTrt, treatment_year==3), aes(x=n, y=mean_change_transform)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('Community Difference') + xlab('N') + ylim(0,1) + xlim(0,60)
mean_n_y4 <- ggplot(data=subset(correSEMdataTrt, treatment_year==4), aes(x=n, y=mean_change_transform)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('Community Difference') + xlab('N') + ylim(0,1) + xlim(0,60)
mean_n_y5 <- ggplot(data=subset(correSEMdataTrt, treatment_year==5), aes(x=n, y=mean_change_transform)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('Community Difference') + xlab('N') + ylim(0,1) + xlim(0,60)
mean_n_y6 <- ggplot(data=subset(correSEMdataTrt, treatment_year==6), aes(x=n, y=mean_change_transform)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('Community Difference') + xlab('N') + ylim(0,1) + xlim(0,60)
mean_n_y7 <- ggplot(data=subset(correSEMdataTrt, treatment_year==7), aes(x=n, y=mean_change_transform)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('Community Difference') + xlab('N') + ylim(0,1) + xlim(0,60)
mean_n_y8 <- ggplot(data=subset(correSEMdataTrt, treatment_year==8), aes(x=n, y=mean_change_transform)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('Community Difference') + xlab('N') + ylim(0,1) + xlim(0,60)
mean_n_y9 <- ggplot(data=subset(correSEMdataTrt, treatment_year==9), aes(x=n, y=mean_change_transform)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('Community Difference') + xlab('N') + ylim(0,1) + xlim(0,60)
mean_n_y10 <- ggplot(data=subset(correSEMdataTrt, treatment_year==10), aes(x=n, y=mean_change_transform)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('Community Difference') + xlab('N') + ylim(0,1) + xlim(0,60)

pushViewport(viewport(layout=grid.layout(10,2)))
print(mean_Rd_y1, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
print(mean_Rd_y2, vp=viewport(layout.pos.row = 2, layout.pos.col = 1))
print(mean_Rd_y3, vp=viewport(layout.pos.row = 3, layout.pos.col = 1))
print(mean_Rd_y4, vp=viewport(layout.pos.row = 4, layout.pos.col = 1))
print(mean_Rd_y5, vp=viewport(layout.pos.row = 5, layout.pos.col = 1))
print(mean_Rd_y6, vp=viewport(layout.pos.row = 6, layout.pos.col = 1))
print(mean_Rd_y7, vp=viewport(layout.pos.row = 7, layout.pos.col = 1))
print(mean_Rd_y8, vp=viewport(layout.pos.row = 8, layout.pos.col = 1))
print(mean_Rd_y9, vp=viewport(layout.pos.row = 9, layout.pos.col = 1))
print(mean_Rd_y10, vp=viewport(layout.pos.row = 10, layout.pos.col = 1))
print(mean_n_y1, vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
print(mean_n_y2, vp=viewport(layout.pos.row = 2, layout.pos.col = 2))
print(mean_n_y3, vp=viewport(layout.pos.row = 3, layout.pos.col = 2))
print(mean_n_y4, vp=viewport(layout.pos.row = 4, layout.pos.col = 2))
print(mean_n_y5, vp=viewport(layout.pos.row = 5, layout.pos.col = 2))
print(mean_n_y6, vp=viewport(layout.pos.row = 6, layout.pos.col = 2))
print(mean_n_y7, vp=viewport(layout.pos.row = 7, layout.pos.col = 2))
print(mean_n_y8, vp=viewport(layout.pos.row = 8, layout.pos.col = 2))
print(mean_n_y9, vp=viewport(layout.pos.row = 9, layout.pos.col = 2))
print(mean_n_y10, vp=viewport(layout.pos.row = 10, layout.pos.col = 2))
#export at 2000 x 4000

###all data-----------
#drought and irrigation are considered negative and positive precip
correModel <-  '
anpp_PC_transform ~ mean_change_transform + n + p + k + CO2 + precip + other_trt
mean_change_transform ~ Sd_transform + Ed_transform + Rd + spd + n + p + k + CO2 + precip + other_trt
Sd_transform ~ n + p + k + CO2 + precip + other_trt
Ed_transform ~ n + p + k + CO2 + precip + other_trt
Rd ~ n + p + k + CO2 + precip + other_trt
spd ~ n + p + k + CO2 + precip + other_trt

#covariances
Sd_transform~~Ed_transform
Sd_transform~~Rd
Sd_transform~~spd
Ed_transform~~Rd
Ed_transform~~spd
Rd~~spd
'
#year 1
correModelFit <- sem(correModel, data=subset(correSEMdataTrt, treatment_year==1), meanstructure=TRUE)
surveyDesign <- svydesign(ids=~site_project_comm, nest=TRUE, data=subset(correSEMdataTrt, treatment_year==1))
survey1Fit <- lavaan.survey(lavaan.fit=correModelFit, survey.design=surveyDesign)
survey1Fit  #gives chi-square
# summary(survey1Fit, fit.measures=TRUE, rsq=TRUE, standardized=TRUE)
stdEst1 <- parameterEstimates(survey1Fit, standardized=TRUE)%>%
  mutate(treatment_year=1)

#year 2
correModel2Fit <- sem(correModel, data=subset(correSEMdataTrt, treatment_year==2), meanstructure=TRUE)
surveyDesign <- svydesign(ids=~site_project_comm, nest=TRUE, data=subset(correSEMdataTrt, treatment_year==2))
survey2Fit <- lavaan.survey(lavaan.fit=correModel2Fit, survey.design=surveyDesign)
survey2Fit  #gives chi-square
# summary(survey2Fit, fit.measures=TRUE, rsq=TRUE, standardized=TRUE)
stdEst2 <- parameterEstimates(survey2Fit, standardized=TRUE)%>%
  mutate(treatment_year=2)

#year 3
correModel3Fit <- sem(correModel, data=subset(correSEMdataTrt, treatment_year==3), meanstructure=TRUE)
surveyDesign <- svydesign(ids=~site_project_comm, nest=TRUE, data=subset(correSEMdataTrt, treatment_year==3))
survey3Fit <- lavaan.survey(lavaan.fit=correModel3Fit, survey.design=surveyDesign)
survey3Fit  #gives chi-square
# summary(survey3Fit, fit.measures=TRUE, rsq=TRUE, standardized=TRUE)
stdEst3 <- parameterEstimates(survey3Fit, standardized=TRUE)%>%
  mutate(treatment_year=3)

#year 4
correModel4Fit <- sem(correModel, data=subset(correSEMdataTrt, treatment_year==4), meanstructure=TRUE)
surveyDesign <- svydesign(ids=~site_project_comm, nest=TRUE, data=subset(correSEMdataTrt, treatment_year==4))
survey4Fit <- lavaan.survey(lavaan.fit=correModel4Fit, survey.design=surveyDesign)
survey4Fit  #gives chi-square
# summary(survey4Fit, fit.measures=TRUE, rsq=TRUE, standardized=TRUE)
stdEst4 <- parameterEstimates(survey4Fit, standardized=TRUE)%>%
  mutate(treatment_year=4)

#year 5
correModel5Fit <- sem(correModel, data=subset(correSEMdataTrt, treatment_year==5), meanstructure=TRUE)
surveyDesign <- svydesign(ids=~site_project_comm, nest=TRUE, data=subset(correSEMdataTrt, treatment_year==5))
survey5Fit <- lavaan.survey(lavaan.fit=correModel5Fit, survey.design=surveyDesign)
survey5Fit  #gives chi-square
# summary(survey5Fit, fit.measures=TRUE, rsq=TRUE, standardized=TRUE)
stdEst5 <- parameterEstimates(survey5Fit, standardized=TRUE)%>%
  mutate(treatment_year=5)

#year 6
correModel6Fit <- sem(correModel, data=subset(correSEMdataTrt, treatment_year==6), meanstructure=TRUE)
surveyDesign <- svydesign(ids=~site_project_comm, nest=TRUE, data=subset(correSEMdataTrt, treatment_year==6))
survey6Fit <- lavaan.survey(lavaan.fit=correModel6Fit, survey.design=surveyDesign)
survey6Fit  #gives chi-square
# summary(survey6Fit, fit.measures=TRUE, rsq=TRUE, standardized=TRUE)
stdEst6 <- parameterEstimates(survey6Fit, standardized=TRUE)%>%
  mutate(treatment_year=6)

#year 7
correModel7Fit <- sem(correModel, data=subset(correSEMdataTrt, treatment_year==7), meanstructure=TRUE)
surveyDesign <- svydesign(ids=~site_project_comm, nest=TRUE, data=subset(correSEMdataTrt, treatment_year==7))
survey7Fit <- lavaan.survey(lavaan.fit=correModel7Fit, survey.design=surveyDesign)
survey7Fit  #gives chi-square
# summary(survey7Fit, fit.measures=TRUE, rsq=TRUE, standardized=TRUE)
stdEst7 <- parameterEstimates(survey7Fit, standardized=TRUE)%>%
  mutate(treatment_year=7)

#year 8
correModel8Fit <- sem(correModel, data=subset(correSEMdataTrt, treatment_year==8), meanstructure=TRUE)
surveyDesign <- svydesign(ids=~site_project_comm, nest=TRUE, data=subset(correSEMdataTrt, treatment_year==8))
survey8Fit <- lavaan.survey(lavaan.fit=correModel8Fit, survey.design=surveyDesign)
survey8Fit  #gives chi-square
# summary(survey8Fit, fit.measures=TRUE, rsq=TRUE, standardized=TRUE)
stdEst8 <- parameterEstimates(survey8Fit, standardized=TRUE)%>%
  mutate(treatment_year=8)

#year 9
correModel9Fit <- sem(correModel, data=subset(correSEMdataTrt, treatment_year==9), meanstructure=TRUE)
surveyDesign <- svydesign(ids=~site_project_comm, nest=TRUE, data=subset(correSEMdataTrt, treatment_year==9))
survey9Fit <- lavaan.survey(lavaan.fit=correModel9Fit, survey.design=surveyDesign)
survey9Fit  #gives chi-square
# summary(survey9Fit, fit.measures=TRUE, rsq=TRUE, standardized=TRUE)
stdEst9 <- parameterEstimates(survey9Fit, standardized=TRUE)%>%
  mutate(treatment_year=9)

#year 10
correModel10Fit <- sem(correModel, data=subset(correSEMdataTrt, treatment_year==10), meanstructure=TRUE)
surveyDesign <- svydesign(ids=~site_project_comm, nest=TRUE, data=subset(correSEMdataTrt, treatment_year==10))
survey10Fit <- lavaan.survey(lavaan.fit=correModel10Fit, survey.design=surveyDesign)
survey10Fit  #gives chi-square
# summary(survey10Fit, fit.measures=TRUE, rsq=TRUE, standardized=TRUE)
stdEst10 <- parameterEstimates(survey10Fit, standardized=TRUE)%>%
  mutate(treatment_year=10)

#bind together output from years 1-10
stdEst <- rbind(stdEst1, stdEst2, stdEst3, stdEst4, stdEst5, stdEst6, stdEst7, stdEst8, stdEst9, stdEst10)%>%
  na.omit()%>%
  mutate(std_alt=ifelse(pvalue>0.1, 0, std.all), std_weighted=std.all/se, est_alt=ifelse(pvalue>0.1, 0, est), est_weighted=est/se, trt_alt=ifelse(rhs=='CO2'|rhs=='precip'|rhs=='k'|rhs=='n'|rhs=='other_trt'|rhs=='p', 'physiology', 'mean_change_transform'))

#all years figs-----------
#std.all includes non-sig effects
#std_alt sets non-sig effects to 0
#est is non-standardized effects
#est_alt sets non-sig non-standardized effects to 0
ggplot(data=subset(stdEst, lhs=='anpp_PC_transform'&rhs!=''&rhs!='anpp_PC_transform'), aes(x=treatment_year, y=std.all)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=std.all-se, ymax=std.all+se), width=0.5) +
  geom_line(size=2) +
  # geom_smooth(method='lm', size=3, se=F) +
  ylab('ANPP Difference Effect Size') +
  xlab('Treatment Year') +
  geom_hline(yintercept=0) +
  facet_wrap(~rhs) +
  theme(strip.text=element_text(size=24))

ggplot(data=subset(stdEst, lhs=='mean_change_transform'&rhs!=''&rhs!='anpp_PC_transform'&rhs!='mean_change_transform'), aes(x=treatment_year, y=std.all)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=std.all-se, ymax=std.all+se), width=0.2) +
  geom_line(size=2) +
  ylab('Community Difference Effect Size') +
  xlab('Treatment Year') +
  geom_hline(yintercept=0) +
  facet_wrap(~rhs) +
  theme(strip.text=element_text(size=24))

ggplot(data=subset(stdEst, lhs!='anpp_PC_transform'&lhs!='mean_change_transform'&rhs!=''&rhs!='anpp_PC_transform'&rhs!='mean_change_transform'&rhs!='Sd_transform'&rhs!='Ed_transform'&rhs!='Rd'&rhs!='spd'), aes(x=treatment_year, y=std.all, color=rhs)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=std.all-se, ymax=std.all+se)) +
  geom_line(size=2) +
  ylab('Effect Size') +
  xlab('Treatment Year') +
  geom_hline(yintercept=0) +
  facet_wrap(~lhs) +
  theme(strip.text=element_text(size=24))

# #bin into first and last 5 years for bar graphs
# stdEstBin <- stdEst%>%
#   mutate(temporal_bin=ifelse(treatment_year>5, '1-5', '6-10'))%>%
#   group_by(lhs, rhs, temporal_bin)%>%
#   summarise(est_mean=mean(est), est_sd=sd(est), est_N=n(), std.all_mean=mean(std.all), std.all_sd=sd(std.all), std.all_N=n())%>%
#   ungroup()
# 
# ggplot(data=subset(stdEstBin, lhs=='anpp_PC_transform'&rhs!=''&rhs!='anpp_PC_transform'), aes(x=temporal_bin, y=std.all_mean, color=rhs)) +
#   geom_bar(stat='identity', position=position_dodge()) +
#   geom_errorbar(aes(ymin=std.all_mean-(std.all_sd/sqrt(std.all_N)), ymax=std.all_mean+(std.all_sd/sqrt(std.all_N))), position=position_dodge(0.9), width=0.2) +
#   ylab('Effect Size') +
#   xlab('Treatment Year')
# 
# ggplot(data=subset(stdEstBin, lhs=='mean_change_transform'&rhs!=''&rhs!='anpp_PC_transform'&rhs!='mean_change_transform'), aes(x=temporal_bin, y=std.all_mean, color=rhs)) +
#   geom_bar(stat='identity', position=position_dodge()) +
#   geom_errorbar(aes(ymin=std.all_mean-(std.all_sd/sqrt(std.all_N)), ymax=std.all_mean+(std.all_sd/sqrt(std.all_N))), position=position_dodge(0.9), width=0.2) +
#   geom_line(size=3) +
#   ylab('Effect Size') +
#   xlab('Treatment Year')
# 
# ggplot(data=subset(stdEstBin, lhs!='anpp_PC_transform'&lhs!='mean_change_transform'&rhs!=''&rhs!='anpp_PC_transform'&rhs!='mean_change_transform'&rhs!='Sd_transform'&rhs!='Ed_transform'&rhs!='Rd'&rhs!='spd'), aes(x=temporal_bin, y=std.all_mean, color=rhs)) +
#   geom_bar(stat='identity', position=position_dodge()) +
#   geom_errorbar(aes(ymin=std.all_mean-(std.all_sd/sqrt(std.all_N)), ymax=std.all_mean+(std.all_sd/sqrt(std.all_N))), position=position_dodge(0.9), width=0.2) +
#   geom_line(size=3) +
#   facet_wrap(~lhs)







###only 10 year datasets-----------
#subset out anything without 10 years
numYears <- correSEMdata%>%
  select(site_code, project_name, community_type, treatment, treatment_year)%>%
  group_by(site_code, project_name, community_type, treatment)%>%
  summarise(num_years=length(treatment_year))

#subset out anything without years 1-10 of data
correSEMdataTrt10 <- correSEMdata%>%
  left_join(trt)%>%
  mutate(n_trt=ifelse(n>0, 1, 0), p_trt=ifelse(p>0, 1, 0), k_trt=ifelse(k>0, 1, 0), CO2_trt=ifelse(CO2>0, 1, 0), irr_trt=ifelse(precip>0, 1, 0), drought_trt=ifelse(precip<0, 1, 0), temp_trt=n_trt+p_trt+k_trt+CO2_trt, other_trt=ifelse((plot_mani-temp_trt)>0, 1, 0))%>%
  left_join(numYears)%>%
  filter(num_years>9&project_name!='TMECE'&project_name!='RaMPs')%>%
  na.omit()

#precip becomes an "other_trt" because only 2 trts manipulate precip (knz irg upland and lowland)

correModel <-  '
anpp_PC_transform ~ mean_change_transform + n + p + k + CO2 + other_trt
mean_change_transform ~ Sd_transform + Ed_transform + Rd + spd + n + p + k + CO2 + other_trt
Sd_transform ~ n + p + k + CO2 + other_trt
Ed_transform ~ n + p + k + CO2 + other_trt
Rd ~ n + p + k + CO2 + other_trt
spd ~ n + p + k + CO2 + other_trt

#covariances
Sd_transform~~Ed_transform
Sd_transform~~Rd
Sd_transform~~spd
Ed_transform~~Rd
Ed_transform~~spd
Rd~~spd
'
#year 1
correModelFit <- sem(correModel, data=subset(correSEMdataTrt10, treatment_year==1), meanstructure=TRUE)
surveyDesign <- svydesign(ids=~site_project_comm, nest=TRUE, data=subset(correSEMdataTrt10, treatment_year==1))
survey1Fit <- lavaan.survey(lavaan.fit=correModelFit, survey.design=surveyDesign)
survey1Fit  #gives chi-square
# summary(survey1Fit, fit.measures=TRUE, rsq=TRUE, standardized=TRUE)
stdEst1 <- parameterEstimates(survey1Fit, standardized=TRUE)%>%
  mutate(treatment_year=1)

#year 2
correModel2Fit <- sem(correModel, data=subset(correSEMdataTrt10, treatment_year==2), meanstructure=TRUE)
surveyDesign <- svydesign(ids=~site_project_comm, nest=TRUE, data=subset(correSEMdataTrt10, treatment_year==2))
survey2Fit <- lavaan.survey(lavaan.fit=correModel2Fit, survey.design=surveyDesign)
survey2Fit  #gives chi-square
# summary(survey2Fit, fit.measures=TRUE, rsq=TRUE, standardized=TRUE)
stdEst2 <- parameterEstimates(survey2Fit, standardized=TRUE)%>%
  mutate(treatment_year=2)

#year 3
correModel3Fit <- sem(correModel, data=subset(correSEMdataTrt10, treatment_year==3), meanstructure=TRUE)
surveyDesign <- svydesign(ids=~site_project_comm, nest=TRUE, data=subset(correSEMdataTrt10, treatment_year==3))
survey3Fit <- lavaan.survey(lavaan.fit=correModel3Fit, survey.design=surveyDesign)
survey3Fit  #gives chi-square
# summary(survey3Fit, fit.measures=TRUE, rsq=TRUE, standardized=TRUE)
stdEst3 <- parameterEstimates(survey3Fit, standardized=TRUE)%>%
  mutate(treatment_year=3)

#year 4
correModel4Fit <- sem(correModel, data=subset(correSEMdataTrt10, treatment_year==4), meanstructure=TRUE)
surveyDesign <- svydesign(ids=~site_project_comm, nest=TRUE, data=subset(correSEMdataTrt10, treatment_year==4))
survey4Fit <- lavaan.survey(lavaan.fit=correModel4Fit, survey.design=surveyDesign)
survey4Fit  #gives chi-square
# summary(survey4Fit, fit.measures=TRUE, rsq=TRUE, standardized=TRUE)
stdEst4 <- parameterEstimates(survey4Fit, standardized=TRUE)%>%
  mutate(treatment_year=4)

#year 5
correModel5Fit <- sem(correModel, data=subset(correSEMdataTrt10, treatment_year==5), meanstructure=TRUE)
surveyDesign <- svydesign(ids=~site_project_comm, nest=TRUE, data=subset(correSEMdataTrt10, treatment_year==5))
survey5Fit <- lavaan.survey(lavaan.fit=correModel5Fit, survey.design=surveyDesign)
survey5Fit  #gives chi-square
# summary(survey5Fit, fit.measures=TRUE, rsq=TRUE, standardized=TRUE)
stdEst5 <- parameterEstimates(survey5Fit, standardized=TRUE)%>%
  mutate(treatment_year=5)

#year 6
correModel6Fit <- sem(correModel, data=subset(correSEMdataTrt10, treatment_year==6), meanstructure=TRUE)
surveyDesign <- svydesign(ids=~site_project_comm, nest=TRUE, data=subset(correSEMdataTrt10, treatment_year==6))
survey6Fit <- lavaan.survey(lavaan.fit=correModel6Fit, survey.design=surveyDesign)
survey6Fit  #gives chi-square
# summary(survey6Fit, fit.measures=TRUE, rsq=TRUE, standardized=TRUE)
stdEst6 <- parameterEstimates(survey6Fit, standardized=TRUE)%>%
  mutate(treatment_year=6)

#year 7
correModel7Fit <- sem(correModel, data=subset(correSEMdataTrt10, treatment_year==7), meanstructure=TRUE)
surveyDesign <- svydesign(ids=~site_project_comm, nest=TRUE, data=subset(correSEMdataTrt10, treatment_year==7))
survey7Fit <- lavaan.survey(lavaan.fit=correModel7Fit, survey.design=surveyDesign)
survey7Fit  #gives chi-square
# summary(survey7Fit, fit.measures=TRUE, rsq=TRUE, standardized=TRUE)
stdEst7 <- parameterEstimates(survey7Fit, standardized=TRUE)%>%
  mutate(treatment_year=7)

#year 8
correModel8Fit <- sem(correModel, data=subset(correSEMdataTrt10, treatment_year==8), meanstructure=TRUE)
surveyDesign <- svydesign(ids=~site_project_comm, nest=TRUE, data=subset(correSEMdataTrt10, treatment_year==8))
survey8Fit <- lavaan.survey(lavaan.fit=correModel8Fit, survey.design=surveyDesign)
survey8Fit  #gives chi-square
# summary(survey8Fit, fit.measures=TRUE, rsq=TRUE, standardized=TRUE)
stdEst8 <- parameterEstimates(survey8Fit, standardized=TRUE)%>%
  mutate(treatment_year=8)

#year 9
correModel9Fit <- sem(correModel, data=subset(correSEMdataTrt10, treatment_year==9), meanstructure=TRUE)
surveyDesign <- svydesign(ids=~site_project_comm, nest=TRUE, data=subset(correSEMdataTrt10, treatment_year==9))
survey9Fit <- lavaan.survey(lavaan.fit=correModel9Fit, survey.design=surveyDesign)
survey9Fit  #gives chi-square
# summary(survey9Fit, fit.measures=TRUE, rsq=TRUE, standardized=TRUE)
stdEst9 <- parameterEstimates(survey9Fit, standardized=TRUE)%>%
  mutate(treatment_year=9)

#year 10
correModel10Fit <- sem(correModel, data=subset(correSEMdataTrt10, treatment_year==10), meanstructure=TRUE)
surveyDesign <- svydesign(ids=~site_project_comm, nest=TRUE, data=subset(correSEMdataTrt10, treatment_year==10))
survey10Fit <- lavaan.survey(lavaan.fit=correModel10Fit, survey.design=surveyDesign)
survey10Fit  #gives chi-square
# summary(survey10Fit, fit.measures=TRUE, rsq=TRUE, standardized=TRUE)
stdEst10 <- parameterEstimates(survey10Fit, standardized=TRUE)%>%
  mutate(treatment_year=10)

#bind together output from years 1-10
stdEst <- rbind(stdEst1, stdEst2, stdEst3, stdEst4, stdEst5, stdEst6, stdEst7, stdEst8, stdEst9, stdEst10)%>%
  na.omit()%>%
  mutate(std_alt=ifelse(pvalue>0.05, 0, std.all))

#10 year figs-----------
ggplot(data=subset(stdEst, lhs=='anpp_PC_transform'&rhs!=''&rhs!='anpp_PC_transform'), aes(x=treatment_year, y=std.all)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=std.all-se, ymax=std.all+se)) +
  geom_line(size=2) +
  ylab('ANPP Difference Effect Size') + xlab('Treatment Year') +
  geom_hline(yintercept=0) +
  facet_wrap(~rhs)

ggplot(data=subset(stdEst, lhs=='mean_change_transform'&rhs!=''&rhs!='anpp_PC_transform'&rhs!='mean_change_transform'), aes(x=treatment_year, y=std.all)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=std.all-se, ymax=std.all+se)) +
  geom_line(size=2) +
  ylab('Community Difference Effect Size') + xlab('Treatment Year') +
  geom_hline(yintercept=0) +
  facet_wrap(~rhs)

ggplot(data=subset(stdEst, lhs!='anpp_PC_transform'&lhs!='mean_change_transform'&rhs!=''&rhs!='anpp_PC_transform'&rhs!='mean_change_transform'&rhs!='Rd'&rhs!='Ed_transform'&rhs!='Sd_transform'&rhs!='spd'), aes(x=treatment_year, y=std.all, color=rhs)) +
  geom_point(size=3) +
  geom_line(size=2) +
  ylab('Effect Size') + xlab('Treatment Year') +
  geom_hline(yintercept=0) +
  facet_wrap(~lhs)
