#filtering to get only non-e002, e001 treatments, later will remove a subset of the treatments from these two projects because they have so many levels and make up a disproportionate amount of the entire dataset
divTrt1 <- div%>%
  filter(pulse==0, plot_mani>0, project_name!="e001"&project_name!="e002")

#calculate average dispersion among just control plots (this is an estimate of average community dissimilarity, to use as a baseline for dissimilarity change between treatment and control plots)
summary(divControls$ctl_dispersion) #mean=0.2811, median=0.2778

#removing a subset of CDR e001 and e002 treatments to prevent the majority of data being from CDR; keeping lowest, highest, and 10 gm-2 N additions (levels most comparable to other studies)
divCDRe001<-div%>%
  filter(site_code=="CDR"&treatment==1|treatment==6|treatment==8|treatment==9,plot_mani>0)
divCDRe002<-div%>%
  filter(site_code=="CDR"&treatment=='1_f_u_n'|treatment=='6_f_u_n'|treatment=='8_f_u_n'|treatment=='9_f_u_n',plot_mani>0)

#combine the two CDR experiments
divTrt<-rbind(divTrt1, divCDRe002, divCDRe001)

##16% of our data is from CDR and 10% is from KNZ
#merge controls and treatments
divCompare <- divControls%>%
  left_join(divTrt)%>%
  #calculate difference in disperion, expH, S, and evenness (note anytime it says change in a header, it is really a difference between treatment and control plots)
  mutate(dispersion_change=dispersion-ctl_dispersion, 
         expH_PC=(expH-ctl_expH)/ctl_expH, 
         S_PC=(S-ctl_S)/ctl_S, 
         SimpEven_change=SimpEven-ctl_SimpEven)%>%
  select(exp_year, treatment_year, treatment, plot_mani, mean_change, dispersion_change, expH_PC,  SimpEven_change, S_PC, site_code, project_name, community_type, calendar_year)


###merging with treatment information
SiteExp<-read.csv("SiteExperimentDetails_Dec2016.csv")%>%
  select(-X)

ForAnalysis<-merge(divCompare, SiteExp, by=c("site_code","project_name","community_type"))

# full dataset
# write.csv(ForAnalysis, "ForBayesianAnalysis_May2017.csv")


###generating treatment categories (resource, non-resource, and interactions)
###4 steps: (1) single resource, (2) single non-resource, (3) 2-way interactions, (4) 3+ way interactions
#step 1: single resource
singleResource <- ForAnalysis%>%
  select(-plot_mani)%>%
  left_join(expInfo)%>%
  #filter pretrt data
  filter(treatment_year!=0)%>%
  #filter just single resource manipulations
  filter(resource_mani==1, plot_mani==1)%>%
  #set CEH Megarich nutrient values to 0 (added to all megaliths, not a treatment)
  mutate(n2=ifelse(site_code=='CEH', 0, n), p2=ifelse(site_code=='CEH', 0, p), k2=ifelse(site_code=='CEH', 0, k))%>%
  #drop lime added, as only one trt does this
  filter(other_trt!='lime added')%>%
  #create categorical treatment type column
  mutate(trt_type=ifelse(n2>0, 'N', ifelse(p2>0, 'P', ifelse(k2>0, 'K', ifelse(precip<0, 'drought', ifelse(precip>0, 'irr', ifelse(CO2>0, 'CO2', 'precip_vari')))))))%>%
  #keep just relevent column names for this analysis
  select(site_code, project_name, community_type, exp_year, treatment_year, calendar_year, treatment, trt_type, mean_change, S_PC, expH_PC, experiment_length, rrich, anpp, MAT, MAP)

#step 2: single non-resource
singleNonresource <- ForAnalysis%>%
  left_join(expInfo)%>%
  #filter pretrt data
  filter(treatment_year!=0)%>%
  #filter just single resource manipulations
  filter(resource_mani==0, plot_mani==1)%>%
  #drop tilled, stone, and fungicide as only one trt each do these; drop 'current pattern' of rainfall, as this is basically a control
  filter(other_trt!='tilled', other_trt!='shallow soil', other_trt!='fungicide added', other_trt!='current pattern')%>%
  #create categorical treatment type column
  mutate(trt_type=ifelse(burn==1, 'burn', ifelse(mow_clip==1, 'mow_clip', ifelse(herb_removal==1, 'herb_rem', ifelse(temp>0, 'temp', ifelse(plant_trt==1, 'plant_mani', 'other'))))))%>%
  #keep just relevent column names for this analysis
  select(site_code, project_name, community_type, exp_year, treatment_year, calendar_year, treatment, trt_type, mean_change, S_PC, expH_PC, experiment_length, rrich, anpp, MAT, MAP)

#step 3: 2-way interactions
twoWay <- ForAnalysis%>%
  left_join(expInfo)%>%
  #filter pretrt data
  filter(treatment_year!=0)%>%
  #filter just single resource manipulations
  filter(plot_mani==2)%>%
  #drop tilled, stone, and fungicide as only one trt each do these; drop 'current pattern' of rainfall, as this is basically a control
  filter(other_trt!='tilled', other_trt!='shallow soil', other_trt!='fungicide added', other_trt!='current pattern')%>%
  #create categorical treatment type column
  mutate(trt_type=ifelse(resource_mani==1&burn==1, 'R*burn', ifelse(resource_mani==1&mow_clip==1, 'R*mow_clip', ifelse(resource_mani==1&herb_removal==1, 'R*herb_rem', ifelse(resource_mani==1&temp>0, 'R*temp', ifelse(resource_mani==1&plant_trt==1, 'R*plant_mani', ifelse(resource_mani==1&other_trt!=0, 'R*other', ifelse(n>0&p>0, 'R*R', ifelse(n>0&CO2>0, 'R*R', ifelse(n>0&precip!=0, 'R*R', ifelse(p>0&k>0, 'R*R', ifelse(CO2>0&precip!=0, 'R*R', 'N*N'))))))))))))%>%
  #drop R*herb_removal (single rep)
  filter(trt_type!='R*herb_rem')%>%
  #keep just relevent column names for this analysis
  select(site_code, project_name, community_type, exp_year, treatment_year, calendar_year, treatment, trt_type, mean_change, S_PC, expH_PC, experiment_length, rrich, anpp, MAT, MAP)

#step 4: 3+ way interactions
threeWay <- ForAnalysis%>%
  left_join(expInfo)%>%
  #filter pretrt data
  filter(treatment_year!=0)%>%
  #filter just single resource manipulations
  filter(plot_mani>2, plot_mani<6)%>%
  #drop tilled, stone, and fungicide as only one trt each do these; drop 'current pattern' of rainfall, as this is basically a control
  filter(other_trt!='tilled', other_trt!='shallow soil', other_trt!='fungicide added', other_trt!='current pattern')%>%
  #create categorical treatment type column
  mutate(trt_type=ifelse(burn==0&mow_clip==0&herb_removal==0&temp==0&plant_trt==0, 'all_resource', ifelse(n==0&p==0&k==0&CO2==0&precip==0, 'all_nonresource', 'both')))%>%
  #drop single all-nonresource treatment (NIN herbdiv 5NF)
  filter(trt_type!='all_nonresource')%>%
  #keep just relevent column names for this analysis
  select(site_code, project_name, community_type, exp_year, treatment_year, calendar_year, treatment, trt_type, mean_change, S_PC, expH_PC, experiment_length, rrich, anpp, MAT, MAP)

#combine for analysis - one big model, 19 trt types
allAnalysis <- rbind(singleResource, singleNonresource, twoWay, threeWay)


#subset out datasets with less than 3 temporal data points;  drops GVN FACE
numPoints <- allAnalysis%>%
  select(site_code, project_name, community_type, treatment, treatment_year)%>%
  unique()%>%
  group_by(site_code, project_name, community_type, treatment)%>%
  summarise(num_datapoints=length(treatment_year))
allAnalysisAllDatasets <- allAnalysis%>%
  left_join(numPoints)%>%
  filter(num_datapoints>2)
# write.csv(allAnalysisAllDatasets, 'ForAnalysis_allAnalysisAllDatasets.csv')


#subset out treatment years 20 or less (i.e., cut off datasets at 20 years)
allAnalysis20yr <- allAnalysisAllDatasets%>%
  filter(treatment_year<21)%>%
  select(-num_datapoints)
numPoints <- allAnalysis20yr%>%
  select(site_code, project_name, community_type, treatment, treatment_year)%>%
  unique()%>%
  group_by(site_code, project_name, community_type, treatment)%>%
  summarise(num_datapoints=length(treatment_year))