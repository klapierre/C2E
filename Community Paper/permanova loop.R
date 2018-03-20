### calculating permanova for each experiment and year
###
### Authors: Kevin Wilcox (wilcoxkr@gmail.com) Andrew Tredennick (atredenn@gmail.com)
### Last updated: March 20 2018

### Set up workspace
setwd("C:\\Users\\wilco\\Dropbox\\")
library(tidyverse)
library(vegan)
library(ggthemes)
corredat1<-read.csv("C2E\\Products\\CommunityChange\\March2018 WG\\data\\SpeciesRelativeAbundance_Oct2017.csv")%>%
  select(-X)%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))%>%
  filter(site_project_comm!="GVN_FACE_0", site_project_comm!="AZI_NitPhos_0", site_project_comm!="JRN_study278_0", site_project_comm!="KNZ_GFP_4F", site_project_comm!="Saskatchewan_CCD_0")

##several studies only have two measurments of a plot. I am dropping those plots
azi<-read.csv("C2E\\Products\\CommunityChange\\March2018 WG\\data\\SpeciesRelativeAbundance_Oct2017.csv")%>%
  select(-X)%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))%>%
  filter(site_code=="AZI")%>%
  filter(plot_id!=11&plot_id!=15&plot_id!=35&plot_id!=37)

jrn<-read.csv("C2E\\Products\\CommunityChange\\March2018 WG\\data\\SpeciesRelativeAbundance_Oct2017.csv")%>%
  select(-X)%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))%>%
  filter(site_project_comm=="JRN_study278_0")%>%
  filter(plot_id!=211&plot_id!=210)

knz<-read.csv("C2E\\Products\\CommunityChange\\March2018 WG\\data\\SpeciesRelativeAbundance_Oct2017.csv")%>%
  select(-X)%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))%>%
  filter(site_project_comm=="KNZ_GFP_4F")%>%
  filter(plot_id!="7_1_1"&plot_id!="7_2_1")

sak<-read.csv("C2E\\Products\\CommunityChange\\March2018 WG\\data\\SpeciesRelativeAbundance_Oct2017.csv")%>%
  select(-X)%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))%>%
  filter(site_project_comm=="Saskatchewan_CCD_0")%>%
  filter(plot_id!=2)

corredat_raw <-rbind(corredat1, azi, jrn, knz, sak)

treatment_info<-read.csv("C2E\\Products\\CommunityChange\\March2018 WG\\data\\ExperimentInformation_Nov2017.csv")%>%
  dplyr::select(site_code, project_name, community_type, treatment,plot_mani)%>%
  unique()%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))

corredat <- corredat_raw %>%
  filter(site_project_comm != "GVN_FACE_0") %>%
  left_join(treatment_info, by=c( "site_code","project_name","community_type", "treatment","site_project_comm"))

site_project_comm_vec <- unique(corredat$site_project_comm)

permanova_out_master <- {}

for(i in 1:length(site_project_comm_vec)) {
  corredat_temp <- filter(corredat, site_project_comm==site_project_comm_vec[i])
  control_temp <- filter(corredat_temp, plot_mani==0)
  treatment_temp <- filter(corredat_temp, plot_mani !=0)
  
  treatment_vec <- unique(treatment_temp$treatment)
  
  for(trt in 1:length(treatment_vec)){
    treatment_temp_2 <- filter(treatment_temp, treatment==treatment_vec[trt])
    control_years_vec <- unique(control_temp$calendar_year)
    treatment_years_vec <- unique(treatment_temp_2$calendar_year)
    
    year_keeper <- intersect(control_years_vec, treatment_years_vec)
    
    one_treatment_temp <- filter(treatment_temp, treatment==treatment_vec[trt]) %>%
      filter(calendar_year %in% year_keeper)
    control_temp_2 <- filter(control_temp, calendar_year %in% year_keeper)

    trt_ctrl_df_temp <- one_treatment_temp %>%
      bind_rows(control_temp_2)

    year_vec <- unique(trt_ctrl_df_temp$calendar_year)
    
    for(yr in 1:length(year_vec)){
      
      trt_ctrl_yr_temp <- trt_ctrl_df_temp %>%
        filter(calendar_year==year_vec[yr]) %>%
        spread(key=genus_species, value=relcov, fill=0)
      
      cover_temp <- trt_ctrl_yr_temp %>%
        dplyr::select(-(site_code:plot_mani))

      env_temp <- trt_ctrl_yr_temp %>%
        dplyr::select(site_code:plot_mani)
      
      permanova_temp <- adonis(cover_temp ~ plot_mani, data=env_temp, permutations=99)

      perm_out_temp <- data.frame(
        site_project_comm = site_project_comm_vec[i],
        treatment = treatment_vec[trt],
        calendar_year = year_vec[yr],
        Pvalue =  permanova_temp$aov.tab$'Pr(>F)'[1]
      )
      
      permanova_out_master <- rbind(permanova_out_master, perm_out_temp)
     
    }
  }
}
       
permanova_out_mod <- permanova_out_master %>%
  mutate(pval_flag = ifelse(Pvalue<.05, 1, 0)) %>%
  mutate(permanova="permanova") %>%
  group_by(site_project_comm, treatment) %>%
  summarise(tot_pval=sum(pval_flag)) %>%
  filter(tot_pval != 0)

#write.csv(permanova_out_mod, file= "C2E\\Products\\CommunityChange\\March2018 WG\\permanova out.csv",row.names=F)





