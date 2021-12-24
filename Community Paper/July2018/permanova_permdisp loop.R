### calculating permanova and permdisp for each experiment and year
###
### Authors: Kevin Wilcox (wilcoxkr@gmail.com) Andrew Tredennick (atredenn@gmail.com) and Meghan Avolio (meghan.avolio@jhu.edu)
### Last updated: Oct 13 2021
### .

### Set up workspace
setwd("C:\\Users\\wilco\\Dropbox\\")
setwd("C:\\Users\\megha\\Dropbox\\")
setwd("C:\\Users\\mavolio2\\Dropbox\\")

library(tidyverse)
library(vegan)
library(ggthemes)

corredat<-read.csv("converge_diverge/datasets/LongForm/SpeciesRelativeAbundance_Nov2019.csv")

#gvn face - only 2 years of data so will only have one point for the dataset, therefore we are removing this dataset from these analyses.
corredat1<-corredat%>%
  select(-X)%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))%>%
  filter(site_project_comm!="GVN_FACE_0", site_project_comm!="AZI_NitPhos_0", site_project_comm!="JRN_study278_0", site_project_comm!="KNZ_GFP_4F", site_project_comm!="Saskatchewan_CCD_0", project_name!="e001", project_name!="e002",site_project_comm!="CHY_EDGE_0", site_project_comm!="HYS_EDGE_0", site_project_comm!="SGS_EDGE_0")

##several studies only have two measurments of a plot. I am dropping those plots
azi<-corredat%>%
  select(-X)%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))%>%
  filter(site_code=="AZI")%>%
  filter(plot_id!=11&plot_id!=15&plot_id!=35&plot_id!=37)

jrn<-corredat%>%
  select(-X)%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))%>%
  filter(site_project_comm=="JRN_study278_0")%>%
  filter(plot_id!=211&plot_id!=210)

knz<-corredat%>%
  select(-X)%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))%>%
  filter(site_project_comm=="KNZ_GFP_4F")%>%
  filter(plot_id!="7_1_1"&plot_id!="7_2_1")

sak<-corredat%>%
  select(-X)%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))%>%
  filter(site_project_comm=="Saskatchewan_CCD_0")%>%
  filter(plot_id!=2)

###remove extra treatments from CDR e001 and e002
cdr <- corredat%>%
  select(-X)%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))%>%
  filter(project_name=="e001"|project_name=="e002")%>%
  filter(treatment==1|treatment==6|treatment==8|treatment==9|treatment=='1_f_u_n'|treatment=='6_f_u_n'|treatment=='8_f_u_n'|treatment=='9_f_u_n')

##remove one of 2 pre-treatment years in edge for CHY, SGS, and HAYS
edge<-corredat%>%
  select(-X)%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))%>%
  filter(site_code=="CHY"|site_code=="SGS"|site_code=="HYS"&project_name=="EDGE")%>%
  filter(calendar_year!=2012)


###final dataset to use
corredat_raw<-rbind(corredat1, azi, jrn, knz, sak, cdr, edge)

treatment_info<-read.csv("converge_diverge/datasets/LongForm/ExperimentInformation_March2019.csv")%>%
  dplyr::select(site_code, project_name, community_type, treatment,plot_mani)%>%
  unique()%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))

corredat <- corredat_raw %>%
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
      
      
      permdisp_temp <- betadisper(vegdist(cover_temp), env_temp$plot_mani, type = "centroid")
      sig_disp <- permutest(permdisp_temp)

      perm_out_temp <- data.frame(
        site_project_comm = site_project_comm_vec[i],
        treatment = treatment_vec[trt],
        calendar_year = year_vec[yr],
        perm_Pvalue =  permanova_temp$aov.tab$'Pr(>F)'[1],
        disp_Pvalue = sig_disp$tab$'Pr(>F)'[1]
        
      )
      
      permanova_out_master <- rbind(permanova_out_master, perm_out_temp)
     
    }
  }
}
       
permanova_out_mod <- permanova_out_master %>%
  mutate(pval_flag_perm = ifelse(perm_Pvalue<.05, 1, 0)) %>%
  mutate(pval_flag_disp = ifelse(disp_Pvalue<.05, 1, 0)) %>%
  mutate(permanova="permanova") %>%
  group_by(site_project_comm, treatment) %>%
  summarise(tot_pval=sum(pval_flag)) %>%
  filter(tot_pval != 0)

#write.csv(permanova_out_mod, file= "C2E\\Products\\CommunityChange\\March2018 WG\\permanova out.csv",row.names=F)

filter(permanova_out_master, site_project_comm =="ASGA_clonal_0")

write.csv(permanova_out_master, file="C2E\\Products\\Testing Hypots\\permanova_permdisp_outputOct2021.csv", row.names=F)

###doing for rare speceis removed, DO THIS FOR 1% AND 5%
noraresp<-corredat%>%
  filter(plot_mani==0)%>%
  group_by(site_project_comm, genus_species)%>%
  summarize(mrelcov=mean(relcov))%>%
  filter(mrelcov>0.05)%>%
  select(-mrelcov)


corredat_norare <- corredat %>%
  right_join(noraresp)

site_project_comm_vec <- unique(corredat_norare$site_project_comm)

permanova_out_master_norare <- {}

for(i in 1:length(site_project_comm_vec)) {
  corredat_temp <- filter(corredat_norare, site_project_comm==site_project_comm_vec[i])
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
      
      
      permdisp_temp <- betadisper(vegdist(cover_temp), env_temp$plot_mani, type = "centroid")
      sig_disp <- permutest(permdisp_temp)
      
      perm_out_temp <- data.frame(
        site_project_comm = site_project_comm_vec[i],
        treatment = treatment_vec[trt],
        calendar_year = year_vec[yr],
        perm_Pvalue =  permanova_temp$aov.tab$'Pr(>F)'[1],
        disp_Pvalue = sig_disp$tab$'Pr(>F)'[1]
        
      )
      
      permanova_out_master_norare <- rbind(permanova_out_master_norare, perm_out_temp)
      
    }
  }
}



write.csv(permanova_out_master_norare, file="C2E\\Products\\CommunityChange\\March2018 WG\\permanova_permdisp_output_norare5_Nov2021.csv", row.names=F)

