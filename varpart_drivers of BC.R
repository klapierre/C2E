### Variance partitioning of drivers of BC dissimilarity
### Author: Kevin Wilcox (wilcoxkr@gmail.com)
###
### Last modified: Oct 14, 2017
################################################################################

##  Clear everything...
rm(list=ls())


####
####  LOAD LIBRARIES ----
####
library(tidyverse)
library(vegan)
library(ggplot2)
library(ggthemes)



####
####  SETUP DIRECTORIES ----
##  KW
#setwd("C:\\Users\\Kevin.Wilcox\\Desktop\\C2E\\BC drivers\\")
#setwd("C:\\Users\\WilcoxKR.WILCOXKR-LAPTOP\\Desktop\\Working groups\\C2E\\BC drivers\\")
#data_dir <- paste0(getwd(),"/data/")

##  AT
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
data_dir <- "../../Google_Drive/C2E/data/"



####
####  READ IN DATA AND FORMAT ----
####
df_raw <- readRDS(paste0(data_dir,"bootstrapped_rac_metrics.RDS")) %>%
  filter(site_project_comm!="CUL_Culardoch_0")

trt_key <- read.csv(paste0(data_dir, "ExperimentInformation_May2017.csv")) %>%
  mutate(site_project_comm = paste(site_code,project_name,community_type,sep="_")) %>%
  dplyr::select(site_project_comm, treatment, plot_mani) %>%
  unique()

df_all <- df_raw %>%
  group_by(site_project_comm, treatment) %>%
  mutate(year_1=min(calendar_year)) %>%
  mutate(treatment_year = calendar_year-year_1+1) %>%
  dplyr::select(-year_1)
  
exp_to_include <- df_all %>%
  group_by(site_project_comm, treatment, calendar_year) %>%
  summarise(avg_s=mean(S,na.rm=T),avg_bc=mean(bc_dissim,na.rm=T)) %>%
  filter(avg_s!=1 & avg_bc!=1) %>%
  dplyr::select(site_project_comm, treatment, calendar_year)

suppressWarnings( # ignore coerce factor to character message
  df <- df_all %>%
    right_join(exp_to_include, by = c("site_project_comm", "treatment", "calendar_year")) %>%
    left_join(trt_key, c("site_project_comm", "treatment"))
)

trt_df <- df %>%
  filter(plot_mani != 0)

cntrl_df <- df %>%
  filter(plot_mani == 0)

site_project_comm_vec <- trt_df %>%
  .$site_project_comm %>%
  unique()

if(nrow(trt_df)+nrow(cntrl_df) != nrow(df)){
  stop("Check data frame splitting; wrong lengths")
}



####
####  LOOP THROUGH AND VARPART IT ----
####
vpart_out <- {} # empty object for storage

##  Progress bar
pb <- txtProgressBar(min=1, max=length(site_project_comm_vec), char="+", style=3, width=65) 
counter <- 1

for(site in 1:length(site_project_comm_vec)){
  df_temp_site <- trt_df %>%
    filter(site_project_comm == site_project_comm_vec[site])
  
  treatment_vec <- df_temp_site %>%
    .$treatment %>%
    unique()
  
  for(trt in 1:length(treatment_vec)){
    df_temp_site_trt <- df_temp_site %>%
      filter(treatment == treatment_vec[trt])
    
    calendar_year_vec <- df_temp_site_trt %>%
      .$calendar_year %>%
      unique()
    
    for(year in 1:length(calendar_year_vec)){
      df_temp_site_trt_year <- df_temp_site_trt %>%
        filter(calendar_year==calendar_year_vec[year]) %>%
        mutate(id=1:nrow(.)) %>%
        gather(key=metric, value=value, appearance:bc_dissim)

      ##  ID metrics that don't vary for removal from this varpart
      metrics_to_remove <- df_temp_site_trt_year %>%
        group_by(metric) %>%
        summarise(var=var(value)) %>%
        filter(var==0) %>%
        .$metric %>%
        unique()
      
      ##  Take out the non-varying metrics
      df_work <- df_temp_site_trt_year %>%
        filter(!metric %in% metrics_to_remove)
      
      ##  Store used metrics as character vector
      metrics_to_include <- unique(df_work$metric) 
      metrics_to_include <- metrics_to_include[which(!metrics_to_include %in% 
                                              c("S", "bc_dissim"))]    
      
      ##  Go wide for matrix-like format (for varpart)
      df_wide <- df_work %>%
        ungroup() %>%
        filter(metric %in% c(metrics_to_include,"bc_dissim")) %>%
        spread(key=metric, value=value) %>%
        dplyr::select(c(metrics_to_include,"bc_dissim"))
      
      cntrl_temp <- cntrl_df %>%
        filter(site_project_comm == site_project_comm_vec[site])
      cntrl_years <- unique(cntrl_temp$calendar_year)
      
      ##  ONLY DO ANALYSIS IF CONTROL SAMPLED IN THIS YEAR
      ##  OTHERWISE, DO NOTHING
      if(calendar_year_vec[year] %in% cntrl_years){
        df_wide_cntrl <- cntrl_df %>%
          filter(site_project_comm == site_project_comm_vec[site]) %>%
          filter(calendar_year==calendar_year_vec[year]) %>%
          mutate(id=1:nrow(.)) %>%
          gather(key = metric, value = value, appearance:bc_dissim) %>%
          filter(metric %in% c(metrics_to_include,"bc_dissim")) %>%
          spread(key=metric, value=value) %>%
          ungroup() %>%
          dplyr::select(c(metrics_to_include,"bc_dissim"))
        
        df_analyze <- rbind(df_wide, df_wide_cntrl)
        if(nrow(df_wide)+nrow(df_wide_cntrl) != nrow(df_analyze)){
          stop("YOU BROKE IT!!! Check lengths of data frames.")
        }
        
        
        ##  Make character vector of individual formula parts
        eqn_work <- paste0("~",metrics_to_include)      
        
        ##  IF all metrics vary, do this
        if(length(eqn_work)==4){
          vpart_temp <- varpart(df_analyze$bc_dissim, 
                                as.formula(eqn_work[1]),
                                as.formula(eqn_work[2]),
                                as.formula(eqn_work[3]),
                                as.formula(eqn_work[4]), 
                                data=df_analyze)
          
          vpart_out_temp <- data.frame(site_project_comm = site_project_comm_vec[site],
                                       treatment = treatment_vec[trt],
                                       calendar_year = calendar_year_vec[year],
                                       treatment_year = df_temp_site_trt_year$treatment_year[1],
                                       metric = c(metrics_to_include,"total_var_expl"),
                                       var_expl_alone = c(vpart_temp$part$indfract$Adj.R.square[1:4],
                                                          vpart_temp$part$fract$Adj.R.square[15]),
                                       var_expl_together = c(vpart_temp$part$fract$Adj.R.square[1:4],
                                                             vpart_temp$part$fract$Adj.R.square[15]))
          
        }
        
        ##  IF only three metrics vary, do this
        if(length(eqn_work)==3){
          vpart_temp <- varpart(df_analyze$bc_dissim, 
                                as.formula(eqn_work[1]),
                                as.formula(eqn_work[2]),
                                as.formula(eqn_work[3]), 
                                data=df_analyze)
          
          vpart_out_temp <- data.frame(site_project_comm = site_project_comm_vec[site],
                                       treatment = treatment_vec[trt],
                                       calendar_year = calendar_year_vec[year],
                                       treatment_year = df_temp_site_trt_year$treatment_year[1],
                                       metric = c(metrics_to_include,"total_var_expl"),
                                       var_expl_alone = c(vpart_temp$part$indfract$Adj.R.square[1:3],
                                                          vpart_temp$part$fract$Adj.R.square[7]),
                                       var_expl_together = c(vpart_temp$part$fract$Adj.R.square[1:3],
                                                             vpart_temp$part$fract$Adj.R.square[7]))
          
        }
        
        ##  IF only two metrics vary, do this
        if(length(eqn_work)==2){
          vpart_temp <- varpart(df_analyze$bc_dissim, 
                                as.formula(eqn_work[1]),
                                as.formula(eqn_work[2]),
                                data=df_analyze)
          
          vpart_out_temp <- data.frame(site_project_comm = site_project_comm_vec[site],
                                       treatment = treatment_vec[trt],
                                       calendar_year = calendar_year_vec[year],
                                       treatment_year = df_temp_site_trt_year$treatment_year[1],
                                       metric = c(metrics_to_include,"total_var_expl"),
                                       var_expl_alone = c(vpart_temp$part$indfract$Adj.R.square[c(1,3)],
                                                          vpart_temp$part$fract$Adj.R.square[3]),
                                       var_expl_together = c(vpart_temp$part$fract$Adj.R.square[c(1:2)],
                                                             vpart_temp$part$fract$Adj.R.square[3]))
          
        }
        
        ##  IF only one metric varies, run lm() instead of varpart()
        if(length(eqn_work)==1){
          eqn_lm <- paste("bc_dissim~",metrics_to_include[1])
          vpart_adj_r2 <- summary(lm(as.formula(eqn_lm),data=df_analyze))$adj.r.squared
          
          vpart_out_temp <- data.frame(site_project_comm = site_project_comm_vec[site],
                                       treatment = treatment_vec[trt],
                                       calendar_year = calendar_year_vec[year],
                                       treatment_year = df_temp_site_trt_year$treatment_year[1],
                                       metric = c(metrics_to_include,"total_var_expl"),
                                       var_expl_alone = vpart_adj_r2,
                                       var_expl_together = vpart_adj_r2)
          
        }
        
        ##  Bind to save
        vpart_out <- rbind(vpart_out, vpart_out_temp)
      } # end if for present year
      
      rm(vpart_out_temp, vpart_temp, df_temp_site_trt_year)
    } # next year
    rm(df_temp_site_trt, calendar_year_vec)
  } # next treatments
  rm(df_temp_site, treatment_vec)
  
  ##  Update progess
  setTxtProgressBar(pb, counter)
  counter <- counter+1
} # next site_project_comm



####
####  SAVE OUTPUT ----
####
write.csv(vpart_out, file=paste0(data_dir,"varpart_output_",Sys.Date(),".csv"))

