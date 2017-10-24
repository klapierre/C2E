### Variance partitioning of drivers of BC dissimilarity
### Author: Kevin Wilcox (wilcoxkr@gmail.com)
###
### Last modified: Oct 14, 2017
rm(list=ls())

library(tidyverse)
library(vegan)
library(ggplot2)
library(ggthemes)

#setwd("C:\\Users\\Kevin.Wilcox\\Desktop\\C2E\\BC drivers\\")
setwd("C:\\Users\\WilcoxKR.WILCOXKR-LAPTOP\\Desktop\\Working groups\\C2E\\BC drivers\\")

df_raw <- readRDS("data/bootstrapped_rac_metrics.RDS") %>%
  filter(site_project_comm!="CUL_Culardoch_0")

df <- df_raw %>%
  group_by(site_project_comm, treatment) %>%
  mutate(year_1=min(calendar_year)) %>%
  mutate(treatment_year = calendar_year-year_1+1) %>%
  dplyr::select(-year_1)


site_project_comm_vec <- df %>%
  .$site_project_comm %>%
  unique()
  
vpart_out <- {}
  
exp_to_include <- df %>%
  group_by(site_project_comm, treatment, calendar_year) %>%
  summarise(avg_s=mean(S,na.rm=T),avg_bc=mean(bc_dissim,na.rm=T)) %>%
  filter(avg_s!=1 & avg_bc!=1) %>%
  dplyr::select(site_project_comm, treatment, calendar_year)
  
df <- df %>%
  right_join(exp_to_include)

for(site in 1:length(site_project_comm_vec)){
  df_temp_site <- df %>%
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

      metrics_to_remove <- df_temp_site_trt_year %>%
        group_by(metric) %>%
        summarise(var=var(value)) %>%
        filter(var==0) %>%
        .$metric %>%
        unique()

      df_work <- df_temp_site_trt_year %>%
        filter(!metric %in% metrics_to_remove)

      metrics_to_include <- unique(df_work$metric) 
      metrics_to_include <- metrics_to_include[which(!metrics_to_include %in% 
                                              c("S", "bc_dissim"))]      
      df_wide <- df_work %>%
        ungroup() %>%
        filter(metric %in% c(metrics_to_include,"bc_dissim")) %>%
        spread(key=metric, value=value) %>%
      dplyr::select(-(site_project_comm:id))
        
      #eqn_work <- paste("df_wide$bc_dissim, ",
      #                  paste(paste0("~",metrics_to_include,","),collapse=" "))      
      eqn_work <- paste0("~",metrics_to_include)      
      
      if(length(eqn_work)==4){
        vpart_temp <- varpart(df_wide$bc_dissim, 
                              as.formula(eqn_work[1]),
                              as.formula(eqn_work[2]),
                              as.formula(eqn_work[3]),
                              as.formula(eqn_work[4]), data=df_wide)

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
      if(length(eqn_work)==3){
        vpart_temp <- varpart(df_wide$bc_dissim, 
                              as.formula(eqn_work[1]),
                              as.formula(eqn_work[2]),
                              as.formula(eqn_work[3]), data=df_wide)

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

      if(length(eqn_work)==2){
        vpart_temp <- varpart(df_wide$bc_dissim, 
                              as.formula(eqn_work[1]),
                              as.formula(eqn_work[2]),
                              data=df_wide)

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

      if(length(eqn_work)==1){
      eqn_lm <- paste("bc_dissim~",metrics_to_include[1])
        vpart_adj_r2 <- summary(lm(as.formula(eqn_lm),data=df_wide))$adj.r.squared

        vpart_out_temp <- data.frame(site_project_comm = site_project_comm_vec[site],
                                     treatment = treatment_vec[trt],
                                     calendar_year = calendar_year_vec[year],
                                     treatment_year = df_temp_site_trt_year$treatment_year[1],
                                     metric = c(metrics_to_include,"total_var_expl"),
                                     var_expl_alone = vpart_adj_r2,
                                     var_expl_together = vpart_adj_r2)
        
      }

      vpart_out <- rbind(vpart_out, vpart_out_temp)
      rm(vpart_out_temp, vpart_temp, df_temp_site_trt_year)
    }
    rm(df_temp_site_trt, calendar_year_vec)
  }
  rm(df_temp_site, treatment_vec)
}

write.csv(vpart_out, file="data//varpart_output_23Oct2017.csv")

ggplot(subset(vpart_out, metric!="total_var_expl"&site_project_comm=="KNZ_pplots_0"),
       aes(x=treatment_year, y=var_expl_together, colour=metric)) +
  geom_line() + facet_wrap(~treatment) + theme_few()

vpart_means <- vpart_out %>%
  group_by(treatment_year, metric) %>%
  summarise(var_expl_alone=mean(var_expl_alone,na.rm=T))

ggplot(subset(vpart_out, metric!="total_var_expl"),
       aes(x=treatment_year, y=var_expl_together, colour=metric)) +
  geom_line() + theme_few()

ggplot(subset(vpart_means, metric!="total_var_expl"),
       aes(x=treatment_year, y=var_expl_alone, fill=metric)) +
  geom_col(position=position_dodge(0.9))


