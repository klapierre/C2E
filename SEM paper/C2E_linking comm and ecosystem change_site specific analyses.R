library(tidyverse)
library(grid)
library(PerformanceAnalytics)
library(lme4)
library(lmerTest)
library(piecewiseSEM)




#source data management code  -- if on desktop, change source code
source('C:\\Users\\lapie\\Desktop\\R files laptop\\C2E\\SEM paper\\C2E_SEM_data processing.R')



#kim's desktop
setwd('C:\\Users\\la pierrek\\Dropbox (Smithsonian)\\working groups\\CoRRE\\converge_diverge\\datasets\\LongForm')

#kim's laptop
setwd('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\CoRRE\\converge_diverge\\datasets\\LongForm')


theme_set(theme_bw())
theme_update(axis.title.x=element_text(size=40, vjust=-0.35, margin=margin(t=15)), axis.text.x=element_text(size=34, color='black'),
             axis.title.y=element_text(size=40, angle=90, vjust=0.5, margin=margin(r=15)), axis.text.y=element_text(size=34, color='black'),
             plot.title = element_text(size=40, vjust=2),
             panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
             legend.title=element_blank(), legend.text=element_text(size=20))


###testing out multiple regressions with random effect of site_project_comm
summary(lmer(anpp_pdiff ~ richness_difference*treatment_year + evenness_diff*treatment_year + rank_difference*treatment_year + species_difference*treatment_year + (1|site_project_comm), data=correSEMdataTrt))

ggplot(correSEMdataTrt, aes(x=treatment_year, y=anpp_pdiff, color=as.character(site_project_comm_trt))) + 
  geom_point() + 
  guides(color = "none") + 
  stat_smooth(method="lm", formula = y ~ x + I(x^2), fill=NA)

#with transformed variables
summary(lmer(anpp_pdiff_transform ~ richness_difference*treatment_year + evenness_diff*treatment_year + rank_difference*treatment_year + species_difference*treatment_year + (1|site_project_comm), data=correSEMdataTrt))

ggplot(correSEMdataTrt, aes(x=treatment_year, y=anpp_pdiff_transform, color=as.character(site_project_comm_trt))) + 
  geom_point() + 
  guides(color = "none") + 
  stat_smooth(method="lm", formula = y ~ x + I(x^2), fill=NA)


###split each site down into long form for just the variables we want
correDataTrt <- correSEMdataTrt%>%
  #filter out projects with less than 2 ANPP datapoints
  filter(!project_name %in% c('GFP', 'gb', 'TIDE', 'CXN'))%>%
  select(site_project_comm_trt, treatment, treatment_year, anpp_pdiff, richness_difference, evenness_diff, rank_difference, species_difference, composition_diff)%>%
  group_by(site_project_comm_trt, treatment, treatment_year)%>%
  gather(key='factor_type', value='value', anpp_pdiff:composition_diff, na.rm=T)%>%
  ungroup()%>%
  #create treatment_year squared variable
  mutate(treatment_year_sq=(treatment_year^2))


###loop across all sites
#generate site list
site_trt_list <- unique(correDataTrt$site_project_comm_trt)


###find the regressions that are significant through time
#create empty dataframe
regressionOutput=data.frame(row.names=1)

for (i in seq_along(site_trt_list)) { 
  
  # run a quadratic regression for each treatment in df for each variable 
  anpp <- lm(value ~ treatment_year + (treatment_year_sq), data=subset(correDataTrt, site_project_comm_trt==site_trt_list[i]&factor_type=='anpp_pdiff'))
  comp <- lm(value ~ treatment_year + (treatment_year_sq), data=subset(correDataTrt, site_project_comm_trt==site_trt_list[i]&factor_type=='composition_diff'))
  evenness <- lm(value ~ treatment_year + (treatment_year_sq), data=subset(correDataTrt, site_project_comm_trt==site_trt_list[i]&factor_type=='evenness_diff'))
  rank <- lm(value ~ treatment_year + (treatment_year_sq), data=subset(correDataTrt, site_project_comm_trt==site_trt_list[i]&factor_type=='rank_difference'))
  richness <- lm(value ~ treatment_year + (treatment_year_sq), data=subset(correDataTrt, site_project_comm_trt==site_trt_list[i]&factor_type=='richness_difference'))
  species <- lm(value ~ treatment_year + (treatment_year_sq), data=subset(correDataTrt, site_project_comm_trt==site_trt_list[i]&factor_type=='species_difference'))
  
  # gather p-values for each regression
  anpp_pval <- as.data.frame(summary(anpp)$coefficients[,4])%>%rename(anpp_pval='summary(anpp)$coefficients[, 4]')
  comp_pval <- as.data.frame(summary(comp)$coefficients[,4])%>%rename(comp_pval='summary(comp)$coefficients[, 4]')
  evenness_pval <- as.data.frame(summary(evenness)$coefficients[,4])%>%rename(evenness_pval='summary(evenness)$coefficients[, 4]')
  rank_pval <- as.data.frame(summary(rank)$coefficients[,4])%>%rename(rank_pval='summary(rank)$coefficients[, 4]')
  richness_pval <- as.data.frame(summary(richness)$coefficients[,4])%>%rename(richness_pval='summary(richness)$coefficients[, 4]')
  species_pval <- as.data.frame(summary(species)$coefficients[,4])%>%rename(species_pval='summary(species)$coefficients[, 4]')
  
  # gather estimates for each regression
  anpp_est <- as.data.frame(summary(anpp)$coefficients[,1])%>%rename(anpp_est='summary(anpp)$coefficients[, 1]')
  comp_est <- as.data.frame(summary(comp)$coefficients[,1])%>%rename(comp_est='summary(comp)$coefficients[, 1]')
  evenness_est <- as.data.frame(summary(evenness)$coefficients[,1])%>%rename(evenness_est='summary(evenness)$coefficients[, 1]')
  rank_est <- as.data.frame(summary(rank)$coefficients[,1])%>%rename(rank_est='summary(rank)$coefficients[, 1]')
  richness_est <- as.data.frame(summary(richness)$coefficients[,1])%>%rename(richness_est='summary(richness)$coefficients[, 1]')
  species_est <- as.data.frame(summary(species)$coefficients[,1])%>%rename(species_est='summary(species)$coefficients[, 1]')
  
  #paste into dataframe
  singleSiteRegressionOutput <- cbind(anpp_pval, anpp_est, comp_pval, comp_est, evenness_pval, evenness_est, rank_pval, rank_est, richness_pval, richness_est, species_pval, species_est)%>%
    mutate(site_project_comm_trt=site_trt_list[i])%>%
    mutate(reg_factor=c('intercept', 'treatment_year', 'treatment_year_sq'))
  
  #merge to all other sites
  regressionOutput <- rbind(regressionOutput, singleSiteRegressionOutput)
  
}

regressionOutput <- regressionOutput%>%
  mutate(anpp_est=ifelse(anpp_pval>0.05, 0, anpp_est), comp_est=ifelse(comp_pval>0.05, 0, comp_est), evenness_est=ifelse(evenness_pval>0.05, 0, evenness_est), rank_est=ifelse(rank_pval>0.05, 0, rank_est), richness_est=ifelse(richness_pval>0.05, 0, richness_est), species_est=ifelse(species_pval>0.05, 0, species_est))

# write.csv(regressionOutput, 'C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\CoRRE\\C2E\\Products\\testing HRF\\corre results\\site_responses_time\\site_responses_time.csv', row.names=F)


### create for loop to produce ggplot2 graphs with only significant lines
#get the values from the regressionOutput dataframe for the lines for each site
#plot the data from correDataTrt and the regression with stat_function from the regressionOutput variables
#paste the five plots together into one "plot" file then print those to the files
for (i in seq_along(site_trt_list)) { 
  
  # create plot for each treatment in df
  anpp_intercept <- regressionOutput$
  anpp_plot <- 
    ggplot(subset(correDataTrt, site_project_comm_trt==site_trt_list[i]), aes(x=treatment_year, y=value)) + 
    geom_point() + 
    stat_function(fun=function(x){}, size=2, xlim=c(0,2)) +
  
  
  
  # save plots as .png
  ggsave(plot, file=paste('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\CoRRE\\C2E\\Products\\testing HRF\\corre results\\site_responses_time\\',
                          site_trt_list[i], ".png", sep=''), scale=3)
  
  # print plots to screen
  print(plot)
}



### create for loop to produce ggplot2 graphs with all lines (including non-significant)
for (i in seq_along(site_trt_list)) { 
  
  # create plot for each treatment in df 
  plot <- 
    ggplot(subset(correDataTrt, site_project_comm_trt==site_trt_list[i]), aes(x=treatment_year, y=value)) + 
    geom_point() + 
    stat_smooth(method="lm", formula = y ~ x + I(x^2), fill=NA) +
    facet_grid(rows=vars(factor_type), scales='free_y')
  
  # save plots as .png
  ggsave(plot, file=paste('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\CoRRE\\C2E\\Products\\testing HRF\\corre results\\site_responses_time\\',
                          site_trt_list[i], ".png", sep=''), scale=3)
  
  # print plots to screen
  print(plot)
}


