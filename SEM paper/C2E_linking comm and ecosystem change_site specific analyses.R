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
  select(site_project_comm_trt, treatment, treatment_year, anpp_pdiff, richness_difference, evenness_diff, rank_difference, species_difference, composition_diff)%>%
  group_by(site_project_comm_trt, treatment, treatment_year)%>%
  gather(key='factor_type', value='value', anpp_pdiff:composition_diff, na.rm=T)%>%
  ungroup()


###loop across all sites
#generate site list
site_trt_list <- unique(correDataTrt$site_project_comm_trt)

# create for loop to produce ggplot2 graphs 
for (i in seq_along(site_trt_list)) { 
  
  # create plot for each county in df 
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

#really we should find the regressions that are significant through time, and then only plot those

