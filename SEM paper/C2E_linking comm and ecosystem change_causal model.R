library(grid)
library(PerformanceAnalytics)
library(sandwich)
library(lmtest)
library(tidyverse)

# ??? ask Laura
# source("stdErrHelperFunctions.R")


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




#######CAUSAL MODEL: use differences between treatment and control plots in each year (horizontal arrows)--------------

#model
compositionModelNitrogen = anpp_pdiff ~
  # dummy variables
  site_project_comm_trt + #intercept
  site_project_comm_trt:treatment_year + #slope
  # variables
  n +
  richness_difference +
  evenness_diff +
  rank_difference +
  species_difference

# fit model
modelTable1 = lm(compositionModelNitrogen, data = subset(correSEMdataTrt, n>0))
# Compute clustered standard errors, and clean up display names for parameters
summary(modelTable1)



#######REGRESSION MODEL: use differences between treatment and control plots in each year (horizontal arrows)--------------

#model
compositionModelNitrogen = anpp_pdiff ~
  richness_difference*n +
  evenness_diff*n +
  rank_difference*n +
  species_difference*n

# fit model
modelTable1 = lm(compositionModelNitrogen, data = subset(correSEMdataTrt, n>0))
# Compute clustered standard errors, and clean up display names for parameters
summary(modelTable1)




#######figures
ggplot(data=subset(correSEMdataTrt, trt_type_2=='added'), aes(x=treatment_year, y=anpp_pdiff, color=trt_type)) +
  geom_line()
