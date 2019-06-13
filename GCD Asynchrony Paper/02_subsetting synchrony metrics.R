### Filter synchrony metrics and sites
###
### Author: Kevin wilcox (kevin.wilcox@uwyo.edu)
### Created: May 14th 2019, last updated: May 14th, 2019

### Set up workspace
setwd("C:\\Users\\wilco\\Desktop\\Working groups\\C2E_May2019\\asynchrony\\data\\")
library(tidyverse)

### Read in old experiment
old_metrics_vec <- read.csv("ecolett paper_site_proj_comm vector.csv") %>%
  .$site_proj_comm
### Read in site information
site_info <- read.csv('SiteExperimentDetails_March2019.csv')

### Read in synchrony metrics and experimental information
synch_metrics <- read.csv("..\\Synchrony metrics response ratios_long form_13May2019.csv") %>%
#  filter(plot_mani <= 1) %>%
  mutate(site_proj_comm = paste(site_code, project_name, community_type, sep="_")) %>%
  filter(site_proj_comm %in% c(as.character(old_metrics_vec), "CHY_EDGE_0","HYS_EDGE_0","KNZ_EDGE_0","SGS_EDGE_0","SEV_EDGE_EB","SEV_EDGE_EG")) %>%
  filter(pulse == 0) %>% ## remove all pulse treatments
  filter(plant_trt == 0) %>% ## remove seeding treatments
  mutate(trt_type2=ifelse(trt_type=="N"|trt_type=="control","N", 
                          ifelse(trt_type=="P", "P", 
                                 ifelse(trt_type=="CO2", "CO2",
                                        ifelse(trt_type=="irr", "Irrigation",
                                               ifelse(trt_type=="temp", "Temperature", 
                                                      ifelse(trt_type=="N*P"|trt_type=="mult_nutrient", "Mult. Nuts.", 
                                                             ifelse(trt_type=="drought", "drought", 
                                                                    ifelse(trt_type=="CO2*temp", "CO2*temp", 
                                                                           ifelse(trt_type=="drought*temp", "drought*temp", 
                                                                                  ifelse(trt_type=="irr*temp", "irr*temp",
                                                                                         ifelse(trt_type=="irr*CO2*temp"|trt_type=="N*CO2*temp"|trt_type=="N*irr*temp"|trt_type=="N*irr*CO2*temp", "mult_res*temp", 
                                                                                                ifelse(trt_type=="irr*herb_removal"|trt_type=="irr*plant_mani"|trt_type=="irr*plant_mani*herb_removal", "irr*NR", 
                                                                                                       ifelse(trt_type=="herb_removal"|trt_type=="till"|trt_type=="mow_clip"|trt_type=="burn"|trt_type=="plant_mani"|trt_type=="stone"|trt_type=="graze"|trt_type=="burn*graze"|trt_type=="fungicide"|trt_type=="plant_mani*herb_removal"|trt_type=="burn*mow_clip", "NR", 
                                                                                                              ifelse(trt_type=="precip_vari", "Precip. Vari.",  
                                                                                                                     ifelse(trt_type=="N*plant_mani"|trt_type=="N*burn"|trt_type=="N*mow_clip"|trt_type=="N*till"|trt_type=="N*stone"|trt_type=="N*burn*graze"|trt_type=="N*burn*mow_clip", "N*NR", 
                                                                                                                            ifelse(trt_type=="N*temp", "N*temp", 
                                                                                                                                   ifelse(trt_type=="N*CO2", "N*CO2",
                                                                                                                                          ifelse(trt_type=="irr*CO2", "irr*CO2",
                                                                                                                                                 ifelse(trt_type=="N*irr", "N*irr",
                                                                                                                                                        ifelse(trt_type=="mult_nutrient*herb_removal"|trt_type=="mult_nutrient*fungicide"|trt_type=="N*P*burn*graze"|trt_type=="N*P*burn"|trt_type=="*P*mow_clip"|trt_type=="N*P*burn*mow_clip"|trt_type=="N*P*mow_clip", "mult_nutrients*NR",
                                                                                                                                                               ifelse(trt_type=="P*mow_clip"|trt_type=="P*burn"|trt_type=="P*burn*graze"|trt_type=="P*burn*mow_clip", "P*NR", 
                                                                                                                                                                      ifelse(trt_type=="precip_vari*temp", "precip_vari*temp",
                                                                                                                                                                             ifelse(trt_type=="light","light",
                                                                                                                                                                             ifelse(trt_type=="N*irr*CO2", "mult_res", 999)))))))))))))))))))))))))

### Vector of treatment types that have at least 5 sites where they are implemented
trt_type2_vector <- c(
  "drought",
  "Irrigation",
  "drought*temp",
  "irr*temp",
  "Temperature",
  "Mult. Nuts.",
  "N",
  "P",
  "Precip. Vari.",
  "N*irr",
  "NR"
)
synch_metrics_sub <- synch_metrics %>%
  filter(trt_type2 %in% trt_type2_vector)

synch_metrics_sub$trt_type2 <- factor(synch_metrics_sub$trt_type2,
                                      levels = c("drought","Irrigation","Precip. Vari.",
                                                 "Temperature","N","P","Mult. Nuts.",
                                                 "drought*temp","irr*temp","N*irr", "NR"))

###
### means and se's for response ratios
###
unique(synch_metrics$trt_type2)
synchrony_rr_summary <- synch_metrics_sub %>%
#  filter(trt_type2 %in% c("CO2", "drought", "irr", "N", "P", "temp")) %>%
  group_by(trt_type2, metric_name) %>%
  summarize(lnRR_mean = mean(lnRR, na.rm=T),
            lnRR_se = sd(lnRR)/sqrt(length(lnRR)))

synch_metrics_water <- synch_metrics_sub %>%
  filter(trt_type2 %in% c("drought","Irrigation","Precip. Vari.", "Temperature",
                          "drought*temp","irr*temp"))

synch_metrics_nuts <- synch_metrics_sub %>%
  filter(trt_type2 %in% c("N","P","Mult. Nuts.", "N*irr"))


### 
### Visualize
###
## boxplots and jitter
ggplot(synch_metrics_sub, aes(x=factor(trt_type2), y=lnRR, col=factor(trt_type2))) +
  geom_boxplot() +
  geom_jitter(shape=21) +
  facet_wrap(~metric_name, scales="free") +
  theme_bw() +
  geom_hline(yintercept=0) + 
  theme(axis.text.x = element_text(angle=60, hjust=1))

## boxplots and jitter -- Water
ggplot(synch_metrics_water, aes(x=factor(trt_type2), y=lnRR, col=factor(trt_type2))) +
  geom_boxplot() +
  geom_jitter(shape=21) +
  facet_wrap(~metric_name, scales="free") +
  theme_bw() +
  geom_hline(yintercept=0) + 
  theme(axis.text.x = element_text(angle=60, hjust=1))

## boxplots and jitter - Nutrients
ggplot(synch_metrics_nuts, aes(x=factor(trt_type2), y=lnRR, col=factor(trt_type2))) +
  geom_boxplot() +
  geom_jitter(shape=21) +
  facet_wrap(~metric_name, scales="free") +
  theme_bw() +
  geom_hline(yintercept=0) + 
  theme(axis.text.x = element_text(angle=60, hjust=1))

## means and se's
ggplot(synchrony_rr_summary, aes(x=trt_type2, y=lnRR_mean, ymin=lnRR_mean-lnRR_se,ymax=lnRR_mean+lnRR_se, col=factor(trt_type2))) +
  geom_point() +
  geom_errorbar(width=0) +
  facet_wrap(~metric_name, scales="free") +
  theme_bw() +
  geom_hline(yintercept=0) +
  theme(axis.text.x = element_text(angle=60, hjust=1))


###
### Look at response ratios across environmental gradients
###
synchrony_rr_long_with_env <- synch_metrics_sub %>%
  full_join(site_info, by=c("site_code", "project_name", "community_type"))


### Nitrogen alone
ggplot(filter(synchrony_rr_long_with_env, trt_type2=="N"), aes(x=MAP, y=lnRR)) +
  geom_hline(yintercept=0) +
  geom_point(size=2) +
  geom_smooth(method='lm') +
  facet_wrap(~metric_name) +
  theme_bw()

ggplot(filter(synchrony_rr_long_with_env, trt_type2=="N"), aes(x=MAT, y=lnRR)) +
  geom_hline(yintercept=0) +
  geom_point(size=2) +
  geom_smooth(method='lm') +
  facet_wrap(~metric_name) +
  theme_bw()

ggplot(filter(synchrony_rr_long_with_env, trt_type2=="N"), aes(x=anpp, y=lnRR)) +
  geom_hline(yintercept=0) +
  geom_point(size=2) +
  geom_smooth(method='lm') +
  facet_wrap(~metric_name) +
  theme_bw()

ggplot(filter(synchrony_rr_long_with_env, trt_type2=="N"), aes(x=rrich, y=lnRR)) +
  geom_hline(yintercept=0) +
  geom_point(size=2) +
  geom_smooth(method='lm') +
  facet_wrap(~metric_name) +
  theme_bw()

ggplot(filter(synchrony_rr_long_with_env, trt_type2=="N"), aes(x=experiment_length, y=lnRR)) +
  geom_hline(yintercept=0) +
  geom_point(size=2) +
  geom_smooth(method='lm') +
  facet_wrap(~metric_name) +
  theme_bw()


###
ggplot(filter(synchrony_rr_long_with_env, trt_type2=="drought"), aes(x=MAP, y=lnRR)) +
  geom_hline(yintercept=0) +
  geom_point(size=2) +
  geom_smooth(method='lm') +
  facet_wrap(~metric_name) +
  theme_bw()

ggplot(filter(synchrony_rr_long_with_env, trt_type2=="drought"), aes(x=MAT, y=lnRR)) +
  geom_hline(yintercept=0) +
  geom_point(size=2) +
  geom_smooth(method='lm') +
  facet_wrap(~metric_name) +
  theme_bw()

ggplot(filter(synchrony_rr_long_with_env, trt_type2=="drought"), aes(x=anpp, y=lnRR)) +
  geom_hline(yintercept=0) +
  geom_point(size=2) +
  geom_smooth(method='lm') +
  facet_wrap(~metric_name) +
  theme_bw()

ggplot(filter(synchrony_rr_long_with_env, trt_type2=="drought"), aes(x=rrich, y=lnRR)) +
  geom_hline(yintercept=0) +
  geom_point(size=2) +
  geom_smooth(method='lm') +
  facet_wrap(~metric_name) +
  theme_bw()

ggplot(filter(synchrony_rr_long_with_env, trt_type2=="drought"), aes(x=experiment_length, y=lnRR)) +
  geom_hline(yintercept=0) +
  geom_point(size=2) +
  geom_smooth(method='lm') +
  facet_wrap(~metric_name) +
  theme_bw()


##basic stats on what we are doing.
df_exp<-synch_metrics%>%
  select(site_proj_comm)%>%
  unique()

df_trt<-synch_metrics%>%
  select(site_proj_comm, treatment, trt_type2)%>%
  unique()

##summary of treatment for table 1
sum_trts<-df_trt%>%
  select(site_proj_comm, trt_type2)%>%
  separate(site_proj_comm, into=c("site_code", "project_name", "community_type"), sep="_")%>%
  group_by(trt_type2)%>%
  summarize(n=length(trt_type2))

sum_sites<-synch_metrics%>%
  separate(site_proj_comm, into=c("site_code", "project_name", "community_type"), sep="_")%>%
  select(site_code, trt_type2)%>%
  unique()%>%
  group_by(trt_type2)%>%
  summarize(n=length(trt_type2))

length(
  synch_metrics %>%
  select(site_proj_comm) %>%
  unique() %>%
    .$site_proj_comm
  )
