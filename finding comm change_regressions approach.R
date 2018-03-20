library(tidyverse)

setwd('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\C2E\\Products\\CommunityChange\\March2018 WG')

multMetrics <- read.csv('CORRE_Mult_Metrics_March2018.csv')
treatments <- read.csv('ExperimentInformation_Nov2017.csv')%>%
  mutate(site_project_comm=as.factor(paste(site_code, project_name, community_type, sep='_')))%>%
  select(-X)


length <- multMetrics%>%
  group_by(site_project_comm, treatment)%>%
  summarise(length=length(treatment_year2))%>%
  ungroup()

multMetricsTrt <- multMetrics%>%
  left_join(treatments)%>%
  left_join(length)%>%
  filter(length>3)%>%
  select(treatment_year, treatment, composition_change, site_project_comm, plot_mani)

  

##test for sig diff between trt-control slopes
spc2<-unique(multMetricsTrt$site_project_comm)
test.lm<-data.frame()
for (i in 1:length(spc2)){
  subset<-multMetricsTrt%>%
    filter(site_project_comm==spc2[i])
  control<-subset%>%
    filter(plot_mani==0)
  treat<-subset%>%
  filter(plot_mani!=0)
trt_list<-unique(treat$treatment)
for (j in 1:length(trt_list)){
  subset2<-treat%>%
    filter(treatment==trt_list[j])
  trt<-trt_list[j]
  ct<-rbind(subset2, control)
  ct.lm<-lm(composition_change~treatment_year*treatment, data=ct)
  output.lm<-data.frame(site_project_comm=unique(subset$site_project_comm),
                        treatment=trt,
                        est_intercept=summary(ct.lm)$coef[3, "Estimate"],
                        val_intercept=summary(ct.lm)$coef[3,"Pr(>|t|)"],
                        est_slope=summary(ct.lm)$coef[4, "Estimate"],
                        val_slope=summary(ct.lm)$coef[4,"Pr(>|t|)"])
  test.lm<-rbind(test.lm, output.lm)
}
}

list <- test.lm%>%
  mutate(sig_intercept=ifelse(val_intercept<0.05, 1, 0), sig_slope=ifelse(val_slope<0.05, 1, 0), sig_any=(sig_intercept+sig_slope), method='regression')%>%
  filter(sig_any>0)%>%
  select(site_project_comm, treatment, sig_intercept, sig_slope, method)

write.csv(list, 'treatments_sig regression.csv')
