library(tidyverse)


dat_mult<-read.csv("~/Dropbox/C2E/Products/CommunityChange/March2018 WG/CORRE_Mult_Metrics_March2018.csv")%>%
  select(-X)

trt<-read.csv("~/Dropbox/converge_diverge/datasets/LongForm/ExperimentInformation_Nov2017.csv")%>%
  select(site_code, project_name, community_type, treatment,plot_mani)%>%
  unique()%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))

dat2 <- dat_mult%>%
  left_join(trt)

##looping through to get at t-tets for each c-t comparision using year as a replicate.
spc<-unique(dat2$site_project_comm)
ttest_output<-data.frame()

for (i in 1:length(spc)){
  
  subset<-dat2%>%
    filter(site_project_comm==spc[i])
  
  control<-subset%>%
    filter(plot_mani==0)
  
  treat_list<-unique(subset(subset, plot_mani>0)$treatment)
  
  for(j in 1:length(treat_list)) {
    spc1<-spc[i]
    trt <- treat_list[j]
    
    subset_trt<-subset%>%
      filter(treatment==treat_list[j])
    
    #dataset of two treatments    
    subset_t12<-rbind(control, subset_trt)
    
    out<-t.test(composition_change~treatment, data=subset_t12)

    output<-data.frame(site_project_comm = spc1,
                         treatment = trt,
                         p.value = out$p.value)
  
    ttest_output <- rbind(ttest_output, output)
  }}

##designate if sig
sig<-ttest_output%>%
  mutate(sig = ifelse(p.value < 0.1, 1, 0))%>%
  filter(sig ==1)%>%
  mutate(ttest = 1)%>%
  select(-p.value, -sig)

write.csv(sig, "~/Dropbox/C2E/Products/CommunityChange/March2018 WG/sig_ttest.csv", row.names= F)
          