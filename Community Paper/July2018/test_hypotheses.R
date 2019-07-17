### Addressing hypotheses in ecosphere paper
###
### Authors:  Meghan Avolio (meghan.avolio@jhu.edu)
### Last updated: Oct 30 2018

setwd("C:\\Users\\megha\\Dropbox\\")
theme_set(theme_bw(12))

library(tidyverse)
library(devtools)
install_github("NCEAS/codyn", ref = "sp_diff_test")
library(codyn)
library(vegan)

###getting ttest differences between rac trt-control and contorl-contorl comparisions

rac_diff <-read.csv("C2E\\Products\\CommunityChange\\March2018 WG\\CORRE_RAC_Diff_Metrics_June2019.csv")

rac_diff_control<- read.csv("C2E\\Products\\CommunityChange\\March2018 WG\\CORRE_RAC_Diff_control_Metrics_Jun2019.csv")


##get average C-T diff for each control replicate
rac_diff_mean <- rac_diff%>%
  group_by(site_project_comm, calendar_year, treatment, treatment2, plot_id)%>%
  summarize(richness_diff = mean(richness_diff),
            evenness_diff = mean(evenness_diff, na.rm = T),
            rank_diff=mean(rank_diff),
            species_diff=mean(species_diff))%>%
  mutate(spc_yr=paste(site_project_comm, calendar_year, sep="_"))


# ##get average C-C diff for each control replicate for each year
rac_diff_control2 <- rac_diff_control%>%
  mutate(spc_yr=paste(site_project_comm, calendar_year, sep="_"))


spc_yr_vec<-unique(rac_diff_mean$spc_yr)

control_means<-data.frame()

for (i in 1:length(spc_yr_vec)){
  subset<-rac_diff_control2%>%
    filter(spc_yr==spc_yr_vec[i])

  control_plots<-unique(subset(rac_diff_mean, spc_yr==spc_yr_vec[i])$plot_id)

  for (i in 1:length(control_plots)){

    subset2<-subset(subset, plot_id==as.character(control_plots[i])|plot_id2==as.character(control_plots[i]))

    average<-subset2%>%
     summarise(C_richness_diff=mean(richness_diff),
                C_evenness_diff=mean(evenness_diff, na.rm = T),
                C_rank_diff = mean(rank_diff),
                C_species_diff=mean(species_diff))%>%
      mutate(site_project_comm = unique(subset$site_project_comm),
             calendar_year = unique(subset$calendar_year),
             plot_id=as.character(control_plots[i]))

    control_means <-rbind(control_means, average)

  }
}


#merge with rank_diff means

rac_diff_c_t<-rac_diff_mean%>%
  left_join(control_means)


##loop through each experiment, year, treatment and ask is if C-T comparisons are different than C-C comparisions with a ttest

#drop problematic experiment, treatment time points
#CAR_salt marsh_SalCus_1999_NPK data is the same for most replicates
#CDR_e001_A_1988_8 evenness is NAN for all trts
#CUL_Culardoch_0_2000_N10burn same as above
#CUL_Culardoch_0_2000_N10burnclip same as above
#CUL_Culardoch_0_2000_N20burn same as above
#CUL_Culardoch_0_2000_N20burnclip same as above
#CUL_Culardoch_0_2000_N50burn same as above
#NANT_wet_Broad_BRC_S_2003_0N1P NAN for controls
#NANT_wet_Broad_BRC_S_2003_1N0P NAN for controls
#NANT_wet_Broad_BRC_S_2003_1N1P NAN for controls
#SERC_CXN_0_2013_t3 evenness NAN for trts
#SERC_TMECE_MX_1998_E richness and sp diff is 0 for all controls
#SERC_TMECE_MX_2001_E same as above
#SERC_TMECE_MX_2002_E same as above
#SERC_TMECE_MX_2004_E same as above
#SERC_TMECE_MX_2006_E same as above
#SERC_TMECE_MX_2008_E same as above
#SERC_TMECE_MX_2009_E same as above
#SERC_TMECE_SC_2010_E same as above
#SERC_TMECE_SC_2011_E same as above
#SERC_TMECE_SP_1997_E problems across the board.
#SERC_TMECE_SP_1999_E evenness is NAN for controls and trts
#SERC_TMECE_SP_2000_E evenness is NAN for controls
#SERC_TMECE_SP_2001_E same as above

rac_diff_c_t_sub<-rac_diff_c_t%>%
  mutate(spc_yr_trt=paste(site_project_comm, calendar_year, treatment2, sep="_"))%>%
  filter(spc_yr_trt!="CAR_salt marsh_SalCus_1999_NPK"&
           spc_yr_trt!="CDR_e001_A_1988_8"&
           spc_yr_trt!="CUL_Culardoch_0_2000_N10burn"&
           spc_yr_trt!="CUL_Culardoch_0_2000_N10burnclip"&
           spc_yr_trt!="CUL_Culardoch_0_2000_N20burn"&
           spc_yr_trt!="CUL_Culardoch_0_2000_N20burnclip"&
           spc_yr_trt!="CUL_Culardoch_0_2000_N50burn"&
           spc_yr_trt!="NANT_wet_Broad_BRC_S_2003_0N1P"&
           spc_yr_trt!="NANT_wet_Broad_BRC_S_2003_1N0P"&
           spc_yr_trt!="NANT_wet_Broad_BRC_S_2003_1N1P"&
           spc_yr_trt!="SERC_CXN_0_2013_t3"&
           spc_yr_trt!="SERC_TMECE_MX_1998_E"&
           spc_yr_trt!="SERC_TMECE_MX_2001_E"&
           spc_yr_trt!="SERC_TMECE_MX_2002_E"&
           spc_yr_trt!="SERC_TMECE_MX_2004_E"&
           spc_yr_trt!="SERC_TMECE_MX_2006_E"&
           spc_yr_trt!="SERC_TMECE_MX_2008_E"&
           spc_yr_trt!="SERC_TMECE_MX_2009_E"&
           spc_yr_trt!="SERC_TMECE_SC_2010_E"&
           spc_yr_trt!="SERC_TMECE_SC_2011_E"&
           spc_yr_trt!="SERC_TMECE_SP_1997_E"&
           spc_yr_trt!="SERC_TMECE_SP_1999_E"&
           spc_yr_trt!="SERC_TMECE_SP_2000_E"&
           spc_yr_trt!="SERC_TMECE_SP_2001_E")


spc_yr_trt_vec<-unique(rac_diff_c_t_sub$spc_yr_trt)

ttests<-data.frame()

for (i in 1:length(spc_yr_trt_vec)){
 subset<-rac_diff_c_t_sub%>%
    filter(spc_yr_trt==spc_yr_trt_vec[i])

 rd<-t.test(subset$richness_diff, subset$C_richness_diff)
 ed<-t.test(subset$evenness_diff, subset$C_evenness_diff)
 rankd<-t.test(subset$rank_diff, subset$C_rank_diff)
 spd<-t.test(subset$species_diff, subset$C_species_diff)

 ttest_temp <- data.frame(
   site_project_comm = unique(subset$site_project_comm),
   treatment = unique(subset$treatment2),
   calendar_year = unique(subset$calendar_year),
   rich_pval =  rd$p.value,
   even_pval =  ed$p.value,
   rank_pval =  rankd$p.value,
   spdiff_pval =  spd$p.value)

 ttests<-rbind(ttests, ttest_temp)

 }

write.csv(ttests, "C2E\\Products\\CommunityChange\\March2018 WG\\RAC_diff_CT_ttests_Jun2019.csv", row.names = F)

####linking RAC differences with compositon/disperison differences

#there are more perm_output b/c did not subset CDR e001/e002
perm_output<-read.csv("C2E\\Products\\CommunityChange\\March2018 WG\\permanova_permdisp_outputJul2019.csv")

mult_diff <- read.csv("C2E\\Products\\CommunityChange\\March2018 WG\\CORRE_Mult_diff_Metrics_Jun2019.csv")%>%
  mutate(treatment = treatment2)

CT_ttests<- read.csv("C2E\\Products\\CommunityChange\\March2018 WG\\RAC_diff_CT_ttests_Jun2019.csv")


#merge perm_output and mult_diff to set up the six scenarios and drop what did not work for ttests.
#1 = no change comp, no change disp
#2 = no change comp, increase disp (T > C)
#3 = no change comp, decrease disp (C > T)
#4 = change comp, no change disp
#5 = change comp, increase disp (T > C )
#6 = change comp, decrease disp (C > T)

scenarios<-perm_output%>%
  right_join(mult_diff)%>%
  mutate(scenario=ifelse(perm_Pvalue>0.0501&disp_Pvalue>0.0501, 1, 
                         ifelse(perm_Pvalue>0.0501&disp_Pvalue<0.0501&greater_disp == "T", 2,
                         ifelse(perm_Pvalue>0.0501&disp_Pvalue<0.0501&greater_disp == "C", 3,
                         ifelse(perm_Pvalue<0.0501&disp_Pvalue>0.0501,4,
                         ifelse(perm_Pvalue<0.0501&disp_Pvalue<0.0501&greater_disp == "T",5,
                         ifelse(perm_Pvalue<0.0501&disp_Pvalue<0.0501&greater_disp == "C",6,999)))))))



##combining to see proportion is different for each RAC metric
RAC_diff_outcomes <- scenarios%>%
  right_join(CT_ttests)%>%
  mutate(rich=ifelse(rich_pval<0.0501, 1, 0),
         even=ifelse(even_pval<0.0601, 1, 0),
         rank=ifelse(rank_pval<0.0501, 1, 0),
         spdiff=ifelse(spdiff_pval<0.501, 1, 0))

num_scen<- RAC_diff_outcomes%>%
  group_by(scenario)%>%
  summarize(n=length(scenario))%>%
  na.omit()

prop_diff<-RAC_diff_outcomes%>%
  na.omit()%>%
  group_by(scenario)%>%
  summarize_at(vars(rich, even, rank, spdiff), funs(sum))%>%
  gather(metric, num, rich:spdiff)%>%
  mutate(prop = ifelse(scenario==1, num/1762, ifelse(scenario==2, num/98, ifelse(scenario==3, num/132, ifelse(scenario==4, num/586, ifelse(scenario==5, num/152, ifelse(scenario==6, num/153, 999)))))))%>%
  mutate(notsig=1-prop)%>%
  select(-num)%>%
  gather(sig, proportion, notsig:prop)

theme_set(theme_bw(12))
ggplot(data=prop_diff, aes(x = metric, y = proportion, fill=sig))+
  geom_bar(stat="identity")+
  scale_fill_manual(name = "", labels=c("No Difference", "C-T Different"), values=c("gray","darkgreen"))+
  scale_x_discrete(limits=c("rich", "even", "rank", "spdiff"))+
  ylab("Proportion")+
  xlab("RAC Difference Metric")+
  facet_wrap(~scenario, ncol=3)+
  geom_hline(yintercept = 0.5)

#people like this figure less
ggplot(data=prop_diff, aes(x = scenario, y = proportion, fill=sig))+
  geom_bar(stat="identity")+
  scale_fill_manual(name = "", labels=c("No Difference", "C-T Different"), values=c("gray","darkgreen"))+
  ylab("Proportion")+
  xlab("Community Difference Scenario")+
  scale_x_continuous(breaks = c(1:6))+
  facet_wrap(~metric, ncol=2)

###redoing this with out the rare species
perm_outputnorare<-read.csv("C2E\\Products\\CommunityChange\\March2018 WG\\permanova_permdisp_output_norare_Jul2019.csv")

mult_diffnorare <- read.csv("C2E\\Products\\CommunityChange\\March2018 WG\\CORRE_Mult_diff_Metrics_norare_Jun2019.csv")%>%
  mutate(treatment = treatment2)

#merge perm_output and mult_diff to set up the six scenarios and drop what did not work for ttests.
#1 = no change comp, no change disp
#2 = no change comp, increase disp (T > C)
#3 = no change comp, decrease disp (C > T)
#4 = change comp, no change disp
#5 = change comp, increase disp (T > C )
#6 = change comp, decrease disp (C > T)

scenariosnorare<-perm_outputnorare%>%
  right_join(mult_diffnorare)%>%
  mutate(scenario=ifelse(perm_Pvalue>0.0501&disp_Pvalue>0.0501, 1, 
                         ifelse(perm_Pvalue>0.0501&disp_Pvalue<0.0501&greater_disp == "T", 2,
                                ifelse(perm_Pvalue>0.0501&disp_Pvalue<0.0501&greater_disp == "C", 3,
                                       ifelse(perm_Pvalue<0.0501&disp_Pvalue>0.0501,4,
                                              ifelse(perm_Pvalue<0.0501&disp_Pvalue<0.0501&greater_disp == "T",5,
                                                     ifelse(perm_Pvalue<0.0501&disp_Pvalue<0.0501&greater_disp == "C",6,999)))))))

num_scenmprare<- scenarios%>%
  group_by(scenario)%>%
  summarize(n=length(scenario))%>%
  na.omit()

####what is the role of species differneces?
abund_diff <- read.csv("C2E\\Products\\CommunityChange\\March2018 WG\\CORRE_Abund_Diff_June2019.csv")%>%
  mutate(treatment = treatment2)

numrep<-abund_diff%>%
  select(site_project_comm, treatment, plot_id2)%>%
  unique()%>%
  group_by(site_project_comm, treatment)%>%
  summarise(rep=length(plot_id2))

scen56<-scenarios%>%
  filter(scenario==5|scenario==6)%>%
  select(site_project_comm, treatment, calendar_year, scenario)

abund_diff56<-abund_diff%>%
  right_join(scen56)%>%
  filter(difference>0.2)%>%
  select(calendar_year, plot_id2, treatment, genus_species, site_project_comm, scenario)%>%
  unique()%>%
  ungroup()%>%
  group_by(calendar_year, treatment, genus_species, site_project_comm, scenario)%>%
  summarize(numplot=length(plot_id2))%>%
  left_join(numrep)%>%
  mutate(propinc=numplot/rep)

ave<-abund_diff56%>%
  group_by(scenario)%>%
  summarize(ave=mean(propinc),
            std=sd(propinc),
            n=length(propinc))%>%
  mutate(se=std/sqrt(n))

s5<-abund_diff56%>%
  filter(scenario==5)
s6<-abund_diff56%>%
  filter(scenario==6)

s5t<-s5$propinc
s6t<-s6$propinc

t.test(s5t, s6t)
