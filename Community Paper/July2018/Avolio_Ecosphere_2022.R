### Addressing hypotheses in Ecosphere paper
###
### This code uses output from RAC Differences.R and permanvoa_permdisp loop.R
###
### Authors:  Meghan Avolio (meghan.avolio@jhu.edu)
### Last updated: Nov 10 2021

setwd("C:\\Users\\mavolio2\\Dropbox\\")
setwd("E:\\Dropbox\\")

library(tidyverse)
library(devtools)
library(codyn)
library(vegan)
library(ggrepel)

theme_set(theme_bw(8))

###step 1. Run RAC Differences overall and for controls only. This is different code and the data is not public.

###Step 2: Do ttests. Getting ttest differences between rac trt-control and contorl-contorl comparisions

rac_diff <-read.csv("C2E\\Products\\Testing Hypots\\CORRE_RAC_Diff_Metrics_Oct2021.csv")

rac_diff_control<- read.csv("C2E\\Products\\Testing Hypots\\CORRE_RAC_Diff_control_Metrics_Oct2021.csv")


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


###get averages
ave_diff<-rac_diff_c_t %>% 
  select(-plot_id) %>% 
  group_by(site_project_comm, calendar_year, treatment, treatment2, spc_yr) %>% 
  summarise_all(funs(mean))


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

write.csv(ttests, "C2E\\Products\\Testing Hypots\\RAC_diff_CT_ttests_Oct2021.csv", row.names = F)

##step 3: run permanova_perm disp loop.R - this data is not public.

###step 4: read in all necessary files and correct for multiple hypothesis testing

perm_output<-read.csv("C2E\\Products\\Testing Hypots\\permanova_permdisp_outputOct2021.csv")%>%
  gather(measure, pval, perm_Pvalue:disp_Pvalue)%>%
  group_by(site_project_comm)%>%
  mutate(adjp=p.adjust(pval, method="BH", n=length(site_project_comm)))%>%
  select(-pval)%>%
  spread(measure, adjp)

num_studies<-perm_output%>%
  select(site_project_comm)%>%
  unique()

mult_diff <- read.csv("C2E\\Products\\Testing Hypots\\CORRE_Mult_diff_Metrics_Oct2021.csv")%>%
  mutate(treatment = treatment2)%>%
  filter(site_project_comm!="NGBER_gb_0"|treatment2!="AMBIENT")

#there are 24 fewer c-t comparisons because we were unable to run the t-test for 24 comparisons.
CT_ttests<- read.csv("C2E\\Products\\Testing Hypots\\RAC_diff_CT_ttests_Oct2021.csv")%>%
  filter(site_project_comm!="NGBER_gb_0"|treatment!="AMBIENT")%>%
  gather(measure, pval, rich_pval:spdiff_pval)%>%
  group_by(site_project_comm)%>%
  mutate(adjp=p.adjust(pval, method="BH", n=length(site_project_comm)))%>%
  select(-pval)%>%
  spread(measure, adjp)

#this is the full dataset of all C-T comparision for which it worked for all measures and p-values are corrected.
fulldataset<-perm_output%>%
  right_join(mult_diff)%>%
  na.omit()%>%
  right_join(CT_ttests)%>%
  na.omit()

site<-fulldataset%>%
  select(site_project_comm)%>%
  unique()

####Step 5: linking RAC differences with composition/disperison differences

#merge perm_output and mult_diff to set up the six scenarios and drop what did not work for ttests.
#1 = no change comp, no change disp
#2 = no change comp, increase disp (T > C)
#3 = no change comp, decrease disp (C > T)
#4 = change comp, no change disp
#5 = change comp, increase disp (T > C )
#6 = change comp, decrease disp (C > T)

scenarios<-fulldataset%>%
  mutate(scenario=ifelse(perm_Pvalue>0.0501&disp_Pvalue>0.0501, 1, 
                         ifelse(perm_Pvalue>0.0501&disp_Pvalue<0.0501&greater_disp == "T", 2,
                         ifelse(perm_Pvalue>0.0501&disp_Pvalue<0.0501&greater_disp == "C", 3,
                         ifelse(perm_Pvalue<0.0501&disp_Pvalue>0.0501,4,
                         ifelse(perm_Pvalue<0.0501&disp_Pvalue<0.0501&greater_disp == "T",5,
                         ifelse(perm_Pvalue<0.0501&disp_Pvalue<0.0501&greater_disp == "C",6,999)))))))%>%
  na.omit

num_scen<- scenarios%>%
  group_by(scenario)%>%
  summarize(n=length(scenario))%>%
  mutate(pct=round(100*(n/2831), 1))

df2 <- num_scen %>% 
  mutate(csum = rev(cumsum(rev(pct))), 
         pos = pct/2 + lead(csum, 1),
         pos = if_else(is.na(pos), pct/2, pos))

###Figure 3 make a pie chart
blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=10, face="bold")
  )

fig3<-ggplot(num_scen, aes(x="",y=pct, fill=as.factor(scenario)))+
  geom_col(width=1, color=1)+
  coord_polar(theta = "y", start=100)+
  scale_fill_brewer(name="Scenario", palette="Pastel2")+
  blank_theme+
  theme(axis.text.x = element_blank())+
  geom_label_repel(data=df2, aes(y = pos, label = paste0(pct, "%")),
                   size = 2.5, nudge_x = 1, show.legend = FALSE)+
  theme(legend.title = element_text(size=8), legend.text = element_text(size=8))
            
ggsave("C:\\Users\\mavolio2\\OneDrive - Johns Hopkins\\Manuscripts\\C2E - testing ecosphere hypotheses\\Revison\\Fig3.jpeg",fig3, units="in", dpi=300, width = 3, height = 3)

##combining to see proportion is different for each RAC metric
RAC_diff_outcomes <- scenarios%>%
  mutate(rich=ifelse(rich_pval<0.0501, 1, 0),
         even=ifelse(even_pval<0.0601, 1, 0),
         rank=ifelse(rank_pval<0.0501, 1, 0),
         spdiff=ifelse(spdiff_pval<0.0501, 1, 0),
         outcome=ifelse(rich==1&even==0&rank==0&spdiff==0, "Richness", ifelse(rich==0&even==1&rank==0&spdiff==0, 'Evenness', ifelse(rich==0&even==0&rank==1&spdiff==0, "Ranks", ifelse(rich==0&even==0&rank==0&spdiff==1, "Species", ifelse(rich==1&even==1&rank==0&spdiff==0, "Rich&Even", ifelse(rich==1&even==0&rank==1&spdiff==0|rich==0&even==1&rank==1&spdiff==0|rich==1&even==1&rank==1&spdiff==0, "Rank&Other", ifelse(rich==1&even==0&rank==0&spdiff==1|rich==0&even==1&rank==0&spdiff==1|rich==1&even==1&rank==0&spdiff==1, "Species&Other", ifelse(rich==1&even==0&rank==1&spdiff==1|rich==0&even==1&rank==1&spdiff==1|rich==1&even==1&rank==1&spdiff==1|rich==0&even==0&rank==1&spdiff==1, "Ranks&Species&Other", "999")))))))))

Check<-RAC_diff_outcomes%>%
  filter(outcome=="999")%>%
  mutate(sum=rich+even+rank+spdiff)

onesigscen<-RAC_diff_outcomes%>%
  mutate(onesig=ifelse(rich|even|rank|spdiff==1, 1, 0))%>%
  replace(is.na(.), 0)%>%
  group_by(scenario, onesig)%>%
  summarize(num=length(onesig))%>%
  mutate(prop = ifelse(scenario==1, num/2194, ifelse(scenario==2, num/80, ifelse(scenario==3, num/99, ifelse(scenario==4, num/300, ifelse(scenario==5, num/58, ifelse(scenario==6, num/100, 999)))))))

sp5<-data.frame("scenario"=5, "outcome"= "Species", "sig"="notsig", "proportion"=1)

prop_diff<-RAC_diff_outcomes%>%
  na.omit()%>%
  filter(outcome!="999")%>%
  group_by(scenario, outcome)%>%
  summarize(num=length(outcome))%>%
  mutate(prop = ifelse(scenario==1, num/2194, ifelse(scenario==2, num/80, ifelse(scenario==3, num/99, ifelse(scenario==4, num/300, ifelse(scenario==5, num/58, ifelse(scenario==6, num/100, 999)))))))%>%
  mutate(notsig=1-prop)%>%
  select(-num)%>%
  gather(sig, proportion, notsig:prop)%>%
  bind_rows(sp5)%>%
  mutate(percent=proportion*100)



theme_set(theme_bw(12))
theme_update(panel.grid.major=element_blank(), panel.grid.minor=element_blank())

#this is the will make each panel in figure 4, just need to change the scenarios on L299 
ggplot(data=subset(prop_diff, scenario==6&sig=="prop"), aes(x = outcome, y = percent))+
  geom_bar(stat="identity")+
  #scale_fill_manual(name = "", labels=c("No Difference", "C-T Different"), values=c("gray90","gray48"))+
  scale_x_discrete(limits=c("Richness", "Evenness", "Ranks", "Species", "Rich&Even", "Rank&Other", "Species&Other", "Ranks&Species&Other"), labels=c("Richness (R)", "Evenness (E)", "Ranks (Ra)", "Species (S)", "R+E", "Ra+R &/or E", "S+R &/or E", "Ra+S+R &/or E"))+
  theme(axis.text.x = element_text(angle = 90, vjust=.5))+
  ylab("% C-T Differences")+
  xlab("RAC Difference Measure")+
  scale_y_continuous(limits = c(0, 50))



###redoing this with out the rare species
perm_outputnorare<-read.csv("C2E\\Products\\CommunityChange\\March2018 WG\\permanova_permdisp_output_norare5_Nov2021.csv")

mult_diffnorare <- read.csv("C2E\\Products\\CommunityChange\\March2018 WG\\CORRE_Mult_diff_Metrics_norare5_Nov2021.csv")%>%
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
                                                     ifelse(perm_Pvalue<0.0501&disp_Pvalue<0.0501&greater_disp == "C",6,999)))))))%>%
  na.omit()

num_scenmprare<- scenariosnorare%>%
  group_by(scenario)%>%
  summarize(n=length(scenario), prop=n/2641)


####what is the role of species differneces?
abund_diff <- read.csv("C2E\\Products\\CommunityChange\\March2018 WG\\CORRE_Abund_Diff_Nov2021.csv")%>%
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
