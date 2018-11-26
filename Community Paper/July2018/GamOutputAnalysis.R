################################################################################
##  GamOutputAnalysis.R This code analyzes the spreadsheet output by Andrew do see when communities changed.
##
##  Author: Meghan Avolio (meghan.avolio@gmail.com)
##  Date: August 6, 2018
##  Updated: Nov 4, 2018
################################################################################


library(tidyverse)
library(gridExtra)
theme_set(theme_bw(20))

#work
setwd("C:\\Users\\megha\\Dropbox\\")
#home
setwd("~/Dropbox/")

#read in and drop sucessional treatments (those that had a big disturbance to start)
gam<-read.csv("C2E/Products/CommunityChange/Summer2018_Results/gam_comparison_table.csv")
gam_exp<-gam%>%
  select(site_proj_comm)%>%
  unique()

#filter the press treatments for the experiments with enough data
trt_touse<-read.csv("C2E/Products/CommunityChange/March2018 WG/treatment interactions_July2018.csv")%>%
  select(site_proj_comm, treatment)%>%
  unique()%>%
  right_join(gam_exp)

gam_touse<-gam%>%
 right_join(trt_touse)

gam_trt<-gam%>%
  select(site_proj_comm, treatment)%>%
  unique()

trts_interactions<-read.csv("C2E/Products/CommunityChange/March2018 WG/treatment interactions_July2018.csv")%>%
  right_join(gam_trt)%>%
  filter(use==1)

#summary table of treatments
trt_summary<-read.csv("C2E/Products/CommunityChange/March2018 WG/treatment interactions_July2018.csv")%>%
  select(site_proj_comm, treatment, trt_type, use, trt_type2)%>%
  unique()%>%
  right_join(gam_exp)

# write.csv(trt_summary, "treatment_summary.csv")

###

#no longer doing this based on those communities that changed.

# #first, how many communities saw sig differences in change between controls and treatments?
# compchange<-gam_touse%>%
#   filter(response_var == "composition_change")%>%
#   group_by(sig_diff_cntrl_trt)%>%
#   summarize(length(sig_diff_cntrl_trt))
# 
# ##did change differ depending on treatment?
# gamtrts<-gam%>%
#   left_join(trts_interactions)%>%
#   filter(response_var == "composition_change", use == 1)%>%
#   group_by(trt_type2) %>%
#   summarise(
#     num_sig = length(which(sig_diff_cntrl_trt == "yes")),
#     num_nonsig = length(which(sig_diff_cntrl_trt == "no"))
#   )%>%
#   filter(trt_type2!="N+irr"&trt_type2!="drought")
# 
# write.csv(gamtrts, "C2E/Products/CommunityChange/Summer2018_Results/chisq_comp_change_newtrts.csv")
# 
# gamtrts_toplot<-gamtrts%>%
#   mutate(sum = num_sig+num_nonsig,
#          psig = num_sig/sum,
#          pnsig = num_nonsig/sum)%>%
#   select(-sum, -num_sig, -num_nonsig)%>%
#   gather(key = sig, value = number, -trt_type2)
# 
# 
# ggplot(gamtrts_toplot, aes(x = trt_type2, y = number, fill = sig)) +
#   geom_col(width = 0.7) +
#   coord_flip() +
#   theme_minimal() +
#   scale_fill_brewer(name = "Treatment vs. Control", labels = c("Not significant", "Significant")) +
#   labs(x = "Global Change Treatment", y = "Proportion of communities") +
#   geom_hline(yintercept = 0.50)+
#   theme(legend.position = "top")+
#   annotate(geom="text", x = 1, y = 0.05, label="n = 7", size=3)+
#   annotate(geom="text", x = 2, y = 0.05, label="n = 12", size=3)+
#   annotate(geom="text", x = 3, y = 0.05, label="n = 51", size=3)+
#   annotate(geom="text", x = 4, y = 0.05, label="n = 34", size=3)+
#   annotate(geom="text", x = 5, y = 0.05, label="n = 11", size=3)
# 
# 
# #second, of those that saw change, what aspect of the community changes
# comchange_sig<-gam%>%
#   filter(response_var == "composition_change" & sig_diff_cntrl_trt == "yes")%>%
#   select(site_proj_comm, treatment)%>%
#   mutate(keep = "yes")
# # 
# # write.csv(comchange_sig, 'C:\\Users\\megha\\Dropbox\\C2E\\Products\\CommunityChange\\Summer2018_Results\\gam_com_sig_change.csv', row.names = F )
# 
# 
# #second, of those that saw change, what aspect of the community changes
# comchange_sig_metrics<-gam%>%
#   filter(response_var != "composition_change" & sig_diff_cntrl_trt == "yes")%>%
#   select(site_proj_comm, treatment, response_var)%>%
#   mutate(keep = "yes")
# 
# overlap<-comchange_sig_metrics%>%
#   select(-keep)%>%
#   left_join(comchange_sig)
# 
# write.csv(comchange_sig_metrics, 'C:\\Users\\megha\\Dropbox\\C2E\\Products\\CommunityChange\\Summer2018_Results\\gam_metrics_sig_change.csv', row.names = F )
# 
# sig_only<-gam%>%
#   right_join(comchange_sig)%>%
#   filter(response_var != "composition_change")
# 
# sig_tally <- sig_only %>%
#   group_by(response_var) %>%
#   summarise(
#     num_sig = length(which(sig_diff_cntrl_trt == "yes")),
#     num_nonsig = length(which(sig_diff_cntrl_trt == "no"))
#   ) %>%
#   mutate(sum = num_sig+num_nonsig,
#          psig = num_sig/sum,
#          pnsig = num_nonsig/sum)%>%
#   select(-sum, -num_sig, -num_nonsig)%>%
#   gather(key = sig, value = value, -response_var) %>%
#   mutate(
#     response_var2 = ifelse(response_var == "evenness_change_abs", "Evenness", response_var),
#     response_var2 = ifelse(response_var == "gains", "Species gains", response_var2),
#     response_var2 = ifelse(response_var == "losses", "Species losses", response_var2),
#     response_var2 = ifelse(response_var == "rank_change", "Rank change", response_var2),
#     response_var2 = ifelse(response_var == "richness_change_abs", "Richness change", response_var2))
# 
# allSERGL<-sig_tally%>%
#   select(-response_var2)%>%
#   spread(sig, value)
# 
# write.csv(allSERGL, 'C:\\Users\\megha\\Dropbox\\C2E\\Products\\CommunityChange\\Summer2018_Results\\gam_com_sig_change_all.csv', row.names = F )
# 
# sergl<-ggplot(sig_tally, aes(x = response_var2, y = value, fill = sig)) +
#   geom_col(width = 0.7) +
#   coord_flip() +
#   theme_minimal() +
#   scale_fill_brewer(name = "Treatment vs. Control", labels = c("Not significant", "Significant")) +
#   labs(x = "Change metric", y = "Proportion of communities") +
#   theme(legend.position = "top")+
#   geom_hline(yintercept = 0.5)
# 
# ##did aspect of community change depend on the treatment?
# gamtrts_metrics<-sig_only%>%
#   left_join(trts_interactions)%>%
#   ungroup()%>%
#   filter(response_var != "composition_change", use == 1)%>%
#   group_by(response_var, trt_type2) %>%
#   summarise(
#     num_sig = length(which(sig_diff_cntrl_trt == "yes")),
#     num_nonsig = length(which(sig_diff_cntrl_trt == "no"))
#   )%>%
#   filter(trt_type2!="N+irr"&trt_type2!="drought")
# 
# write.csv(gamtrts_metrics, "C:\\Users\\megha\\Dropbox\\C2E\\Products\\CommunityChange\\Summer2018_Results\\chisq_metrics_change.csv")
# 
# tograph_metrics_trt<-gamtrts_metrics%>%
#   mutate(sum = num_sig + num_nonsig,
#          psig = num_sig/sum,
#          pnonsig = num_nonsig/sum)%>%
#   select(-num_sig, -num_nonsig, -sum)%>%
#   gather(key = sig, value = value, -trt_type2, -response_var) %>%
#   mutate(
#     response_var2 = ifelse(response_var == "evenness_change_abs", "Evenness", ifelse(response_var == "gains", "Species gains", ifelse(response_var == "losses", "Species losses", ifelse(response_var == "rank_change", "Rank change", "Richness change")))))
# 
# 
# trt_sergl<-ggplot(subset(tograph_metrics_trt, trt_type2!="P"), aes(x = response_var2, y = value, fill = sig)) +
#   geom_col(width = 0.7) +
#   coord_flip() +
#   theme_minimal() +
#   scale_fill_brewer(name = "Treatment vs. Control", labels = c("Not significant", "Significant")) +
#   labs(x = "Change metric", y = "Proportion of communities") +
#   theme(legend.position = "top")+
#   geom_hline(yintercept=0.50)+
#   facet_wrap(~trt_type2, ncol = 2)


#Regradless of whether the community changed, how much to the metrics change?
metrics_sig<-gam_touse%>%
  filter(response_var != "composition_change")

metrics_sig_tally <- metrics_sig %>%
  group_by(response_var) %>%
  summarise(
    num_sig = length(which(sig_diff_cntrl_trt == "yes")),
    num_nonsig = length(which(sig_diff_cntrl_trt == "no"))
  )


#overall diff in metrics of change - NO
prop.test(x=as.matrix(metrics_sig_tally[c('num_sig', 'num_nonsig')]), alternative='two.sided')

sig_tograph<-metrics_sig_tally%>%
  mutate(sum = num_sig+num_nonsig,
         psig = num_sig/sum,
         pnsig = num_nonsig/sum)%>%
  select(-sum, -num_sig, -num_nonsig)%>%
  gather(key = sig, value = value, -response_var) %>%
  mutate(
    response_var2 = ifelse(response_var == "evenness_change_abs", "Evenness", response_var),
    response_var2 = ifelse(response_var == "gains", "Species gains", response_var2),
    response_var2 = ifelse(response_var == "losses", "Species losses", response_var2),
    response_var2 = ifelse(response_var == "rank_change", "Rank change", response_var2),
    response_var2 = ifelse(response_var == "richness_change_abs", "Richness change", response_var2))


sergl<-ggplot(sig_tograph, aes(x = response_var2, y = value, fill = sig)) +
  geom_col(width = 0.7) +
  coord_flip() +
  theme_minimal() +
  scale_fill_brewer(name = "Treatment vs. Control", labels = c("Not significant", "Significant")) +
  scale_x_discrete(limits=c("Richness change","Evenness","Rank change","Species gains","Species losses"), labels = c("Richness change","Evenness Change","Rank Change","Species Gains","Species Losses"))+
  labs(x = "Change metric", y = "Proportion of communities") +
  theme(legend.position = "top")+
  geom_hline(yintercept = 0.5)

###how does this differ by GCD?

gamtrts_metrics<-metrics_sig%>%
  left_join(trts_interactions)%>%
  ungroup()%>%
  filter(response_var != "composition_change", use == 1)%>%
  group_by(response_var, trt_type2) %>%
  summarise(
    num_sig = length(which(sig_diff_cntrl_trt == "yes")),
    num_nonsig = length(which(sig_diff_cntrl_trt == "no"))
  )

CO2 <- gamtrts_metrics%>%filter(trt_type2=='CO2')
prop.test(x=as.matrix(CO2[c('num_sig', 'num_nonsig')]), alternative='two.sided')

#N
N <- gamtrts_metrics%>%filter(trt_type2=='N')
prop.test(x=as.matrix(N[c('num_sig', 'num_nonsig')]), alternative='two.sided')

#P
P <- gamtrts_metrics%>%filter(trt_type2=='P')
prop.test(x=as.matrix(P[c('num_sig', 'num_nonsig')]), alternative='two.sided')

#irr
irr <- gamtrts_metrics%>%filter(trt_type2=='Irrigation')
prop.test(x=as.matrix(irr[c('num_sig', 'num_nonsig')]), alternative='two.sided')

#irr+temp
irrT <- gamtrts_metrics%>%filter(trt_type2=='Irr + Temp')
prop.test(x=as.matrix(irrT[c('num_sig', 'num_nonsig')]), alternative='two.sided')

#mult nuts
multnuts <- gamtrts_metrics%>%filter(trt_type2=='Mult. Nuts.')
prop.test(x=as.matrix(multnuts[c('num_sig', 'num_nonsig')]), alternative='two.sided')

#mult nuts
temp <- gamtrts_metrics%>%filter(trt_type2=='Temperature')
prop.test(x=as.matrix(temp[c('num_sig', 'num_nonsig')]), alternative='two.sided')



tograph_metrics_trt<-gamtrts_metrics%>%
  mutate(sum = num_sig + num_nonsig,
         psig = num_sig/sum,
         pnonsig = num_nonsig/sum)%>%
  select(-num_sig, -num_nonsig, -sum)%>%
  gather(key = sig, value = value, -trt_type2, -response_var) %>%
  mutate(
    response_var2 = ifelse(response_var == "evenness_change_abs", "Evenness", ifelse(response_var == "gains", "Species gains", ifelse(response_var == "losses", "Species losses", ifelse(response_var == "rank_change", "Rank change", "Richness change")))))


trt_sergl<-ggplot(subset(tograph_metrics_trt, trt_type2!="P"), aes(x = response_var2, y = value, fill = sig)) +
  geom_col(width = 0.7) +
  coord_flip() +
  theme_minimal() +
  scale_fill_brewer(name = "Treatment vs. Control", labels = c("Not significant", "Significant")) +
  labs(x = "Change metric", y = "Proportion of communities") +
  theme(legend.position = "top")+
  geom_hline(yintercept=0.50)+
  facet_wrap(~trt_type2, ncol = 2)
