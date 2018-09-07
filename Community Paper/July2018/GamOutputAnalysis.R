################################################################################
##  GamOutputAnalysis.R This code analyzes the spreadsheet output by Andrew do see when communities changed.
##
##  Author: Meghan Avolio (meghan.avolio@gmail.com)
##  Date: August 6, 2018
################################################################################


library(tidyverse)
theme_set(theme_bw(20))

#works
gam<-read.csv("C:\\Users\\megha\\Dropbox\\C2E\\Products\\CommunityChange\\Summer2018_Results\\gam_comparison_table.csv")
trts_interactions<-read.csv("C:\\Users\\megha\\Dropbox\\C2E\\Products\\CommunityChange\\March2018 WG\\treatment interactions_July2018.csv")%>%
  select(-site_project_comm)%>%
  mutate(site_proj_comm = paste(site_code, project_name, community_type, sep = "_"))

#mac
gam<-read.csv("~/Dropbox/C2E/Products/CommunityChange/Summer2018_Results/gam_comparison_table.csv")

trts_interactions<-read.csv("~/Dropbox/C2E/Products/CommunityChange/March2018 WG/treatment interactions_July2018.csv")%>%
  select(-site_project_comm)%>%
  mutate(site_proj_comm = paste(site_code, project_name, community_type, sep = "_"))

#first, how many communities saw sig differences in change between controls and treatments?
compchange<-gam%>%
  filter(response_var == "composition_change")%>%
  group_by(sig_diff_cntrl_trt)%>%
  summarize(length(sig_diff_cntrl_trt))

##did change differ depending on treatment?
gamtrts<-gam%>%
  left_join(trts_interactions)%>%
  filter(response_var == "composition_change", use == 1)%>%
  group_by(trt_type) %>%
  summarise(
    num_sig = length(which(sig_diff_cntrl_trt == "yes")),
    num_nonsig = length(which(sig_diff_cntrl_trt == "no"))
  )%>%
  filter(trt_type!="N+irr"&trt_type!="drought")

write.csv(gamtrts, "C:\\Users\\megha\\Dropbox\\C2E\\Products\\CommunityChange\\Summer2018_Results\\chisq_comp_change.csv")

gamtrts_toplot<-gamtrts%>%
  mutate(sum = num_sig+num_nonsig,
         psig = num_sig/sum,
         pnsig = num_nonsig/sum)%>%
  select(-sum, -num_sig, -num_nonsig)%>%
  gather(key = sig, value = number, -trt_type)


  ggplot(gamtrts_toplot, aes(x = trt_type, y = number, fill = sig)) +
  geom_col(width = 0.7) +
  coord_flip() +
  theme_minimal() +
  scale_fill_brewer(type = "qual", name = "Treatment vs. Control", labels = c("Not significant", "Significant")) +
  labs(x = "Global Change Treatment", y = "Proportion of communities") +
  geom_hline(yintercept = 0.50)+
  theme(legend.position = "top")


#second, of those that saw change, what aspect of the community changes
comchange_sig<-gam%>%
  filter(response_var == "composition_change" & sig_diff_cntrl_trt == "yes")%>%
  select(site_proj_comm, treatment)%>%
  mutate(keep = "yes")

write.csv(comchange_sig, 'C:\\Users\\megha\\Dropbox\\C2E\\Products\\CommunityChange\\Summer2018_Results\\gam_com_sig_change.csv', row.names = F )

sig_only<-gam%>%
  right_join(comchange_sig)%>%
  filter(response_var != "composition_change")

sig_tally <- sig_only %>%
  group_by(response_var) %>%
  summarise(
    num_sig = length(which(sig_diff_cntrl_trt == "yes")),
    num_nonsig = length(which(sig_diff_cntrl_trt == "no"))
  ) %>%
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


ggplot(sig_tally, aes(x = response_var2, y = value, fill = sig)) +
  geom_col(width = 0.7) +
  coord_flip() +
  theme_minimal() +
  scale_fill_brewer(type = "qual", name = "Treatment vs. Control", labels = c("Not significant", "Significant")) +
  labs(x = "Change metric", y = "Proportion of communities") +
  theme(legend.position = "top")+
  geom_hline(yintercept = 0.5)

##did aspect of community change depend on the treatment?
gamtrts_metrics<-sig_only%>%
  left_join(trts_interactions)%>%
  ungroup()%>%
  filter(response_var != "composition_change", use == 1)%>%
  group_by(response_var, trt_type) %>%
  summarise(
    num_sig = length(which(sig_diff_cntrl_trt == "yes")),
    num_nonsig = length(which(sig_diff_cntrl_trt == "no"))
  )%>%
  filter(trt_type!="N+irr"&trt_type!="drought"&trt_type!="N+P+K")

write.csv(gamtrts_metrics, "C:\\Users\\megha\\Dropbox\\C2E\\Products\\CommunityChange\\Summer2018_Results\\chisq_metrics_change.csv")

tograph_metrics_trt<-gamtrts_metrics%>%
  mutate(sum = num_sig + num_nonsig,
         psig = num_sig/sum,
         pnonsig = num_nonsig/sum)%>%
  select(-num_sig, -num_nonsig, -sum)%>%
  gather(key = sig, value = value, -trt_type, -response_var) %>%
  mutate(
    response_var2 = ifelse(response_var == "evenness_change_abs", "Evenness", ifelse(response_var == "gains", "Species gains", ifelse(response_var == "losses", "Species losses", ifelse(response_var == "rank_change", "Rank change", "Richness change")))))


ggplot(subset(tograph_metrics_trt, trt_type!="P"), aes(x = response_var2, y = value, fill = sig)) +
  geom_col(width = 0.7) +
  coord_flip() +
  theme_minimal() +
  scale_fill_brewer(type = "qual", name = "Treatment vs. Control", labels = c("Not significant", "Significant")) +
  labs(x = "Change metric", y = "Proportion of communities") +
  theme(legend.position = "top")+
  geom_hline(yintercept=0.50)+
  facet_wrap(~trt_type, ncol = 2)
