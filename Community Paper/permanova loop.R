### calculating permanova for each experiment and year
###
### Authors: Kevin Wilcox (wilcoxkr@gmail.com) Andrew Tredennick (atredenn@gmail.com)
### Last updated: March 20 2018

corredat1<-read.csv("~/Dropbox/converge_diverge/datasets/LongForm/SpeciesRelativeAbundance_Oct2017.csv")%>%
  select(-X)%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))%>%
  filter(site_project_comm!="GVN_FACE_0", site_project_comm!="AZI_NitPhos_0", site_project_comm!="JRN_study278_0", site_project_comm!="KNZ_GFP_4F", site_project_comm!="Saskatchewan_CCD_0")

##several studies only have two measurments of a plot. I am dropping those plots
azi<-read.csv("~/Dropbox/converge_diverge/datasets/LongForm/SpeciesRelativeAbundance_Oct2017.csv")%>%
  select(-X)%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))%>%
  filter(site_code=="AZI")%>%
  filter(plot_id!=11&plot_id!=15&plot_id!=35&plot_id!=37)

jrn<-read.csv("~/Dropbox/converge_diverge/datasets/LongForm/SpeciesRelativeAbundance_Oct2017.csv")%>%
  select(-X)%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))%>%
  filter(site_project_comm=="JRN_study278_0")%>%
  filter(plot_id!=211&plot_id!=210)

knz<-read.csv("~/Dropbox/converge_diverge/datasets/LongForm/SpeciesRelativeAbundance_Oct2017.csv")%>%
  select(-X)%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))%>%
  filter(site_project_comm=="KNZ_GFP_4F")%>%
  filter(plot_id!="7_1_1"&plot_id!="7_2_1")

sak<-read.csv("~/Dropbox/converge_diverge/datasets/LongForm/SpeciesRelativeAbundance_Oct2017.csv")%>%
  select(-X)%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))%>%
  filter(site_project_comm=="Saskatchewan_CCD_0")%>%
  filter(plot_id!=2)

corredat<-rbind(corredat1, azi, jrn, knz, sak)

trt<-read.csv("~/Dropbox/converge_diverge/datasets/LongForm/ExperimentInformation_Nov2017.csv")%>%
  select(site_code, project_name, community_type, treatment,plot_mani)%>%
  unique()%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))


control<-ranks%>%
  filter(C_T=="Control")

treat_list<-unique(subset(ranks, C_T=="Treatment")$treatment)

for(i in 1:length(treat_list)) {
  subset_trt<-ranks%>%
    filter(treatment==treat_list[i])
  
  #dataset of two treatments    
  subset_t12<-rbind(control, subset_trt)
  
  result <- subset_t12 %>%
    group_by(time) %>%
    do({
      y <- unique(.$treatment)###assumption this is a length 2 list
      df1 <- filter(., treatment==y[[1]])
      df2 <- filter(., treatment==y[[2]])
      sf1 <- stepfun(df1$relrank, c(0, df1$cumabund))
      sf2 <- stepfun(df2$relrank, c(0, df2$cumabund))
      r <- sort(unique(c(0, df1$relrank, df2$relrank)))
      h <- abs(sf1(r) - sf2(r))
      w <- c(diff(r), 0)
      data.frame(CC=sum(w*h))#do has to output a dataframe
    })