library(tidyverse)
library(gridExtra)
library(reldist)
library(grid)
library(gtable)
library(codyn)
library(vegan)
library(Kendall)


#read in data and subset only experiemnts that are 10 or more years

dat<-read.csv("C:\\Users\\megha\\Dropbox\\converge_diverge\\datasets\\LongForm\\CORRE_RAC_Metrics_Oct2017_compareyears.csv")%>%
  select(-X)

longset<-dat%>%
  select(site_project_comm, treatment_year)%>%
  group_by(site_project_comm)%>%
  summarize(len=max(treatment_year))%>%
  filter(len>9)%>%
  select(-len)

datsub<-merge(dat, longset, by="site_project_comm")

##LOOP THIS
#subset to get a single treatment for a site_proj_comm
test_c<-datsub%>%
  filter(site_project_comm=="KNZ_pplots_0"&treatment=="N1P0")

#do vp and pull out what we want
vp<-varpart(test_np$mean_change, ~even, ~gain, ~loss, ~MRSc, data=test_np)

standdev<-sd(test_np$mean_change)

adjR.temp <- data.frame(metric=c("even","gain","loss","MRSc","resid"),
                        adj.r2=vp[["part"]][["indfract"]][c(1:4,16),"Adj.R.square"],
                        sd_meanchange=standdev)

# looking at this
 plot(vp)
# plot(test_c$calendar_year,test_c$mean_change)
# points(test_np$calendar_year, test_np$mean_change, col="red", pch=19)
