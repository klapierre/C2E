library(tidyverse)
library(randomForest)
# library(MultivariateRandomForest)

source('C:\\Users\\lapie\\Desktop\\R files laptop\\C2E\\Community Paper\\calculate glass delta.R')

setwd('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\C2E\\Products\\CommunityChange\\March2018 WG')

theme_set(theme_bw())
theme_update(axis.title.x=element_text(size=20, vjust=-0.35), axis.text.x=element_text(size=16),
             axis.title.y=element_text(size=20, angle=90, vjust=0.5), axis.text.y=element_text(size=16),
             plot.title = element_text(size=24, vjust=2),
             panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
             legend.title=element_blank(), legend.text=element_text(size=20))

###read in data
#treatments
trt <- read.csv('ExperimentInformation_Nov2017.csv')%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep='_'))%>%
  select(-X, -plot_mani)

#site information
site <- read.csv('SiteExperimentDetails_Dec2016.csv')%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep='_'))%>%
  select(-X)

#spp abundances (just to link plot_id to treatments)
abund <- read.csv('SpeciesRelativeAbundance_Oct2017.csv')%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep='_'))%>%
  select(site_project_comm, treatment, plot_id, treatment_year)%>%
  unique()

#community change metrics
metrics <- read.csv('CORRE_RAC_Metrics_March2018_trtyr.csv')%>%
  select(-X)%>%
  left_join(abund)%>%
  left_join(trt)%>%
  left_join(site)%>%
  select(site_project_comm, site_code, project_name, community_type, treatment, treatment_year, plot_id, richness_change, evenness_change, rank_change, gains, losses, MAP, MAT, rrich, anpp, n, p, k, CO2, precip, temp)

#glass' delta
glass_metrics <- change_glass_d%>%
  select(-plot_mani)%>%
  mutate(treatment=treatment.x)%>%
  left_join(trt)%>%
  left_join(site)%>%
  select(site_project_comm, site_code, project_name, community_type, treatment, treatment_year, abs_richness_glass, abs_evenness_glass, rank_glass, gains_glass, losses_glass, MAP, MAT, rrich, anpp, n, p, k, CO2, precip, temp)%>%
  group_by(site_project_comm, site_code, project_name, community_type, treatment, MAP, MAT, rrich, anpp, n, p, k, CO2, precip, temp)%>%
  summarise(abs_richness_mean=mean(abs_richness_glass), abs_evenness_mean=mean(abs_evenness_glass), rank_mean=mean(rank_glass), gains_mean=mean(gains_glass), losses_mean=mean(losses_glass))%>%
  ungroup()



###random forest
set.seed(71)
#richness
glass_richness <- glass_metrics%>%
  select(MAP, MAT, rrich, anpp, n, p, k, CO2, precip, temp, abs_richness_mean)
glass_richness <- glass_richness[complete.cases(glass_richness), ]
richnessRandomForest <- randomForest(abs_richness_mean ~ ., data=glass_richness, importance=TRUE,
                        proximity=TRUE)
print(richnessRandomForest)
richnessImportance <- as.data.frame(round(importance(richnessRandomForest), 2))%>%
  rownames_to_column()%>%
  mutate(metric='abs_richness')
#evenness
glass_evenness <- glass_metrics%>%
  select(MAP, MAT, rrich, anpp, n, p, k, CO2, precip, temp, abs_evenness_mean)
glass_evenness <- glass_evenness[complete.cases(glass_evenness), ]
evennessRandomForest <- randomForest(abs_evenness_mean ~ ., data=glass_evenness, importance=TRUE,
                                     proximity=TRUE)
print(evennessRandomForest)
evennessImportance <- as.data.frame(round(importance(evennessRandomForest), 2))%>%
  rownames_to_column()%>%
  mutate(metric='evenness')
#rank
glass_rank <- glass_metrics%>%
  select(MAP, MAT, rrich, anpp, n, p, k, CO2, precip, temp, rank_mean)
glass_rank <- glass_rank[complete.cases(glass_rank), ]
rankRandomForest <- randomForest(rank_mean ~ ., data=glass_rank, importance=TRUE,
                                     proximity=TRUE)
print(rankRandomForest)
rankImportance <- as.data.frame(round(importance(rankRandomForest), 2))%>%
  rownames_to_column()%>%
  mutate(metric='rank')
#gains
glass_gains <- glass_metrics%>%
  select(MAP, MAT, rrich, anpp, n, p, k, CO2, precip, temp, gains_mean)%>%
  filter(!is.infinite(gains_mean))
glass_gains <- glass_gains[complete.cases(glass_gains), ]
gainsRandomForest <- randomForest(gains_mean ~ ., data=glass_gains, importance=TRUE,
                                     proximity=TRUE)
print(gainsRandomForest)
gainsImportance <- as.data.frame(round(importance(gainsRandomForest), 2))%>%
  rownames_to_column()%>%
  mutate(metric='gains')
#losses
glass_loss <- glass_metrics%>%
  select(MAP, MAT, rrich, anpp, n, p, k, CO2, precip, temp, losses_mean)%>%
  filter(!is.infinite(losses_mean))
glass_loss <- glass_loss[complete.cases(glass_loss), ]
lossRandomForest <- randomForest(losses_mean ~ ., data=glass_loss, importance=TRUE,
                                     proximity=TRUE)
print(lossRandomForest)
lossImportance <- as.data.frame(round(importance(lossRandomForest), 2))%>%
  rownames_to_column()%>%
  mutate(metric='loss')
#all together
allImportance <- rbind(richnessImportance, evennessImportance, rankImportance, gainsImportance, lossImportance)%>%
  mutate(factor=rowname)%>%
  select(-rowname)
names(allImportance)[names(allImportance) == '%IncMSE'] <- 'importance'

#figure
ggplot(data=allImportance, aes(x=factor, y=importance)) +
  geom_bar(stat='identity') +
  scale_x_discrete(limits=c("MAP","MAT","anpp", "rrich", "n", "p", "k", "CO2", "precip", "temp")) +
  facet_wrap(~metric)


#what drives compositional change?
comp <- read.csv('CORRE_Mult_Metrics_March2018.csv')%>%
  select(-X)%>%
  group_by(site_project_comm, treatment)%>%
  summarise(comp_mean=mean(composition_change))%>%
  ungroup()

glassComp <- glass_metrics%>%
  left_join(comp)%>%
  filter(!is.infinite(gains_mean))%>%
  filter(!is.infinite(losses_mean))%>%
  select(abs_richness_mean, abs_evenness_mean, rank_mean, gains_mean, losses_mean, MAT, MAP, rrich, anpp, n, p, k, CO2, precip, temp, comp_mean)
  
glassComp <- glassComp[complete.cases(glassComp), ]

#random forest
compRandomForest <- randomForest(comp_mean ~ ., data=glassComp, importance=TRUE,
                                 proximity=TRUE)
print(compRandomForest)
compImportance <- as.data.frame(round(importance(compRandomForest), 2))%>%
  rownames_to_column()%>%
  mutate(factor=rowname)%>%
  select(-rowname)
names(compImportance)[names(compImportance) == '%IncMSE'] <- 'importance'


#figure
ggplot(data=compImportance, aes(x=factor, y=importance)) +
  geom_bar(stat='identity') +
  scale_x_discrete(limits=c("MAP","MAT","anpp", "rrich", "n", "p", "k", "CO2", "precip", "temp", "abs_richness_mean", "abs_evenness_mean", "rank_mean", "gains_mean", "losses_mean"))


#without the site-level predictors
glassComp2 <- glass_metrics%>%
  left_join(comp)%>%
  filter(!is.infinite(gains_mean))%>%
  filter(!is.infinite(losses_mean))%>%
  select(abs_richness_mean, abs_evenness_mean, rank_mean, gains_mean, losses_mean, comp_mean)

glassComp2 <- glassComp[complete.cases(glassComp), ]
compRandomForest2 <- randomForest(comp_mean ~ ., data=glassComp2, importance=TRUE,
                                 proximity=TRUE)
print(compRandomForest2)
compImportance2 <- as.data.frame(round(importance(compRandomForest2), 2))%>%
  rownames_to_column()%>%
  mutate(factor=rowname)%>%
  select(-rowname)
names(compImportance2)[names(compImportance2) == '%IncMSE'] <- 'importance'

#figure
ggplot(data=compImportance2, aes(x=factor, y=importance)) +
  geom_bar(stat='identity') +
  scale_x_discrete(limits=c("abs_richness_mean", "abs_evenness_mean", "rank_mean", "gains_mean", "losses_mean"))


# #### example for randomForest
# ##data(iris)
# set.seed(71)
# iris.rf <- randomForest(Species ~ ., data=iris, importance=TRUE,
#                         proximity=TRUE)
# print(iris.rf)
# ## Look at variable importance:
# round(importance(iris.rf), 2)
# ## Do MDS on 1 - proximity:
# iris.mds <- cmdscale(1 - iris.rf$proximity, eig=TRUE)
# op <- par(pty="s")
# pairs(cbind(iris[,1:4], iris.mds$points), cex=0.6, gap=0,
#       col=c("red", "green", "blue")[as.numeric(iris$Species)],
#       main="Iris Data: Predictors and MDS of Proximity Based on RandomForest")
# par(op)
# print(iris.mds$GOF)
