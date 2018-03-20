### Bootstrapping ###
require(tidyr)
require(dplyr)
require(vegan)
require(mvtnorm)
require(plyr)

setwd("C:\\Users\\WilcoxKR.WILCOXKR-LAPTOP\\Documents\\GitHub\\C2E\\")

relabun <- read.csv("CORRE_relative_abundance.csv")
head(relabun)
colnames(relabun)
pplots <- relabun%>%
  filter(project_name=="pplots"&(treatment=="N1P0"|treatment=="N2P0"))
pplots <- pplots[,c(-1,-8)]


# Calculating richness and evenness ---------------------------------------

pplots.wide <- spread(pplots, key=genus_species, value=relcov, fill=0)
pplots.rich <- specnumber(pplots.wide[,8:ncol(pplots.wide)])
pplots.invD <- diversity(pplots.wide[,8:ncol(pplots.wide)],"inv")
pplots.even <- pplots.invD/pplots.rich

control <- ddply(subset(pplots.2, treatment=="N1P0"),
                 .(calendar_year), summarize,
                 even.ctrl = mean(pplots.even),
                 S.ctrl = mean(pplots.rich,na.rm=TRUE))

pplots.2 <- cbind(pplots.wide[,1:7], cbind(pplots.rich, pplots.even))
pplots.3 <- merge(pplots.2, control, by=c("calendar_year"),all=TRUE)
pplots.3$even.lnRR <- log(pplots.3$pplots.even/pplots.3$even.ctrl)
pplots.3$S.lnRR <- log(pplots.3$pplots.rich/pplots.3$S.ctrl)

pplots.trt <- subset(pplots.3, treatment=="N2P0")
head(pplots.trt)

div.trt <- pplots.trt[,c("calendar_year","plot_id","even.lnRR","S.lnRR")]
head(div.trt)

# CAlculating BC distances ------------------------------------------------

years <- levels(as.factor(pplots.wide$calendar_year))
trt.dist <- list()
for(i in 1:length(years)){
  data.temp <- subset(pplots.wide, calendar_year==years[i])
  sp.temp <- data.temp[,8:ncol(data.temp)]
  bc.temp <- vegdist(sp.temp, method="bray")  
  bc.small <- as.matrix(bc.temp)[7:12,1:6]
  dist.temp <- rowMeans(bc.small)
  plotid.temp <- as.numeric(as.character(subset(data.temp, treatment=="N2P0")$plot_id))
  trt.dist.temp <- data.frame(calendar_year=years[i],plot_id=plotid.temp, dist=dist.temp)
  trt.dist <- rbind(trt.dist,trt.dist.temp)
      }
head(trt.dist)

# Combine all into data frame -------------------------------------------

all.trt <- merge(trt.dist,div.trt, by=c("calendar_year","plot_id"),all=TRUE)
head(all.trt)
head(pplots.trt)

# Generate means for each year --------------------------------------------

means <- ddply(all.trt, .(calendar_year), summarize,
                     mean.even = mean(even.lnRR),
                     mean.s = mean(S.lnRR),
                     mean.dist = mean(dist))


# Loop through years and combine covariance matrics for each year ---------

years <- levels(as.factor(all.trt$calendar_year))
cov.list <- list()
for(i in 1:length(years)){
  data.temp <- subset(all.trt, calendar_year==years[i])
  cov.temp <- cov(as.matrix(data.temp[,3:5]))
  cov.list[[years[i]]] <- cov.temp
}
cov.list[[1]]



# Generate variables using covariance matrix ------------------------------
i=1
years <- levels(as.factor(all.trt$calendar_year))
random.list <- list()
for(i in 1:length(years)){
  data.temp <- as.numeric(subset(means, calendar_year==years[i])[,-1])
  covariance.temp <- as.matrix(cov.list[[i]])
  random.temp <- rmvnorm(1000, data.temp, covariance.temp)
  random.list[[years[i]]] <- random.temp
}
draws <- sample(1:1000, size=100, replace=TRUE)  

i=12
varpart.df <- list()
for(i in 1:length(years)){
  data.temp <- random.list[[i]]
  short.temp <- as.data.frame(data.temp[draws,])
  colnames(short.temp) <- c("even","s","dist")
  varpart.temp <- varpart(short.temp$dist, ~s, ~even, data=short.temp)
  adjR.temp <- data.frame(year=years[i],
                          metric=c("even","shared","rich","resid"),
                          adj.r2=varpart.temp[["part"]][["indfract"]][1:4,"Adj.R.squared"])
  varpart.df <- rbind(varpart.df, adjR.temp)
}

require(ggplot2)  
ggplot()+
  geom_line(data=subset(varpart.df,metric%in%c("even","rich")), aes(x=year, y=adj.r2, colour=metric, group=metric))+
  geom_point(data=means,aes(x=calendar_year,y=mean.dist))
means


