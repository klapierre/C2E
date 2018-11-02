################################################################################
##  test_proportional_differences.R: Chi-squared test for differences among responses by treatment for compositional responses and its component parts.
##
##  Author: Kimberly La Pierre
##  Date created: August 7, 2018
################################################################################


library(tidyverse)


#kim
setwd('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\C2E\\Products\\CommunityChange\\Summer2018_Results')

#meghan
setwd('C:\\Users\\megha\\Dropbox\\C2E\\Products\\CommunityChange\\Summer2018_Results')

#comp change data
compChange <- read.csv('chisq_comp_change.csv')

###treatment type - compositional responses
prop.test(x=as.matrix(compChange[,c('num_sig', 'num_nonsig')]), alternative='two.sided')

FUN = function(i,j){     
  chisq.test(matrix(c(compChange[i,2], compChange[i,3],
                      compChange[j,2], compChange[j,3]),
                    nrow=2,
                    byrow=TRUE))$ p.value
}

pairwise.table(FUN,
               rownames(compChange),
               p.adjust.method="none")

###SERGL change data
SERGL <- read.csv('chisq_metrics_change.csv')

#across GCD treatments
SERGLall <- SERGL%>%
  select(-X)%>%
  group_by(response_var)%>%
  summarise(num_sig2=sum(num_sig), num_nonsig2=sum(num_nonsig))
prop.test(x=as.matrix(SERGLall[c('num_sig2', 'num_nonsig2')]), alternative='two.sided')

#evenness
evenness <- SERGL%>%filter(response_var=='evenness_change_abs')
prop.test(x=as.matrix(evenness[c('num_sig', 'num_nonsig')]), alternative='two.sided')

#gains
gains <- SERGL%>%filter(response_var=='gains')
prop.test(x=as.matrix(gains[,c('num_sig', 'num_nonsig')]), alternative='two.sided')

#losses
losses <- SERGL%>%filter(response_var=='losses')
prop.test(x=as.matrix(losses[,c('num_sig', 'num_nonsig')]), alternative='two.sided')

#rank_change
rank_change <- SERGL%>%filter(response_var=='rank_change')
prop.test(x=as.matrix(rank_change[,c('num_sig', 'num_nonsig')]), alternative='two.sided')

#richness_change_abs
richness_change_abs <- SERGL%>%filter(response_var=='richness_change_abs')
prop.test(x=as.matrix(richness_change_abs[,c('num_sig', 'num_nonsig')]), alternative='two.sided')


###looking at difference by GCD
#CO2
CO2 <- SERGL%>%filter(trt_type=='CO2')
prop.test(x=as.matrix(CO2[c('num_sig', 'num_nonsig')]), alternative='two.sided')

#N
N <- SERGL%>%filter(trt_type=='N')
prop.test(x=as.matrix(N[c('num_sig', 'num_nonsig')]), alternative='two.sided')

#irr
irr <- SERGL%>%filter(trt_type=='irr')
prop.test(x=as.matrix(irr[c('num_sig', 'num_nonsig')]), alternative='two.sided')

#N+P
NP <- SERGL%>%filter(trt_type=='N+P')
prop.test(x=as.matrix(NP[c('num_sig', 'num_nonsig')]), alternative='two.sided')
