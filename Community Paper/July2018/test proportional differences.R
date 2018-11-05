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
setwd('~/Dropbox/C2E/Products/CommunityChange/Summer2018_Results')

#comp change data - DO trts differ for all the data no only those that saw change
compChange <- read.csv('chisq_comp_change_newtrts.csv')

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

###SERGL change data  #only for sig changes
SERGL_all <- read.csv('gam_com_sig_change_all.csv')

#across GCD treatments
prop.test(x=as.matrix(SERGL_all[c('pnsig', 'psig')]), alternative='two.sided')

#no longer doing this.
# #evenness
# evenness <- SERGL_all%>%filter(response_var=='evenness_change_abs')
# prop.test(x=as.matrix(evenness[c('num_sig', 'num_nonsig')]), alternative='two.sided')
# 
# #gains
# gains <- SERGL%>%filter(response_var=='gains')
# prop.test(x=as.matrix(gains[,c('num_sig', 'num_nonsig')]), alternative='two.sided')
# 
# #losses
# losses <- SERGL%>%filter(response_var=='losses')
# prop.test(x=as.matrix(losses[,c('num_sig', 'num_nonsig')]), alternative='two.sided')
# 
# #rank_change
# rank_change <- SERGL%>%filter(response_var=='rank_change')
# prop.test(x=as.matrix(rank_change[,c('num_sig', 'num_nonsig')]), alternative='two.sided')
# 
# #richness_change_abs
# richness_change_abs <- SERGL%>%filter(response_var=='richness_change_abs')
# prop.test(x=as.matrix(richness_change_abs[,c('num_sig', 'num_nonsig')]), alternative='two.sided')
# 
###SERGL change data  #only for sig changes for select treatments
SERGL <- read.csv('chisq_metrics_change.csv')

###looking at difference by GCD
#CO2
CO2 <- SERGL%>%filter(trt_type2=='CO2')
prop.test(x=as.matrix(CO2[c('num_sig', 'num_nonsig')]), alternative='two.sided')

#N
N <- SERGL%>%filter(trt_type2=='N')
prop.test(x=as.matrix(N[c('num_sig', 'num_nonsig')]), alternative='two.sided')

#irr
irr <- SERGL%>%filter(trt_type2=='Irr')
prop.test(x=as.matrix(irr[c('num_sig', 'num_nonsig')]), alternative='two.sided')

#mult nuts
multnuts <- SERGL%>%filter(trt_type2=='Mult. Nuts.')
prop.test(x=as.matrix(multnuts[c('num_sig', 'num_nonsig')]), alternative='two.sided')
