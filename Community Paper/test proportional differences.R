################################################################################
##  test_proportional_differences.R: Chi-squared test for differences among responses by treatment for compositional responses and its component parts.
##
##  Author: Kimberly La Pierre
##  Date created: August 7, 2018
################################################################################


library(tidyverse)


#kim
setwd('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\C2E\\Products\\CommunityChange\\Summer2018_Results')

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

#evenness
evenness <- SERGL%>%filter(response_var=='evenness_change_abs')
prop.test(x=as.matrix(evenness[,c('num_sig', 'num_nonsig')]), alternative='two.sided')

#gains
gains <- SERGL%>%filter(response_var=='gains')
prop.test(x=as.matrix(gains[,c('num_sig', 'num_nonsig')]), alternative='two.sided')

#losses
losses <- SERGL%>%filter(response_var=='losses')
prop.test(x=as.matrix(losses[,c('num_sig', 'num_nonsig')]), alternative='two.sided')
