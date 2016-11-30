library(tidyr)
library(dplyr)
library(codyn)
library(vegan)
library(Kendall)


#' Generate a random integer partition through the Chinese
#' Restaurant Process (CRP).

#STEP1 - Run the function to generate communities

rCRP = function(n, theta, alpha = 0, kappa = NULL, m = NULL, zeros = TRUE) {
  if (!is.null(kappa) & !is.null(m)) {
    if (m == round(m) & m > 0 & kappa > 0) {
      alpha <- -kappa
      theta <- kappa * m
    } else {
      stop("Parameter m must be a positive integer and kappa must be non-negative.")
    }
  } else {
    if (alpha < 0 | 1 <= alpha | -alpha < theta) { 
      stop("Without kappa or m, parameters must satisfy 0 <= alpha < 1 & theta > -alpha.")
    }
  }
  if (!is.null(kappa) & !is.null(m) & zeros) {
    # Sample so as to preserve zero abundances, given integer m
    result <- rep(0, m)
    extant <- which(result != 0)
    j <- m - length(extant)
    p <- rep(NA, m)
    idx <- sample(m, 1)
    for (k in 1:n) {
      result[[idx]] <- result[[idx]] + 1
      if (result[[idx]] == 1) {
        # respond to new class
        extant <- c(extant, idx)
        j <- j + 1
        p[-extant] <- (theta + j * alpha) / (m - j)
        p[[idx]] <- 1 - alpha
      } else {
        # step up the sampling prob for extant class
        p[[idx]] <- p[[idx]] + 1
      }
      # sample a class for the next object
      idx <- sample(m, 1, prob = p / (theta + k))
    }
  } else if (zeros) {
    stop("Preserving zeros not implemented for non-null kappa and m")
  } else {
    # Initialize with a single instance of the first class.
    result <- c(1)
    k <- length(result)
    # Iterate according to sample size
    for (j in 2:n) {
      if (runif(1, 0, 1) < (theta + k * alpha) / (theta + j)) {
        # Add a new class.
        result <- c(result, 1)
        k <- k + 1
      } else {
        # Add to an existing class,
        # with probability related to current class abundance
        i <- sample(k, size = 1, prob = result - alpha)
        result[[i]] <- result[[i]] + 1
      }
    }
  }
  return(result)
}


####
#loop to get dissimilarity for each community

community=data.frame(row.names=1)
time<-as.data.frame(seq(1:10))
names(time)[1]<-paste("time")

for(i in 1:length(time$time)) {
  #simulate community
  spool<-as.data.frame(rCRP(n=500, kappa=1, m=50))
  names(spool)[1]<-paste("abundance")
  spool$timestep<-time$time[i]
  spool$species<-seq(1:50)
  spool$rank=rank(-spool$abundance, ties.method="average")
  
  community=rbind(spool, community)  
}

turnover<-turnover(df=community, time.var="timestep", species.var="species", abundance.var="abundance", metric="total")

immigration<-turnover(df=community, time.var="timestep", species.var="species", abundance.var="abundance", metric="appearance")

loss<-turnover(df=community, time.var="timestep", species.var="species", abundance.var="abundance", metric="disappearance")

mrs<-rank_shift(df=community, time.var = "timestep", species.var = "species", abundance.var="abundance")%>%
  separate(year_pair, c("year1","year2"), sep="-", remove=F)%>%
  select(-year_pair, -year1)

m<-merge(turnover, immigration, by="timestep")
m1<-merge(m, loss, by="timestep")
codyn_metrics<-merge(m1, mrs, by.x="timestep", by.y="year2")


####so kendall rank correlation
##1 code all species by rank in year 1

year1<-community%>%
  filter(timestep==1)%>%
  mutate(rankyr1=rank)%>%
  select(-abundance, -timestep, -rank)

##get matrix to correlate
community_rank<-merge(year1, community, by="species")%>%
  select(-abundance)%>%
  mutate(timestep2=paste("t", timestep, sep=""))%>%
  select(-timestep)%>%
  spread(timestep2, rank)

###
tao<-c(
as.numeric(Kendall(community_rank$rankyr1, community_rank$t2)[1]),
as.numeric(Kendall(community_rank$t2, community_rank$t3)[1]),
as.numeric(Kendall(community_rank$t3, community_rank$t4)[1]),
as.numeric(Kendall(community_rank$t4, community_rank$t5)[1]),
as.numeric(Kendall(community_rank$t5, community_rank$t6)[1]),
as.numeric(Kendall(community_rank$t6, community_rank$t7)[1]),
as.numeric(Kendall(community_rank$t7, community_rank$t8)[1]),
as.numeric(Kendall(community_rank$t8, community_rank$t9)[1]),
as.numeric(Kendall(community_rank$t9, community_rank$t10)[1]))

codyn_metrics$Cor<-tao

pairs(codyn_metrics[2:6])

test

