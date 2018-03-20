### calculating permanova for each experiment and year
###
### Authors: Kevin Wilcox (wilcoxkr@gmail.com) Andrew Tredennick (atredenn@gmail.com)
### Last updated: March 20 2018



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