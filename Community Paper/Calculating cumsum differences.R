?cumsum
rm(list=ls())
df <- data.frame(group=c(rep("a",6),rep("b",6)), a=1:6, b=1:12)

### pull # of comparisons
cumsum(df$a)

g_csum <- df %>%
  group_by(group) %>%
  mutate(cumsum=cumsum(b))

library(tidyverse)
library(ggplot2)

rr_full <- read.csv("C:\\Users\\wilco\\Dropbox\\C2E\\Products\\CommunityChange\\March2018 WG\\CORRE_RAC_LogRR_March2018_trtyr.csv")

site_project_comm_vec <- unique(rr_full$site_project_comm)

rr_test <- filter(rr_full, site_project_comm==site_project_comm_vec[1])

rr_cumsum <- rr_test %>%
  
