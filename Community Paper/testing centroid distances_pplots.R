### Testing distance between centroids with Pplots
##
### Author: Kevin Wilcox (wilcoxkr@gmail.com)
### Last updated March 19th, 2018

### set up workspace
library(vegan)
library(devtools)
library(ggplot2)
library(ggthemes)
library(tidyverse)
library(gridExtra)

install_github("mavolio/codyn",ref="RACs_cleaner")

library(codyn)

### read in pplots data
data(pplots)

### calculate distance between centroids using pplots data with all treatments
pplots_out_full <- centroid_difference(pplots, 
                          time.var="year", 
                          species.var="species", 
                          abundance.var="relative_cover", 
                          treatment.var="treatment", 
                          replicate.var="plot")

### calculate distance between centroids using pplots data with one treatment removed
pplots_small <- subset(pplots, treatment != "N2P3")
pplots_small$treatment <- factor(pplots_small$treatment)

test_small <- centroid_difference(df_small, 
                          time.var="year", 
                          species.var="species", 
                          abundance.var="relative_cover", 
                          treatment.var="treatment", 
                          replicate.var="plot")

### create data frames with output (data type issues)
df_small_2 <- data.frame(year=test_small$year,
                         treatment=test_small$treatment,
                         treatment2=test_small$treatment2,
                         dist_small=test_small$centroid_distance_diff)

df_full_2 <- data.frame(year=test_full$year,
                        treatment=test_full$treatment,
                         treatment2=test_full$treatment2,
                         dist_full=test_full$centroid_distance_diff)

### Combine small and full output
df_both <- df_small_2 %>%
  full_join(df_full_2, by=c("year", "treatment","treatment2"))

### Plot distances
a <- ggplot(test_full, aes(x=year_pair, y=centroid_distance_change, col=treatment)) +
  geom_point()+
  theme_classic()

b <- ggplot(test_short, aes(x=year_pair, y=centroid_distance_change, col=treatment)) +
  geom_point()+
  theme_classic()
 
grid.arrange(a,b)

### for 1:1 line
full_sub <- subset(test_full, year_pair=="2004-2005")
short_sub <- subset(test_short, year_pair=="2004-2005")

test_df <- full_sub %>%
  cbind(short_sub$centroid_distance_change)

ggplot(df_both, aes(x=dist_full, y=dist_small))+
  geom_point() +
  geom_abline(slope=1, intercept=0) +
  geom_smooth(method="lm",se=F)

### PCoA ###
### create wide format data frames
pplots_wide <- pplots %>%
  spread(key=species, value=relative_cover, fill=0)
pplots_small_wide <- pplots_small %>%
  spread(key=species, value=relative_cover, fill=0)

### Dissimilarity matrices
dist_full <- vegdist(pplots_wide[,5:ncol(pplots_wide)], method="bray")
dist_small <- vegdist(pplots_small_wide[,5:ncol(pplots_small_wide)], method="bray")

### run PCoAs on full and small datasets
pcoa_out_full <- cmdscale(dist_full) %>%
  cbind(pplots_wide[,1:4]) %>%
  rename("axis_1_full" = "1") %>%
  rename("axis_2_full" = "2")

pcoa_out_small <- cmdscale(dist_small) %>%
  cbind(pplots_small_wide[,1:4]) %>%
  rename("axis_1_small" = "1") %>%
  rename("axis_2_small" = "2")

### Combine
pcoa_out_both <- pcoa_out_full %>%
  full_join(pcoa_out_small, by=c("treatment","plot","block","year"))

### plot 1:1 lines
ggplot(pcoa_out_both, aes(x=axis_1_small,y=axis_1_full)) +
  geom_point()+
  theme_classic() +
  geom_abline(slope=1, intercept=0) +
  geom_smooth(method="lm",se=F)

ggplot(pcoa_out_both, aes(x=axis_2_small,y=axis_2_full)) +
  geom_point()+
  theme_classic() +
  geom_abline(slope=1, intercept=0) +
  geom_smooth(method="lm",se=F)

########## CDR data ################## 
cdr <- read.csv("C:\\Users\\wilco\\Dropbox\\C2E\\Products\\CommunityChange\\March2018 WG\\SpeciesRelativeAbundance_Oct2017.csv") %>%
  filter(site_code=="CDR" & project_name=="e001" & community_type=="B")

cdr_small <- cdr %>%
  filter(treatment!=8)

diff_out_full <- centroid_difference(cdr, 
                                  time.var="calendar_year", 
                                  species.var="genus_species", 
                                  abundance.var="relcov", 
                                  treatment.var="treatment", 
                                  replicate.var="plot_id")

diff_out_small <- centroid_difference(cdr_small, 
                                  time.var="calendar_year", 
                                  species.var="genus_species", 
                                  abundance.var="relcov", 
                                  treatment.var="treatment", 
                                  replicate.var="plot_id")

diff_out_small_2 <- data.frame(year=diff_out_small$calendar_year,
                         treatment=diff_out_small$treatment,
                         treatment2=diff_out_small$treatment2,
                         dist_small=diff_out_small$centroid_distance_diff)

diff_out_full_2 <- data.frame(year=diff_out_full$calendar_year,
                        treatment=diff_out_full$treatment,
                        treatment2=diff_out_full$treatment2,
                        dist_full=diff_out_full$centroid_distance_diff)

diff_out_both <- diff_out_small_2 %>%
  full_join(diff_out_full_2, by=c("year", "treatment","treatment2"))

ggplot(diff_out_both, aes(x=dist_full, y=dist_small))+
  geom_point() +
  geom_abline(slope=1, intercept=0) +
  geom_smooth(method="lm",se=F)


