library(tidyverse)
library(grid)


#source data management code  -- if on desktop, change source code
source('C:\\Users\\lapie\\Desktop\\R files laptop\\C2E\\SEM paper\\C2E_SEM_data processing.R')


#kim's desktop
# setwd('C:\\Users\\la pierrek\\Dropbox (Smithsonian)\\working groups\\CoRRE\\converge_diverge\\datasets\\LongForm')

#kim's laptop
setwd('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\CoRRE\\converge_diverge\\datasets\\LongForm')


theme_set(theme_bw())
theme_update(axis.title.x=element_text(size=40, vjust=-0.35, margin=margin(t=15)), axis.text.x=element_text(size=34, color='black'),
             axis.title.y=element_text(size=40, angle=90, vjust=0.5, margin=margin(r=15)), axis.text.y=element_text(size=34, color='black'),
             plot.title = element_text(size=40, vjust=2),
             panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
             legend.title=element_blank(), legend.text=element_text(size=20))






##anpp difference through time by trt type
#all yrs
ggplot(data=correSEMdataTrt, aes(x=treatment_year, y=anpp_pdiff)) +
  geom_smooth(method='lm', formula=y~x+I(x^2), se=F, aes(color=trt_type)) +
  geom_smooth(method='lm', formula=y~x+I(x^2), size=3, color='black') +
  xlab('Treatment Year') + ylab('ANPP (%) Difference')
#10 yrs
ggplot(data=subset(correSEMdataTrt, treatment_year<11), aes(x=treatment_year, y=anpp_pdiff)) +
  geom_smooth(method='lm', formula=y~x+I(x^2), se=F, aes(color=trt_type)) +
  geom_smooth(method='lm', formula=y~x+I(x^2), size=3, color='black') +
  xlab('Treatment Year') + ylab('ANPP (%) Difference')
#export at 2000 x 800



###anpp difference correlations
anpp_mean_y1 <- ggplot(data=subset(correSEMdataTrt, treatment_year==1), aes(x=composition_diff, y=anpp_pdiff_transform)) +
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('ANPP Difference') + xlab('Community Difference') + ylim(0,0.7) + xlim(0,1)
anpp_mean_y2 <- ggplot(data=subset(correSEMdataTrt, treatment_year==2), aes(x=composition_diff, y=anpp_pdiff_transform)) +   geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('ANPP Difference') + xlab('Community Difference') + ylim(0,0.7) + xlim(0,1)
anpp_mean_y3 <- ggplot(data=subset(correSEMdataTrt, treatment_year==3), aes(x=composition_diff, y=anpp_pdiff_transform)) +   geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('ANPP Difference') + xlab('Community Difference') + ylim(0,0.7) + xlim(0,1)
anpp_mean_y4 <- ggplot(data=subset(correSEMdataTrt, treatment_year==4), aes(x=composition_diff, y=anpp_pdiff_transform)) +   geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('ANPP Difference') + xlab('Community Difference') + ylim(0,0.7) + xlim(0,1)
anpp_mean_y5 <- ggplot(data=subset(correSEMdataTrt, treatment_year==5), aes(x=composition_diff, y=anpp_pdiff_transform)) +   geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('ANPP Difference') + xlab('Community Difference') + ylim(0,0.7) + xlim(0,1)
anpp_mean_y6 <- ggplot(data=subset(correSEMdataTrt, treatment_year==6), aes(x=composition_diff, y=anpp_pdiff_transform)) +   geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('ANPP Difference') + xlab('Community Difference') + ylim(0,0.7) + xlim(0,1)
anpp_mean_y7 <- ggplot(data=subset(correSEMdataTrt, treatment_year==7), aes(x=composition_diff, y=anpp_pdiff_transform)) +   geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('ANPP Difference') + xlab('Community Difference') + ylim(0,0.7) + xlim(0,1)
anpp_mean_y8 <- ggplot(data=subset(correSEMdataTrt, treatment_year==8), aes(x=composition_diff, y=anpp_pdiff_transform)) +   geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('ANPP Difference') + xlab('Community Difference') + ylim(0,0.7) + xlim(0,1)
anpp_mean_y9 <- ggplot(data=subset(correSEMdataTrt, treatment_year==9), aes(x=composition_diff, y=anpp_pdiff_transform)) +   geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('ANPP Difference') + xlab('Community Difference') + ylim(0,0.7) + xlim(0,1)
anpp_mean_y10 <- ggplot(data=subset(correSEMdataTrt, treatment_year==10), aes(x=composition_diff, y=anpp_pdiff_transform)) +   geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('ANPP Difference') + xlab('Community Difference') + ylim(0,0.7) + xlim(0,1)
#with N
anpp_n_y1 <- ggplot(data=subset(correSEMdataTrt, treatment_year==1), aes(x=n, y=anpp_pdiff_transform)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('ANPP Difference') + xlab('N') + ylim(0,0.7) + xlim(0,60)
anpp_n_y2 <- ggplot(data=subset(correSEMdataTrt, treatment_year==2), aes(x=n, y=anpp_pdiff_transform)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('ANPP Difference') + xlab('N') + ylim(0,0.7) + xlim(0,60)
anpp_n_y3 <- ggplot(data=subset(correSEMdataTrt, treatment_year==3), aes(x=n, y=anpp_pdiff_transform)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('ANPP Difference') + xlab('N') + ylim(0,0.7) + xlim(0,60)
anpp_n_y4 <- ggplot(data=subset(correSEMdataTrt, treatment_year==4), aes(x=n, y=anpp_pdiff_transform)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('ANPP Difference') + xlab('N') + ylim(0,0.7) + xlim(0,60)
anpp_n_y5 <- ggplot(data=subset(correSEMdataTrt, treatment_year==5), aes(x=n, y=anpp_pdiff_transform)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('ANPP Difference') + xlab('N') + ylim(0,0.7) + xlim(0,60)
anpp_n_y6 <- ggplot(data=subset(correSEMdataTrt, treatment_year==6), aes(x=n, y=anpp_pdiff_transform)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('ANPP Difference') + xlab('N') + ylim(0,0.7) + xlim(0,60)
anpp_n_y7 <- ggplot(data=subset(correSEMdataTrt, treatment_year==7), aes(x=n, y=anpp_pdiff_transform)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('ANPP Difference') + xlab('N') + ylim(0,0.7) + xlim(0,60)
anpp_n_y8 <- ggplot(data=subset(correSEMdataTrt, treatment_year==8), aes(x=n, y=anpp_pdiff_transform)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('ANPP Difference') + xlab('N') + ylim(0,0.7) + xlim(0,60)
anpp_n_y9 <- ggplot(data=subset(correSEMdataTrt, treatment_year==9), aes(x=n, y=anpp_pdiff_transform)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('ANPP Difference') + xlab('N') + ylim(0,0.7) + xlim(0,60)
anpp_n_y10 <- ggplot(data=subset(correSEMdataTrt, treatment_year==10), aes(x=n, y=anpp_pdiff_transform)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('ANPP Difference') + xlab('N') + ylim(0,0.7) + xlim(0,60)

pushViewport(viewport(layout=grid.layout(9,2)))
print(anpp_mean_y2, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
print(anpp_mean_y3, vp=viewport(layout.pos.row = 2, layout.pos.col = 1))
print(anpp_mean_y4, vp=viewport(layout.pos.row = 3, layout.pos.col = 1))
print(anpp_mean_y5, vp=viewport(layout.pos.row = 4, layout.pos.col = 1))
print(anpp_mean_y6, vp=viewport(layout.pos.row = 5, layout.pos.col = 1))
print(anpp_mean_y7, vp=viewport(layout.pos.row = 6, layout.pos.col = 1))
print(anpp_mean_y8, vp=viewport(layout.pos.row = 7, layout.pos.col = 1))
print(anpp_mean_y9, vp=viewport(layout.pos.row = 8, layout.pos.col = 1))
print(anpp_mean_y10, vp=viewport(layout.pos.row = 9, layout.pos.col = 1))
print(anpp_n_y2, vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
print(anpp_n_y3, vp=viewport(layout.pos.row = 2, layout.pos.col = 2))
print(anpp_n_y4, vp=viewport(layout.pos.row = 3, layout.pos.col = 2))
print(anpp_n_y5, vp=viewport(layout.pos.row = 4, layout.pos.col = 2))
print(anpp_n_y6, vp=viewport(layout.pos.row = 5, layout.pos.col = 2))
print(anpp_n_y7, vp=viewport(layout.pos.row = 6, layout.pos.col = 2))
print(anpp_n_y8, vp=viewport(layout.pos.row = 7, layout.pos.col = 2))
print(anpp_n_y9, vp=viewport(layout.pos.row = 8, layout.pos.col = 2))
print(anpp_n_y10, vp=viewport(layout.pos.row = 9, layout.pos.col = 2))
#export at 2000 x 4000


###Community Difference correlations
#with mean_change
mean_rank_change_transform_y1 <- ggplot(data=subset(correSEMdataTrt, treatment_year==1), aes(x=rank_difference, y=composition_diff)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('Community Difference') + xlab('rank_difference') + ylim(0,1) + xlim(0,0.41)
mean_rank_change_transform_y2 <- ggplot(data=subset(correSEMdataTrt, treatment_year==2), aes(x=rank_difference, y=composition_diff)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('Community Difference') + xlab('rank_difference') + ylim(0,1) + xlim(0,0.41)
mean_rank_change_transform_y3 <- ggplot(data=subset(correSEMdataTrt, treatment_year==3), aes(x=rank_difference, y=composition_diff)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('Community Difference') + xlab('rank_difference') + ylim(0,1) + xlim(0,0.41)
mean_rank_change_transform_y4 <- ggplot(data=subset(correSEMdataTrt, treatment_year==4), aes(x=rank_difference, y=composition_diff)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('Community Difference') + xlab('rank_difference') + ylim(0,1) + xlim(0,0.41)
mean_rank_change_transform_y5 <- ggplot(data=subset(correSEMdataTrt, treatment_year==5), aes(x=rank_difference, y=composition_diff)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('Community Difference') + xlab('rank_difference') + ylim(0,1) + xlim(0,0.41)
mean_rank_change_transform_y6 <- ggplot(data=subset(correSEMdataTrt, treatment_year==6), aes(x=rank_difference, y=composition_diff)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('Community Difference') + xlab('rank_difference') + ylim(0,1) + xlim(0,0.41)
mean_rank_change_transform_y7 <- ggplot(data=subset(correSEMdataTrt, treatment_year==7), aes(x=rank_difference, y=composition_diff)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('Community Difference') + xlab('rank_difference') + ylim(0,1) + xlim(0,0.41)
mean_rank_change_transform_y8 <- ggplot(data=subset(correSEMdataTrt, treatment_year==8), aes(x=rank_difference, y=composition_diff)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('Community Difference') + xlab('rank_difference') + ylim(0,1) + xlim(0,0.41)
mean_rank_change_transform_y9 <- ggplot(data=subset(correSEMdataTrt, treatment_year==9), aes(x=rank_difference, y=composition_diff)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('Community Difference') + xlab('rank_difference') + ylim(0,1) + xlim(0,0.41)
mean_rank_change_transform_y10 <- ggplot(data=subset(correSEMdataTrt, treatment_year==10), aes(x=rank_difference, y=composition_diff)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('Community Difference') + xlab('rank_difference') + ylim(0,1) + xlim(0,0.41)
#with N
mean_n_y1 <- ggplot(data=subset(correSEMdataTrt, treatment_year==1), aes(x=n, y=composition_diff)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('Community Difference') + xlab('N') + ylim(0,1) + xlim(0,60)
mean_n_y2 <- ggplot(data=subset(correSEMdataTrt, treatment_year==2), aes(x=n, y=composition_diff)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('Community Difference') + xlab('N') + ylim(0,1) + xlim(0,60)
mean_n_y3 <- ggplot(data=subset(correSEMdataTrt, treatment_year==3), aes(x=n, y=composition_diff)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('Community Difference') + xlab('N') + ylim(0,1) + xlim(0,60)
mean_n_y4 <- ggplot(data=subset(correSEMdataTrt, treatment_year==4), aes(x=n, y=composition_diff)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('Community Difference') + xlab('N') + ylim(0,1) + xlim(0,60)
mean_n_y5 <- ggplot(data=subset(correSEMdataTrt, treatment_year==5), aes(x=n, y=composition_diff)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('Community Difference') + xlab('N') + ylim(0,1) + xlim(0,60)
mean_n_y6 <- ggplot(data=subset(correSEMdataTrt, treatment_year==6), aes(x=n, y=composition_diff)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('Community Difference') + xlab('N') + ylim(0,1) + xlim(0,60)
mean_n_y7 <- ggplot(data=subset(correSEMdataTrt, treatment_year==7), aes(x=n, y=composition_diff)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('Community Difference') + xlab('N') + ylim(0,1) + xlim(0,60)
mean_n_y8 <- ggplot(data=subset(correSEMdataTrt, treatment_year==8), aes(x=n, y=composition_diff)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('Community Difference') + xlab('N') + ylim(0,1) + xlim(0,60)
mean_n_y9 <- ggplot(data=subset(correSEMdataTrt, treatment_year==9), aes(x=n, y=composition_diff)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('Community Difference') + xlab('N') + ylim(0,1) + xlim(0,60)
mean_n_y10 <- ggplot(data=subset(correSEMdataTrt, treatment_year==10), aes(x=n, y=composition_diff)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  ylab('Community Difference') + xlab('N') + ylim(0,1) + xlim(0,60)

pushViewport(viewport(layout=grid.layout(9,2)))
print(mean_rank_change_transform_y2, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
print(mean_rank_change_transform_y3, vp=viewport(layout.pos.row = 2, layout.pos.col = 1))
print(mean_rank_change_transform_y4, vp=viewport(layout.pos.row = 3, layout.pos.col = 1))
print(mean_rank_change_transform_y5, vp=viewport(layout.pos.row = 4, layout.pos.col = 1))
print(mean_rank_change_transform_y6, vp=viewport(layout.pos.row = 5, layout.pos.col = 1))
print(mean_rank_change_transform_y7, vp=viewport(layout.pos.row = 6, layout.pos.col = 1))
print(mean_rank_change_transform_y8, vp=viewport(layout.pos.row = 7, layout.pos.col = 1))
print(mean_rank_change_transform_y9, vp=viewport(layout.pos.row = 8, layout.pos.col = 1))
print(mean_rank_change_transform_y10, vp=viewport(layout.pos.row = 9, layout.pos.col = 1))
print(mean_n_y2, vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
print(mean_n_y3, vp=viewport(layout.pos.row = 2, layout.pos.col = 2))
print(mean_n_y4, vp=viewport(layout.pos.row = 3, layout.pos.col = 2))
print(mean_n_y5, vp=viewport(layout.pos.row = 4, layout.pos.col = 2))
print(mean_n_y6, vp=viewport(layout.pos.row = 5, layout.pos.col = 2))
print(mean_n_y7, vp=viewport(layout.pos.row = 6, layout.pos.col = 2))
print(mean_n_y8, vp=viewport(layout.pos.row = 7, layout.pos.col = 2))
print(mean_n_y9, vp=viewport(layout.pos.row = 8, layout.pos.col = 2))
print(mean_n_y10, vp=viewport(layout.pos.row = 9, layout.pos.col = 2))
#export at 2000 x 4000