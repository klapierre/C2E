library(tidyverse)
library(grid)
library(PerformanceAnalytics)
library(lme4)
library(lmerTest)
library(piecewiseSEM)




#source data management code  -- if on desktop, change source code
source('C:\\Users\\lapie\\Desktop\\R files laptop\\C2E\\SEM paper\\C2E_SEM_data processing.R')



#kim's desktop
setwd('C:\\Users\\la pierrek\\Dropbox (Smithsonian)\\working groups\\CoRRE\\converge_diverge\\datasets\\LongForm')

#kim's laptop
setwd('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\CoRRE\\converge_diverge\\datasets\\LongForm')


theme_set(theme_bw())
theme_update(axis.title.x=element_text(size=40, vjust=-0.35, margin=margin(t=15)), axis.text.x=element_text(size=34, color='black'),
             axis.title.y=element_text(size=40, angle=90, vjust=0.5, margin=margin(r=15)), axis.text.y=element_text(size=34, color='black'),
             plot.title = element_text(size=40, vjust=2),
             panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
             legend.title=element_blank(), legend.text=element_text(size=20))


###difference metrics
#all years
summary(compositionModel <- psem(
  lmer(anpp_pdiff ~ composition_diff_transform + (1|site_project_comm), data=correSEMdataTrt),
  lmer(composition_diff_transform ~ richness_difference + evenness_diff + rank_difference + species_difference + (1|site_project_comm), data=correSEMdataTrt))
  )

###testing our multiple regressions with random effect of site_project_comm
summary(lmer(anpp_pdiff ~ richness_difference*treatment_year + evenness_diff*treatment_year + rank_difference*treatment_year + species_difference*treatment_year + (1|site_project_comm), data=correSEMdataTrt))

#with transformed variables
summary(lmer(anpp_pdiff_transform ~ richness_difference*treatment_year + evenness_diff*treatment_year + rank_difference*treatment_year + species_difference*treatment_year + (1|site_project_comm), data=correSEMdataTrt))

ggplot(correSEMdataTrt, aes(x=treatment_year, y=anpp_pdiff_transform, color=as.character(site_project_comm_trt))) + 
  geom_point() + 
  guides(color = "none") + 
  stat_smooth(method="lm", formula = y ~ x + I(x^2), fill=NA)



###all data, split into single years
#year 1
summary(compositionModel1 <- psem(
  lm(anpp_pdiff_transform ~ composition_diff_transform + richness_difference + evenness_diff + rank_difference + species_difference, data=subset(correSEMdataTrt, treatment_year==1)),
  lm(composition_diff_transform ~ richness_difference + evenness_diff + rank_difference + species_difference, data=subset(correSEMdataTrt, treatment_year==1)),
  richness_difference %~~% evenness_diff,
  richness_difference %~~% rank_difference,
  richness_difference %~~% species_difference,
  evenness_diff %~~% rank_difference,
  evenness_diff %~~% species_difference,
  rank_difference %~~% species_difference
  )
)
coefs1 <- coefs(compositionModel1, standardize = "scale", standardize.type = "latent.linear", intercepts = FALSE)%>%
  select(Response, Predictor, Estimate, Std.Error, DF, Crit.Value, P.Value, Std.Estimate)%>%
  mutate(treatment_year=1)

#year 2
summary(compositionModel2 <- psem(
  lm(anpp_pdiff ~ richness_difference + evenness_diff + rank_difference + species_difference, data=subset(correSEMdataTrt, treatment_year==2)))
)
coefs2 <- coefs(compositionModel2, standardize = "scale", standardize.type = "latent.linear", intercepts = FALSE)%>%
  select(Response, Predictor, Estimate, Std.Error, DF, Crit.Value, P.Value, Std.Estimate)%>%
  mutate(treatment_year=2)

#year 3
summary(compositionModel3 <- psem(
  lm(anpp_pdiff ~ richness_difference + evenness_diff + rank_difference + species_difference, data=subset(correSEMdataTrt, treatment_year==3)))
)
coefs3 <- coefs(compositionModel3, standardize = "scale", standardize.type = "latent.linear", intercepts = FALSE)%>%
  select(Response, Predictor, Estimate, Std.Error, DF, Crit.Value, P.Value, Std.Estimate)%>%
  mutate(treatment_year=3)

#year 4
summary(compositionModel4 <- psem(
  lm(anpp_pdiff ~ richness_difference + evenness_diff + rank_difference + species_difference, data=subset(correSEMdataTrt, treatment_year==4)))
)
coefs4 <- coefs(compositionModel4, standardize = "scale", standardize.type = "latent.linear", intercepts = FALSE)%>%
  select(Response, Predictor, Estimate, Std.Error, DF, Crit.Value, P.Value, Std.Estimate)%>%
  mutate(treatment_year=4)

#year 5
summary(compositionModel5 <- psem(
  lm(anpp_pdiff ~ richness_difference + evenness_diff + rank_difference + species_difference, data=subset(correSEMdataTrt, treatment_year==5)))
)
coefs5 <- coefs(compositionModel5, standardize = "scale", standardize.type = "latent.linear", intercepts = FALSE)%>%
  select(Response, Predictor, Estimate, Std.Error, DF, Crit.Value, P.Value, Std.Estimate)%>%
  mutate(treatment_year=5)

#year 6
summary(compositionModel6 <- psem(
  lm(anpp_pdiff ~ richness_difference + evenness_diff + rank_difference + species_difference, data=subset(correSEMdataTrt, treatment_year==6)))
)
coefs6 <- coefs(compositionModel6, standardize = "scale", standardize.type = "latent.linear", intercepts = FALSE)%>%
  select(Response, Predictor, Estimate, Std.Error, DF, Crit.Value, P.Value, Std.Estimate)%>%
  mutate(treatment_year=6)

#year 7
summary(compositionModel7 <- psem(
  lm(anpp_pdiff ~ richness_difference + evenness_diff + rank_difference + species_difference, data=subset(correSEMdataTrt, treatment_year==7)))
)
coefs7 <- coefs(compositionModel7, standardize = "scale", standardize.type = "latent.linear", intercepts = FALSE)%>%
  select(Response, Predictor, Estimate, Std.Error, DF, Crit.Value, P.Value, Std.Estimate)%>%
  mutate(treatment_year=7)

#year 8
summary(compositionModel8 <- psem(
  lm(anpp_pdiff ~ richness_difference + evenness_diff + rank_difference + species_difference, data=subset(correSEMdataTrt, treatment_year==8)))
)
coefs8 <- coefs(compositionModel8, standardize = "scale", standardize.type = "latent.linear", intercepts = FALSE)%>%
  select(Response, Predictor, Estimate, Std.Error, DF, Crit.Value, P.Value, Std.Estimate)%>%
  mutate(treatment_year=8)

#year 9
summary(compositionModel9 <- psem(
  lm(anpp_pdiff ~ richness_difference + evenness_diff + rank_difference + species_difference, data=subset(correSEMdataTrt, treatment_year==9)))
)
coefs9 <- coefs(compositionModel9, standardize = "scale", standardize.type = "latent.linear", intercepts = FALSE)%>%
  select(Response, Predictor, Estimate, Std.Error, DF, Crit.Value, P.Value, Std.Estimate)%>%
  mutate(treatment_year=9)

#year 10
summary(compositionModel10 <- psem(
  lm(anpp_pdiff ~ richness_difference + evenness_diff + rank_difference + species_difference, data=subset(correSEMdataTrt, treatment_year==10)))
)
coefs10 <- coefs(compositionModel10, standardize = "scale", standardize.type = "latent.linear", intercepts = FALSE)%>%
  select(Response, Predictor, Estimate, Std.Error, DF, Crit.Value, P.Value, Std.Estimate)%>%
  mutate(treatment_year=10)

#combine the coeffecients
metricsModelCoef <- rbind(coefs1, coefs2, coefs3, coefs4, coefs5, coefs6, coefs7, coefs8, coefs9, coefs10)
metricsSum <- metricsModelCoef%>%
  group_by(treatment_year)%>%
  summarise(total_std_estimate=sum(abs(Std.Estimate)))%>%
  ungroup()

###figure of drivers of ANPP percent difference
ggplot(data=metricsModelCoef, aes(x=treatment_year, y=Std.Estimate)) +
  geom_point(size=5) +
  geom_smooth() +
  facet_wrap(~Predictor)
