setwd("/Users/egrman/Documents/RESEARCH/LTER convergence divergence/C2E working group march 2018/c2ecommunitypaper")

setwd("~/Dropbox/")

dat=read.csv("C2E/Products/CommunityChange/Summer2018_Results/sig_diff_by_trts_May2019.csv")


library(ggplot2)
library(reshape2)
library(car)


melt.dat=melt(dat, id=c("site_project_comm", "treatment", "response_var", "site_code", "project_name", "community_type", "trt_type2"))
counts.samplesize=dcast(melt.dat, response_var + trt_type2 ~ variable, fun=length)
names(counts.samplesize)[names(counts.samplesize)=="sigdiff"]="number.studies"
counts.signif=dcast(melt.dat, response_var + trt_type2 ~ variable, fun=sum)
counts=merge(counts.samplesize, counts.signif, by=c("response_var", "trt_type2"))
counts$perc.change=counts$sigdiff/counts$number.studies
ggplot(aes(response_var, perc.change), data=counts) + coord_flip() + geom_bar(stat="identity") + facet_wrap(~trt_type2)

#are some GCD more likely to change some aspects of the community than other aspects?: NO

#binomial distributions can suffer from overdispersion too! values of this (function below) should be about 1 if a model is not overdispersed; values above 1.5 are problematic. overdispersion doesn't affect parameter estimates but it does affect uncertainty around them. If you do have overdispersion, then have to use quasibinomial, which isn't a real distribution so you can't get AIC from them. 

dispersion=function(model) {
  sum(residuals(model, type="pearson")^2)/(length(model$y) - length(model$coefficients))
}

mod=glm(sigdiff/number.studies~response_var * trt_type2, data=counts, weights=number.studies, family=binomial)
mod2=glm(sigdiff/number.studies~response_var + trt_type2, data=counts, weights=number.studies, family=binomial)
anova(mod, mod2, test="Chi")
dispersion(mod)

mod=glm(sigdiff/number.studies~response_var * trt_type2, data=counts, weights=number.studies, family=quasibinomial)
dispersion(mod) #infinite??

#are some aspects of community change more likely to occur than others?: NO
mod=glm(sigdiff/number.studies~response_var, data=counts, weights=number.studies, family=quasibinomial); Anova(mod); Anova(mod, test="F") #p=0.7 with either Chisq or F test

#are some GCD more likely to cause ANY change in the community?: NO
anychange=dcast(melt.dat, site_project_comm + treatment + site_code + project_name + community_type + trt_type2 ~ variable, fun=max) #stack overflow says it's OK to ignore the warning: https://stackoverflow.com/questions/24282550/no-non-missing-arguments-warning-when-using-min-or-max-in-reshape2
counts.anychange=dcast(anychange, trt_type2~., value.var="sigdiff", fun=length)
names(counts.anychange)[names(counts.anychange)=="."]="number.studies"
counts.signif.anychange=dcast(anychange, trt_type2~., value.var="sigdiff", fun=sum)
names(counts.signif.anychange)[names(counts.signif.anychange)=="."]="sigdiff"
counts.anychange=merge(counts.anychange, counts.signif.anychange, by="trt_type2")
counts.anychange$perc.change= counts.anychange$sigdiff/counts.anychange$number.studies
ggplot(aes(trt_type2, perc.change), data=counts.anychange) + coord_flip() + geom_bar(stat="identity")

mod3=glm(sigdiff/number.studies~trt_type2, data=counts.anychange, weights=number.studies, family=binomial(link="logit"))
dispersion(mod3)

mod3=glm(sigdiff/number.studies~trt_type2, data=counts.anychange, weights=number.studies, family=quasibinomial); Anova(mod3); Anova(mod3, test="F")

#Even though the interaction isn't significant, do GCD differ in their effects on each aspect of community change (one at a time)?: NO (marginally signifiant CO2)
mod=glm(sigdiff/number.studies~trt_type2, data=counts[counts$response_var=="evenness_change_abs",], weights=number.studies, family=binomial(link="logit"))
mod2=glm(sigdiff/number.studies~1, data=counts[counts$response_var=="evenness_change_abs",], weights=number.studies, family=binomial(link="logit"))
anova(mod, mod2, test="Chi")

mod=glm(sigdiff/number.studies~trt_type2, data=counts[counts$response_var=="gains",], weights=number.studies, family=binomial(link="logit"))
mod2=glm(sigdiff/number.studies~1, data=counts[counts$response_var=="gains",], weights=number.studies, family=binomial(link="logit"))
anova(mod, mod2, test="Chi")

mod=glm(sigdiff/number.studies~trt_type2, data=counts[counts$response_var=="losses",], weights=number.studies, family=binomial(link="logit"))
mod2=glm(sigdiff/number.studies~1, data=counts[counts$response_var=="losses",], weights=number.studies, family=binomial(link="logit"))
anova(mod, mod2, test="Chi")

mod=glm(sigdiff/number.studies~trt_type2, data=counts[counts$response_var=="rank_change",], weights=number.studies, family=binomial(link="logit"))
mod2=glm(sigdiff/number.studies~1, data=counts[counts$response_var=="rank_change",], weights=number.studies, family=binomial(link="logit"))
anova(mod, mod2, test="Chi")
summary(mod)

mod=glm(sigdiff/number.studies~trt_type2, data=counts[counts$response_var=="richness_change_abs",], weights=number.studies, family=binomial(link="logit"))
mod2=glm(sigdiff/number.studies~1, data=counts[counts$response_var=="richness_change_abs",], weights=number.studies, family=binomial(link="logit"))
anova(mod, mod2, test="Chi")
dispersion(mod)


