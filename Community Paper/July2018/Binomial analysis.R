setwd("/Users/egrman/Documents/RESEARCH/LTER convergence divergence/C2E working group march 2018/c2ecommunitypaper")

dat=read.csv("sig_diff_by_trts.csv")

dat<-sig_diff_by_trt<-gamtrts_metrics_em

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
mod=glm(sigdiff/number.studies~response_var * trt_type2, data=counts, weights=number.studies, family=binomial(link="logit"))
mod2=update(mod, ~.-response_var:trt_type2)
anova(mod, mod2, test="Chi") #p=0.98 so number of significants for response vars don't depend on treatment

#are some aspects of community change more likely to occur than others?: NO
mod=glm(sigdiff/number.studies~response_var, data=counts, weights=number.studies, family=binomial(link="logit"))
mod2=glm(sigdiff/number.studies~1, data=counts, weights=number.studies, family=binomial(link="logit"))
anova(mod, mod2, test="Chi") #p=0.7

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
mod4=glm(sigdiff/number.studies~1, data=counts.anychange, weights=number.studies, family=binomial(link="logit"))
anova(mod3, mod4, test="Chi") #p=0.2

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



