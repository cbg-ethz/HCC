library(survival)
library(survminer)
library(ggplot2)

source("helpresfns.R")

#survival analysis
memb2check<-bnres$memb #plug in memberships by other algorithms here as well
HCC_surv_df<-readRDS("HCCinputs/info.rds") #load clinical data
pats<-HCC_surv_df$patientID

#for non-adjusted model
curru<-checkduplicated(memb2check,pats) #check patients with duplicated biopsies (same cluster or not, if not - exclude both, if same include one)
survobj<-HCC_surv_df[curru,]
survobj<-cbind(survobj,as.factor(memb2check[curru]))
nrow(survobj) #how many samples were included
colnames(survobj)[ncol(survobj)]<-"cluster"
surv_object <- Surv(time = survobj$survival_time, event = survobj$death)
fit2 <- survfit(surv_object ~ cluster, data = survobj)
ggsurvplot(fit2, data = survobj)

fitmod <- coxph(surv_object ~  cluster, survobj)
#here is the result for non-adjusted model
print(fitmod, digits=3)
#optional hazard ratio graph
#ggforest(fitmod,fontsize=1.2)


#for adjusted model, exclude patient 41, group '0'
curru<-checkduplicated(memb2check,pats,addex=c(41))
survobj<-HCC_surv_df[curru,]
survobj<-cbind(survobj,as.factor(memb2check[curru]))
nrow(survobj) #how many samples were included
colnames(survobj)[ncol(survobj)]<-"cluster"
surv_object <- Surv(time = survobj$survival_time, event = survobj$death)

fitbase <- coxph(surv_object ~  bclc, survobj)
print(fitbase, digits=3)
fitbn <- update(fitbase, . ~ . + cluster)
print(fitbn, digits=3)
#here is the result for adjusted model
anova(fitbase, fitbn)

