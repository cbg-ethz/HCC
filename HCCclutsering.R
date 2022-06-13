library(bnclustOmics)

#load HCC omics data
HCCdata<-readRDS("HCCinputs/HCCdata.rds")
names(HCCdata)
dim(HCCdata$M) #24M nodes
dim(HCCdata$CN) #292CN nodes
dim(HCCdata$T) #188 T nodes
dim(HCCdata$P) #116 P nodes
dim(HCCdata$PP)#158 PP nodes

#load penalization matrix
HCCpm<-readRDS("HCCinputs/HCCpm.rds")
dim(HCCpm) #all 778 nodes
#histogram of penalization factors of all possible edges in the search space
hist(pmHCC)


#load blacklist matrix
HCCbl<-readRDS("HCCinputs/HCCbl.rds")

#bnInfo object, needed only if different IDs aer used for different omics types
#can be created from scratch with bnInfo function: see ?bnInfo
#mappings<-readRDS("HCCinputs/mappings.rds")
#namesHCC<-bnInfo(HCCdata,c("b","o","c","c","c"),c("M","CN","T","P","PP"),mappings=mappings,attachtype=FALSE)
namesHCC<-readRDS("HCCinputs/HCCinfo.rds")


#clustering for k=3
#will take a while ~24 hours to finish
#faster with epmatrix=FALSE (5 hours), but we need epmatrix for downstream analysis
#epmatrix=TRUE -> posterioirs of single edges will be estimated
bnres<-bnclustOmics(HCCdata,namesHCC,HCCbl,HCCpm,epmatrix = TRUE,
                    kclust=3,chixi=0,seed=100,maxEM=10,startpoint="mclustPCA",
                    baseprob=0.4,hardlim=6,deltahl=5,commonspace=TRUE)

#load result, same as above
bnres<-readRDS("HCCresults/res_main.rds")
bnresn<-readRDS("/Users/polinasuter/Downloads/jun223TT.rds")


compareDAGs(bnres$DAGs[[3]],bnresn$DAGs[[3]])

bnres$likel
bnresn$likel



#other k
#it is important that epmatrix=FALSE for below runs
#sampling from posterior distribution (epmatrix=TRUE) takes longer
#than finding MAP model, for choosing k we do not need epmatrix
#all runs bnres1 - bners5 take 20hrs together
bnres1<-bnclustOmics(HCCdata,namesHCC,HCCbl,HCCpm,epmatrix = FALSE,
                     kclust=2,baseprob=0.6)

bnres2<-bnclustOmics(HCCdata,namesHCC,HCCbl,HCCpm,epmatrix = FALSE,
                     kclust=2, baseprob=0.6)

bnres4<-bnclustOmics(HCCdata,namesHCC,HCCbl,HCCpm,epmatrix = FALSE,
                     kclust=4,baseprob=0.3)

bnres5<-bnclustOmics(HCCdata,namesHCC,HCCbl,HCCpm,epmatrix = FALSE,
                     kclust=5,baseprob=0.3)

#or load the result from the folder
bnres1<-readRDS("HCCresults/res1.rds")
bnres2<-readRDS("HCCresults/res2.rds")
bnres4<-readRDS("HCCresults/res4.rds")
bnres5<-readRDS("HCCresults/res5.rds")

bnlist<-list(bnres1,bnres2,bnres,bnres4,bnres5)
chooseK(bnlist,fun="AIC")
chooseK(bnlist,fun="BIC")

#we choose bnres (k=3) for the downstream analysis

#load clinical data
HCC_surv_df<-readRDS("HCCinputs/info.rds")
source("helpresfns.R")
nonadjustedCox(bnres$memb,HCC_surv_df)
adjustedCox(bnres$memb,HCC_surv_df)


#load interactions from databases STRING; Omnipath: kinase-substrate and transcription factor-target
DBlist<-readRDS(file="HCCinputs/DBlist.rds")
#annotate edges from discovered networks
#in the resulting data frame pcl denotes posterior probability of an edge in respective cluster
intconstot<-annotateEdges(bnres,namesHCC,sump=1.2,minp=0.4,minkp=0.9,dblist=DBlist)

#plot neibouthoods of nodes
library(igraph)
library(plotrix)
#TP53-M
plotNode(intconstot,"TP53",rmult=11)
#CTNNB-M
plotNode(intconstot,"CTNNB1",rmult=11)

#common hubs: 2T, 9M, 9P
#set minkp to higher than one to exclude cluster-specific connections
intconstotC<-annotateEdges(bnres,namesHCC,sump=1.2,minp=0.4,minkp=1.1,maxkp=0.9,dblist=DBlist)
comhubs<-names(sort(table(c(intconstotC$from,intconstotC$to)),decreasing=TRUE)[1:20])
comhubs
intconstotC<-intconstotC[which(intconstotC$from%in%comhubs | intconstotC$to%in%comhubs),]
head(intconstotC)

#differential hubs: 17PP, 3P
#set minp to 1 to exclude connections, present with non-zero posterior in more than one cluster
intconstotD<-annotateEdges(bnres,namesHCC,sump=1.2,minp=1,minkp=0.9,maxkp=0.9,dblist=DBlist)
diffhubs<-names(sort(table(c(intconstotD$from,intconstotD$to)),decreasing=TRUE)[1:20])
diffhubs
intconstotD<-intconstotD[which(intconstotD$from%in%diffhubs | intconstotD$to%in%diffhubs),]
head(intconstotD)





