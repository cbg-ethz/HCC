library(bnclustOmics)

#load HCC omics data
HCCdata<-readRDS("/Users/polinasuter/Downloads/HCC/submit/HCCinputs/HCCdata.rds")
names(HCCdata)
dim(HCCdata$M) #24 M nodes
dim(HCCdata$CN) #292 CN nodes
dim(HCCdata$T) #188 T nodes
dim(HCCdata$P) # 116 P nodes
dim(HCCdata$PP)#158 PP nodes

#load penalization matrix
HCCpm<-readRDS("/Users/polinasuter/Downloads/HCC/submit/HCCinputs/HCCpm.rds")
dim(HCCpm) #all 778 nodes
length(which(HCCpm==2))
length(which(HCCpm==1))
length(which(HCCpm<2))


#load blacklist matrix
HCCbl<-readRDS("/Users/polinasuter/Downloads/HCC/submit/HCCinputs/HCCbl.rds")

#bnInfo object, needed only if different IDs aer used for different omics types
#can be created from scratch with bnInfo function: see ?bnInfo
namesHCC<-readRDS("HCCinputs/HCCinfo.rds")

#clustering for k=3
#will take a while ~24 hours to finish
#epmatrix=TRUE -> posterioirs of single edges will be estimated
bnres<-bnclustOmics(HCCdata,namesHCC,HCCbl,HCCpm,epmatrix = TRUE,
                    kclust=3,chixi=0,seed=100,maxEM=10,startpoint="mclustPCA",
                    baseprob=0.4,hardlim=6,deltahl=5,commonspace=TRUE)

#consensus graphs and edge annotation
#load result, same as above
bnres<-readRDS("HCCresults/res_main.rds")

#load p-values from DGE
DElist<-readRDS(file="/Users/polinasuter/Downloads/HCC/submit/HCCresults/DElist.rds")
#load interactions from databases STRING; Omnipath: kinase-substrate and transcription factor-target
DBlist<-readRDS(file="/Users/polinasuter/Downloads/HCC/submit/HCCinputs/DBlist.rds")
#annotate edges from discovered networks
#in the resulting data frame pcl denotes posterior probability of an edge in respective cluster
intconstot<-annotateEdges(bnres,namesHCC,sump=1.2,minp=0.5,minkp=0.9,dblist=DBlist)

#interactions in G_2 with posterioir>0.9
intconstot[which(intconstot$type1=="M" & intconstot$gene1=="TP53" & intconstot$pcl2>0.9),]

#interactions of GLUL-T
intconstot[which((intconstot$gene1=="GLUL" & intconstot$type1=="T") | (intconstot$gene2=="GLUL"  & intconstot$type2=="T")),]

#plot neibouthoods of nodes
library(igraph)
library(plotrix)
#TP53-M
localint<-intconstot[which(intconstot$gene1=="TP53" & intconstot$type1=="M"),]
plot1nodecircle(localint,"TP53",rmult=10,thresh=0.4)
#CTNNB-M
localint<-intconstot[which(intconstot$gene1=="CTNNB1" & intconstot$type1=="M"),]
plot1nodecircle(localint,"CTNNB1",rmult=10,thresh=0.4)

#common hubs: 2T, 9M, 9P
intconstotC<-annotateEdges(bnres,namesHCC,sump=1.2,minp=0.4,minkp=1.1,maxkp=0.9,dblist=DBlist)
comhubs<-names(sort(table(c(intconstotC$from,intconstotC$to)),decreasing=TRUE)[1:20])
comhubs
intconstotC<-intconstotC[which(intconstotC$from%in%comhubs | intconstotC$to%in%comhubs),]
head(intconstotC)

#differential hubs: 17PP, 3P
intconstotD<-annotateEdges(bnres,namesHCC,sump=1.2,minp=1,minkp=0.9,maxkp=0.9,dblist=DBlist)
diffhubs<-names(sort(table(c(intconstotD$from,intconstotD$to)),decreasing=TRUE)[1:20])
diffhubs
intconstotD<-intconstotD[which(intconstotD$from%in%diffhubs | intconstotD$to%in%diffhubs),]
head(intconstotD)


#other k
bnres1<-bnclustOmics(HCCdata,namesHCC,HCCbl,HCCpm,epmatrix = FALSE,
                    kclust=3,chixi=0,seed=100,maxEM=10,startpoint="mclustPCA",
                    baseprob=0.4,hardlim=6,deltahl=5,commonspace=TRUE)

bnres2<-bnclustOmics(HCCdata,namesHCC,HCCbl,HCCpm,epmatrix = FALSE,
                    kclust=3,chixi=0,seed=100,maxEM=10,startpoint="mclustPCA",
                    baseprob=0.4,hardlim=6,deltahl=5,commonspace=TRUE)

bnres4<-bnclustOmics(HCCdata,namesHCC,HCCbl,HCCpm,epmatrix = FALSE,
                    kclust=3,chixi=0,seed=100,maxEM=10,startpoint="mclustPCA",
                    baseprob=0.4,hardlim=6,deltahl=5,commonspace=TRUE)

bnres5<-bnclustOmics(HCCdata,namesHCC,HCCbl,HCCpm,epmatrix = FALSE,
                    kclust=3,chixi=0,seed=100,maxEM=10,startpoint="mclustPCA",
                    baseprob=0.4,hardlim=6,deltahl=5,commonspace=TRUE)



