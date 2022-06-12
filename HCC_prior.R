######################################################################
#Constructing graphical prior (penalization matrix) for bnClustOmics #
#using the interactions list from STRING database                    #
######################################################################


#loading the object of class bnInfo
#for construction instructions see HCCclustering.R
namesHCC<-readRDS("HCCinputs/HCCinfo.rds")

#we will need the saved penalization matrix to compare it to the
#constructed penalization matrix
HCCpm<-readRDS("HCCinputs/HCCpm.rds")

#loading the list of interactions saved from the STRING web-site
#https://string-db.org/
stringintlist<-read.csv2("DBdata/HCCint3.tsv",sep="\t")

#loading the the mapping between STRING IDs and IDs used in the analysis
#such mapping can be downloaded from the STRING web-site
#https://string-db.org/
stringmapping<-read.csv2("DBdata/HCCmap3.tsv",sep="\t")

#interactions from BIOGRID human database
#can be downloaded from the web-site
biogridhuman<-read.delim("DBdata/biogrid_human.txt", header = TRUE, sep = "\t")


stringchange<-stringmapping[(which(stringmapping$queryItem!=stringmapping$preferredName)),c(1:4)]
stringchange<-stringchange[-which(duplicated(stringchange$preferredName)),]
rownames(stringchange)<-stringchange$preferredName
changename1<-which(stringintlist$X.node1%in%stringchange$preferredName)
changename2<-which(stringintlist$node2%in%stringchange$preferredName)
stringintlist$X.node1[changename1]<-stringchange[stringintlist$X.node1[changename1] ,"queryItem"]
stringintlist$node2[changename2]<-stringchange[stringintlist$node2[changename2] ,"queryItem"]
stringintlist<-stringintlist[,c(1,2,13)]
colnames(stringintlist)<-c("gene1","gene2","score")
stringintlist$score<-as.numeric(stringintlist$score)

#histogram of interaction scores
hist(as.numeric(stringintlist$combined_score))



#fill the penalization matrix with values according to interactions from the STRING database
#pfbase base interaction factor, will be used for edges, for which no interaction was found in STRING
#intpf starting value for interactions found in the database
#intsame value for penalization factor foe IDs encoding for the same gene, e.g.
#transcript TP3 and protein TP53
#the penalization factor for edges found in the STRING database will be calculated as
#max(1,intpf-2*interaction_score)
#may take 2 minutes or so
pmHCC<-penInit(namesHCC,pfbase=2,intpf=2,intsame=2,intlist=stringintlist,usescore=TRUE)

#now we define gene products from different omics types that should not be panalized
#when both gene products correspond to the same gene X
#we simply update existing penalization matrix using intsame=1
pmHCC<-penUpdateInter(pmHCC,namesHCC,"CN","T",NULL,pfbase=2,intpf=2,intsame=1)
pmHCC<-penUpdateInter(pmHCC,namesHCC,"CN","P",NULL,pfbase=2,intpf=2,intsame=1)
pmHCC<-penUpdateInter(pmHCC,namesHCC,"CN","PP",NULL,pfbase=2,intpf=2,intsame=1)
pmHCC<-penUpdateInter(pmHCC,namesHCC,"M","T",NULL,pfbase=2,intpf=2,intsame=1)
pmHCC<-penUpdateInter(pmHCC,namesHCC,"M","P",NULL,pfbase=2,intpf=2,intsame=1)
pmHCC<-penUpdateInter(pmHCC,namesHCC,"M","PP",NULL,pfbase=2,intpf=2,intsame=1)
pmHCC<-penUpdateInter(pmHCC,namesHCC,"T","P",NULL,pfbase=2,intpf=2,intsame=1)
pmHCC<-penUpdateInter(pmHCC,namesHCC,"T","PP",NULL,pfbase=2,intpf=2,intsame=1)


#check that our ID mappings have worked
#GRHL2 and ENSG00000083307 are IDs of the same gene, hence penalization factor should be 1 (no penalization)
#SQLE and Q14534 are IDs of the same gene, hence penalization factor should be 1 (no penalization)
pmHCC["GRHL2","ENSG00000083307"]
pmHCC["GRHL2","ENSG00000104549"]
pmHCC["SQLE","Q14534"]

#add interactions from omnipath
#first, we download interactions
library(OmnipathR)
org<-9606
uniqueprots<-function(phlist,u=TRUE) {
  if(u) {
    return(unique(sub('_.*',"",phlist)))
  } else {
    return(sub('_.*',"",phlist))
  }
}
#get kinase-substrate interactions
kinaseint <-import_kinaseextra_interactions(resources=c("PhosphoPoint",
                                              "PhosphoSite"), organism = org)
#get protein-protein interactions
KSint<-as.data.frame(kinaseint)[which(kinaseint$target%in%uniqueprots(namesHCC$allnamesonebn) &
                                        kinaseint$source%in%uniqueprots(namesHCC$allnamesonebn)),]
nrow(KSint)
KSint<-KSint[,c(3,4)]
colnames(KSint)[c(1,2)]<-c("gene1","gene2")

pmHCC<-penUpdateIntra(pmHCC,namesHCC,"PP",KSint,pfbase=2,intpf=1,intsame=1,bi=FALSE)
pmHCC<-penUpdateIntra(pmHCC,namesHCC,"P",KSint,pfbase=2,intpf=1,intsame=2,bi=FALSE)
pmHCC<-penUpdateInter(pmHCC,namesHCC,"P","PP",KSint,pfbase=2,intpf=1,intsame=1,bi=FALSE)
pmHCC<-penUpdateInter(pmHCC,namesHCC,"PP","P",KSint,pfbase=2,intpf=1,intsame=1,bi=FALSE)

hist(pmHCC)

#check if the constructed object corresponds to the saved object
nrow(which(pmHCC!=HCCpm,arr.ind=TRUE))

#compare STRING and BIOGRID priors
#download the latest version of BIOGRID interactions in human
stringintlist$biogrid<-FALSE
for(i in 1:nrow(stringintlist)){
  intgene<-which(biogridhuman$Official.Symbol.Interactor.A==stringintlist$gene1[i])
  if(length(intgene)>0) {
    ints<-biogridhuman$Official.Symbol.Interactor.B[intgene]
    if(stringintlist$gene2[i]%in%ints)
      stringintlist$biogrid[i]<-TRUE
  }

  intgene<-which(biogridhuman$Official.Symbol.Interactor.B==stringintlist$gene1[i])
  if(length(intgene)>0) {
    ints<-biogridhuman$Official.Symbol.Interactor.A[intgene]
    if(stringintlist$gene2[i]%in%ints)
      stringintlist$biogrid[i]<-TRUE
  }
}

#plot the scores of interactions from STRING that are also present in BIOGRID
#to interactions from STRING that are absent in BIOGRID
boxplot(stringintlist$score[which(stringintlist$biogrid==TRUE)],
        stringintlist$score[which(stringintlist$biogrid==FALSE)],
        names=c("in BIOGRID", "absent"),ylab=c("STRING score"))





