#new
initDAG<-function(n, d,lB=0.1, uB=1,type=c("cont", "bin", "mixed"),
                  wmpct=0, wmmin=1,wmmax=3, shdpct=0, wm=NULL,
                  mixedpar=list(nbin,avpar,avmixed,lB,uB)) {
  res<-list()
  res$type=type
  if(res$type=="mixed") res$mixedpar<-mixedpar
  res$lB<-lB
  res$uB<-uB
  res$wmpct<-wmpct
  res$shdpct<-shdpct
  res$wmmin<-wmmin
  res$wmmax<-wmmax
  res$wm<-wm
  
  return(res)
}

#not finished new
initsimBN<-function(nrep=10, DAGpar, algorithm=c("pc","ges","mcmcMAP","mcmcsample","iterativesample"), 
                pcparam=list(alphas=c(0.01,0.05,0.1,0.2,0.3,0.4),test="gaussCItest"),
                gesparam=list(lambdas=c(0.5,1,1.5,2,3,5,7), bootstrap=FALSE),
                mcmcpar=list(alpha=0.1, p=c(0.95,0.99,0.7,0.5,0.4),plus1it=10, hardlimit = 14, optimize=TRUE,
                             MAP=TRUE),
                itermcmcpar=list(alpha=0.1,pbarrier=0.5,hardlimit=12,plus1it=10,accum=FALSE, MAP=TRUE,
                                 p=c(0.95,0.99,0.7,0.5,0.4))) {
  
  res$DAGpar<-DAGpar
  allalgonames<-c("pc", "ges", "mcmcMAP", "mcmcsample","iterativesample")
  if(!(algorithm%in%allalgonames)) {
    stop(paste("unknown algorithm ", algorithm, "\t", sep=""))
  }
  res$algorithm<-algorithm
  
  
  
}



simBN<-function(paramDAG,paramsimBN)

initSimBN<-function(DAGparam=list(nrep=10, DAGtype, n, d, par1, par2, lB=0.1, uB=1, labels=NULL, 
                                  ordered=TRUE, startseed=100, ss, type=c("cont", "bin", "mixed"),
                                  nbin=0, binmean=-1, cna=FALSE), 
                    DAGdatalist=NULL, DAGdatapath=NULL,
                    algorithms=c("pc", "ges", "orderMCMC.MAP.base", "orderMCMC.sample.base", 
                                 "orderMCMC.MAP.extended","orderMCMC.sample.extended",
                                 "orderMCMC.MAP.optimized", "orderMCMC.sample.optimized",
                                 "iterativeMCMC.MAP","iterativeMCMC.sample"), out=c("ROC","score","time"),
                    usralgoresult=NULL,usralgoname=NULL,saveDAGdata=TRUE,saveDAGdatapath=NULL,
                    saverespath=NULL, resname="simulationresults",
                    pcparam=list(alphas=c(0.01,0.05,0.1,0.2,0.3,0.4),test="gaussCItest"),
                    gesparam=list(lambdas=c(0.5,1,1.5,2,3,5,7), bootstrap=FALSE),
                    mcmcpar=list(alpha=0.1, p=c(0.95,0.99,0.7,0.5,0.4),plus1it=10, hardlimit = 14),
                    itermcmcpar=list(alpha=0.1,pbarrier=0.5,hardlimit=12,plus1it=10,accum=FALSE,
                                       p=c(0.95,0.99,0.7,0.5,0.4)))  {
  allalgonames<-c("pc", "ges", "orderMCMC.MAP.base", "orderMCMC.sample.base", 
                  "orderMCMC.MAP.extended","orderMCMC.sample.extended",
                  "orderMCMC.MAP.optimized", "orderMCMC.sample.optimized",
                  "iterativeMCMC.MAP","iterativeMCMC.sample")
 if(length(which(algorithms%in%allalgonames))<length(algorithms)) {
   stop(paste("unknown algorithm ", setdiff(algorithms,allalgonames), "\t", sep=""))
 }


 param<-list()
 param$binmean<-DAGparam$binmean
 param$resname<-resname
 if(!is.null(DAGdatalist)) {
   param$nrep<-length(DAGdatalist)
   param$DAGparam<-NULL
   param$DAGsource<-"listx"
 } else if (!is.null(DAGdatapath)) {
   fillist<-list.files(path=DAGdatapath)
   param$nrep<-length(filelist)
   param$DAGparam<-NULL
   param$DAGsource<-"pathx"
 } else {
   param$nrep<-DAGparam$nrep
   param$DAGparam<-DAGparam
   param$DAGsource<-"genx"
 }
 
 if(!is.null(usralgoresult)) {
   if(is.null(usralgoname)) {
     param$usralgoname<-"useralgo" 
   } else {
     param$usralgoname<-"useralgo" 
   }
   param$usralgoresult<-usralgoresult
 }
 param$mcmcnames<-c("orderMCMC.MAP.base", "orderMCMC.sample.base", 
                    "orderMCMC.MAP.extended","orderMCMC.sample.extended")
 param$itermcmcnames<-c("iterativeMCMC.MAP","iterativeMCMC.sample",
                        "orderMCMC.MAP.optimized", "orderMCMC.sample.optimized")
 param$algoparamlist<-vector()
 
 if(DAGparam$type=="cont") {
  scoretype<-"bge"
 } else if ((DAGparam$type=="bin")) {
   scoretype<-"bde"
 } else {
   scoretype<-"mixed"
 }
 whichmcmc<-which(algorithms%in%param$mcmcnames)
 if(length(whichmcmc)>0) {
   param$algoparamlist[whichmcmc]<-"mcmcpar"
   param$mcmcpar<-mcmcpar
   param$mcmcpar$score<-scoretype
 }
 whichitermcmc<-which(algorithms%in%param$itermcmcnames)
 if(length(whichitermcmc)>0) {
   param$algoparamlist[whichitermcmc]<-"itermcmcpar"
   param$itermcmcpar<-itermcmcpar
   param$itermcmcpar$score<-scoretype
 }
 whichpc<-which(algorithms%in%c("pc"))
 if(length(whichpc)>0){
   param$algoparamlist[whichpc]<-"pcparam"
   param$pcparam<-pcparam 
 }
 whichges<-which(algorithms%in%c("ges"))
 if(length(whichges)>0){
   param$algoparamlist[whichges]<-"gesparam"
   param$gesparam<-gesparam 
 }
 param$algorithms<-algorithms
 param$out<-out
 param$saveDAGdata<-saveDAGdata
 param$saveDAGdatapath<-saveDAGdatapath
 param$saverespath<-saverespath
 return(param)
}




#algosettings, GES lambda=1.5, 
#base MAP, plus1 MAP, optimized MAP, optimized sample, ges
#ROCdf, cl1, cl2, cl3
