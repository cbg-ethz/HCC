#Run benchmarking of structure learning algorithms with multiple replicates

#This function can be used to compare the performance of algorithms for structure learning of Bayesian networks.
#The user can specify paramters of the simulation: DAG characteristics, sample size, size of the network etc. The performance
#can be evaluated according to different graph destance metrics as well as runtime.
#make bootstrap, add algos, publish!

#'@param param an object of class simparam which contains all the
#'@param parallel logical, if TRUE the replicates are parallelized
#'@param nnodes integer, used when the parallel equals TRUE, indicates the number of nodes used for parrallelization

#'@return an object of class BNsimulation
#'@export
simBN<-function(param, parallel=TRUE, nnodes=4){

  if(parallel) {

    loopsize<-param$nrep%/%nnodes
    afterloop<-param$nrep%%nnodes
    if(afterloop>0) {loopsize<-loopsize+1}
    dfs<-list()

    for(i in 1:loopsize) {
      if(i==loopsize && afterloop>0) {
        ind<-c(1:afterloop+(i-1)*nnodes)
        cl <- makeCluster(afterloop)
      } else {
        ind<-c(1:nnodes+(i-1)*nnodes)
        cl <- makeCluster(nnodes)
      }

      outputClApply <- clusterApply(cl, ind, simBNcore, param=param)
      stopCluster(cl)

      dfs[[i]]<-makeDF(outputClApply)

      if(!is.null(param$saverespath)) {
        saveRDS(dfs,file=paste(param$saverespath, "/temp/", param$resname,i, ".rds", sep=""))
      }
      alldfs<-makeDF(dfs)
    }
  } else {
    simx<-list()
    for(i in 1:param$nrep) {
      simx[[i]]<-simBNcore(i,param)
    }
    simx<-makeDF(simx)
  }
  return(simx)
}
simBNcore<-function(ind,param){
  DAGdata<-getDAGdata(ind, param)
  DAGinfo<-getDAGinfo(DAGdata$DAG, param, ind, nrow(DAGdata$data))
  rocx<-NULL
  scorex<-NULL
  timex<-NULL
  countalgo<-1
  for(i in param$algorithms) {
    algoparams<-param[[param$algoparamlist[countalgo]]]
    if(param$DAGparam$type=="mixed") {
      resx<-learnBN.mixed(DAGdata$data, i, algoparams, param$DAGparam$n, param$DAGparam$nbin)
    } else {
      resx<-learnBN(DAGdata$data, i, algoparams)
    }

    for(j in param$out) {
      switch(j,
             "ROC"={
               rocx<-rbind(rocx,compareBNdf(resx,DAGdata$DAG,DAGinfo))
             },
             "score"={
               scorex<-rbind(scorex,getscoredf(resx,DAGdata$DAG, DAGinfo))
             },
             "time"={
               timex<-rbind(timex,gettimedf(resx,DAGdata$DAG, DAGinfo))
             }
      )
    }
    countalgo<-countalgo+1

  }

  rocx<-cbind(rocx, DAGinfo[rep(1,nrow(rocx)),])

  if(!is.null(scorex)) {
    scorex<-cbind(scorex, DAGinfo[rep(1,nrow(scorex)),])
  }
  if(!is.null(timex)) {
    timex<-cbind(timex, DAGinfo[rep(1,nrow(timex)),])
  }

  res<-list()
  res$ROC<-rocx
  res$time<-timex
  res$score<-scorex

  return(res)
}
simiterBN<-function(nrep, DAGparam=list, algolabels=c("mcmc1", "mcmc2"),
                          itermcmcpar1=list(alpha=0.1,pbarrier=0.5,hardlimit=12,plus1it=10,accum=FALSE,
                                            p=c(0.95,0.99,0.7,0.5,0.4), searchspace="hybrid"),
                          itermcmcpar2=list(alpha=0.1,pbarrier=0.5,hardlimit=12,plus1it=10,accum=FALSE,
                                            p=c(0.95,0.99,0.7,0.5,0.4), searchspace="hybrid")) {
  if(parallel) {

    loopsize<-param$nrep%/%nnodes
    afterloop<-param$nrep%%nnodes
    if(afterloop>0) {loopsize<-loopsize+1}
    dfs<-list()

    for(i in 1:loopsize) {
      if(i==loopsize && afterloop>0) {
        ind<-c(1:afterloop+(i-1)*nnodes)
        cl <- makeCluster(afterloop)
      } else {
        ind<-c(1:nnodes+(i-1)*nnodes)
        cl <- makeCluster(nnodes)
      }

      outputClApply <- clusterApply(cl, ind, SampMapcore, param=list(itermcmcpar1, itermcmcpar2))
      stopCluster(cl)

      dfs[[i]]<-makeDF(outputClApply)

      if(!is.null(param$saverespath)) {
        saveRDS(dfs,file=paste(param$saverespath, "/temp/", param$resname,i, ".rds", sep=""))
      }
      alldfs<-makeDF(dfs)
    }
  } else {
    simx<-list()
    for(i in 1:param$nrep) {
      simx[[i]]<-simBNcore(i,param)
    }
    simx<-makeDF(simx)
  }
  return(simx)

}
simiterBNcore<-function(ind,param){

  DAGdata<-getDAGdata(ind, param)
  DAGinfo<-getDAGinfo(DAGdata$DAG, param, ind, nrow(DAGdata$data))
  rocx<-NULL
  scorex<-NULL
  timex<-NULL
  countalgo<-1

  for(i in 1:2) {
    algoparams<-param[[i]]
    resx<-learnBN(DAGdata$data, i, algoparams)
    for(j in param$out) {
      switch(j,
             "ROC"={
               rocx<-rbind(rocx,compareBNdf(resx,DAGdata$DAG,DAGinfo))
             },
             "score"={
               scorex<-rbind(scorex,getscoredf(resx,DAGdata$DAG, DAGinfo))
             },
             "time"={
               timex<-rbind(timex,gettimedf(resx,DAGdata$DAG, DAGinfo))
             }
      )
    }
    countalgo<-countalgo+1

  }

}





