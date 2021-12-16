makeDF<-function(simres) {
  alldfs<-list()
  if(!is.null(simres[[1]]$ROC)) {
    rocdf<-lapply(simres, function(x)x$ROC)
    rocdf<-do.call("rbind", rocdf)
    alldfs$roc<-rocdf
  }
  if(!is.null(simres[[1]]$time)) {
    timedf<-lapply(simres, function(x)x$time)
    timedf<-do.call("rbind", timedf)
    alldfs$time<-timedf
  }
  if(!is.null(simres[[1]]$score)) {
    scoredf<-lapply(simres, function(x)x$score)
    scoredf<-do.call("rbind", scoredf)
    alldfs$score<-scoredf
  }
  return(alldfs)
}
makeDFclust<-function(simres) {
  alldfs<-list()
  if(!is.null(simres[[1]]$roc)) {
    rocdf<-lapply(simres, function(x)x$ROC)
    rocdf<-do.call("rbind", rocdf)
    alldfs$roc<-rocdf
  }
  if(!is.null(simres[[1]]$lambda)) {
    lambdadf<-lapply(simres, function(x)x$lambda)
    lambdadf<-do.call("rbind", lambdadf)
    alldfs$lambda<-lambdadf
  }
  if(!is.null(simres[[1]]$accuracy)) {
    accuracydf<-lapply(simres, function(x)x$accuracy)
    accuracydf<-do.call("rbind", accuracydf)
    alldfs$accuracy<-accuracydf
  }
  return(alldfs)
}
