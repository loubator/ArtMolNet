addWeights<-function(NETall){
  # add Weights
  NETall[,"Weights"]<-runif(nrow(NETall),10,1000)
  NETall[,"Weights"]<-ifelse(NETall$interaction_directed_signed%in%"ACTIVATE",  NETall[,"Weights"],- NETall[,"Weights"])
  return(NETall)
}
