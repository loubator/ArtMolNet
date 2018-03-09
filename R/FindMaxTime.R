
FindMaxTime<-function(NETall, NameProj, timeMaxi=100, NodeToEval=NULL){
  setwd(gsub("/tmp","",getwd()))
  setwd(paste(getwd(),"/tmp/",sep=""))
  getwd()
  ### Write the cfg file
  ### find minimal max_time for steady states 
  i=0
  MAXpre=10
  while(MAXpre<timeMaxi){
    i=i+1
    print(paste("N simulation =",i))
  # time parameters
  nameSimIndiv<-paste("findMaxTime",i,sep = "")
  dest_file <- paste(getwd(),"/tmp/" ,nameSimIndiv,"_",NameProj,".cfg", sep="")
  
#  all_spec<-rownames(NetMat)
  all_spec<-union(NETall$source_hgnc,NETall$target_hgnc)
  # rates
  upRate<-paste("$u_", all_spec,"=1;",sep = "")
  write(upRate, file = dest_file, append = TRUE)
  dnRate<-paste("$d_", all_spec,"=1;",sep = "")
  write(dnRate, file = dest_file, append = TRUE)
  
  # istates random = dont need to write it
  
  #isinternal is which node you want to look at
  # may be selected on relevance, or number of interactions
  
  if(is.null(NodeToEval)){
    NodeToEval<-tail(sort(table(c(as.character(NETall$source_hgnc),as.character(NETall$target_hgnc)))))
    NodeToEval<-names(NodeToEval)
  }
  
  isinternal<- paste(NodeToEval, ".is_internal=", 0, ";", sep="")
  write(isinternal, file = dest_file, append = TRUE)
  isinternal<- paste(all_spec[!all_spec%in%isinternal], ".is_internal=", 1, ";", sep="")
  write(isinternal, file = dest_file, append = TRUE)
  
  # ====== Max time of the simulation ====== #
  # choose duration for the 1st simulation
  if(i==1) {timeMaxi=300}
  #
  runTime <- paste("max_time = ", timeMaxi+0.10, ";", sep="")
  write(runTime, file = dest_file, append = TRUE) #1
#  paste("time_tick =", 0.01,";",sep = "")
  
  write( paste("time_tick = 0.1;
discrete_time = 0;
use_physrandgen = FALSE;
seed_pseudorandom = 100;
sample_count = 1000;
thread_count = 4;
statdist_traj_count = 100;
statdist_cluster_threshold = 0.9;"),file = dest_file, append = TRUE)
    
    # run MaBoSS
    # check dir for .bnd
    

system(paste("./MBSS_FormatTable.pl",NameProj,".bnd ",dest_file, sep=""))
    
#  timeMaxi <- testOnSteadyState(i, NameProj, CNOlist, modelCut, timeMaxi)
    SlopeRSE <- testOnSteadyState(nameSimIndiv, timeMaxi)
    # returnes slopeValues & rseValues
setwd(gsub("/tmp","",getwd()))

    slopeValues<-SlopeRSE$slopeValues
    rseValues<-SloperRSE$rseValues
    
    MAX<-SlopeRSE$slopeValues[SlopeRSE$slopeValues==max(SlopeRSE$slopeValues)]
    #  if (addTime == TRUE) {
    MAXpre<-timeMaxi
    
    timeMaxi <- if(MAX>1) timeMaxi+100 else if(1>=MAX&MAX>0.0001) timeMaxi+10 else timeMaxi
    
    #  }
    if(i==10){
      print("stopped at N? 10")
      file.copy(dest_file, gsub(nameSimIndiv,"Base",dest_file))
      file.remove(dest_file)
      break
    } else if(MAXpre==timeMaxi){
      file.copy(dest_file, gsub(nameSimIndiv,"Base",dest_file))
      file.remove(dest_file)
    } else {
      file.remove(dest_file)
   }
  }
   return(timeMaxi)
}
