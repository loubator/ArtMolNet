DrawPath<-function(TotAttractors){
  
  # if state transition too big, reduce it
  # but can impare the detection of last species before attraction
  if(min(do.call("rbind",lapply(TotAttractors$LAttractors,dim))[,1])>20){
    TotAttractors$LAttractors<-lapply(TotAttractors$LAttractors,function(x){
      x<-rbind(x[seq(1,nrow(x)-20,length.out = 20),],x[(nrow(x)-20):nrow(x),])
      return(x)
    })
  }
  
  # true attractors
  # all unique possible states
  LAttractors<-lapply(TotAttractors$LAttractors,function(x){
    apply(x,1,function(y){paste(y,collapse = "")})
  })
  
  # give a number to each state
  TotiStateUnique<-unique(as.character(unlist(LAttractors)))
  TotiStateUnique<-data.frame(TotiStateUnique,NUM=seq(length(TotiStateUnique)))
  nrow(TotiStateUnique)
  
  # compute all the path
  # put all path the same number of steps by repeting the last state for those who finished early
  MAX<-max(unlist(lapply(TotAttractors$LAttractors,nrow)))
    
  Path<-lapply(TotAttractors$LAttractors,function(x){
    PathAt<-as.character(apply(x,1,function(y){paste(y,collapse = "")}))
    PathAt<-as.numeric(sapply(PathAt,function(h){which(TotiStateUnique$TotiStateUnique%in%h)}))
    Attract<-names(table(PathAt))[table(PathAt)%in%max(table(PathAt))]
    
    Path1<-c(PathAt,rep(PathAt[length(PathAt)], MAX-length(PathAt))) # repeat the last state
    GeneTransition=c(rownames(x),rep(tail(rownames(x),1),MAX-length(rownames(x)) ))
    list(Path1=Path1, GeneTransition=GeneTransition, Attract=Attract)
  })
  
  return(Path)
  
  #  GeneTransition<-lapply(Path,function(x)x[["GeneTransition"]])
  Path<-lapply(Path,function(x)x[["Path1"]])
  
  # store all path in data frame
  Path<-as.data.frame(Path)
  Path[,ncol(Path)+1]<-seq(nrow(Path))
  Path<-as.matrix(Path)
  
  ##### do the representation of the trajectories
  Latt<-length(list.files(getwd(),pattern =  "transition_graph"))
  NAME<-paste("transition_graph",Latt+1,sep = "")
  
  pdf(paste(getwd(),"/",NAME,".pdf",sep = ""))
  par(mar=c(5, 4, 4, 2) + 0.1)
  matplot(Path[,-ncol(Path)],type="l", 
          main="Trajectories",xlab="Time steps", ylab="State trasitions",axes=F)
  arrows(x0 = 0,y0 = -max(Path)/10,x1 = nrow(Path),y1 = -max(Path)/10,length = 0.2,xpd=T)
  axis(2,las=2)
  
  matplot(Path,type="l",xlim = c(13,18), 
          main="Trajectories",xlab="Time steps", ylab="State trasitions",axes=F)
  dev.off()
  
}
