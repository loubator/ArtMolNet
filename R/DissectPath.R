DissectPath<-function(Path,NTrans=3){
  
  # load function for 3D
  Plot3DAttract<-function(TopTrans){
    TopT<-TopTrans[TopTrans$OverallFreq>0,]
    
    kd<-with(TopT,kde2d(seq(nrow(TopT)),OverallFreq))
    
    zlim <- range(kd$y)
    zlen <- zlim[2] - zlim[1] + 1
    colorlut <- terrain.colors(10) # height color lookup table
    col <- cm.colors(20)[1 + round(19*(kd$z - min(kd$z))/diff(range(kd$z)))]
    
    open3d()
    surface3d(x = kd$x,y = kd$y,z = -kd$z, col=col)
    aspect3d(1,1,1)
    decorate3d()
    
    Latt<-length(list.files(getwd(),pattern =  "BestAttractLog"))
    NAME<-paste("BestAttractLog",Latt+1,sep = "")
    movie3d(spin3d(axis = c(0.5,0,1), rpm = 10), duration=6, movie = NAME,dir = getwd())
  }
  
  GeneTransition<-lapply(Path,function(x)x[["GeneTransition"]])
  Attract<-lapply(Path,function(x)as.numeric(x[["Attract"]]))
  Path<-lapply(Path,function(x)x[["Path1"]])
  
  # store all path in data frame
  Path<-as.data.frame(Path)
  Path[,ncol(Path)+1]<-seq(nrow(Path))
  Path<-as.matrix(Path)
  
  # compute list of vectors of all possible transitions from A to B
  ROWPath<-apply(Path,1,function(x)paste(x,collapse = ""))
  AllPath<-apply(Path,1,function(x){
    rbind(Path[grep(paste(x,collapse = ""),ROWPath)-1,-ncol(Path)],x[-ncol(Path)])
  })
  allTrans<-lapply(AllPath,function(x){
    apply(x,2,function(y){paste(y,collapse = " ")})
  })
  allTrans<-unlist(allTrans)
  allTrans<-allTrans[grep(" ",allTrans)]
  
  # do a table, estimate frequencies by transition
  allTrans<-data.frame(allTrans,WhichTotPath = names(allTrans))
  Tab<-table(allTrans)
  Tab<-data.frame(Tab)
  TabTrans<-table(allTrans$allTrans)
  Tab$Freq<-TabTrans[Tab$allTrans]
  
  Tab<-cbind(t(matrix(as.numeric(unlist(strsplit( as.character(Tab$allTrans)," " ) )),nrow = 2)),Tab$Freq)
  colnames(Tab)<-c("From","To","OverallFreq")
  
  # the frequent transitions and the attractors:
  
  # attractors sorted by frequency
  #AttractNames<-unique(unlist(Attract))
  AttractNames<-as.numeric(names(sort(table(unlist(Attract)),decreasing = T)))
  
  TopTrans<-Tab[Tab[,2]%in%AttractNames,]
  
  for(i in seq(NTrans)){
    TopTrans<-rbind(TopTrans, Tab[Tab[,2]%in%TopTrans[,1],]) # identify and select thoses just before
  }
  
  TopTrans<-TopTrans[!duplicated(paste(TopTrans[,1],TopTrans[,2],TopTrans[,3])),] #remove duplicates
  TopTrans<-TopTrans[!TopTrans[,1]==TopTrans[,2],] # remove the transition within a single path
  
  TopTrans<-as.data.frame(TopTrans)
  TopTrans$OverallFreq<-log(TopTrans$OverallFreq) # only none zero relations
  
  # plot
  Latt<-length(list.files(getwd(),pattern =  "CordPlotTransitionStates"))
  NAME<-paste("CordPlotTransitionStates",Latt+1,sep = "")
  pdf(paste(getwd(),"/",NAME,".pdf",sep = ""))
  par(mar=c(5, 4, 4, 2) + 0.1)
  chordDiagram(TopTrans,directional = 1,direction.type = "arrows",col=ifelse(TopTrans$To%in%AttractNames,3,4))
  dev.off()
  
  Plot3DAttract(TopTrans)
  
  #identify which species are involved in last paths
  GeneTransition<-as.data.frame(GeneTransition)
  TopGeneTransition<-apply(Path[,-ncol(Path)],2,function(x){
    NCOL<-grep(paste(x,collapse = ""),apply(Path,2,function(h)paste(h,collapse = "")))
    Attractor<-AttractNames[AttractNames%in%x][1]
    
    # positions of unique states and attractors, but some >1
    States<-which(!duplicated(x))
    PosAttract<-which(x%in%Attractor)[1]
    
    if(grep(PosAttract,States)>NTrans){
      Last<-gsub("_[0-9].*","",GeneTransition[States[seq(grep(PosAttract,States)-NTrans,grep(PosAttract,States))],NCOL])
      Attractor<-x[States[seq(grep(PosAttract,States)-NTrans,grep(PosAttract,States))]]
    } else {
      Last<-gsub("_[0-9].*","",GeneTransition[States[seq(0,grep(PosAttract,States))],NCOL])
      Attractor<-x[States[seq(0,grep(PosAttract,States))]]
    }
    #  Attractor<-Path[States[seq(grep(PosAttract,States)-3,grep(PosAttract,States))],NCOL]
    #  Last<-gsub("_[0-9].*","",GeneTransition[seq(grep(Attractor,x)[1]-3,grep(Attractor,x)[1]),NCOL])
    #  Attractor<-Path[seq(grep(Attractor,x)[1]-3,grep(Attractor,x)[1]),NCOL]
    list(Last=Last,Attractor=Attractor)
  })
  
  TopGeneTransition<-as.data.frame(do.call("rbind",TopGeneTransition))
  Transition<-as.data.frame(do.call("rbind",TopGeneTransition$Last))
  Transition<-cbind(Transition,as.data.frame(do.call("rbind",TopGeneTransition$Attractor)))
  colnames(Transition)<-paste("V",seq(ncol(Transition)),sep = "")
  
  SpeciesImportance<-head(sort(table(as.character(unlist(
    Transition[,which(!sapply(Transition[1,],is.numeric))]))),
    decreasing = T),4) # or NTrans?
  
  # careful, 2 best attractors maybe too low in some cases
  BestTransition<-Transition[Transition[,ncol(Transition)]%in%AttractNames[1:2],]
  SpeciesImportanceAttract<-head(sort(table(as.character(unlist(
    BestTransition[,which(!sapply(BestTransition[1,],is.numeric))]))),
    decreasing = T),4)
  
  return(list(Transition=Transition,SpeciesImportance=SpeciesImportance,SpeciesImportanceAttract=SpeciesImportanceAttract))
}
