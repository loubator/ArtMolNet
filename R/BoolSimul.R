
BoolSimul<-function(NETall, Boolean=F, ASYNC=T, NumiStates=20, GOF=NULL, LOF=NULL){
  
  if(!Boolean){Sigmoid=T}else{Sigmoid=F}
  
  if(Boolean){
    Update<-function(NETall,Ospecies,iStates){
    
    TabAND<-NETall[NETall$target_hgnc%in%Ospecies,]

    GREP<-TabAND$LogicRule%in%"OR"
    LOGOR<-iStates[TabAND[GREP,1]]
    if(length(LOGOR)>0){
      # Boolean
      LOGOR<-ifelse(TabAND[GREP,"interaction_directed_signed"]%in%"ACTIVATE",LOGOR,F)
      LOGOR<-paste("(", paste(LOGOR,collapse = "|"),")",sep = "")
    }
    
    # the species with whom is the AND connection is the other species with a AND annotation!
    GREP<-TabAND$LogicRule%in%"AND"
    
    LOGAND<-iStates[TabAND[grep("AND",TabAND$LogicRule),1]]
    if(length(LOGAND)>0){
      # Boolean
      LOGAND<-ifelse(TabAND[GREP,"interaction_directed_signed"]%in%"ACTIVATE",LOGAND,F)
      LOGAND<-paste("(", paste(LOGAND,collapse = "&"),")",sep = "")
    }
    
    Rule<-paste(LOGOR,ifelse(length(LOGOR)>0&length(LOGAND)>0,"&",""),LOGAND,sep = "")

    # Boolean
    return(eval(parse(text=Rule)))
  }
  } else if(Sigmoid){
    sigmoid<-function(z){1/(1+exp(-z))}
    Update<-function(NETall,Ospecies,iStates){
      
      TabAND<-NETall[NETall$target_hgnc%in%Ospecies,]
      TOT<-sum(TabAND[,"Weights"]*iStates[TabAND$source_hgnc])
      
      return(round(sigmoid(TOT)))
    }
  }
  
  # put connectors
  if(length(grep("LogicRule",colnames(NETall)))==0){
    NETall[,"LogicRule"]<-"OR" # begin with OR connectors only
    NETall[NETall$interaction_directed_signed%in%"INHIBIT","LogicRule"]<-"AND"
  }
  
  # initial states
  Species<-union(NETall$target_hgnc,NETall$target_hgnc)
  Nspecies<-length(Species)
  
  if(Nspecies<20){
    if(Boolean){
      PossiStates<-expand.grid(rep(list(c(T,F) ),Nspecies))
    } else if(Sigmoid){
      PossiStates<-expand.grid(rep(list(c(0,1)), Nspecies))
    }
   names(PossiStates)<-Species
  } else{
    # other way, by dividing the matrix
    if(Boolean){
      SubPossiStates<-matrix(sample(rep(c(T,F),Nspecies*NumiStates*10),Nspecies*NumiStates),ncol = Nspecies)
    } else if(Sigmoid){
      SubPossiStates<-matrix(sample(rep(c(0,1),Nspecies*NumiStates*10),Nspecies*NumiStates),ncol = Nspecies)
    }
     
    SubPossiStates<-SubPossiStates[!duplicated(apply(SubPossiStates,1,function(x)paste(x,collapse = ""))),]
    PossiStates<-SubPossiStates
   }
  
  # do the table of possible attractors and multiply it by the number of species to
  # perform the simulation with random initial species to update
  ListAttractors <- split(PossiStates, seq(nrow(PossiStates)))
  
  ListAttractors<-lapply(ListAttractors,function(x){
    names(x)<-Species
    return(x)
    })
  
  iStates<-ListAttractors[[9]]
  
########### simulations
    Simulation<-function(iStates, NETall, ASYNC=T){
      if(!ASYNC){SYNC=T}else{SYNC=F}
      Species<-union(NETall$target_hgnc,NETall$target_hgnc)
      Nspecies<-length(Species)
      
    TotiState<-as.data.frame(t(iStates),row.names = "")
    NewSpecies<-sample(names(iStates),1)
    
    iStates<-unlist(iStates)
      
    # random selection of the first one / select how many? tot, 1/2, 1?
    MonitorSpecies<-NewSpecies

      if(ASYNC){
        
        CHECK<-TRUE
        while(ifelse(nrow(TotiState)<(Nspecies*2),TRUE,CHECK)){
     
          for(Yspecies in NewSpecies){
            
            if(Yspecies%in%GOF){
              iStates[Yspecies]<-1
            } else if(Yspecies%in%LOF){
              iStates[Yspecies]<-0
            } else {
              iStates[Yspecies]<-Update(NETall,Ospecies = Yspecies,iStates)
            }
            
            TotiState<-rbind(TotiState,iStates)
            rownames(TotiState)[nrow(TotiState)]<-paste(Yspecies,length(grep(Yspecies,rownames(TotiState)))+1,sep = "_")
          }
          
        NewSpecies<-unique(NETall[NETall$source_hgnc%in%NewSpecies,"target_hgnc"])

        # if the last of monitorspecies are the same as NewSpecies, change species
        if(all(NewSpecies%in%tail(MonitorSpecies,100))){
          NewSpecies<-sample(Species,1)
        }

        MonitorSpecies<-c(MonitorSpecies,NewSpecies)

        # stoping rule
        ROWAttractors<-apply(TotiState,1,function(y)paste(y,collapse = ""))
        if(length(unique(tail(ROWAttractors,Nspecies*2)))==1){
          CHECK<-!ROWAttractors[length(ROWAttractors)]%in%unique(tail(ROWAttractors[-length(ROWAttractors)],Nspecies*2))
        }else{
          CHECK<-TRUE 
        }
        
        } # end of while loop
      } else if(SYNC){
        
      # that may be only for synchronous update
          CHECK<-TRUE
          while(ifelse(nrow(TotiState)<log2(Nspecies),TRUE,CHECK)){
            
        Updated<-sapply(Species,function(Ospecies){ Update(NETall,Ospecies,iStates) })
        
        iStates[Species]<-Updated
        
        if(length(c(LOF,GOF))>0){
          iStates[GOF]<-1
          iStates[LOF]<-0
        }
        TotiState<-rbind(TotiState,iStates)

      # stoping rule
      ROWAttractors<-apply(TotiState,1,function(y)paste(y,collapse = ""))
      
      #unique(tail(ROWAttractors,log2(Nspecies)))
      
      if(length(unique(tail(ROWAttractors,log2(Nspecies))))<3){
        CHECK<-!ROWAttractors[length(ROWAttractors)]%in%unique(tail(ROWAttractors[-length(ROWAttractors)],log2(Nspecies) ))
      }else {
        CHECK<-TRUE 
        }
      if(nrow(TotiState)>1000){break}
        } # #end of while loop
      } # end of SYNC
      
    return(TotiState)
  }  # end of simulation function
  
   # Initiate cluster
   no_cores <- detectCores() - 1
   cl <- makeCluster(no_cores)
   clusterExport(cl,c("Update","sigmoid","Simulation","iStates","NETall", "ASYNC"), envir = environment())
   
   # run simulations
   LAttractors<-parLapply(cl,ListAttractors,function(iStates){
     Simulation(iStates,NETall,ASYNC)
   })
   stopCluster(cl)
   
  # identify same initial states
  COORDiStates<-t(sapply(seq(nrow(PossiStates)),function(x)seq(x,length(LAttractors),nrow(PossiStates))))
  
  # add stable states to possible istate table
  TotAttractors<-lapply(LAttractors,function(y){
    y[nrow(y),]
  })

  TotAttractors<-data.frame(matrix(unlist(TotAttractors), nrow=length(TotAttractors), byrow=T))
  rownames(TotAttractors)<-names(LAttractors)
  colnames(TotAttractors)<-Species
  
  return(list(TotAttractors=TotAttractors, PossiStates=PossiStates, LAttractors=LAttractors))
  
}
