
NetBuilding <- function(NET,nblayers=3,FamilyGene=ImmunGenes){
  ################
  # simplified
  # keep all genes at each layers and then filter those with feedbacks
  
  NETall<-NET
  i=0
  Nspecies<-c(10,20)
  FUNC=TRUE
  while( length(unique(tail(round(Nspecies,-1),n = 2)))!=1 )  {
    
    i=i+1
    #  A<-NETall$source_hgnc[3]
    if(FUNC){
      LAYER<-toupper(letters[seq(nblayers+1)])
      N<-list()
      N[[LAYER[1]]]<-unique(NETall$source_hgnc)
      y=1
      for(y in seq(nblayers) ){
        N[[y+1]]<-OMNI$target_hgnc[OMNI$source_hgnc%in% N[[y]]]
        N[[y+1]]<-N[[y+1]][!N[[y+1]]%in%NETall$source_hgnc]
        FB<-OMNI[OMNI$source_hgnc%in%N[[y+1]] & OMNI$target_hgnc%in%union(NETall$source_hgnc,NETall$target_hgnc),]
        dim(FB)
        if(nrow(FB)>0){
          NETall<-rbind(NETall,FB)
        }
      }
      NETall<-NETall[NETall$source_hgnc%in%ImmunGenes|NETall$target_hgnc%in%ImmunGenes,] # keep immunogenes only
      NETall<-rbind(NETall,NET) # restaure initial genes
      
      
    } else {
      for(A in unique(NETall$source_hgnc)){
        names<-toupper(letters[seq(nblayers+1)])
        N<-list()
        N[[names[1]]]<-A
        # y=1
        for(y in seq(nblayers)){
          N[[names[y+1]]]<-OMNI$target_hgnc[OMNI$source_hgnc%in% N[[names[y]]]]
          N[[names[y+1]]]<-N[[names[y+1]]][!N[[names[y+1]]]%in%NETall$source_hgnc]
          #B<-OMNI[OMNI$source_hgnc%in%N[[names[y+1]]],]
          FB<-OMNI[OMNI$source_hgnc%in%N[[names[y+1]]] & OMNI$target_hgnc%in%union(NETall$source_hgnc,NETall$target_hgnc),]
          if(nrow(FB)>0){
            NETall<-rbind(NETall,FB)
          }
        }
        
        NETall<-NETall[NETall$source_hgnc%in%ImmunGenes|NETall$target_hgnc%in%ImmunGenes,] # keep immunogenes only
        NETall<-rbind(NETall,NET) # restaure initial genes
      }
    }
    NETall<-NETall[!duplicated(paste(NETall$source_hgnc,NETall$target_hgnc)),] # no dups
    
    Nspecies[i]<-length(unique(NETall$source_hgnc))[1]
    print(paste( "Round",i,"=",length(unique(NETall$source_hgnc))[1], "species"))
    if(i==10){break}
  }
  
  #write.csv2(NETall,"C:/Users/L_VERLINGUE/Desktop/ModelK/Roma/networkexpansion3_PDL1_PD1_CTLA4_flash.csv")
  #write.csv2(NETall,"C:/Users/L_VERLINGUE/Desktop/ModelK/Roma/networkexpansion3_PDL1_PD1_CTLA4_flash_restricted.csv")
  
  Nspecies<-c(1,2)
  while( length(unique(tail(round(Nspecies),n = 2)))!=1 )  {
    i=i+1
    
    # then do a restriciton to fb edges
    TS<-NETall$target_hgnc%in%NETall$source_hgnc # node without output
    ST<-NETall$source_hgnc%in%NETall$target_hgnc # node without input
    table(TS&ST)
    
    while(length(unique(TS&ST))!=1){
      TS<-NETall$target_hgnc%in%NETall$source_hgnc # node without output
      ST<-NETall$source_hgnc%in%NETall$target_hgnc # node without input
      #data.frame(NETall[,c(1:3)],TS,ST)
      NETall<-NETall[TS&ST,]
      # dim(NETall)
      if(dim(NETall)[1]<10) break
    }
    dim(NETall)
    print(paste("post restriction to FB N=",length(unique(NETall$source_hgnc)) ))
    
    # are nodes only "passengers": discard them
    TAB<-sort(table(c(as.character(NETall$source_hgnc),as.character(NETall$target_hgnc))))
    #TAB[names(TAB)%in%unique(NET$source_hgnc)]
    PASSENGERS<-names(TAB[TAB==2])
    
    # function to store the passenger inputs and outputs
    x<-PASSENGERS[3]
    PASSV<-sapply(PASSENGERS,function(x){
      TAR<-NETall[NETall$source_hgnc%in%x,"target_hgnc"]
      ITAR<-NETall[NETall$source_hgnc%in%x,"interaction_directed_signed"]
      SOU<-NETall[NETall$target_hgnc%in%x,"source_hgnc"]
      ISOU<-NETall[NETall$target_hgnc%in%x,"interaction_directed_signed"]
      
      V<-c(ITAR,ISOU)
      names(V)<-c(TAR,SOU)
      return(list(V))
    })
    
    KEEP<-as.numeric(which((lapply(PASSV,length)==2)))
    PASSV<-PASSV[KEEP]
    rm(KEEP)
    
    # remove passengers from the network table
    #table((NETall$source_hgnc%in%names(PASSV)|NETall$target_hgnc%in%names(PASSV)))
    NETall<-NETall[!(NETall$source_hgnc%in%names(PASSV)|NETall$target_hgnc%in%names(PASSV)),]
    dim(NETall)
    #x<-PASSV[[1]]
    if(length(PASSV)>0){
      for(y in seq(length(PASSV))){
        x<-PASSV[[y]]
        if( length(unique(paste(names(x),x,sep = "_")))==1 ){
          next
        }else if(length(unique(x))==1){ # if only activ or 2 inhib
          NETall[nrow(NETall)+1,"target_hgnc"]<-names(x)[1]
          NETall[nrow(NETall),"source_hgnc"]<-names(x)[2]
          NETall[nrow(NETall),"interaction_directed_signed"]<-"ACTIVATE"
        } else {
          NETall[nrow(NETall)+1,"target_hgnc"]<-names(x)[1]
          NETall[nrow(NETall),"source_hgnc"]<-names(x)[2]
          NETall[nrow(NETall),"interaction_directed_signed"]<-"INHIBIT"
        }
      }
      
    }
    dim(NETall)
    print(paste("post passenger filtering N =",length(unique(NETall$source_hgnc)) ))
    
    NETall<-NETall[!duplicated(paste(NETall$source_hgnc,NETall$target_hgnc)),] # no dups
    
    Nspecies[i]<-length(unique(NETall$source_hgnc))[1]
    print(paste( "Round",i,"=",length(unique(NETall$source_hgnc))[1], "species"))
    if(i==10){break}
  }
  
  #write.csv2(NETall,"C:/Users/L_VERLINGUE/Desktop/ModelK/Roma/networkexpansion3_PDL1_PD1_CTLA4_all_flash_multi_reduced2.csv")
  write.csv2(NETall,paste(getwd(),"/model/NETall.csv",sep = ""))
  dim(NETall)
  
  return(NETall)
  
}
