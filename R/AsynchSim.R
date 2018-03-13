NETall<-data.frame(source_hgnc=c("G0","G1","G2","G0"),
           interaction_directed_signed=c("ACTIVATE","INHIBIT","ACTIVATE","ACTIVATE"),
           target_hgnc=c("G1","G2","G0","G2"))

#DUPS<-NETall$target_hgnc[duplicated(NETall$target_hgnc)]
#NETall[NETall$target_hgnc%in%DUPS,"LogicRule"]<-"OR"
NETall[,"LogicRule"]<-"OR" # begin with OR connectors only
NETall[NETall$interaction_directed_signed%in%"INHIBIT","LogicRule"]<-"AND"

Species<-unique(NETall$source_hgnc)
Nspecies<-length(Species)

# asynch

#iStates<-round(runif(Nspecies, 0,1))
#names(iStates)<-Species
#iStates
PossiStates<-expand.grid(rep(list(0:1), length(Species)))
names(PossiStates)<-Species

###################
  # do all in 1 loop
  # order the vector

iStates<-as.vector(PossiStates[2,])
iStates
ListAttractors<-list()
for(i in seq(nrow(PossiStates))){
  iStates<-as.vector(PossiStates[i,])
  TotiState<-(as.data.frame(iStates))
  while(!identical(paste(TotiState[nrow(TotiState),],collapse = ""),
                   paste(TotiState[nrow(TotiState)-2,],collapse = ""))|nrow(TotiState)<2){
    RandSpecie=sample(Species,1) # first one
    Species<-relevel(as.factor(Species),as.character(RandSpecie))
    Ospecies="G0"
    for(Ospecies in levels(Species)){
      SourceSpecies<-NETall[NETall$target_hgnc%in%Ospecies,"source_hgnc"]
      Interact<-NETall[NETall$target_hgnc%in%Ospecies,"interaction_directed_signed"]
      data.frame(Ospecies,
                 SourceSpecies,
                 Interact,iStates[SourceSpecies])
      
      # with AND connectors
      TabAND<-NETall[NETall$target_hgnc%in%Ospecies,]
      TabAND<-TabAND[iStates[SourceSpecies]!=0,]
      if(nrow(TabAND)>1){
        TabAND<-TabAND[order(TabAND$LogicRule,decreasing = T),] # put the AND last
        FUNC<-ifelse(TabAND$interaction_directed_signed=="INHIBIT","!","") #build the rule: inhibition
        Rule<-paste(ifelse(TabAND$LogicRule%in%"AND","&","|"),FUNC,sep = "") #build the rule: connectors
        Rule<-paste(Rule,ifelse(iStates[SourceSpecies]==0,FALSE,TRUE),sep = "",collapse = "")# activity of the node
        Rule<-substring(Rule,2,nchar(Rule))
        eval(parse(text=Rule))
        iStates[Ospecies]<-ifelse(eval(parse(text=Rule)),1,0)
      } else {
        FUNC<-ifelse(Interact=="INHIBIT",FALSE,TRUE)[iStates[SourceSpecies]!=0]
        if(length(FUNC)>0){
          iStates[Ospecies]<-ifelse(all(FUNC),1,0)
        } else {
          iStates[Ospecies]<-0
        }
      }
      
      #
      TotiState<-rbind(TotiState,iStates)
    }
    if(nrow(TotiState)>30){break}
  }
  ListAttractors[[i]]<-  TotiState
  
}
ListAttractors
