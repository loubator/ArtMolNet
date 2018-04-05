
PlotAttractors<-function(TotAttractors, PDF=F){
  PossiStates<-TotAttractors$PossiStates
  TotAttractors<-TotAttractors$TotAttractors
  
  ROWAttractors<-apply(TotAttractors,1,function(y)paste(y,collapse = ""))
  if(length(setdiff(as.numeric(table(ROWAttractors)),1))!=0){
    
    setwd(gsub("/tmp","",getwd()))
    setwd(paste(getwd(),"/tmp/",sep = ""))
    getwd()
    
    Latt<-length(list.files(getwd(),pattern =  "attractors"))
    NAME<-paste("attractors",Latt+1,sep = "")
    
    # the uniques
    TotAttractorsUnique<-TotAttractors[!duplicated(ROWAttractors),]
    rownames(TotAttractorsUnique)<-ROWAttractors[!duplicated(ROWAttractors)]
    
    # freq per uniques: carefull order!
    FreqAttractors<-table(apply(TotAttractors,1,function(y)paste(y,collapse = "")))
    
    # calculate the frequence and order
    FreqAttractors<-as.data.frame(FreqAttractors)
    
    rownames(FreqAttractors)<-FreqAttractors$Var1
    FreqAttractors[,"NumAttract"]<-seq(nrow(FreqAttractors))
    FreqAttractors[,"Percent"]<-100*FreqAttractors$Freq/sum(FreqAttractors$Freq)
    
    ## best attractors
    NamesAttract<-rownames(FreqAttractors[FreqAttractors$Freq>=quantile(FreqAttractors$Freq,0.85),])
    PositionBestAttract<-FreqAttractors[NamesAttract,]
    
    # order by frequency
    FreqAttractors<-FreqAttractors[rownames(TotAttractorsUnique),]
    TotAttractorsUnique<-TotAttractorsUnique[order(FreqAttractors$Freq,decreasing = F),]
    FreqAttractors<-FreqAttractors[order(FreqAttractors$Freq,decreasing = F),]
    
    # finish best attractor ordering
    BestAttract<-TotAttractorsUnique[apply(TotAttractorsUnique,1,function(y)paste(y,collapse = ""))%in%NamesAttract,]
    #BestAttract<-as.data.frame(BestAttract)
    #BestAttract[,"Percent"]<- as.numeric(FreqAttractors[rownames(BestAttract),"Percent"])
    dim(BestAttract)
    
    Percent<-as.numeric(FreqAttractors[rownames(BestAttract),"Percent"])
    if(PDF){
      pdf(paste(getwd(),"/",NAME,".pdf",sep = ""))
    }
  
    par(mar=c(7, 4, 10, 2) + 0.1)
    image(t(as.matrix(BestAttract)),axes=T,col=c(4,2))
    mtext(colnames(BestAttract),side = 3,line = 1, at = seq(par()$xaxp[1],par()$xaxp[2],length.out = ncol(BestAttract)),las=2,cex=0.5)
    #text(0.5,seq(0,nrow(BestAttract))/(nrow(BestAttract)-1 ), 
    text(0.5,seq(par()$yaxp[1],par()$yaxp[2],length.out = nrow(BestAttract))/diff(par()$yaxp[1:2]), 
         paste(round(Percent,3),"%") ) #FreqAttractors[rownames(BestAttract),"Percent"]
    par(new=T)
    par(mar=c(1, 4, 4, 2) + 0.1)
    plot(1,1,type="n",axes=F,xlab="",ylab="")
    legend("bottomleft",legend = c("Node = 0","Node = 1"),fill = c(4,2), title ="Activity")
    
    # map with arrows
    FreqAttractorsArr<-FreqAttractors[apply(TotAttractors,1,function(y)paste(y,collapse = "")),]
    
    par(mar=c(5, 4, 4, 2) + 0.1)
    plot(1,1,type='n',xlim=c(0,nrow(TotAttractorsUnique) ),ylim=c(-1,5),xlab="",ylab="Best attractors",axes=F)
    X0<-(seq(nrow(PossiStates))*nrow(TotAttractorsUnique))/(nrow(PossiStates))
    arrows(x0 = X0, x1 = FreqAttractorsArr$NumAttract,y0 = 4,y1 = 3, length = 0.1,lwd=0.4)
    
    PositionBestAttract<-FreqAttractorsArr[NamesAttract,]
    
    if(length(c(grep("TRUE", rownames(PositionBestAttract)),
                grep("FALSE",rownames(PositionBestAttract))))>0){
      BestAttract1<-gsub("FALSE",0,rownames(PositionBestAttract))
      BestAttract1<-gsub("TRUE",1,BestAttract1)
      BestAttract1<-paste(BestAttract1,collapse = "")
      text(PositionBestAttract$NumAttract,2.9,labels = BestAttract1,pos = 2,srt=90,cex=0.3)
      
    } else{
      text(PositionBestAttract$NumAttract,2.9,labels = rownames(PositionBestAttract),pos = 2,srt=90,cex=0.3)
    }
    
    text(nrow(TotAttractorsUnique)/2,4.5,labels = "Randomly selected initial nodes")
    
    # plot density not scaled.... to be improved
    par(new=T)
    par(mar=c(5, 4, 4, 2) + 0.1)
    DENS<-density(FreqAttractorsArr$Freq)
    plot(rev(DENS$x),DENS$y,type="l",
         xlim=c(quantile(DENS$x,0.05),max(DENS$x)),ylim=c(min(DENS$y),max(DENS$y)*4),
         axes=F,main="",ylab="",xlab="") #,xlim=c(0,max(FreqAttractorsArr$Freq))
    mtext(text = "density", side = 2,line = 3,at = mean(density(FreqAttractorsArr$NumAttract)$y))
    if(PDF){
      dev.off()
    }
    return(BestAttract)
  }
}
