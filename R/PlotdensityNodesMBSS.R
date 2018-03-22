
# ====== Function to plot the density of the nodes along the simulation ====== #
PlotdensityNodesMBSS <- function(nameSimIndiv) {
  setwd(gsub("/tmp","",getwd()))
  setwd(paste(getwd(),"/tmp/",sep = ""))
  
  nameFolder <-list.files(getwd(),pattern = nameSimIndiv)
  nameFolder<-nameFolder[-grep("cfg",nameFolder)]
  nameFile <-paste(nameFolder, "_probtraj_table.csv", sep="")
  
  #  load table
  dataMB <- read.table(paste(getwd(),nameFolder,nameFile,sep="/"), header = TRUE)
  
  dataMB <- dataMB[,str_detect(colnames(dataMB), "^[Prob|Time]")]
  colnames(dataMB) <- str_replace(colnames(dataMB), "^Prob.", "")
  names <- colnames(dataMB)
  
  # ====== List with indices of each node in the dataMB data frame ====== #
  #  species <- model$namesSpecies
  species <- unique(unlist(str_split(colnames(dataMB),"\\.")))
  species <- species[!species%in%c("","Time","TH","nil") ]
  
  indexRO <- list()
  for (node in species){
    indexRO[[node]] <- c()
    for (aTransState in names){
      if ((str_detect(aTransState, node)==TRUE)){
        indexRO[[node]] <- c(indexRO[[node]], which(names==aTransState))
      }
    }
  }
  
  qtty <- c()
  for(node in species) {
    indices <- indexRO[[node]]
    if(is.null(dim(dataMB[,indices])) == TRUE){
      # ====== the node/species does not get any activation during the simulation ====== #
      qtty[[node]] <- dataMB[,indices]
    } else {
      # ====== sum of the probabilites for each time tick ====== #
      qtty[[node]] <- apply(dataMB[,indices], 1, sum)
    }
  }
  # plot
  timeSeq <- dataMB[,which(colnames(dataMB) == "Time")]
  pdf(paste(getwd(),"/",nameFolder,nameFile,".pdf", sep=""))
  par(mfrow=c(1,1))
  vecCol <- rainbow(length(species)-1)
  plot(timeSeq, unlist(qtty[[1]]), type="l", col=1, ylim = c(0,1),
       main = paste(nameFolder1,"_", sep=""))
  for (i in 2:length(species)){
    spec <- species[i]
    lines(timeSeq, unlist(qtty[[spec]]), col=vecCol[i-1])
    abline(origines[i],slopeValues[i], col=vecCol[i-1], lty="dotdash")
  }
  abline(v=c(timeSeq[length(timeSeq)-20]), col="black")
  legend("right", legend = names(qtty), col = c(1,vecCol), lwd = 1)
  dev.off()
}
