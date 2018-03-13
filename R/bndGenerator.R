bndGenerator <- function(NETall, NameProj, treatmt){
  
  
    ### Does the transcription from a given topology in model to logical rules
  ### Creates a new file
  
  # ====== Creation of a new file ====== #
  fileName <- paste(getwd(),"/tmp/",NameProj,".bnd", sep = "")
  write("", file = fileName)
  
  # ====== Extraction of informations from the topology ====== #
  nameReac<-as.character(apply(NETall[,1:3],1,function(x){
    paste(x,collapse = "_")
    }))
  
  nameSpecies <- union(NETall$source_hgnc,NETall$target_hgnc)
  
  #  nameSpecies <- c(model$namesSpecies)
#1  nameSpecies <- rownames(NetMat)
#  nameReac <- c(model$reacID)
#1  nameReac <- colnames(NetMat)
#  stim <- colnames(CNOlist@stimuli) # are inputs
  NETsif<-str_split(nameReac,"_",simplify = T)
  ST<-NETsif[,1]%in%NETsif[,ncol(NETsif)] # node without input
  stim <- NETsif[!ST,1]
  
  TS<-NETsif[,ncol(NETsif)]%in%NETsif[,1] # node without output
  
  # ====== Write the Boolean rules according to the topology tested ====== # 
  if (length(nameReac) == 0) {
    for (spec in setdiff(nameSpecies,stim)) {
      write(paste("Node", spec, "{", sep = " "), file = fileName, append = TRUE)
      write("  logic = (0);", file = fileName, append=TRUE)
      write(paste("  rate_up = @logic ? $u_",spec," : 0;",sep = ""),
            file = fileName, append = TRUE)
      write(paste("  rate_down = @logic ? 0 : $d_",spec,";",sep = ""),
            file = fileName, append = TRUE)
      write("}", file = fileName, append = TRUE)
    }
    for (spec in stim) {
      write(paste("Node", spec, "{", sep = " "), file = fileName, append = TRUE)
      write(paste("  logic = (",spec,");", sep=""), file = fileName, append=TRUE)
      write(paste("  rate_up = @logic ? $u_",spec," : 0;",sep = ""),
            file = fileName, append = TRUE)
      write(paste("  rate_down = @logic ? 0 : $d_",spec,";",sep = ""),
            file = fileName, append = TRUE)
      write("}", file = fileName, append = TRUE)
    }
  }
  
  # ====== Case when there is only one reaction in the topology ====== #
  else if (length(nameReac) == 1) {
#    target <- str_split(model$reacID[1], "=", simplify = TRUE)[2]
#    reg <- str_split(model$reacID[1], "=", simplify = TRUE)[1]
    
    target <- str_split(colnames(NetMat)[1], "_", simplify = TRUE)
    target<-target[length(target)]
    reg <- str_split(colnames(NetMat)[1], "_", simplify = TRUE)[1]

    write(paste("Node", target, "{", sep = " "), file = fileName, append = TRUE)
    write("  logic = ", file = fileName, append=TRUE)
    
    
    # ====== Write the Boolean rule for the regulated node ====== #
    # == target regulated by only one node == #
    if (str_detect(reg, "[+]") == FALSE){
      write(paste("(",reg,");", sep=""), file = fileName, append = TRUE)
      write(paste("  rate_up = @logic ? $u_",target," : 0;",sep = ""),
            file = fileName, append = TRUE)
      write(paste("  rate_down = @logic ? 0 : $d_",target,";",sep = ""),
            file = fileName, append = TRUE)
      write("}", file = fileName, append = TRUE)
      
    } else {
      # == Rule with an AND gate == #
      
      reg <- str_split(reg,"[+]", simplify = TRUE)
      firstElt = TRUE
      write("(", file = fileName, append = TRUE)
      for (aReg in reg) {
        ##### Add of AND condition in case of necessity #####
        if (firstElt != TRUE){
          write(" & ", file = fileName, append = TRUE)
        } else {
          firstElt = FALSE
        }
        write(aReg, file = fileName, append = TRUE)
      }
      write(");", file = fileName, append = TRUE)
      write(paste("  rate_up = @logic ? $u_",target," : 0;",sep = ""),
            file = fileName, append = TRUE)
      write(paste("  rate_down = @logic ? 0 : $d_",target,";",sep = ""),
            file = fileName, append = TRUE)
      write("}", file = fileName, append = TRUE)
    }
    
    
    # ====== Write the rule for all the other nodes : stimuli and not targeted nodes ====== #
    nameSpecies <- setdiff(nameSpecies,target)
    
    for (spec in stim){
      write(paste("Node", spec, "{", sep = " "), file = fileName, append = TRUE)
      write(paste("  logic = (", spec, ");", sep = ""), file = fileName, append=TRUE)
      write(paste("  rate_up = @logic ? $u_",spec," : 0;",sep = ""),
            file = fileName, append = TRUE)
      write(paste("  rate_down = @logic ? 0 : $d_",spec,";",sep = ""),
            file = fileName, append = TRUE)
      write("}", file = fileName, append = TRUE)
    }
    
    for (spec in setdiff(nameSpecies,union(target,stim))) {
      write(paste("Node", spec, "{", sep = " "), file = fileName, append = TRUE)
      write(paste("  logic = (0);", sep = ""), file = fileName, append=TRUE)
      write(paste("  rate_up = @logic ? $u_",spec," : 0;",sep = ""),
            file = fileName, append = TRUE)
      write(paste("  rate_down = @logic ? 0 : $d_",spec,";",sep = ""),
            file = fileName, append = TRUE)
      write("}", file = fileName, append = TRUE)
    }
  }
  
  
  # ====== Case when there are more than 1 reaction ====== #
  else {
    
    target<-nameSpecies[1]
    for (target in nameSpecies){
      write(paste("Node", target, "{", sep = " "), file = fileName, append = TRUE)
      write("  logic = ", file = fileName, append=TRUE)
      
      # ====== Build the boolean rules ====== #
      bRule = FALSE # == indicates if the node is regulated by another or not
      firstReg = TRUE
      
     SOU<-unlist(lapply(str_split(grep(target,nameReac,value = T),"_"),function(x)x[1]))
      TAR<-unlist(lapply(str_split(grep(target,nameReac,value = T),"_"),function(x)x[3]))
      SOU1<-SOU[-setdiff(grep(target,SOU),grep(target,TAR))] # keep auto-regulated nodes
      INFL<-unlist(lapply(str_split(grep(target,nameReac,value = T),"_"),function(x)x[2]))
      INFL<-INFL[-setdiff(grep(target,SOU),grep(target,TAR))]
      
      SOU1<-paste(ifelse(INFL=="INHIBIT","!",""),SOU1,sep = "")  
      SOU1
      
      if(length(grep("INHIB",INFL))>0){
      # order by inhibition or activation
      ORD <- order(INFL,decreasing = T)
      INFL<-INFL[ORD]
      SOU1 <- SOU1[ORD]
      
      # group inhibitors
      SOUI<-SOU1[grep("INHIB",INFL)]
      SOUA<-SOU1[grep("ACTIV",INFL)]
      SOUI<-paste(SOUI,collapse = " & ")
      if(length(SOUA)>0){
        TOT<-paste("(", SOUI," & ", SOUA ,")",sep = "")
        TOT<-paste(TOT,collapse = "|")
      } else if(length(SOUA)==0){
        TOT<-paste("(", SOUI,")",sep = "")
      }
      } else{
        TOT<-paste("(", paste(SOU1,collapse = " | "),")",sep = "")
      }
      
      write(TOT, file = fileName, append = TRUE)

      write(";", file = fileName, append = TRUE)
      write(paste("  rate_up = @logic ? $u_",target," : 0;",sep = ""), file = fileName, append = TRUE)
      write(paste("  rate_down = @logic ? 0 : $d_",target,";",sep = ""), file = fileName, append = TRUE)
      write("}", file = fileName, append = TRUE)
    }
  }
  # ====== End of the creation of the Boolean rules ====== #
  
  # ====== Combine the lines in order to write the file as it is required by MaBoSS ====== #
  myFile <- readLines(fileName)
  write("", file=fileName)
  for (wd in myFile){
    if (str_detect(wd, "^$") != TRUE){
      if (str_detect(wd, "^Node") == TRUE){ ### Beginning of the rule
        write(wd, file = fileName, append = TRUE)
      } else if (str_detect(wd, "^  logic") == TRUE){ ### Write the boolean rule
        lgcl <- wd
      } else if (str_detect(wd, "^  rate_up") == TRUE){ ### Write the rules regarding rates
        write(lgcl, file = fileName, append = TRUE)
        write(wd, file = fileName, append = TRUE)
      } else if (str_detect(wd, "^  rate_down") == TRUE){ ### kind of useless because written in the previous step
        write(wd, file = fileName, append = TRUE)
      } else if (str_detect(wd, "[}]") == TRUE){ ### End of the rule
        write(paste(wd,"\n", sep = ""), file = fileName, append = TRUE)
      } else {
        lgcl <- paste(lgcl, wd, sep="")
      }
    }
  }
}
