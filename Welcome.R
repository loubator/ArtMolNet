
#######
# welcome to AMoNet
# Artificial Molecular Network creator, simulation and optimisation
# dedicated to anti-cancer treatments
# Loic Verlingue, MD PhDc, DITEP, U900

# load pack
library(reshape2)
#library(recount)
library(glmnet)
library(stringr)
library(mgcv)
library(rgl)
library(circlize)
library(parallel)
library(MASS)


# + install MaBoSS and MBSS_FormatTable.pl in the R/ direction

########
# set directory
setwd("C:/Users/L_VERLINGUE/Desktop/ModelK/Rpack/ArtMolNet/")
getwd()

#######
# create folders
dir.create(file.path(getwd(), "model"), showWarnings = FALSE)
dir.create(file.path(getwd(), "tmp"), showWarnings = FALSE)

########
# load data net and FamilyGenes
list.files(paste(getwd(),"/data/",sep = ""))
OMNI<-read.csv2(paste(getwd(),"/data/OMNI.csv",sep = ""), stringsAsFactors = F)

#######
# name your project
NameProj<-"Angio"

#######
# select your genes manually
# may propose a selection in a list
GENESman<-c("VEGFA", "EGFR")

#######
# select the treatment targets manually
treatmt<-c("VEGFA", "EGFR")
treatmt<-paste(treatmt,"_TTT",sep = "")

#### 
# build NET
NET<-OMNI[OMNI$source_hgnc%in%union(GENESman,gsub("_TTT","",treatmt))|OMNI$target_hgnc%in%union(GENESman,gsub("_TTT","",treatmt)),]
# and go to NetBuilding()
