#######
# welcome to ArtMolNet
# Artificial Molecular Network building, simulation and optimisation
# dedicated to anti-cancer treatments
# Loic Verlingue, MD PhDc, DITEP, U900

# load packages
library(reshape2)
#library(recount)
library(glmnet)
library(stringr)

# instal MaBoSS in the /tmp/ direction

#######
# create folders
getwd()
setwd("~/R/")
dir.create(file.path(getwd(), "model"), showWarnings = FALSE)
dir.create(file.path(getwd(), "tmp"), showWarnings = FALSE)

########
# load data net and FamilyGenes
list.files(paste(getwd(),"/data/",sep = ""))
OMNI<-read.csv2(paste(getwd(),"/data/OMNI.csv",sep = ""), stringsAsFactors = F)
ImmunGenes<-read.csv2(paste(getwd(),"/data/ImmunGenes.csv",sep = ""), stringsAsFactors = F)
ImmunGenes<-unique(ImmunGenes$gene)

#######
# name your project
NameProj<-"Immuno"

#######
# select your genes manually
# may propose a selection in a list
GENESman<-c("PDCD1", "CD274", "CTLA4") #"EGFR

#######
# select the treatment targets manually
treatmt<-c("CD274", "CTLA4") #"EGFR"
treatmt<-paste(treatmt,"_TTT",sep = "")

#######
# define the gene family you want to keep
ImmunGenes

#### 
# build NET
NET<-OMNI[OMNI$source_hgnc%in%union(GENESman,gsub("_TTT","",treatmt))|OMNI$target_hgnc%in%union(GENESman,gsub("_TTT","",treatmt)),]
# and go to NetBuilding()
