# A single fonction to perform the entire workflow

# choose your paramaters
# in Welcome

# Cleaning the net
NETall<-NetBuilding(NET = NET,nblayers = 3,FamilyGene = ImmunGenes)

#######
# MaBoSS simulation
# build .bnd
bndGenerator(NETall = NETall,NameProj = NameProj)

# built .cfg in multiple steps
timeMaxi<-FindMaxTime(NETall, NameProj, timeMaxi=100, NodeToEval=NULL)
