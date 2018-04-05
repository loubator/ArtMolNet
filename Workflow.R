# choose your paramaters
# in Welcome

# Cleaning the net
NETall<-NetBuilding(NET = NET,nblayers = 3,FamilyGene =  NULL)

# add weights
NETall<-addWeights(NETall)

#######
# MaBoSS simulation
# build .bnd
bndGenerator(NETall = NETall,NameProj = NameProj)

# built .cfg in multiple steps
timeMaxi<-FindMaxTime(NETall, NameProj, timeMaxi=100, NodeToEval=NULL)

# plot your results
PlotdensityNodesMBSS(nameSimIndiv)
