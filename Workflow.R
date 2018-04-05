# choose your paramaters
# in Welcome

# Cleaning the net
NETall<-NetBuilding(NET = NET,nblayers = 3,FamilyGene =  NULL)

# add weights
NETall<-addWeights(NETall)

# do simulations to defined attractors
# speed up? Do proba on it # simplify to optimize ?
# add rates to calculate speed of transition
# add random number selection in the update <-> MaboSS ?
TotAttractors<-BoolSimul(NETall,ASYNC = T, Boolean = F, NumiStates = 100)

# plot the attractors (in /tmp/attractor.pdf) and write table with best attractors
BestAttract<-PlotAttractors(TotAttractors,PDF = T)

# plot species activities
SpeciesActivity<-plotProbaState(Species = union(NETall$source_hgnc,NETall$target_hgnc),
               TotAttractors = TotAttractors, Stochastic = F) 


#######
# MaBoSS simulation
# build .bnd
bndGenerator(NETall = NETall,NameProj = NameProj)

# built .cfg in multiple steps
timeMaxi<-FindMaxTime(NETall, NameProj, timeMaxi=100, NodeToEval=NULL)

# plot your results
PlotdensityNodesMBSS(nameSimIndiv)
