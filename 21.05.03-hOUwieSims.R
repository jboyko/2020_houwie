source("~/2020_hOUwie/hOUwie.R")
source("~/2020_hOUwie/Utils.R")
# source("/space_2/jamesboyko/2020_hOUwie/hOUwie.R")
# source("/space_2/jamesboyko/2020_hOUwie/Utils.R")

require(OUwie)
require(corHMM)
require(parallel)
require(phytools)

nTip <- 100
phy <- sim.bdtree(b = 1, d = 0.5, stop = "taxa", n = nTip) 
phy <- drop.extinct(phy)
phy$edge.length <- phy$edge.length/max(branching.times(phy))

fit.CIDOUM.cor <- equateStateMatPars(getRateCatMat(2), c(1,2))
fit.CIDOUM.cor <- getFullMat(list(fit.CIDOUM.cor, fit.CIDOUM.cor), fit.CIDOUM.cor)
fit.CIDOUM.cor <- equateStateMatPars(fit.CIDOUM.cor, c(1,2,3))
fit.CIDOUM.ou <- getOUParamStructure("OUM", "three.point", FALSE, FALSE, dim(fit.CIDOUM.cor)[1])
fit.CIDOUM.ou[3,] <- c(3,3,4,4)

pars = c(1, 10, 5, 5, 10)  # mk, alpha, sigma, theta1, theta2

data.houwie <- generateData(phy, fit.CIDOUM.cor, fit.CIDOUM.ou, pars)


cols<-setNames(c("gold","purple","red","black"),
               c("1","2","3","4"))
plotDataSet(data.houwie)



#Rsaves <- getRsaves("CIDOUM", nTip, nMap)[1:82]
Rsaves <- getRsaves("CIDOUM", nMap=100)[1:82]
DataTables <- AICTables <- list()
for(j in 1:length(Rsaves)){
  load(Rsaves[j])
  AICTables[[j]] <- getModelTable(out$model.fits$model.fits, "BIC")$dBIC
  # # AICTables[[j]] <- out$model.fits$model.table$dAICc
  # TruePars <- getMapWeightedPars(out$simulating.data, "Root")
  DataTables[[j]] <- out$simulating.data
}
AICTables <- do.call(rbind, AICTables)
colnames(AICTables) <-  c("CDBMS", "CDOUM", "CIDBMS", "CIDOUM", "HYBBMS", "HYBOUM")
AICTables[21,]
cols<-setNames(c("gold","purple","red","black"),
               c("1","2","3","4"))
plotDataSet(DataTables[[21]]); legend("bottomleft", legend = c("1A", "2A", "1B", "2B"), pch=16, col = cols)
load(Rsaves[21])
names(out$model.fits$model.fits) <- c("CDBMS", "CDOUM", "CIDBMS", "CIDOUM", "HYBBMS", "HYBOUM")
out$model.fits$model.fits
# $pars
# [1]  1  5  5  5 10
datum <- out$model.fits$model.fits[[3]]
hOUwie(phy = datum$phy, data = datum$data, rate.cat = 2, nSim = 100, index.cor = datum$index.cor, index.ou = datum$index.ou, p = c(0.1413891,275.79977,3.20614,4.64587))

datum <- out$model.fits$model.fits[[4]]
hOUwie(phy = datum$phy, data = datum$data, rate.cat = 2, nSim = 100, index.cor = datum$index.cor, index.ou = datum$index.ou, p = c(0.1877104,0.000657168,4.37957450,1.54974124,5.81677889))

datum <- out$model.fits$model.fits[[6]]
hOUwie(phy = datum$phy, data = datum$data, rate.cat = 2, nSim = 100, index.cor = datum$index.cor, index.ou = datum$index.ou, p = c(0.3131504, 0.5033651,3.8857972, 7.7850543,7.2854634, 12.8029686,12.8029686))

plotSimmap(out$simulating.data$simmap[[1]])

lam.table <- cbind(AICTables, lambda = 0)
for(i in 1:length(DataTables)){
  phy <- DataTables[[i]]$simmap
  data <- DataTables[[i]]$data[,3]
  names(data) <- DataTables[[i]]$data[,1]
  lambda.res <- fitContinuous(phy[[1]], data, model = "lambda")
  lam.table[i,7] <- lambda.res$opt$lambda
}





