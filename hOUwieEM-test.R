source("~/2020_hOUwie/hOUwieEM.R")
source("~/2020_hOUwie/Utils.R")
require(OUwie)
require(corHMM)
data(tworegime)
phy <- tree
phy$node.label <- NULL
# data <- trait
discrete_model <- "ARD"
continuous_model <- "OUMV"
rate.cat <- 1
shift.point <- 0.5
mserr <- "none"
dual = FALSE
collapse <- TRUE
root.station <- FALSE
get.root.theta <- FALSE
lb.disc = NULL
ub.disc = NULL
lb.cont = NULL
ub.cont = NULL
opts = NULL

#  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# BMS
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fit.cor <- equateStateMatPars(getRateCatMat(2), c(1,2))
fit.ou <- getOUParamStructure("BMS", "three.point", FALSE, FALSE, dim(fit.cor)[1])
pars = c(0.5, 1, 10, 8)  # mk, alpha, sigma, theta1, theta2
data.houwie <- generateData(phy, fit.cor, fit.ou, pars)
cols<-setNames(c("gold","red", "purple", "black"),
               c("1","2","3","4"))
plotDataSet(data.houwie); legend("bottomleft", legend = c("1A","2A","1B","2B"), pch=16, col = cols)

fitBM1 <- hOUwie.EM(phy, data.houwie$data, "ER", "BM1", rate.cat, shift.point, mserr, dual, collapse, root.station, get.root.theta, lb.disc, ub.disc, lb.cont, ub.cont, opts, 100)
fitOU1 <- hOUwie.EM(phy, data.houwie$data, "ER", "OU1", rate.cat, shift.point, mserr, dual, collapse, root.station, get.root.theta, lb.disc, ub.disc, lb.cont, ub.cont, opts, 100)
fitBMS <- hOUwie.EM(phy, data.houwie$data, "ER", "BMS", rate.cat, shift.point, mserr, dual, collapse, root.station, get.root.theta, lb.disc, ub.disc, lb.cont, ub.cont, opts, 100)
fitOUM <- hOUwie.EM(phy, data.houwie$data, "ER", "OUM", rate.cat, shift.point, mserr, dual, collapse, root.station, get.root.theta, lb.disc, ub.disc, lb.cont, ub.cont, opts, 100)

#  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# OUM
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fit.cor <- equateStateMatPars(getRateCatMat(2), c(1,2))
fit.ou <- getOUParamStructure("OUM", "three.point", FALSE, FALSE, dim(fit.cor)[1])
pars = c(0.5, 2, 2, 3, 8)  # mk, alpha, sigma, theta1, theta2
data.houwie <- generateData(phy, fit.cor, fit.ou, pars)
cols<-setNames(c("gold","red", "purple", "black"),
               c("1","2","3","4"))
plotDataSet(data.houwie); legend("bottomleft", legend = c("1A","2A","1B","2B"), pch=16, col = cols)

fitBM1 <- hOUwie.EM(phy, data.houwie$data, "ER", "BM1", rate.cat, shift.point, mserr, dual, collapse, root.station, get.root.theta, lb.disc, ub.disc, lb.cont, ub.cont, opts, 500)
fitOU1 <- hOUwie.EM(phy, data.houwie$data, "ER", "OU1", rate.cat, shift.point, mserr, dual, collapse, root.station, get.root.theta, lb.disc, ub.disc, lb.cont, ub.cont, opts, 500)
fitBMS <- hOUwie.EM(phy, data.houwie$data, "ER", "BMS", rate.cat, shift.point, mserr, dual, collapse, root.station, get.root.theta, lb.disc, ub.disc, lb.cont, ub.cont, opts, 500)
fitOUM <- hOUwie.EM(phy, data.houwie$data, "ER", "OUM", rate.cat, shift.point, mserr, dual, collapse, root.station, get.root.theta, lb.disc, ub.disc, lb.cont, ub.cont, opts, 500)


