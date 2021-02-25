require(corHMM)
require(OUwie)
require(parallel)
require(geiger)
require(proftools)

phy <- sim.bdtree(b = 1, d = 0.5, stop = "taxa", n = 250)
phy <- drop.extinct(phy)
phy$edge.length <- phy$edge.length/max(branching.times(phy))
root.p = c(1, 0)
p.mk <- c(2, 2)
alpha = c(1, 1)
sig2= c(0.5, 0.5)
theta = c(3, 8)
theta0 = 3
rate.cat = 1
model.cor = "ER"
model.ou = "OUM"

singleRun <- function(phy, p.mk, root.p, alpha, sig2, theta, theta0, rate.cat, model.cor, model.ou){
  Q = matrix(c(-p.mk[1],p.mk[2],p.mk[1],-p.mk[2]), 2, 2)
  data <- OUwie:::hOUwie.sim(phy, Q, root.p, alpha, sig2, theta0, theta)
  no.changes <- sum(unlist(lapply(data[[2]]$maps, function(x) length(x)))-1)
  data <- data[[1]]
  hOUwie.dat <- OUwie:::organizeHOUwieDat(data)
  nObs <- length(hOUwie.dat$ObservedTraits)
  model.set.final <- corHMM:::rate.cat.set.corHMM.JDB(phy=phy,data=hOUwie.dat$data.cor,rate.cat=rate.cat, ntraits = nObs, model = model.cor)
  phy <- reorder(phy, "pruningwise")
  index.ou <- OUwie:::getParamStructure(model.ou, "three.point", FALSE, FALSE, dim(model.set.final$Q)[2])

  p = c(p.mk[1], alpha[1], sig2[1], theta[1], theta[2])
  o1 <- OUwie:::hOUwie.dev(p = log(p), phy = phy, data.cor = OUwie:::organizeHOUwieDat(data)$data.cor, data.ou = OUwie:::organizeHOUwieDat(data)$data.ou, liks = model.set.final$liks, Q=model.set.final$Q, rate=model.set.final$rate, root.p="yang", rate.cat=rate.cat, index.ou=index.ou, model.ou=model.ou, nSim=100, nCores=1, algorithm = "three.point")
  cat("\n")
  p = c(p.mk[1]/2, alpha[1], sig2[1], theta[1], theta[2])
  o01 <- OUwie:::hOUwie.dev(p = log(p), phy = phy, data.cor = OUwie:::organizeHOUwieDat(data)$data.cor, data.ou = OUwie:::organizeHOUwieDat(data)$data.ou, liks = model.set.final$liks, Q=model.set.final$Q, rate=model.set.final$rate, root.p="yang", rate.cat=rate.cat, index.ou=index.ou, model.ou=model.ou, nSim=100, nCores=1, algorithm = "three.point")
  cat("\n")
  p = c(p.mk[1]*2, alpha[1], sig2[1], theta[1], theta[2])
  o10 <- OUwie:::hOUwie.dev(p = log(p), phy = phy, data.cor = OUwie:::organizeHOUwieDat(data)$data.cor, data.ou = OUwie:::organizeHOUwieDat(data)$data.ou, liks = model.set.final$liks, Q=model.set.final$Q, rate=model.set.final$rate, root.p="yang", rate.cat=rate.cat, index.ou=index.ou, model.ou=model.ou, nSim=100, nCores=1, algorithm = "three.point")
  cat("\n")
  return(c(o1, o01, o10, no.changes))
}

A <- singleRun(phy, p.mk, root.p, alpha, sig2, theta, theta0, rate.cat, model.cor, model.ou)

B <- mclapply(1:100, function(x) singleRun(phy, p.mk, root.p, alpha, sig2, theta, theta0, rate.cat, model.cor, model.ou), mc.cores = 20)
tmp <- do.call(rbind, B)
length(which(apply(tmp[,1:3], 1, function(x) which.min(x)) == 2))/100
colMeans(tmp)


