#simplest case
require(OUwie)
require(phylolm)
require(geiger)
data(tworegime)
#Plot the tree and the internal nodes to highlight the selective regimes:
select.reg<-character(length(tree$node.label))
select.reg[tree$node.label == 1] <- "black"
select.reg[tree$node.label == 2] <- "red"
tree$edge.length <- tree$edge.length/max(branching.times(tree))
StdErr <- c(0.0, 0.02, 0.05, 0.1, 0.2, 0.5, 1)
LikRes <- matrix(NA, 3, length(StdErr), dimnames = list(c("Phylolm", "Invert", "3-point"), StdErr))

for(i in 1:length(StdErr)){
  PhyloLMDat <- trait[,3]
  trait[,4] <- StdErr[i]
  names(PhyloLMDat) <- trait[,1]
  alpha=c(2, 2)
  sigma.sq=c(0.1, 0.1)
  theta=c(3, 3)
  PLM <- OU1d.loglik(trait=PhyloLMDat, phy=tree, model="OUfixedRoot", parameters=list(ancestral.state=theta[1], alpha=alpha[1],sigma2=sigma.sq[1], optimal.value=theta[1], sigma2_error=StdErr[i]^2))
  INV <- OUwie.fixed(tree,trait,model=c("OU1"), simmap.tree=FALSE, scaleHeight=FALSE, mserr = "known", clade=NULL, alpha=alpha,sigma.sq=sigma.sq,theta=theta, algorithm = "invert")$loglik
  TPT <- OUwie.fixed(tree,trait,model=c("OU1"), simmap.tree=FALSE, scaleHeight=FALSE, mserr = "known", clade=NULL, alpha=alpha,sigma.sq=sigma.sq,theta=theta, algorithm = "three.point")$loglik
  LikRes[,i] <- c(PLM, INV, TPT)
}
LikRes



require(OUwie)
require(phylolm)
require(geiger)
data(tworegime)
#Plot the tree and the internal nodes to highlight the selective regimes:
select.reg<-character(length(tree$node.label))
select.reg[tree$node.label == 1] <- "black"
select.reg[tree$node.label == 2] <- "red"
tree$edge.length <- tree$edge.length/max(branching.times(tree))
StdErr <- c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1)
LikRes <- matrix(NA, 3, length(StdErr), dimnames = list(c("Geiger", "Invert", "3-point"), StdErr))
for(i in 1:length(StdErr)){
  PhyloLMDat <- trait[,3]
  StdErrG <- abs(rnorm(length(tree$tip.label), mean = StdErr[i], sd = StdErr[i]/10))
  trait[,4] <- StdErrG
  names(PhyloLMDat) <- trait[,1]
  names(StdErrG) <- trait[,1]
  fitCon <- fitContinuous(tree, PhyloLMDat, SE = StdErrG, model = "OU")
  if(fitCon$opt$alpha > 1e-10){
    alpha=c(fitCon$opt$alpha, fitCon$opt$alpha)
  }else{
    alpha=c(1e-10, 1e-10)
  }
  sigma.sq=c(fitCon$opt$sigsq, fitCon$opt$sigsq)
  theta=c(fitCon$opt$z0, fitCon$opt$z0)
  INV <- OUwie.fixed(tree,trait,model=c("OU1"), simmap.tree=FALSE, scaleHeight=FALSE, mserr = "known", clade=NULL, alpha=alpha,sigma.sq=sigma.sq,theta=theta, algorithm = "invert")$loglik
  GEI <- fitCon$opt$lnL
  TPT <- OUwie.fixed(tree,trait,model=c("OU1"), simmap.tree=FALSE, scaleHeight=FALSE, mserr = "known", clade=NULL, alpha=alpha,sigma.sq=sigma.sq,theta=theta, algorithm = "three.point")$loglik
  LikRes[,i] <- c(GEI, INV, TPT)
}




require(OUwie)
require(phylolm)
require(geiger)
data(tworegime)
#Plot the tree and the internal nodes to highlight the selective regimes:
select.reg<-character(length(tree$node.label))
select.reg[tree$node.label == 1] <- "black"
select.reg[tree$node.label == 2] <- "red"
tree$edge.length <- tree$edge.length/max(branching.times(tree))
StdErr <- c(0.00, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2)
LikRes <- matrix(NA, 2, length(StdErr), dimnames = list(c("Invert", "3-point"), StdErr))
for(i in 1:length(StdErr)){
  trait[,4] <- abs(rnorm(length(tree$tip.label), mean = StdErr[i], sd = StdErr[i]/10))
  alpha=c(1, 2)
  sigma.sq=c(1, 2)
  theta=c(5, 10)
  INV <- OUwie.fixed(tree,trait,model=c("OUMVA"), simmap.tree=FALSE, scaleHeight=FALSE, mserr = "known", clade=NULL, alpha=alpha,sigma.sq=sigma.sq,theta=theta, algorithm = "invert")$loglik
  TPT <- OUwie.fixed(tree,trait,model=c("OUMVA"), simmap.tree=FALSE, scaleHeight=FALSE, mserr = "known", clade=NULL, alpha=alpha,sigma.sq=sigma.sq,theta=theta, algorithm = "three.point")$loglik
  LikRes[,i] <- c(INV, TPT)
}

debug(vcvPhylo)
vcvPhylo(tree, anc.nodes = FALSE)

-142.7562

prm = list(alpha = p$alpha, sigma2_error = p$sigma2_error * 2 * p$alpha/p$sigma2)
errEdge = rep(p$sigma2_error, n)
errEdge = errEdge * exp(-2 * alpha * D[des[externalEdge]])

## preparing for general use of "parameter" for branch length transformation
sigma2_error = p$sigma2_error * 2 * p$alpha/p$sigma2 
# note that sigma2_error = true_sigma2_error/sigma2
# sigma2_error = p$sigma2_error / (p$sigma2 / 2 * p$alpha)
errEdge = rep(p$sigma2_error,n)
errEdge = errEdge*exp(-2*alpha*D[des[externalEdge]]) # adjust measurement errors for OU models





require(corHMM)
require(OUwie)
require(parallel)
require(geiger)
require(proftools)

phy <- sim.bdtree(b = 1, d = 0.5, stop = "taxa", n = 250)
phy <- drop.extinct(phy)
root.p = c(0.5, 0.5)
p.mk <- c(0.5, 0.5)
alpha = c(1, 2)
sig2= c(0.1, 0.1)
theta = c(3, 8)
Q = matrix(c(-p.mk[1],p.mk[2],p.mk[1],-p.mk[2]), 2, 2)
theta0 = 5
rate.cat = 1
model.cor = "ER"
model.ou = "OUMA"

data <- OUwie:::hOUwie.sim(phy, Q, root.p, alpha, sig2, theta0, theta)[[1]]
p = c(0.1, 0.01, 0.1, 1, 3, 8)
hOUwie.dat <- OUwie:::organizeHOUwieDat(data)
nObs <- length(hOUwie.dat$ObservedTraits)
model.set.final <- corHMM:::rate.cat.set.corHMM.JDB(phy=phy,data=hOUwie.dat$data.cor,rate.cat=rate.cat, ntraits = nObs, model = model.cor)
phy <- reorder(phy, "pruningwise")
index.ou <- OUwie:::getParamStructure(model.ou, "three.point", FALSE, FALSE, dim(model.set.final$Q)[2])
# phy$edge.length <- phy$edge.length/max(branching.times(phy))

OUwie:::hOUwie.dev(p = log(p), phy = phy, data.cor = OUwie:::organizeHOUwieDat(data)$data.cor, data.ou = OUwie:::organizeHOUwieDat(data)$data.ou, liks = model.set.final$liks, Q=model.set.final$Q, rate=model.set.final$rate, root.p=root.p, rate.cat=rate.cat, index.ou=index.ou, model.ou=model.ou, nSim=1, nCores=1, algorithm = "three.point")

INV <- unlist(mclapply(simmap, function(x) OUwie.fixed(x, data.ou, model = model.ou, simmap.tree = TRUE, scaleHeight = FALSE, clade = NULL, alpha = alpha, sigma.sq = sigma.sq, theta = theta, check.identify = FALSE, algorithm = "invert", tip.paths = tip.paths, quiet = TRUE)$loglik, mc.cores = nCores))

TPT <- unlist(mclapply(simmap, function(x) OUwie.fixed(x, data.ou, model = model.ou, simmap.tree = TRUE, scaleHeight = FALSE, clade = NULL, alpha = alpha, sigma.sq = sigma.sq, theta = theta, check.identify = FALSE, algorithm = "three.point", tip.paths = tip.paths, quiet = TRUE)$loglik, mc.cores = nCores))


