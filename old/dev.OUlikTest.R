require(corHMM)
require(OUwie)
require(parallel)
require(geiger)
require(proftools)

phy <- sim.bdtree(b = 1, d = 0.5, stop = "taxa", n = 250)
phy <- drop.extinct(phy)
phy$edge.length <- phy$edge.length/max(branching.times(phy))
root.p = c(0.5, 0.5)
p.mk <- c(0.1)
alpha = c(0.01, 0.01)
sig2= c(0.1, 0.1)
theta = c(3, 8)
Q = matrix(c(-p.mk[1],p.mk[1],p.mk[1],-p.mk[1]), 2, 2)
theta0 = 5
rate.cat = 1
model.cor = "ER"
model.ou = "OUM"
root.age = NULL
algorithm = "three.point"

data <- OUwie:::hOUwie.sim(phy, Q, root.p, alpha, sig2, theta0, theta)[[1]]
p = c(0.1, 1, 0.1, 3, 8)
p.ou <- c(1, 0.1, 3, 8)
hOUwie.dat <- OUwie:::organizeHOUwieDat(data)
nObs <- length(hOUwie.dat$ObservedTraits)
model.set.final <- corHMM:::rate.cat.set.corHMM.JDB(phy=phy,data=hOUwie.dat$data.cor,rate.cat=rate.cat, ntraits = nObs, model = model.cor)
phy <- reorder(phy, "pruningwise")
index.ou <- OUwie:::getParamStructure(model.ou, algorithm, FALSE, FALSE, dim(model.set.final$Q)[2])

Rate.mat <- matrix(1, 3, dim(model.set.final$rate)[2])
Rate.mat[] <- c(p.ou, 1e-10)[index.ou]

n <- max(phy$edge[,1])
ntips <- length(phy$tip.label)
age.table <- OUwie:::MakeAgeTable(phy, root.age=root.age)
edges <- cbind(c(1:(n-1)),phy$edge, age.table)
Tmax <- max(age.table)
# automatic rescaling
edges[,4:5]<-edges[,4:5]/Tmax
root.age <-  1
edges <- edges[sort.list(edges[,3]),]
edges <- edges[sort.list(edges[,1]),]
if(algorithm == "three.point"){
  x <- hOUwie.dat$data.ou[,3]
  names(x) <- hOUwie.dat$data.ou[,1]
}else{
  x <- as.matrix(data[,3])
}
Q <- model.set.final$Q
Q[] <- c(p.mk, 0)[model.set.final$rate]
diag(Q) <- -rowSums(Q)
tot.states <- as.factor(c(1,2))

simmap.tree <- makeSimmap(phy, OUwie:::organizeHOUwieDat(data)$data.cor, Q, 1)[[1]]
# simmap.tree$edge.length <- simmap.tree$edge.length/Tmax
# simmap.tree$maps <- lapply(simmap.tree$maps, function(x) x/Tmax)
pars <- matrix(c(Rate.mat[3,], Rate.mat[2,], Rate.mat[1,]), dim(Rate.mat)[2], 3, dimnames = list(levels(tot.states), c("opt", "sig", "alp")))
root.state <- names(simmap.tree$maps[which(phy$edge[,1] == min(phy$edge[,1]))][[1]][1])
# W <- OUwie:::weight.mat(simmap.tree, edges, Rate.mat, root.state=NULL, simmap.tree=TRUE, root.age=root.age, scaleHeight=FALSE, assume.station=TRUE, shift.point=0.5)
W <- OUwie:::weight.mat(simmap.tree, edges, Rate.mat, root.state = as.numeric(root.state), simmap.tree = TRUE, root.age = NULL, scaleHeight = FALSE, assume.station = TRUE, shift.point = 0.5)

expected.vals <- colSums(t(W) * pars[,1])
names(expected.vals) <- phy$tip.label

transformPhy(simmap.tree, simmap.tree$maps, pars)

OUwie.fixed(simmap.tree, OUwie:::organizeHOUwieDat(data)$data.ou, model.ou, TRUE, alpha = alpha, sigma.sq = sig2, theta = theta, algorithm = "three.point")
OUwie:::dev.loglik.ouwie(simmap.tree, x, edges, Rate.mat)
OUwie:::getOULik(simmap.tree, x, expected.vals, simmap.tree$maps, pars)
