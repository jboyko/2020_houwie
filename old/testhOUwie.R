require(corHMM)
require(OUwie)
require(parallel)
data(tworegime)

trans.rt=c(0.01607089, 0.01707089)
alpha=c(0.5632459,0.1726052)
sigma.sq=c(0.1064417,0.3461386)
theta=c(1.678196,0.4185894)

hOUwie.dat <- organizeHOUwieDat(data = trait)
p <- c(trans.rt, alpha, sigma.sq, theta)

#Plot the tree and the internal nodes to highlight the selective regimes:
dat <- data.frame(sp = tree$tip.label, d = trait[,2])
MK <- corHMM(phy = tree, data = dat, rate.cat = 1, model = "ARD", p = p[1:2], get.tip.states = TRUE)

phy <- MK$phy
data <- MK$data
model <- MK$solution
model[is.na(model)] <- 0
diag(model) <- -rowSums(model)
states <- MK$states
tip.states <- MK$tip.states
phy$node.label <- NULL
## run get simmap (can be plotted using phytools)
simmap <- makeSimmap(tree=phy, tip.states=tip.states, states=states, model=model, nSim=1, nCores=1)

OURes <- OUwie.fixed(simmap[[1]],trait,model=c("OUMVA"), simmap.tree=TRUE, scaleHeight=FALSE, clade=NULL, alpha=alpha,sigma.sq=sigma.sq,theta=theta, algorithm="invert")



