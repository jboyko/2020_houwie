
# test data
require(OUwie)
require(phylolm)
data(tworegime)
#Plot the tree and the internal nodes to highlight the selective regimes:
select.reg<-character(length(tree$node.label))
select.reg[tree$node.label == 1] <- "black"
select.reg[tree$node.label == 2] <- "red"
plot(tree)
nodelabels(pch=21, bg=select.reg)
PhyloLMDat <- trait[,3]
names(PhyloLMDat) <- trait[,1]

# OUwieRes <- OUwie(tree,trait,model=c("OUM"),root.station=FALSE)
# alpha=c(OUwieRes$solution[1,])
# sigma.sq=c(OUwieRes$solution[2,])
# theta=c(OUwieRes$theta[1,1], OUwieRes$theta[1,1])
alpha=c(0.01, 0.01)
sigma.sq=c(0.02, 0.02)
theta=c(0, 0)

OUwie.fixed(tree,trait,model=c("OUM"), simmap.tree=FALSE, scaleHeight=FALSE,
            clade=NULL, alpha=alpha,sigma.sq=sigma.sq,theta=theta)
OU1d.loglik(trait=PhyloLMDat, phy=tree, model="OUfixedRoot", parameters=list(ancestral.state=theta[1], alpha=alpha[1],sigma2=sigma.sq[1], optimal.value=theta[1]))

## define our variables
# y = the dependent variable
y = PhyloLMDat

# X = the independent variable(s)
# the expected values of this model's OU process (one expectation per tip)
# calculate the expected values using the weight matrix and optima
# need phy, edges, rate.mat
n <- max(tree$edge[, 1])
edges <- cbind(c(1:(n - 1)), tree$edge, OUwie:::MakeAgeTable(tree,root.age = NULL))
edges <- edges[sort.list(edges[, 3]), ]
edges <- cbind(edges, matrix(c(1,0), dim(edges)[1], 2, byrow = TRUE))
OUwie:::varcov.ou(phy = tree, edges = edges, Rate.mat = OUwieRes$solution)
W <- OUwie:::weight.mat(phy = tree, edges = edges, Rate.mat = OUwieRes$solution, root.state = 1)
X <- W %*% c(OUwieRes$theta[1,1], 0)

## transform the tree
# based on the OU model being proposed change the branch lengths of the tree
# this step is necessary to use the three point algorithm and avoid the VCV inversion
times <- branching.times(tree)
distFromRoot <- exp(-2 * alpha * times) * (1 - exp(-2 * alpha * (max(times) - times)))
d1 = distFromRoot[anc - n]
d2 = numeric(N)
d2[externalEdge] = exp(-2 * alpha * D[des[externalEdge]]) * (1 - exp(-2 * alpha * (Tmax - D[des[externalEdge]])))d2[!externalEdge] = distFromRoot[des[!externalEdge] - n]

## calculate those pesky terms that slow computation
# use the three point algorithm on the transformed tree, X, and y to avoid matrix inversion
getThreePoint(tree, X, y)

## calculate the log likelihood for these parameters
# the equation is still a mystery to me and is ripped from phylolm code

debug(OUwie:::weight.mat)
debug(OUwie)
debug(OU1d.loglik)
