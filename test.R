# testing
#
require(corHMM)
require(OUwie)
require(parallel)
source("~/2020_hOUwie/hOUwie.R")

data(tworegime)
data <- trait
phy <- tree
phy$node.label <- NULL
root.p = c(0.5, 0.5)
p.mk <- c(0.5, 0.5)
alpha = c(0.1, 0.1)
sigma.sq= c(0.5, 0.5)
theta = c(3, 8)
Q = matrix(c(-p.mk[1],p.mk[2],p.mk[1],-p.mk[2]), 2, 2)
theta0 = 3

full.data <- hOUwie.sim(phy, Q, root.p, alpha, sigma.sq, theta0, theta)
phy <- tree
phy$node.label <- NULL
dat <- full.data$data

index.cor <- corHMM:::getFullMat(list(corHMM:::getStateMat4Dat(data[,c(1,2)], "ER")$rate.mat, corHMM:::getStateMat4Dat(data[,c(1,2)], "ER")$rate.mat), corHMM:::getStateMat4Dat(data[,c(1,2)], "ER")$rate.mat)
# index.cor <- corHMM:::getStateMat4Dat(data[,c(1,2)], "ER")$rate.mat
index.ou <- getOUParamStructure("BMS", "three.point", FALSE, FALSE, dim(index.cor)[2])

test <- hOUwie(phy = phy, data = dat, rate.cat = 2, 
               index.cor = index.cor, root.p="yang", lb.cor=1e-5, ub.cor=1,
               index.ou = index.ou, root.station=FALSE, get.root.theta=FALSE, lb.ou=c(1e-5,1e-5,1e-1), ub.ou=c(2, 10, 10),  
               nSim = 10, weighted = FALSE, opts = list("algorithm"="NLOPT_LN_SBPLX", "maxeval"="10", "ftol_rel"=.Machine$double.eps^0.5))

test

undebug(hOUwie)
debug(organizeHOUwiePars)

source("~/2020_hOUwie/hOUwie.R")
test
# 