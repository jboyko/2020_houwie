# fit the OU likelihood to the stochastic mapping
optFunc_CD <- function(p, tree, data, Q){
  out <- getAllJointProbs(tree, data, 1, 1.1, Q, rep(p[1],2), rep(p[2], 2), p[3:4], FALSE)
  print(c(log(colSums(exp(out))), p))
  return(-log(sum(exp(out[,3]))))
}

optFunc_CID <- function(p, tree, data, Q){
  out <- getAllJointProbs(tree, data, 1, 1.1, Q, rep(p[1],2), rep(p[2], 2), rep(p[3], 2), FALSE)
  print(c(log(colSums(exp(out))), p))
  return(-log(sum(exp(out[,3]))))
}

optFunc_CIDx <- function(p, tree, data, Q){
  out <- getAllJointProbs(tree, data, 2, time_slice, Q, rep(p[1],4), rep(p[2], 4), p[c(3,3,4,4)], FALSE)
  print(c(log(colSums(exp(out))), p))
  return(-log(sum(exp(out[,3]))))
}

setwd("~/2020_houwie/")
require(geiger)
require(corHMM)
require(OUwie)
require(partitions)
require(expm)
require(MASS)
require(phytools)
source("hOUwieNode.R")


balenced_tree <- read.tree(text = "((S1:0.5,S2:0.5):0.5,(S3:0.5,S4:0.5):0.5);")
pectinate_tree <- read.tree(text = "(((S1:0.33, S2:0.33):0.33,S3:0.66):0.34,S4:1);")
plot(balenced_tree)

plot(pectinate_tree)

nTip <- 4
possible_discrete <- sapply(1:16, function(x) as.numeric(intToBits(x)[1:4]))
possible_continuous <- (possible_discrete * 5) + 5
possible_combinations <- expand.grid(Disc=1:16, Cont=1:16)


time_slice <- 1.1
rate.cat <- 1
rate <- 1
alpha = c(4,4,4,4)
sigma.sq = c(1,1,1,1)
theta0 = 5
theta = c(5,5,10,10)
root <- c(1,0)
Q1 <- matrix(c(-rate[1], rate[1],rate[1],-rate[1]),2,2)
Q2 <- equateStateMatPars(getFullMat(list(getRateCatMat(2), getRateCatMat(2)), getRateCatMat(2)), 1:6)
Q2[Q2 > 0] <- rate
diag(Q2) <- -rowSums(Q2)

data <- data.frame(sp = c("S1", "S2", "S3", "S4"), reg = possible_datasets[,3], x = c(5, 5, 10, 10))
# out_CD <- getAllJointProbs(balenced_tree, data, 1, time_slice, Q1, alpha[1:2], sigma.sq[1:2], theta[2:1])
# out_CID <- getAllJointProbs(balenced_tree, data, 1, time_slice, Q1, alpha[1:2], sigma.sq[1:2], c(7.5,7.5))
# out_CIDx <- getAllJointProbs(balenced_tree, data, 2, time_slice, Q2, alpha, sigma.sq, theta[c(1,1,2,2)])
# balenced_out <- getAllJointProbs(balenced_tree, data, rate.cat, time_slice, Q, alpha, sigma.sq, theta)
# pectinate_out <- getAllJointProbs(pectinate_tree, data, rate.cat, time_slice, Q, alpha, sigma.sq, theta)

opts <- list("algorithm"="NLOPT_LN_SBPLX", "maxeval"="1000000", "ftol_rel"=.Machine$double.eps^0.1)

optFunc_CD(c(2,1,10,5), pectinate_tree, data, Q1)
optFunc_CID(c(2,1,7.5), pectinate_tree, data, Q1)
optFunc_CIDx(c(2,1,5,10), pectinate_tree, data, Q2)

out_CD = nloptr(x0=c(2,1,10,5), eval_f=optFunc_CD, lb=c(1e-5, 1e-5, 1, 1), ub=c(10,10,50,50), opts=opts, 
             tree=balenced_tree, data=data, Q = Q1)

out_CID = nloptr(x0=c(2,1,7.5), eval_f=optFunc_CID, lb=c(1e-5, 1e-5, 1), ub=c(10,10,50), opts=opts, 
                tree=balenced_tree, data=data, Q = Q1)

out_CIDx = nloptr(x0=c(2,1,10,5), eval_f=optFunc_CIDx, lb=c(1e-5, 1e-5, 1, 1), ub=c(10,10,50,50), opts=opts, 
                tree=balenced_tree, data=data, Q = Q2)





