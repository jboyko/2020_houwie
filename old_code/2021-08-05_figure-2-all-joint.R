# fit the OU likelihood to the stochastic mapping
setwd("~/2020_houwie/")
require(geiger)
require(corHMM)
require(OUwie)
require(partitions)
require(expm)
require(MASS)
require(phytools)
library(plyr)
library(ggplot2)
source("hOUwieNode.R")

tree <- read.tree(text = "((S1:0.5,S2:0.5):0.5,(S3:0.5,S4:0.5):0.5);")
rate <- 1
alpha = 2
sigma.sq = 1
theta = c(5,10)
Q1 <- matrix(c(-rate[1], rate[1],rate[1],-rate[1]),2,2)
Q2 <- equateStateMatPars(getFullMat(list(getRateCatMat(2), getRateCatMat(2)), getRateCatMat(2)), 1:6)
Q2[Q2 > 0] <- rate
diag(Q2) <- -rowSums(Q2)

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### CD model sampling 
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 

plot(tree)
data <- data.frame(sp = tree$tip.label, reg = c(1,1,2,2), x = c(5, 5, 10, 10))
out_CD <- getAllJointProbs(tree, data, 1, 0.5, Q1, alpha[c(1,1)], sigma.sq[c(1,1)], theta, FALSE)
log(colSums(exp(out_CD)))

data <- as.data.frame(out_CD)
data$disc_prob <- exp(data$disc)
data$disc_weight <- data$disc_prob/sum(data$disc_prob)
results <- data.frame()
for (nmodels in sequence(nrow(data)-1)) {
  for (rep in sequence(20)) {
    for (type in sequence(2)) {
      if(type==1) {
        focal_rows <- sample(sequence(nrow(data)), size=nmodels, prob=data$disc_weight, replace=FALSE)
      } else {
        focal_rows <- sample(sequence(nrow(data)), size=nmodels, prob=NULL, replace=FALSE)
      }
      results <- plyr::rbind.fill(results, data.frame(nmodel=nmodels, type=ifelse(type==1,"discrete-weighted", "flat-weighted"), total_lnl=log(sum(exp(data$total[focal_rows])))))
    }
  }
}
ggplot(results, aes(x=nmodel, y=total_lnl)) + geom_point(aes(colour=factor(type), alpha=0.3)) + facet_wrap(~type)
plot(x = data$disc, y = data$total, xlab = "disc_lnl", ylab = "total_lnl", pch = 16)

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### CID model sampling
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 

plot(tree)
data <- data.frame(sp = tree$tip.label, reg = c(1,1,2,2), x = c(7.5,7.5,7.5,7.5))
out_CID <- getAllJointProbs(tree, data, 1, 1.1, Q1, alpha[c(1,1)], sigma.sq[c(1,1)], c(7.5,7.5), FALSE)
log(colSums(exp(out_CID)))

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### CID model sampling
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 

plot(tree)
data <- data.frame(sp = tree$tip.label, reg = c(1,1,2,2), x = c(5,10,5,10))
out_CIDx <- getAllJointProbs(tree, data, 2, 1.1, Q2, alpha[c(1,1,1,1)], sigma.sq[c(1,1,1,1)], theta[c(1,2,1,2)], FALSE)
log(colSums(exp(out_CIDx)))




