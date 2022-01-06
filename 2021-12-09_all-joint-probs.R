setwd("~/2020_houwie/")
require(geiger)
require(parallel)
require(corHMM)
require(OUwie)
require(partitions)
require(expm)
require(MASS)
require(phytools)
require(ggtree)
require(ggplot2)
require(reshape2)
require(gridExtra)
require(ggplotify)
source("hOUwieNode.R")
source("Utils.R")


# balenced_tree <- read.tree(text = "((S1:0.5,S2:0.5):0.5,(S3:0.5,S4:0.5):0.5);")
pectinate_tree <- read.tree(text = "(((S1:0.33, S2:0.33):0.33,S3:0.66):0.34,S4:1);")

plot(pectinate_tree)

# generate all possible discrete datasets
nTip <- 4
possible_discrete <- sapply(1:16, function(x) as.numeric(intToBits(x)[1:4]))
possible_continuous <- (possible_discrete * 5) + 5
possible_combinations <- expand.grid(Disc=1:16, Cont=1:16)

# parameters for estimation (the MLE is an alpha that goes to infinity and sigma goes to 0)
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
Q2[is.na(Q2)] <- 0
Q2[Q2 > 0] <- rate
diag(Q2) <- -rowSums(Q2)

# a function for examing joint probs for all possible datasets
getJoint <- function(poss_disc){
  data <- data.frame(sp = c("S1", "S2", "S3", "S4"), reg = poss_disc, x = c(5,10)[poss_disc + 1])
  out_CD <- getAllJointProbs(pectinate_tree, data, 1, time_slice, Q1, alpha[1:2], sigma.sq[1:2], theta[2:3])
  return(out_CD)
}

# the results of every possible discrete character, but we'll only be interested in number 5
full_joint_res <- apply(possible_discrete[,1:14], 2, getJoint)


fitAllJoint <- function(p, phy, dat, rate.cat){
  if(rate.cat == 1){
    Q <- matrix(c(-p[1], p[1],p[1],-p[1]),2,2)
    alpha = p[c(2,2)]
    sigma2 = p[c(3,3)]
    theta = p[c(4,5)]
  }
  if(rate.cat == 2){
    Q <- equateStateMatPars(getFullMat(list(getRateCatMat(2), getRateCatMat(2)), getRateCatMat(2)), 1:6)
    Q[is.na(Q)] <- 0
    Q[Q > 0] <- p[1]
    diag(Q) <- -rowSums(Q)
    alpha = p[c(2,2,2,2)]
    sigma2 = p[c(3,3,3,3)]
    theta = p[c(4,4,5,5)]
  }
  res <- getAllJointProbs(phy, dat, rate.cat, 1.1, Q, alpha, sigma2, theta)
  llik <- log(sum(exp(res$joint_probability_table[,3])))
  cat("\r LnLik=", llik, " p=", round(p,3), "        ")
  return(-llik)
}

getOptimOut <- function(optim){
  pars <- optim$solution
  names(pars) <- c("rate", "alpha", "sigma2", "theta_1", "theta_2")
  lnlik <- -optim$objective
  k <- length(pars)
  aic <- (-2 * lnlik) + (2 * k)
  return(c(lnlik=lnlik, k=k, aic=aic, pars))
}

getTable <- function(cd_optim, cid_optim){
  return(rbind(CD=getOptimOut(cd_optim), CID=getOptimOut(cid_optim)))
}

# data is consistent with character dependence
# we will use map 7 as our simuating map
regime_mapping_cd <- full_joint_res[[5]]$simmap_list[[7]]
plot(regime_mapping_cd)
dat.ou <- OUwie.sim(regime_mapping_cd, simmap.tree = TRUE, alpha = c(10,10), sigma.sq = c(5,5), theta0 = 5, theta = c(5,10))
dat_cd <- data.frame(sp = dat.ou$Genus_species, reg = c(2,1,2,1), x = dat.ou$X)
p = c(1, 10, 5, 5, 10)

dat_cd <- data.frame(sp = dat.ou$Genus_species, reg = c(2,1,2,1), x = dat.ou$X)
p = c(1, 10, 5, 5, 10)
opts <- list("algorithm"="NLOPT_LN_SBPLX", "maxeval"="1000", "ftol_rel"=.Machine$double.eps^0.1)
cd_optimized_cd_dat <- nloptr(x0 = p, eval_f = fitAllJoint, lb=rep(0, length(p)), ub=rep(50, length(p)), opts=opts, phy = pectinate_tree, dat = dat_cd, rate.cat = 1)
cid_optimized_cd_dat <- nloptr(x0 = p, eval_f = fitAllJoint, lb=rep(0, length(p)), ub=rep(50, length(p)), opts=opts, phy = pectinate_tree, dat = dat_cd, rate.cat = 2)

cont_model_cd <- getOUParamStructure("OUM", "three.point", FALSE, FALSE, 2)
cont_model_cid <- getOUParamStructure("OUM", "three.point", FALSE, FALSE, 4)
cont_model_cid[3,] <- c(3,3,4,4)
discrete_model_cd <- equateStateMatPars(getRateCatMat(2), 1:2)
discrete_model_cid <- getFullMat(list(discrete_model_cd, discrete_model_cd), getRateCatMat(2))
discrete_model_cid <- equateStateMatPars(discrete_model_cid, c(1:4))

cd_res <- hOUwie(pectinate_tree, dat_cd, rate.cat = 1, discrete_model_cd, cont_model_cd, nSim = 10, time_slice = 0.5, optimizer = "nlopt_ln", opts=opts, sample_nodes = TRUE)
cid_res <- hOUwie(pectinate_tree, dat_cd, rate.cat = 2, discrete_model_cid, cont_model_cid, nSim = 200, time_slice = 0.5, optimizer = "nlopt_ln", opts=opts, sample_nodes = TRUE)






# data is consistent with character independence
regime_mapping_cid <- full_joint_res[[3]]$simmap_list[[5]]
plot(regime_mapping_cid)
dat.ou <- OUwie.sim(regime_mapping_cid, simmap.tree = TRUE, alpha = c(10,10), sigma.sq = c(5,5), theta0 = 5, theta = c(5,10))
dat_cid <- data.frame(sp = dat.ou$Genus_species, reg = c(2,1,2,1), x = dat.ou$X)

p = c(1, 10, 5, 5, 10)
opts <- list("algorithm"="NLOPT_LN_SBPLX", "maxeval"="1000", "ftol_rel"=.Machine$double.eps^0.01)
cd_optimized_cid_dat <- nloptr(x0 = p, eval_f = fitAllJoint, lb=rep(0, length(p)), ub=rep(50, length(p)), opts=opts, phy = pectinate_tree, dat = dat_cid, rate.cat = 1)
cid_optimized_cid_dat <- nloptr(x0 = p, eval_f = fitAllJoint, lb=rep(0, length(p)), ub=rep(50, length(p)), opts=opts, phy = pectinate_tree, dat = dat_cid, rate.cat = 2)


getTable(cd_optimized_cd_dat, cid_optimized_cd_dat)
getTable(cd_optimized_cid_dat, cid_optimized_cid_dat)

# par(mfrow=c(1,2))
cols<-setNames(c("#e41a1c", "#377eb8"), c("1","2"))
xmax <- 1.4
xmin <- 1.12
cd_scalars <- dat_cd$x/max(dat_cd$x)
x_plot <- xmin + ((xmax - xmin) * cd_scalars)
pdf(file = "figures/raw/cd_data.pdf", height = 5, width = 8)
plotSimmap(regime_mapping_cd, color=cols, fsize = 1e-10, lwd = 20, outline = TRUE, xlim = c(0, xmax))
tiplabels(pch=c(15,16,15,16), offset = 0.06, cex = c(3,3.5,3,3.5))
lines(c(xmin, x_plot[1]), c(1,1), lwd=20)
lines(c(xmin, x_plot[2]), c(2,2), lwd=20)
lines(c(xmin, x_plot[3]), c(3,3), lwd=20)
lines(c(xmin, x_plot[4]), c(4,4), lwd=20)
dev.off()
# cid
cols<-setNames(c("#4daf4a", "#984ea3"), c("1","2"))
cid_scalars <- dat_cid$x/max(dat_cid$x)
x_plot <- xmin + ((xmax - xmin) * cid_scalars)
pdf(file = "figures/raw/cid_data.pdf", height = 5, width = 8)
plotSimmap(regime_mapping_cid, color=cols, fsize = 1e-10, lwd = 20, outline = TRUE, xlim = c(0, xmax))
tiplabels(pch=c(15,16,15,16), offset = 0.06, cex = c(3,3.5,3,3.5))
lines(c(xmin, x_plot[1]), c(1,1), lwd=20)
lines(c(xmin, x_plot[2]), c(2,2), lwd=20)
lines(c(xmin, x_plot[3]), c(3,3), lwd=20)
lines(c(xmin, x_plot[4]), c(4,4), lwd=20)
dev.off()

g <- as.grob(~plotSimmap(regime_mapping, color=cols, fsize = 1e-10, lwd = 5, ylim=c(-0.25,4.5), outline = TRUE))




