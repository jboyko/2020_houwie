# check OUwie basic against OUwie fixed
require(corHMM)
require(OUwie)
source("~/OUwie/R/hOUwie.R")

data(tworegime)

data <- trait
phy <- tree
phy$node.label <- NULL
root.p = c(0.5, 0.5)
p.mk <- c(0.5, 0.5)
alpha = c(5, 5)
sigma.sq= c(0.1, 0.1)
theta = c(3, 8)
Q = matrix(c(-p.mk[1],p.mk[2],p.mk[1],-p.mk[2]), 2, 2)
theta0 = 5

full.data <- OUwie:::hOUwie.sim(phy, Q, root.p, alpha, sigma.sq, theta0, theta)
phy <- full.data$simmap[[1]]
dat <- full.data$data

# BM1
alpha = c(1e-10, 1e-10)
sigma.sq= c(0.1, 0.1)
theta = c(8, 8)
fixed.lik <- OUwie.fixed(phy, dat, model="BM1", simmap.tree=TRUE, scaleHeight=FALSE, clade=NULL, sigma.sq=sigma.sq, theta=theta, check.identify=FALSE, algorithm="three.point", tip.paths=NULL, mserr="none", quiet=TRUE)$loglik[1]
basic.lik <- OUwie.basic(phy, dat, simmap.tree=TRUE, scaleHeight=FALSE, alpha=alpha, sigma.sq=sigma.sq, theta=theta, algorithm="three.point", tip.paths=NULL, mserr="none")
identical(round(fixed.lik, 5), round(basic.lik, 5))


# BMS
alpha = c(1e-10, 1e-10)
sigma.sq= c(0.1, 0.5)
theta = c(8, 8)
fixed.lik <- OUwie.fixed(phy, dat, model="BMS", simmap.tree=TRUE, scaleHeight=FALSE, clade=NULL, sigma.sq=sigma.sq, theta=theta, check.identify=FALSE, algorithm="three.point", tip.paths=NULL, mserr="none", quiet=TRUE)$loglik[1]
basic.lik <- OUwie.basic(phy, dat, simmap.tree=TRUE, scaleHeight=FALSE, alpha=alpha, sigma.sq=sigma.sq, theta=theta, algorithm="three.point", tip.paths=NULL, mserr="none")
identical(round(fixed.lik, 5), round(basic.lik, 5))

# OUM
alpha = c(5, 5)
sigma.sq= c(0.1, 0.1)
theta = c(3, 8)
fixed.lik <- OUwie.fixed(phy, dat, model="OUM", simmap.tree=TRUE, scaleHeight=FALSE, clade=NULL, alpha=alpha, sigma.sq=sigma.sq, theta=theta, check.identify=FALSE, algorithm="three.point", tip.paths=NULL, mserr="none", quiet=TRUE)$loglik[1]
basic.lik <- OUwie.basic(phy, dat, simmap.tree=TRUE, scaleHeight=FALSE, alpha=alpha, sigma.sq=sigma.sq, theta=theta, algorithm="three.point", tip.paths=NULL, mserr="none")
identical(round(fixed.lik, 5), round(basic.lik, 5))

# OUMA
alpha = c(5, 1)
sigma.sq= c(0.1, 0.1)
theta = c(3, 8)
fixed.lik <- OUwie.fixed(phy, dat, model="OUMA", simmap.tree=TRUE, scaleHeight=FALSE, clade=NULL, alpha=alpha, sigma.sq=sigma.sq, theta=theta, check.identify=FALSE, algorithm="three.point", tip.paths=NULL, mserr="none", quiet=TRUE)$loglik[1]
basic.lik <- OUwie.basic(phy, dat, simmap.tree=TRUE, scaleHeight=FALSE, alpha=alpha, sigma.sq=sigma.sq, theta=theta, algorithm="three.point", tip.paths=NULL, mserr="none")
identical(round(fixed.lik, 5), round(basic.lik, 5))

# OUMV
alpha = c(1, 1)
sigma.sq= c(0.1, 0.5)
theta = c(3, 8)
fixed.lik <- OUwie.fixed(phy, dat, model="OUMV", simmap.tree=TRUE, scaleHeight=FALSE, clade=NULL, alpha=alpha, sigma.sq=sigma.sq, theta=theta, check.identify=FALSE, algorithm="three.point", tip.paths=NULL, mserr="none", quiet=TRUE)$loglik[1]
basic.lik <- OUwie.basic(phy, dat, simmap.tree=TRUE, scaleHeight=FALSE, alpha=alpha, sigma.sq=sigma.sq, theta=theta, algorithm="three.point", tip.paths=NULL, mserr="none")
identical(round(fixed.lik, 5), round(basic.lik, 5))

# OUMVA
alpha = c(1, 5)
sigma.sq= c(0.1, 0.5)
theta = c(3, 8)
fixed.lik <- OUwie.fixed(phy, dat, model="OUMVA", simmap.tree=TRUE, scaleHeight=FALSE, clade=NULL, alpha=alpha, sigma.sq=sigma.sq, theta=theta, check.identify=FALSE, algorithm="three.point", tip.paths=NULL, mserr="none", quiet=TRUE)$loglik[1]
basic.lik <- OUwie.basic(phy, dat, simmap.tree=TRUE, scaleHeight=FALSE, alpha=alpha, sigma.sq=sigma.sq, theta=theta, algorithm="three.point", tip.paths=NULL, mserr="none")
identical(round(fixed.lik, 5), round(basic.lik, 5))

