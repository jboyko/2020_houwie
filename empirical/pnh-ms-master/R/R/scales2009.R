setwd("~/2020_hOUwie/empirical/pnh-ms-master/R/R")

library(treeplyr)
library(devtools)
library(geiger)
library(phytools)
library(foreach)
library(doParallel)

scales <- read.csv("../data/scales2009.csv")
tetTree <- read.tree("../data/tetrapods.tre")
tetTree$tip.label[tetTree$tip.label == "Eumeces_schneideri"] <- "Eumeces_fasciatus"
tetTree <- multi2di(tetTree)
tetTree$edge.length[tetTree$edge.length==0] <- .Machine$double.eps
td <- make.treedata(tetTree, scales)
rm(tetTree)

td <- reorder(td, "postorder")
phenogram(td$phy, td[['FG.frac']], spread.labels=FALSE)

## Let's make the hypotheses from Scales 2009
#H.FM <- identifyBranches(td$phy, 7)
H.FM <- list(k=3, ntheta=3, sb=c(34, 36, 40), t2=c(3, 3, 2), loc=c(0,0,0)) 
H.PE <- list(k=5, ntheta=3, sb=c(18, 24, 29, 25, 37), t2=c(3, 3, 2, 2, 2), loc=c(0,0,0,0,0))
H.FMPE <- list(k=8, ntheta=6, sb=c(29, 25, 24, 18, 34, 37, 36, 40), t2=c(2, 2, 6, 6, 5, 3, 2, 4), loc=rep(0,8))
H.phryno <- list(k=1, ntheta=2, sb=c(24), t2=2, loc=0)
par(mfrow=c(2,3))
plotBayoupars(H.FM, td$phy, col=setNames(c("yellow", "purple", "darkgreen"), 1:3), lwd=3)
plotBayoupars(H.PE, td$phy, col=setNames(c("purple", "yellow",  "darkgreen"), 1:3), lwd=3)
plotBayoupars(H.FMPE, td$phy, col=setNames(c("red", "orange", "yellow","purple", "blue", "darkgreen"), 1:6), lwd=3)

phenogram(pars2simmap(H.FM, tree)$tr, dat, spread.labels=FALSE, fsize=0.5)
phenogram(pars2simmap(H.PE, tree)$tr, dat, spread.labels=FALSE, fsize=0.5)
phenogram(pars2simmap(H.FMPE, tree)$tr, dat, spread.labels=FALSE, fsize=0.5)


td <- reorder(td, "postorder")
tree <- td$phy
#tree$edge.length <- tree$edge.length/max(branching.times(tree$edge.length))
dat <- td[['FG.frac']]


prior.rj <- make.prior(tree, 
                       dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", 
                                  dk="cdpois", dtheta="dnorm"),
                       param=list(dalpha=list(scale=0.1), dsig2=list(scale=0.1),
                                  dk=list(lambda=3, kmax=10), dsb=list(bmax=1, prob=1), 
                                  dtheta=list(mean=mean(dat), sd=1.5*sd(dat)))
)

prior.FM <- make.prior(tree, 
                       dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", 
                                  dk="fixed", dtheta="dnorm", dsb="fixed"),
                       param=list(dalpha=list(scale=0.1), dsig2=list(scale=0.1),
                                  dk="fixed", dsb="fixed",
                                  dtheta=list(mean=mean(dat), sd=1.5*sd(dat))), 
                       fixed=list(k=H.FM$k, ntheta=H.FM$ntheta, sb=H.FM$sb, t2=H.FM$t2)
)

prior.PE <- make.prior(tree, 
                       dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", 
                                  dk="fixed", dtheta="dnorm", dsb="fixed"),
                       param=list(dalpha=list(scale=0.1), dsig2=list(scale=0.1),
                                  dk="fixed", dsb="fixed", 
                                  dtheta=list(mean=mean(dat), sd=1.5*sd(dat))), 
                       fixed=list(k=H.PE$k, ntheta=H.PE$ntheta, sb=H.PE$sb, t2=H.PE$t2)
)

prior.FMPE <- make.prior(tree, 
                          dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", 
                                     dk="fixed", dtheta="dnorm", dsb="fixed"),
                          param=list(dalpha=list(scale=0.1), dsig2=list(scale=0.1),
                                     dk="fixed", dsb="fixed",
                                     dtheta=list(mean=mean(dat), sd=1.5*sd(dat))), 
                          fixed=list(k=H.FMPE$k, ntheta=H.FMPE$ntheta, sb=H.FMPE$sb, t2=H.FMPE$t2)
)

prior.Clade <- make.prior(tree, 
                            dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", 
                                       dk="fixed", dtheta="dnorm", dsb="fixed"),
                            param=list(dalpha=list(scale=0.1), dsig2=list(scale=0.1),
                                       dk="fixed", dsb="fixed",
                                       dtheta=list(mean=mean(dat), sd=1.5*sd(dat))), 
                            fixed=list(k=4, ntheta=5, sb=c(25, 24, 37, 36), t2=2:5)
)

prior.phryno <-  make.prior(tree, 
                            dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", 
                                       dk="fixed", dtheta="dnorm", dsb="fixed"),
                            param=list(dalpha=list(scale=0.1), dsig2=list(scale=0.1),
                                       dk="fixed", dsb="fixed",
                                       dtheta=list(mean=mean(dat), sd=1.5*sd(dat))), 
                            fixed=list(k=1, ntheta=2, sb=c(24), t2=2)
)


mcmc.rj <- bayou.makeMCMC(tree, dat, SE=0.01, prior=prior.rj, 
                         new.dir="../output/scales/", outname="rj_r001", plot.freq=NULL, samp=50) # Set up the MCMC
mcmc.FM <- bayou.makeMCMC(tree, dat, SE=0.01, prior=prior.FM, 
                          new.dir="../output/scales/", outname="FM_r001", plot.freq=NULL, samp=50) 
mcmc.PE <- bayou.makeMCMC(tree, dat, SE=0.01, prior=prior.PE, 
                          new.dir="../output/scales/", outname="PE_r001", plot.freq=NULL, samp=50) 
mcmc.FMPE <- bayou.makeMCMC(tree, dat, SE=0.01, prior=prior.FMPE, 
                          new.dir="../output/scales/", outname="FMPE_r001", plot.freq=NULL, samp=50) 
mcmc.Clade <- bayou.makeMCMC(tree, dat, SE=0.01, prior=prior.Clade, 
                            new.dir="../output/scales/", outname="Clade_r001", plot.freq=NULL, samp=50) 
mcmc.phryno <- bayou.makeMCMC(tree, dat, SE=0.01, prior=prior.phryno, 
                              new.dir="../output/scales/", outname="phryno_r001", plot.freq=NULL, samp=50) 

mcmc.FM$run(1000000)
mcmc.PE$run(1000000)
mcmc.FMPE$run(1000000)
mcmc.rj$run(1000000) # Run the MCMC
mcmc.Clade$run(1000000)
mcmc.phryno$run(1000000)

chain.rj <- mcmc.rj$load()
chain.FM <- mcmc.FM$load()
chain.PE <- mcmc.PE$load()
chain.FMPE <- mcmc.FMPE$load()
chain.Clade <- mcmc.Clade$load()
chain.phryno <- mcmc.phryno$load()

chain.rj <- set.burnin(chain.rj, 0.3)
chain.FM <- set.burnin(chain.FM, 0.3)
chain.PE <- set.burnin(chain.PE, 0.3)
chain.FMPE <- set.burnin(chain.FMPE, 0.3)
chain.Clade <- set.burnin(chain.Clade, 0.3)
chain.phryno <- set.burnin(chain.phryno, 0.3)

Bk <- qbeta(seq(0,1, length.out=10), 0.3,1)
require(foreach)
registerDoParallel(cores=10)
ss.rj <- mcmc.rj$steppingstone(1000000, chain.rj, Bk=Bk, burnin=0.3, plot=FALSE)
ss.FM <- mcmc.FM$steppingstone(1000000, chain.FM, Bk=Bk, burnin=0.3)
ss.PE <- mcmc.PE$steppingstone(1000000, chain.PE, Bk=Bk, burnin=0.1)
ss.FMPE <- mcmc.FMPE$steppingstone(1000000, chain.FMPE, Bk=Bk, burnin=0.3)
ss.Clade <- mcmc.Clade$steppingstone(1000000, chain.Clade, Bk=Bk, burnin=0.3)
ss.phryno <- mcmc.phryno$steppingstone(1000000, chain.phryno, Bk=Bk, burnin=0.3)

saveRDS(ss.rj, file="ss.rj.rds")
saveRDS(ss.rj, file="ss.FM.rds")
saveRDS(ss.rj, file="ss.PE.rds")
saveRDS(ss.rj, file="ss.FMPE.rds")
saveRDS(ss.rj, file="ss.Clade.rds")
saveRDS(ss.rj, file="ss.phryno.rds")

saveRDS(chain.rj, file="chain.rj.rds")
saveRDS(chain.FM, file="chain.FM.rds")
saveRDS(chain.PE, file="chain.PE.rds")
saveRDS(chain.FMPE, file="chain.FMPE.rds")
saveRDS(chain.Clade, file="chain.Clade.rds")
saveRDS(chain.phryno, file="chain.phryno.rds")

