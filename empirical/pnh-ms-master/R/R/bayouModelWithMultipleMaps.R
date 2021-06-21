setwd("~/repos/pnh-ms/R/R")

library(bayou)
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


pars <- list(alpha=0.15, sig2=0.001, k=1, ntheta=2, theta=c(0.55, 0.68), sb=c(25), loc=rep(0, 1), t2=2)

simpars1 <- bayou::priorSim(prior.phryno, td$phy, nsim=1)
simpars1$pars[[1]]$alpha <- 0.15
simpars1$pars[[1]]$sig2 <- 0.001
simpars1$pars[[1]]$theta <- c(0.6, 0.3)
simdat1 <- bayou::dataSim(simpars1$pars[[1]], model="OU", td$phy)

simpars2 <- bayou::priorSim(prior.PE, td$phy, nsim=1)
simpars2$pars[[1]]$alpha <- 0.15
simpars2$pars[[1]]$sig2 <- 0.001
simpars2$pars[[1]]$theta <- c(0.5, 0.7, 0.3)
simdat2 <- bayou::dataSim(simpars2$pars[[1]], model="OU", td$phy)

## Now this is the trait-based answer for PE
pars$thPE1 <- 0.55
pars$thPE2 <- 0.7
pars$thPE3 <- 0.3
fixedpars <- list()
fixedpars$k <- H.PE$k
fixedpars$ntheta <- H.PE$ntheta
fixedpars$sb <- H.PE$sb
fixedpars$loc <- H.PE$loc
fixedpars$t2 <- H.PE$t2

pars$w <- 0.5

PE.lik <- function(pars, cache, X, model="Custom"){
  n <- cache$n
  X <- cache$dat
  pred <- cache$pred
  EX.map1 <- bayou:::C_weightmatrix(cache, pars)$E
  pars_PE <- list(alpha=pars$alpha, sig2=pars$sig2, k=fixedpars$k, ntheta=fixedpars$ntheta, theta=c(pars$thPE1, pars$thPE2, pars$thPE3), sb=fixedpars$sb, loc=fixedpars$loc, t2=fixedpars$t2)
  EX.map2 <- bayou:::C_weightmatrix(cache, pars_PE)$E  
  X.c <- X - pars$w*(EX.map1) - (1-pars$w)*(EX.map2)
  ## This part adds the endothermy parameter to the theta for Mammal and Bird branches
  #dpars <- pars
  #dpars$theta[dpars$t2[which(dpars$sb %in% c(857, 1719))]] <- dpars$theta[dpars$t2[which(dpars$sb %in% c(857, 1719))]]+dpars$endo
  ### The part below mostly does not change
  transf.phy <- bayou:::C_transf_branch_lengths(cache, 1, X.c, pars$alpha)
  transf.phy$edge.length[cache$externalEdge] <- transf.phy$edge[cache$externalEdge] + cache$SE[cache$phy$edge[cache$externalEdge, 2]]^2*(2*pars$alpha)/pars$sig2
  comp <- bayou:::C_threepoint(list(n=n, N=cache$N, anc=cache$phy$edge[, 1], des=cache$phy$edge[, 2], diagMatrix=transf.phy$diagMatrix, P=X.c, root=transf.phy$root.edge, len=transf.phy$edge.length))
  if(pars$alpha==0){
    inv.yVy <- comp$PP
    detV <- comp$logd
  } else {
    inv.yVy <- comp$PP*(2*pars$alpha)/(pars$sig2)
    detV <- comp$logd+n*log(pars$sig2/(2*pars$alpha))
  }
  llh <- -0.5*(n*log(2*pi)+detV+inv.yVy)
  #llh <- llh + gs.lik(c(pars$pred.sig2, pars$pred.root), root=ROOT.GIVEN) #$impute
  return(list(loglik=llh, theta=pars$theta,resid=X.c, comp=comp, transf.phy=transf.phy))
}

PE.monitor = function(i, lik, pr, pars, accept, accept.type, j){
  names <- c("gen", "lnL", "prior", "alpha","sig2", "rtheta", "rthPE","w", "k")
  string <- "%-8i%-8.2f%-8.2f%-8.2f%-8.2f%-8.2f%-8.2f%-8.2f%-8i"
  acceptratios <- tapply(accept, accept.type, mean)
  names <- c(names, names(acceptratios))
  if(j==0){
    cat(sprintf("%-7.7s", names), "\n", sep=" ")                           
  }
  cat(sprintf(string, i, lik, pr, pars$alpha, pars$sig2, pars$theta[1], pars$thPE1, pars$w, pars$k), sprintf("%-8.2f", acceptratios),"\n", sep="")
}


model.weightedPE <- list(moves = list(alpha=".multiplierProposal", sig2=".multiplierProposal", 
                                      k=".splitmergePrior", theta=".adjustTheta", thPE1=".slidingWindowProposal", 
                                      thPE2=".slidingWindowProposal", 
                                      thPE3=".slidingWindowProposal", 
                                      w= ".multiplierProposal", slide=".slide"
                                      ),
                        control.weights = list(alpha=5, sig2=3, theta=15, thPE1=3,
                                               thPE2=3, thPE3=3,
                                              k=10, w=3, slide=1),
                        D = list(alpha=0.5, sig2= 0.5, k=1, theta=1, thPE1=1, 
                                 thPE2=1, thPE3=1, w=1, slide=1),
                        parorder = c("alpha", "sig2", "w","thPE1", "thPE2", "thPE3", "k", "ntheta",  "theta", "sb", "loc", "t2"),
                        rjpars = c("theta"),
                        shiftpars = c("sb", "loc", "t2"),
                        monitor.fn = PE.monitor,
                        lik.fn = PE.lik)


#pars <- pars[model.weightedPE$parorder]
#custom.lik(pars, cache, X)$loglik

#ct <- bayou:::.buildControl(pars, prior.wPE, move.weights=model.weightedPE$control.weights)
#pp <- bayou:::.splitmergePrior(pars, cache, d=1, ct=ct, prior=prior.wPE, move = "k")$pars
#bayou:::.vectorSlidingWindow(cache, pp, d=1, move="thPE1", prior=prior.wPE, ct=ct)$pars

prior.wPE <- make.prior(td$phy, plot.prior = FALSE, 
                    dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", 
                               dk="cdpois", dtheta="dnorm", dthPE1="dnorm",
                               dthPE2="dnorm",
                               dthPE3="dnorm",
                               dw="dbeta", dsb="dsb"
                              ), 
                    param=list(dalpha=list(scale=0.1), dsig2=list(scale=0.1), 
                               dk=list(lambda=0.5, kmax=10), dtheta=list(mean=0.5, sd=0.25),
                               dthPE1=list(mean=0.5, sd=0.25), 
                               dthPE2=list(mean=0.5, sd=0.25), 
                               dthPE3=list(mean=0.5, sd=0.25), 
                               dw=list(shape1=0.8, shape2=0.8), dsb=list(bmax=1, prob=1)
                              )
                    #fixed=list(sb=startpar$sb, k=startpar$k)
)


mcmc.wPE <- bayou.makeMCMC(td$phy, td[['FG.frac']], SE=0, prior=prior.wPE, model=model.weightedPE, startpar=pars,
                          new.dir="../output/scales/", outname="wPE_r001", plot.freq=NULL, samp=10, ticker.freq=1000) # Set up the MCMC

mcmc.wPE.sim1 <- bayou.makeMCMC(td$phy, simdat1$dat, SE=0, prior=prior.wPE, model=model.weightedPE, startpar=pars,
                           new.dir="../output/scales/", outname="wPEsim1_r001", plot.freq=NULL, samp=10, ticker.freq=1000) # Set up the MCMC

mcmc.wPE.sim2 <- bayou.makeMCMC(td$phy, simdat2$dat, SE=0, prior=prior.wPE, model=model.weightedPE, startpar=pars,
                                new.dir="../output/scales/", outname="wPEsim2_r001", plot.freq=NULL, samp=10, ticker.freq=1000) # Set up the MCMC

mcmc.wPE$run(200000)
mcmc.wPE.sim1$run(200000)
mcmc.wPE.sim2$run(200000)
chain.wPE <- mcmc.wPE$load()
chain.wPE.sim1 <- mcmc.wPE.sim1$load()
chain.wPE.sim2 <- mcmc.wPE.sim2$load()

summary(chain.wPE)


plotSimmap.mcmc(chain.wPE, burnin=0.3)
plotBayoupars(H.PE, td$phy)
postburn <- floor(0.3*length(chain.wPE$gen)):length(chain.wPE$gen)

pdf("~/repos/pnh-ms/R/output/scales/ScalesReanalysis.pdf", height=8, width=12)
par(mfrow=c(2,2), mar=c(5,5,0.75, 0.25))
phenogram(pars2simmap(H.PE, td$phy)$tr, dat, colors = setNames(viridis::viridis(3), c(2,3, 1)), fsize=1, ylab="FG fraction", xlab="Time (my)")

priorDens <- density(rbeta(length(postburn)*90, 0.8, 0.8), bw=0.005)
plot(priorDens$x, priorDens$y*2000, type="n", lty=2, ylim=c(0, 7000), xlab="Weight to PE Hypothesis", ylab="Density", yaxt="n")
lines(c(0, priorDens$x[priorDens$x>0 & priorDens$x<1], 1), c(0, priorDens$y[priorDens$x>0 & priorDens$x<1]*1500, 0), lty=2)
cols <- sapply(viridis::viridis(3), function(x) bayou:::makeTransparent(x, 50))
hist(1-chain.wPE$w[postburn], col=cols[1], add=TRUE, breaks=seq(0,1,0.025))
hist(1-chain.wPE.sim1$w[postburn], col=cols[2], add=TRUE, breaks=seq(0,1,0.025))
hist(1-chain.wPE.sim2$w[postburn], col=cols[3], add=TRUE, breaks=seq(0,1,0.025))
legend(0, 7000, legend = c("Prior","Scales data", "Simulated data, phrynos only", "Simulated data under PE Hypothesis"),
       pch=rep(22,4) , pt.bg=c("black", cols),lty=c(2,NA,NA,NA),  pt.cex=c(0,2,2,2))

L1 <- Lposterior(chain.wPE,tree, burnin=0.3)
L2 <- Lposterior(chain.wPE.sim2, tree, burnin=.3)
plot(L1$pp, L2$pp, ylab="Branch posterior prob. for simulated data (PE Hyp)", xlab="Branch posterior prob. for empirical data", 
     xlim=c(0,1), ylim=c(0,1), pch=21, bg="gray90", cex=1.5)
text(0.85, 0.1, "Phrynosoma")
abline(0,1, lty=2)
plotSimmap.mcmc(chain.wPE, burnin=0.3, pp.cutoff=0.5,lwd=2, pal=colorRampPalette(viridis::viridis(3)[c(1, 2)]))
dev.off()



####################
###Try again with Phryno only simulated data


Phryno.lik <- function(pars, cache, X, model="Custom"){
  n <- cache$n
  X <- cache$dat
  pred <- cache$pred
  EX.map1 <- bayou:::C_weightmatrix(cache, pars)$E
  pars_phryno <- list(alpha=pars$alpha, sig2=pars$sig2, k=H.phryno$k, ntheta=H.phryno$ntheta, theta=c(pars$thPhryno1, pars$thPhryno2), sb=H.phryno$sb, loc=H.phryno$loc, t2=H.phryno$t2)
  EX.map2 <- bayou:::C_weightmatrix(cache, pars_phryno)$E  
  X.c <- X - pars$w*(EX.map1) - (1-pars$w)*(EX.map2)
  ## This part adds the endothermy parameter to the theta for Mammal and Bird branches
  #dpars <- pars
  #dpars$theta[dpars$t2[which(dpars$sb %in% c(857, 1719))]] <- dpars$theta[dpars$t2[which(dpars$sb %in% c(857, 1719))]]+dpars$endo
  ### The part below mostly does not change
  transf.phy <- bayou:::C_transf_branch_lengths(cache, 1, X.c, pars$alpha)
  transf.phy$edge.length[cache$externalEdge] <- transf.phy$edge[cache$externalEdge] + cache$SE[cache$phy$edge[cache$externalEdge, 2]]^2*(2*pars$alpha)/pars$sig2
  comp <- bayou:::C_threepoint(list(n=n, N=cache$N, anc=cache$phy$edge[, 1], des=cache$phy$edge[, 2], diagMatrix=transf.phy$diagMatrix, P=X.c, root=transf.phy$root.edge, len=transf.phy$edge.length))
  if(pars$alpha==0){
    inv.yVy <- comp$PP
    detV <- comp$logd
  } else {
    inv.yVy <- comp$PP*(2*pars$alpha)/(pars$sig2)
    detV <- comp$logd+n*log(pars$sig2/(2*pars$alpha))
  }
  llh <- -0.5*(n*log(2*pi)+detV+inv.yVy)
  #llh <- llh + gs.lik(c(pars$pred.sig2, pars$pred.root), root=ROOT.GIVEN) #$impute
  return(list(loglik=llh, theta=pars$theta,resid=X.c, comp=comp, transf.phy=transf.phy))
}

Phryno.monitor = function(i, lik, pr, pars, accept, accept.type, j){
  names <- c("gen", "lnL", "prior", "alpha","sig2", "rtheta", "rthPh","w", "k")
  string <- "%-8i%-8.2f%-8.2f%-8.2f%-8.2f%-8.2f%-8.2f%-8.2f%-8i"
  acceptratios <- tapply(accept, accept.type, mean)
  names <- c(names, names(acceptratios))
  if(j==0){
    cat(sprintf("%-7.7s", names), "\n", sep=" ")                           
  }
  cat(sprintf(string, i, lik, pr, pars$alpha, pars$sig2, pars$theta[1], pars$thPhryno1, pars$w, pars$k), sprintf("%-8.2f", acceptratios),"\n", sep="")
}

simpars1 <- bayou::priorSim(prior.phryno, td$phy, nsim=1)
simpars1$pars[[1]]$alpha <- 0.15
simpars1$pars[[1]]$sig2 <- 0.001
simpars1$pars[[1]]$theta <- c(0.6, 0.3)
simdat1 <- bayou::dataSim(simpars1$pars[[1]], model="OU", td$phy)

simpars2 <- bayou::priorSim(prior.PE, td$phy, nsim=1)
simpars2$pars[[1]]$alpha <- 0.15
simpars2$pars[[1]]$sig2 <- 0.001
simpars2$pars[[1]]$theta <- c(0.5, 0.7, 0.3)
simdat2 <- bayou::dataSim(simpars2$pars[[1]], model="OU", td$phy)

model.weightedPhryno <- list(moves = list(alpha=".multiplierProposal", sig2=".multiplierProposal", 
                                      k=".splitmergePrior", theta=".adjustTheta", thPhryno1=".slidingWindowProposal", 
                                      thPhryno2=".slidingWindowProposal", 
                                      w= ".multiplierProposal", slide=".slide"
                                      ),
          control.weights = list(alpha=5, sig2=3, theta=15, thPhryno1=3,
                       thPhryno2=3,
                       k=10, w=3, slide=1),
          D = list(alpha=0.5, sig2= 0.5, k=1, theta=1, thPhryno1=1, 
                       thPhryno2=1, w=1, slide=1),
          parorder = c("alpha", "sig2", "w","thPhryno1", "thPhryno2", "k", "ntheta",  "theta", "sb", "loc", "t2"),
          rjpars = c("theta"),
          shiftpars = c("sb", "loc", "t2"),
          monitor.fn = Phryno.monitor,
          lik.fn = Phryno.lik
)

prior.wPhryno <- make.prior(td$phy, plot.prior = FALSE, 
                        dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", 
                                   dk="cdpois", dtheta="dnorm", dthPhryno1="dnorm",
                                   dthPhryno2="dnorm",
                                   dw="dbeta", dsb="dsb"
                        ), 
                        param=list(dalpha=list(scale=0.1), dsig2=list(scale=0.1), 
                                   dk=list(lambda=0.5, kmax=10), dtheta=list(mean=0.5, sd=0.25),
                                   dthPhryno1=list(mean=0.5, sd=0.25), 
                                   dthPhryno2=list(mean=0.5, sd=0.25), 
                                   dw=list(shape1=0.8, shape2=0.8), dsb=list(bmax=1, prob=1)
                        )
                        #fixed=list(sb=startpar$sb, k=startpar$k)
)

pars <- list(alpha=0.15, sig2=0.001, k=1, ntheta=2, theta=c(0.55, 0.68), sb=c(25), loc=rep(0, 1), t2=2)
pars$thPhryno1 <- 0.55
pars$thPhryno2 <- 0.33
pars$w <- 0.5


mcmc.wPhryno <- bayou.makeMCMC(td$phy, td[['FG.frac']], SE=0, prior=prior.wPhryno, model=model.weightedPhryno, startpar=pars,
                           new.dir="../output/scales/", outname="wPhryno_r001", plot.freq=NULL, samp=10, ticker.freq=1000) # Set up the MCMC

mcmc.wPhryno.sim1 <- bayou.makeMCMC(td$phy, simdat1$dat, SE=0, prior=prior.wPhryno, model=model.weightedPhryno, startpar=pars,
                                new.dir="../output/scales/", outname="wPhrynosim1_r001", plot.freq=NULL, samp=10, ticker.freq=1000) # Set up the MCMC

mcmc.wPhryno.sim2 <- bayou.makeMCMC(td$phy, simdat2$dat, SE=0, prior=prior.wPhryno, model=model.weightedPhryno, startpar=pars,
                                new.dir="../output/scales/", outname="wPhrynosim2_r001", plot.freq=NULL, samp=10, ticker.freq=1000) # Set up the MCMC


mcmc.wPhryno$run(200000)
mcmc.wPhryno.sim1$run(200000)
mcmc.wPhryno.sim2$run(200000)

chain.wPhryno <- mcmc.wPhryno$load()
chain.wPhryno.sim1 <- mcmc.wPhryno.sim1$load()
chain.wPhryno.sim2 <- mcmc.wPhryno.sim2$load()

postburn <- floor(0.3*length(chain.wPhryno$gen)):length(chain.wPhryno$gen)
plot(density(1-chain.wPhryno$w[postburn]))
lines(density(1-chain.wPhryno.sim1$w[postburn]), col="red")
lines(density(1-chain.wPhryno.sim2$w[postburn]), col="blue")


pdf("~/repos/pnh-ms/R/output/scales/PhrynosOnlyAnalysis.pdf", height=8, width=12)
par(mfrow=c(2,2), mar=c(5,5,0.75, 0.25))
phenogram(pars2simmap(H.phryno, td$phy)$tr, simdat1$dat, colors = setNames(viridis::viridis(2), c(2,1)), fsize=1, ylab="FG fraction", xlab="Time (my)")

priorDens <- density(rbeta(length(postburn)*90, 0.8, 0.8), bw=0.005)
plot(priorDens$x, priorDens$y*2000, type="n", lty=2, ylim=c(0, 5000), xlab="Weight to Phryno Hypothesis", ylab="Density", yaxt="n")
lines(c(0, priorDens$x[priorDens$x>0 & priorDens$x<1], 1), c(0, priorDens$y[priorDens$x>0 & priorDens$x<1]*1500, 0), lty=2)
cols <- sapply(viridis::viridis(3), function(x) bayou:::makeTransparent(x, 50))
hist(1-chain.wPhryno$w[postburn], col=cols[1], add=TRUE, breaks=seq(0,1,0.025))
hist(1-chain.wPhryno.sim1$w[postburn], col=cols[2], add=TRUE, breaks=seq(0,1,0.025))
hist(1-chain.wPhryno.sim2$w[postburn], col=cols[3], add=TRUE, breaks=seq(0,1,0.025))
legend(0, 5000, legend = c("Prior","Scales data", "Simulated data, phrynos only", "Simulated data under PE Hypothesis"),
       pch=rep(22,4) , pt.bg=c("black", cols),lty=c(2,NA,NA,NA),  pt.cex=c(0,2,2,2))

L1 <- Lposterior(chain.wPhryno,tree, burnin=0.3)
L2 <- Lposterior(chain.wPhryno.sim2, tree, burnin=.3)
plot(L1$pp, L2$pp, ylab="Branch posterior prob. for simulated data (PE Hyp)", xlab="Branch posterior prob. for empirical data", 
     xlim=c(0,1), ylim=c(0,1), pch=21, bg="gray90", cex=1.5)
#text(0.85, 0.1, "Phrynosoma")
abline(0,1, lty=2)
plotSimmap.mcmc(chain.wPhryno, burnin=0.3, pp.cutoff=0.1,lwd=2, pal=colorRampPalette(viridis::viridis(2)[c(1, 2)]))
dev.off()


####################
###Try again with all 3. 
library(DirichletReg)



.dirichletMove <- function(cache=NULL, pars, d=NULL, move, ct=NULL, prior=NULL){
  prop <- MCMCpack::rdirichlet(1, pars[[move]]*d)
  if(any(prop==0)){
    lnHastingsRatio <- -Inf
    prop[prop==0] <- .Machine$double.eps
    prop <- prop/sum(prop)
  } else {
    lnHastingsRatio <- DirichletReg::ddirichlet(matrix(pars[[move]], nrow=1), prop, log=TRUE)  - DirichletReg::ddirichlet(matrix(prop, nrow=1), pars[[move]], log=TRUE)
  }
  pars.new <- pars
  pars.new[[move]] <- prop
  return(list(pars=pars.new, hr=lnHastingsRatio))
}

ddirichlet <- function(x, alpha, log=TRUE, sum.up=FALSE){
  DirichletReg::ddirichlet(matrix(x, nrow=1), alpha, log=log, sum.up=sum.up)
}
prop.W <- NULL
W <- NULL
pars$w <- matrix(rep(1,3)/3, nrow=1)
#for(i in 1:1000){
#  tmp <- .dirichletMove(pars=pars, move="w")
#  
#  if(runif(1) < exp(prior.wH3(tmp$pars) - prior.wH3(pars) + tmp$hr)){
#    pars <- tmp$pars
#  }
#  W <- rbind(W, pars$w)
#  prop.W <- rbind(prop.W, tmp$pars$w)
#}
H3.lik <- function(pars, cache, X, model="Custom"){
  n <- cache$n
  X <- cache$dat
  pred <- cache$pred
  EX.map1 <- bayou:::C_weightmatrix(cache, pars)$E
  pars_phryno <- list(alpha=pars$alpha, sig2=pars$sig2, k=H.phryno$k, ntheta=H.phryno$ntheta, theta=c(pars$thPhryno1, pars$thPhryno2), sb=H.phryno$sb, loc=H.phryno$loc, t2=H.phryno$t2)
  EX.map2 <- bayou:::C_weightmatrix(cache, pars_phryno)$E  
  pars_PE <- list(alpha=pars$alpha, sig2=pars$sig2, k=fixedpars$k, ntheta=fixedpars$ntheta, theta=c(pars$thPE1, pars$thPE2, pars$thPE3), sb=fixedpars$sb, loc=fixedpars$loc, t2=fixedpars$t2)
  EX.map3 <- bayou:::C_weightmatrix(cache, pars_PE)$E 
  X.c <- X - pars$w[1]*(EX.map1) - pars$w[2]*(EX.map2) - pars$w[3]*EX.map3
  ## This part adds the endothermy parameter to the theta for Mammal and Bird branches
  #dpars <- pars1
  #dpars$theta[dpars$t2[which(dpars$sb %in% c(857, 1719))]] <- dpars$theta[dpars$t2[which(dpars$sb %in% c(857, 1719))]]+dpars$endo
  ### The part below mostly does not change
  transf.phy <- bayou:::C_transf_branch_lengths(cache, 1, X.c, pars$alpha)
  transf.phy$edge.length[cache$externalEdge] <- transf.phy$edge[cache$externalEdge] + cache$SE[cache$phy$edge[cache$externalEdge, 2]]^2*(2*pars$alpha)/pars$sig2
  comp <- bayou:::C_threepoint(list(n=n, N=cache$N, anc=cache$phy$edge[, 1], des=cache$phy$edge[, 2], diagMatrix=transf.phy$diagMatrix, P=X.c, root=transf.phy$root.edge, len=transf.phy$edge.length))
  if(pars$alpha==0){
    inv.yVy <- comp$PP
    detV <- comp$logd
  } else {
    inv.yVy <- comp$PP*(2*pars$alpha)/(pars$sig2)
    detV <- comp$logd+n*log(pars$sig2/(2*pars$alpha))
  }
  llh <- -0.5*(n*log(2*pi)+detV+inv.yVy)
  #llh <- llh + gs.lik(c(pars$pred.sig2, pars$pred.root), root=ROOT.GIVEN) #$impute
  return(list(loglik=llh, theta=pars$theta,resid=X.c, comp=comp, transf.phy=transf.phy))
}

H3.monitor = function(i, lik, pr, pars, accept, accept.type, j){
  names <- c("gen", "lnL", "prior", "alpha","sig2", "rtheta", "rthPh","w", "k")
  string <- "%-8i%-8.2f%-8.2f%-8.2f%-8.2f%-8.2f%-8.2f%-8.2f%-8i"
  acceptratios <- tapply(accept, accept.type, mean)
  names <- c(names, names(acceptratios))
  if(j==0){
    cat(sprintf("%-7.7s", names), "\n", sep=" ")                           
  }
  cat(sprintf(string, i, lik, pr, pars$alpha, pars$sig2, pars$theta[1], pars$thPhryno1, pars$w[1], pars$k), sprintf("%-8.2f", acceptratios),"\n", sep="")
}

simpars1 <- bayou::priorSim(prior.phryno, td$phy, nsim=1)
simpars1$pars[[1]]$alpha <- 0.15
simpars1$pars[[1]]$sig2 <- 0.001
simpars1$pars[[1]]$theta <- c(0.6, 0.3)
simdat1 <- bayou::dataSim(simpars1$pars[[1]], model="OU", td$phy)

simpars2 <- bayou::priorSim(prior.PE, td$phy, nsim=1)
simpars2$pars[[1]]$alpha <- 0.15
simpars2$pars[[1]]$sig2 <- 0.001
simpars2$pars[[1]]$theta <- c(0.5, 0.7, 0.3)
simdat2 <- bayou::dataSim(simpars2$pars[[1]], model="OU", td$phy)

model.weightedH3 <- list(moves = list(alpha=".multiplierProposal", sig2=".multiplierProposal", 
                                          k=".splitmergePrior", theta=".adjustTheta", thPhryno1=".slidingWindowProposal", 
                                          thPhryno2=".slidingWindowProposal", thPE1=".slidingWindowProposal",
                                          thPE2=".slidingWindowProposal", thPE3=".slidingWindowProposal",
                                          w= ".dirichletMove", slide=".slide"
),
control.weights = list(alpha=5, sig2=3, theta=15, thPhryno1=3,
                       thPhryno2=3, thPE1=3, thPE2=3, thPE3=3,
                       k=10, w=3, slide=1),
D = list(alpha=0.5, sig2= 0.5, k=1, theta=1, thPhryno1=1, 
         thPhryno2=1, thPE1=1, thPE2=1, thPE3=1, w=100, slide=1),
parorder = c("alpha", "sig2", "w","thPhryno1", "thPhryno2","thPE1", "thPE2", "thPE3", "k", "ntheta",  "theta", "sb", "loc", "t2"),
rjpars = c("theta"),
shiftpars = c("sb", "loc", "t2"),
monitor.fn = H3.monitor,
lik.fn = H3.lik
)

prior.wH3 <- make.prior(td$phy, plot.prior = FALSE, 
                            dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", 
                                       dk="cdpois", dtheta="dnorm", dthPhryno1="dnorm",
                                       dthPhryno2="dnorm", dthPE1="dnorm", dthPE2="dnorm",
                                       dthPE3="dnorm",
                                       dw="ddirichlet", dsb="dsb"
                            ), 
                            param=list(dalpha=list(scale=0.1), dsig2=list(scale=0.1), 
                                       dk=list(lambda=0.5, kmax=10), dtheta=list(mean=0.5, sd=0.25),
                                       dthPhryno1=list(mean=0.5, sd=0.25), 
                                       dthPhryno2=list(mean=0.5, sd=0.25), 
                                       dthPE1=list(mean=0.5, sd=0.25), 
                                       dthPE2=list(mean=0.5, sd=0.25),
                                       dthPE3=list(mean=0.5, sd=0.25),
                                       dw=list(alpha=c(1,1,1)/3), dsb=list(bmax=1, prob=1)
                            )
                            #fixed=list(sb=startpar$sb, k=startpar$k)
)

pars <- list(alpha=0.15, sig2=0.001, k=1, ntheta=2, theta=c(0.55, 0.68), sb=c(25), loc=rep(0, 1), t2=2)
pars$thPhryno1 <- 0.55
pars$thPhryno2 <- 0.33
pars$thPE1 <- 0.55
pars$thPE2 <- 0.7
pars$thPE3 <- 0.3
pars$w <- matrix(c(1,1,1)/3, nrow=1)


mcmc.wH3 <- bayou.makeMCMC(td$phy, td[['FG.frac']], SE=0, prior=prior.wH3, model=model.weightedH3, startpar=pars,
                               new.dir="../output/scales/", outname="wH3_r001", plot.freq=NULL, samp=10, ticker.freq=1000) # Set up the MCMC

mcmc.wH3.sim1 <- bayou.makeMCMC(td$phy, simdat1$dat, SE=0, prior=prior.wH3, model=model.weightedH3, startpar=pars,
                                    new.dir="../output/scales/", outname="wH3sim1_r001", plot.freq=NULL, samp=10, ticker.freq=1000) # Set up the MCMC

mcmc.wH3.sim2 <- bayou.makeMCMC(td$phy, simdat2$dat, SE=0, prior=prior.wH3, model=model.weightedH3, startpar=pars,
                                    new.dir="../output/scales/", outname="wH3sim2_r001", plot.freq=NULL, samp=10, ticker.freq=1000) # Set up the MCMC


mcmc.wH3$run(200000)
mcmc.wH3.sim1$run(200000)
mcmc.wH3.sim2$run(200000)

chain.wH3 <- mcmc.wH3$load()
chain.wH3.sim1 <- mcmc.wH3.sim1$load()
chain.wH3.sim2 <- mcmc.wH3.sim2$load()

pdf("~/repos/pnh-ms/R/output/scales/H3weights.pdf", height=10, width=6)
par(mfrow=c(3,1))
postburn <- floor(0.3*length(chain.wH3$gen)):length(chain.wH3$gen)
wtmp <- do.call(rbind, chain.wH3$w[postburn])
cols <- sapply(viridis::viridis(3), function(x) bayou:::makeTransparent(x, 50))
hist(wtmp[,1], col=cols[1], breaks=seq(0,1,0.025), xlim=c(0,1), ylim=c(0,4000), main="Scales Data", xlab="Weight")
hist(wtmp[,2], col=cols[2], add=TRUE, breaks=seq(0,1,0.025))
hist(wtmp[,3], col=cols[3], add=TRUE, breaks=seq(0,1,0.025))
legend(0,4000, legend=c("rjMCMC", "Phrynos only", "PE"), pch=22, pt.bg=cols, pt.cex=2)

wtmp <- do.call(rbind, chain.wH3.sim1$w[postburn])
cols <- sapply(viridis::viridis(3), function(x) bayou:::makeTransparent(x, 50))
hist(wtmp[,1], col=cols[1], breaks=seq(0,1,0.025), xlim=c(0,1), ylim=c(0,4000), main="Simulated u/Phrynos Only", xlab="Weight")
hist(wtmp[,2], col=cols[2], add=TRUE, breaks=seq(0,1,0.025))
hist(wtmp[,3], col=cols[3], add=TRUE, breaks=seq(0,1,0.025))

wtmp <- do.call(rbind, chain.wH3.sim2$w[postburn])
cols <- sapply(viridis::viridis(3), function(x) bayou:::makeTransparent(x, 50))
hist(wtmp[,1], col=cols[1], breaks=seq(0,1,0.025), xlim=c(0,1), ylim=c(0,4000), main="Simulated u/PE", xlab="Weight")
hist(wtmp[,2], col=cols[2], add=TRUE, breaks=seq(0,1,0.025))
hist(wtmp[,3], col=cols[3], add=TRUE, breaks=seq(0,1,0.025))


pdf("~/repos/pnh-ms/R/output/scales/PhrynosOnlyAnalysis.pdf", height=8, width=12)
par(mfrow=c(2,2), mar=c(5,5,0.75, 0.25))
phenogram(pars2simmap(H.wH3, td$phy)$tr, simdat1$dat, colors = setNames(viridis::viridis(2), c(2,1)), fsize=1, ylab="FG fraction", xlab="Time (my)")

#priorDens <- density(rbeta(length(postburn)*90, 0.8, 0.8), bw=0.005)
#plot(priorDens$x, priorDens$y*2000, type="n", lty=2, ylim=c(0, 5000), xlab="Weight to Phryno Hypothesis", ylab="Density", yaxt="n")
lines(c(0, priorDens$x[priorDens$x>0 & priorDens$x<1], 1), c(0, priorDens$y[priorDens$x>0 & priorDens$x<1]*1500, 0), lty=2)
cols <- sapply(viridis::viridis(3), function(x) bayou:::makeTransparent(x, 50))
hist(1-chain.wPhryno$w[postburn], col=cols[1], add=TRUE, breaks=seq(0,1,0.025))
hist(1-chain.wPhryno.sim1$w[postburn], col=cols[2], add=TRUE, breaks=seq(0,1,0.025))
hist(1-chain.wPhryno.sim2$w[postburn], col=cols[3], add=TRUE, breaks=seq(0,1,0.025))
legend(0, 5000, legend = c("Prior","Scales data", "Simulated data, phrynos only", "Simulated data under PE Hypothesis"),
       pch=rep(22,4) , pt.bg=c("black", cols),lty=c(2,NA,NA,NA),  pt.cex=c(0,2,2,2))

L1 <- Lposterior(chain.wPhryno,tree, burnin=0.3)
L2 <- Lposterior(chain.wPhryno.sim2, tree, burnin=.3)
plot(L1$pp, L2$pp, ylab="Branch posterior prob. for simulated data (PE Hyp)", xlab="Branch posterior prob. for empirical data", 
     xlim=c(0,1), ylim=c(0,1), pch=21, bg="gray90", cex=1.5)
#text(0.85, 0.1, "Phrynosoma")
abline(0,1, lty=2)
plotSimmap.mcmc(chain.wPhryno, burnin=0.3, pp.cutoff=0.1,lwd=2, pal=colorRampPalette(viridis::viridis(2)[c(1, 2)]))
dev.off()