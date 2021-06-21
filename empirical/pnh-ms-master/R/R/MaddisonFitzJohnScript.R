require(geiger)
require(bayou)
require(diversitree)
D <- NULL
for(i in 1:100){
  tree <- sim.bdtree(b=1, d=0, stop="taxa", n=20)
  tree <- reorder(tree, "postorder")
  tree$edge.length <- tree$edge.length/max(branching.times(tree))
  dat <- cbind(x=rep(0, 20), y=rep(0, 20))
  sb=sample(1:nrow(tree$edge), 1, replace=FALSE)
  pars <- list("sb"=sb, "loc"=runif(1, 0, tree$edge.length[sb]))

  pars$k <- 1; pars$ntheta=2;pars$t2 <-2
  dat <- data.frame(x=bayou:::.tipregime(pars, tree)-1,y=bayou:::.tipregime(pars, tree)-1)

  lik.u <- make.musse.multitrait(tree, dat)
  lik.u <- constrain(lik.u, lambda0~1, lambdax~1, lambday~1, lambdaxy~1, mu0~0, mux~0, muy~0, muxy~0)
  lik.c <- constrain(lik.u, qx01.y~ qx01.0, qx10.y~qx10.0,qy01.x~qy01.0,qy10.x~qy10.0)

  lik.h <- constrain(lik.u, qx01.y ~ 10000, qy01.x ~ 10000, qx10.0 ~ 0, qy10.0 ~ 0, qx10.y ~ 0, qy10.x ~ 0)

  hResult <- find.mle(lik.h, c(0.25, 0.25))
  cResult <- find.mle(lik.c, c(0.5, 0.5, 0.5, 0.5))

  diff <- cResult$lnLik - hResult$lnLik

  predDiff <- (log(tree$edge.length[pars$sb])-log(sum(tree$edge.length)))

  D <- rbind(D, "Emp"=diff, "Pred"=predDiff)

}









uRes <- NULL
for(i in 1:20){
  uResult <- find.mle(lik.u, runif(8, 0, 1), 
                    method="nlminb", control=list(upper=c(6,6,6,6,6,6,6,6),
                                                  lower=c(0,0,0,0,0,0,0,0)))
  uRes <- rbind(uRes, c("lnLik"=uResult$lnLik, uResult$par.full))
}
max(uRes[,1])
cResult$lnLik
max(uRes[,1])-cResult$lnLik
#uResult$lnLik

pagel <- phytools::fitPagel(tree, setNames(dat$x, tree$tip.label),
                            setNames(dat$y, tree$tip.label))

pagel
2*(log(tree$edge.length[pars$sb])-log(sum(tree$edge.length)))


