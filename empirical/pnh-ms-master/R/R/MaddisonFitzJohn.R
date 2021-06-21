require(geiger)
require(bayou)
require(diversitree)
D <- NULL
Q <- NULL
L <- NULL
ntip <- 20
for(i in 1:100){
  tree <- sim.bdtree(b=1, d=0, stop="taxa", n=ntip)
  tree <- reorder(tree, "postorder")
  tree$edge.length <- tree$edge.length/max(branching.times(tree))
  dat <- cbind(x=rep(0, ntip), y=rep(0, ntip))
  tipbranch <- which(tree$edge[,2] <= ntip)
  sb=sample((1:nrow(tree$edge))[-tipbranch], 1, replace=FALSE)
  pars <- list("sb"=sb, "loc"=runif(1, 0, tree$edge.length[sb]))

  pars$k <- 1; pars$ntheta=2;pars$t2 <-2
  dat <- data.frame(x=bayou:::.tipregime(pars, tree)-1,y=bayou:::.tipregime(pars, tree)-1)

  lik.u <- set.defaults(make.musse.multitrait(tree, dat), root.p=c(1,0,0,0), root=ROOT.GIVEN)
  lik.u <- constrain(lik.u, lambda0~1, lambdax~1, lambday~1, lambdaxy~1, mu0~0, mux~0, muy~0, muxy~0)
  lik.q <- constrain(lik.u, qx01.y~ qx01.0, qx10.y~qx10.0,qy01.x~qy01.0,qy10.x~qy10.0)
  lik.c <- constrain(lik.u, qx01.y~ qx01.0, qx10.y~qx10.0,qy01.x~qy01.0,qy10.x~qy10.0, qx10.0 ~ 0, qy10.0~0)
  lik.h <- constrain(lik.u, qx01.y ~ 10000, qy01.x ~ 10000, qx10.0 ~ 0, qy10.0 ~ 0, qx10.y ~ 0, qy10.x ~ 0, qx01.0 ~ qy01.0)
  hRes <- list(lnLik = -Inf)
  cRes <- list(lnLik = -Inf)
  uRes <- list(lnLik = -Inf)
  qRes <- list(lnLik = -Inf)
  for(j in 1:10){
    hResult <- find.mle(lik.h, method="optim", runif(1, 0, 0.5))
    cResult <- find.mle(lik.c, method="optim", runif(2, 0, 0.5))
    qResult <- find.mle(lik.q, runif(4, 0, 0.5), method="nlminb", lower=rep(0,4), upper=rep(1000, 4))
    uResult <- find.mle(lik.u, runif(8, 0, 0.5), method="nlminb", lower=rep(0,8), upper=rep(1000, 8))
    if(hRes$lnLik < hResult$lnLik) hRes <- hResult
    if(cRes$lnLik < cResult$lnLik) cRes <- cResult
    if(uRes$lnLik < uResult$lnLik) uRes <- uResult
    if(qRes$lnLik < qResult$lnLik) qRes <- qResult
  }
  diff <- cRes$lnLik - hRes$lnLik
  diff2 <- qRes$lnLik - uRes$lnLik

  predDiff <- (log(tree$edge.length[pars$sb])-log(sum(tree$edge.length)))
  Q <- rbind(Q, c(hRes$par, cRes$par))
  D <- rbind(D, c("Emp"=diff, "Pred"=predDiff, "Emp2"=diff2))
  L <- rbind(L, c(cRes$lnLik, hRes$lnLik, uRes$lnLik, qRes$lnLik))
  print(i)
}

ramp <- viridis::viridis(150)[c(1,5,10,15,20,25,30,35,40,45, 61:150)]
qs <- apply(Q, 1, max)
qi <- floor(qs/(max(qs))*100) 
ptcol <- sapply(ramp[qi], function(x) bayou:::makeTransparent(x, 175))#ifelse(qs < 2/sum(tree$edge.length), bayou:::makeTransparent("black", 100), bayou:::makeTransparent("white", 100)) #
lncol <-  sapply(ramp[qi], function(x) bayou:::makeTransparent(x, 200))


pdf("~/repos/pnh-ms/R/output/maddfitz/darwinsscenario_unconstrained.pdf", height=8, width=8)
par(mar=c(6,6,1,1))
plot(D[,3], D[,2], type="n", xlim=c(-10, 0), ylim=c(-10, 0), xlab="Empirical difference in likelihood", ylab="Predicted difference in likelihood", cex.lab=1.5)
abline(0,1, lty=2)
points(D[,3], D[,2],pch=21, col=1, bg=ptcol, cex=2)

legend_image <- as.raster(matrix(ramp, ncol=1))
text(-9.25, 0.1, labels="max(Q)")
text(x=-9, y = -0.05+c(-0.14,-0.60, -1.06, -1.52, -1.98), labels = round(seq(0.09, 1.25, l=5),2))
text(x=-9.5, y = -0.1+c(-0.14,-0.60, -1.06, -1.52, -1.98), labels = "-", cex=1.5)
rasterImage(legend_image, -10, -2.1, -9.5,-0.1)
invT <- -0.1 - (c(0.1362242)/max(qs)*100)/50
lines(c(-10, -8), c(invT, invT), lty=3)
text(-7.75, invT, labels="1/T")
dev.off()

pdf("~/repos/pnh-ms/R/output/maddfitz/darwinsscenario_constrained.pdf", height=8, width=8)
plot(D[,1], D[,2], type="n", xlim=c(-10, 0), ylim=c(-10, 0), xlab="Empirical difference in likelihood", ylab="Predicted difference in likelihood", cex.lab=1.5)
abline(0,1, lty=2)
points(D[,1], D[,2],pch=21, col=1, bg=ptcol, cex=2)
legend_image <- as.raster(matrix(ramp, ncol=1))
text(-9.25, 0.1, labels="max(Q)")
text(x=-9, y = -0.05+c(-0.14,-0.60, -1.06, -1.52, -1.98), labels = round(seq(0.09, 1.25, l=5),2))
text(x=-9.5, y = -0.1+c(-0.14,-0.60, -1.06, -1.52, -1.98), labels = "-", cex=1.5)
rasterImage(legend_image, -10, -2.1, -9.5,-0.1)
invT <- -0.1 - (c(0.1362242)/max(qs)*100)/50
lines(c(-10, -8), c(invT, invT), lty=3)
text(-7.75, invT, labels="1/T")


dev.off()


plot(apply(Q, 1, max), D[,1]-D[,2])

tmp <- lm(abs(D[,1]-D[,2])~apply(Q, 1, max))
summary(tmp)


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


Trs <- NULL
for(i in 1:1000){
  tree <- sim.bdtree(b=1, d=0, stop="taxa", n=ntip)
  tree <- reorder(tree, "postorder")
  tree$edge.length <- tree$edge.length/max(branching.times(tree))
  Trs <- c(Trs, sum(tree$edge.length))
}

