## The Felsenstein Problem
require(bayou)
require(treeplyr)
require(geiger)

felsTree <- function(n, stemL=1, lambda=1, b=c(1,1), d=c(0,0), seeds=NULL){
  if(length(lambda) <2){
    lambda <- c(lambda, lambda)
  }
  if(is.null(seeds)){
    seeds <- sample(1:10000, 2)
  }
  phy1 <- reorder(drop.extinct(sim.bdtree(b=b[1], d=d[1], stop="taxa", n=n, seed=seeds[1])), "postorder")
  phy2 <- reorder(drop.extinct(sim.bdtree(b=b[2], d=d[2], stop="taxa", n=n, seed=seeds[2])), "postorder")
  phy1 <- rescale(phy1, model="lambda", lambda[1])
  phy2 <- rescale(phy2, model="lambda", lambda[2])
  phy1$edge.length <- phy1$edge.length/max(branching.times(phy1))
  phy2$edge.length <- phy2$edge.length/max(branching.times(phy2))
  o <- list(phy1=phy1, phy2=phy2)
  phy2$tip.label <- paste("s", (n+1):(2*n), sep="")
  phy <- list(edge=NULL, Nnode=(n-1)*2 + 1, tip.label=c(phy1$tip.label, phy2$tip.label), edge.length=NULL)
  phy1$edge[phy1$edge > n] <- phy1$edge[phy1$edge > n]+(n+1)
  phy2$edge[phy2$edge > n] <- phy2$edge[phy2$edge > n]+2*n
  phy2$edge[phy2$edge <= n] <- phy2$edge[phy2$edge <=n]+n

  phy$edge <- rbind(phy1$edge, phy2$edge, c((2*n+1), (2*n+n+1)), c((2*n+1), (2*n+2)))
  phy$edge.length <- c(phy1$edge.length, phy2$edge.length, rep(stemL, 2))
  class(phy) <- class(o$phy1)
  phy <- reorder(phy, "postorder")
  plot(phy)
  return(phy)
}

phyI <- felsTree(10, stemL=0.5, lambda=1, seeds=1:2)
phyII <- felsTree(10, stemL=0.5, lambda=1, seeds=1:2)
phyIII <- felsTree(10, stemL=0.5, lambda=1, seeds=1:2)
phyIV <- felsTree(50, stemL=0.5, lambda=c(1, 0.75), b=c(100, 0.1), d=c(0, 0), seeds=1:2)
phyIV <- drop.tip(phyIV, 61:100)
plot(phyIV)

##Phy I - Felsenstein
#identifyBranches(phyI, 1)
pars <- list(alpha=1, sig2=1, k=1, ntheta=2, theta=c(0, 4), sb=38, loc=0.25, t2=2)
set.seed(8)
dat1 <- dataSim(pars, model="OU", phyI, phenogram=FALSE)
set.seed(2)
dat2 <- dataSim(pars, model="OU", phyI, phenogram=FALSE)
trI <- pars2simmap(pars, phyI)
cols <- setNames(viridis::viridis(4)[c(1,3)], 1:2)
cols2 <- setNames(viridis::viridis(4)[c(2,4)], 1:2)
rangey <- c(min(c(dat1$dat, dat2$dat)), max(c(dat1$dat, dat2$dat)))
#pr <- make.prior(phyI, dists=list(dk="cdpois"), param=list(dk=list(lambda=3,kmax=10)), model="OU")
#mcmc <- bayou.makeMCMC(phyI, dat1$dat, prior=pr)
#mcmc$run(50000)
#chain <- mcmc$load()
cmap <- contMap(phyI, dat2$dat, legend=FALSE)
cmap <- setMap(cmap, colors=c(viridis::viridis(5)[1:4], rep(viridis::viridis(10)[9:10], 1)))

cmap2 <- contMap(phyI, dat1$dat, legend=FALSE)
cmap2 <- setMap(cmap2, colors=c(viridis::viridis(5)[1:4], rep(viridis::viridis(10)[9:10], 1)))

names(cmap$tree$maps[[37]]) <- rep(352, length(cmap$tree$maps[[37]]))
names(cmap$tree$maps[[38]])[1:length(cmap$tree$maps[[38]])/2] <- c(floor(seq(352,425,length.out=length(cmap$tree$maps[[38]])/2)))

names(cmap2$tree$maps[[37]]) <- rep(352, length(cmap2$tree$maps[[37]]))
names(cmap2$tree$maps[[38]])[1:length(cmap2$tree$maps[[38]])/2] <- c(floor(seq(272,400,length.out=length(cmap2$tree$maps[[38]])/2)))


dyrange.c <- c(-1.5, 22.5); lwd.c <- 2; xrange.c <- c(-0.25, 1.4)

pdf("~/repos/pnh-ms/R/output/Figure1.pdf", width=11, height=4)
par(mfrow=c(1,8))
par(mar=c(5, 2, 0, 0))
#plot(c(0,2), c(-1, 5), type="n", yaxt="n", xaxt="n", bty="n", xlab="", ylab="", xlim=c(0,1.5), ylim=rangey)
#phenogram(trI$tree, dat2$dat, ftype="off", axes=list(time=c(0,2), trait=c(-1,5)), add=TRUE, colors=cols, lwd=2)
plot(cmap, legend=FALSE, ftype="off", ylim=yrange.c,  lwd=lwd.c, outline=FALSE, mar=c(5, 2, 0, 0))
text(0, 21.5, labels="I.", cex=2)
text(1, 21.5, labels="Trait Y (C)")
text(0.5, 0, labels="IC", cex=2)
.Lpp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
points(.Lpp$xx[22]/2, .Lpp$yy[22],pch="|", cex=2)

par(mar=c(5, 0, 0, 2))
#plot(c(2,0), c(-1, 5), type="n", yaxt="n", xaxt="n", bty="n", xlab="", ylab="", xlim=c(1.5,0), ylim=rangey)
#phenogram(trI$tree, dat1$dat, axes=list(time=c(2,0), trait=c(-1, 5)), ftype="off", add=TRUE, colors=cols2, lwd=2)
plot(cmap2, legend=FALSE, ftype="off", direction="leftwards", ylim=yrange.c, lwd=lwd.c, outline=FALSE, mar=c(5, 0, 0, 2))
text(0.5, 21.5, labels="Trait X (C)")
.Lpp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
points(.Lpp$xx[22]+.Lpp$xx[22]/4, .Lpp$yy[22],pch="|", cex=2)

##Phy II - OU models
par(mar=c(5, 2, 0, 0))
#plot(c(0,2), c(-1, 5), type="n", yaxt="n", xaxt="n", bty="n", xlab="", ylab="", xlim=c(0,1.5), ylim=rangey)
#phenogram(trI$tree, dat2$dat, ftype="off", axes=list(time=c(0,2), trait=c(-1,5)), add=TRUE, colors=cols, lwd=2)
plot(cmap, legend=FALSE, ftype="off", direction="rightwards",ylim=yrange.c, lwd=lwd.c, outline=FALSE, mar=c(5, 2, 0, 0))
text(0.1, 21.5, labels="II.", cex=2)
text(1, 21.5, labels="Trait Y (C)")
text(0.5, 0, labels="OU", cex=2)
.Lpp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
points(.Lpp$xx[22]/2, .Lpp$yy[22],pch="|", cex=2)

par(mar=c(5, 0, 0, 2))
plotRegimes(trI$tree, pal=function(x) rev(cols2), show.tip.label=FALSE, direction="leftwards", y.lim=c(-1.5, 22.5), lwd=2)
tiplabels(pch=21, bg=cols2[bayou:::.tipregime(pars, trI$tree)], cex=1.5, col=cols2[bayou:::.tipregime(pars, trI$tree)])
text(0.5, 21.5, labels="Trait X (D)")
.Lpp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
points(.Lpp$xx[22]+.Lpp$xx[22]/4, .Lpp$yy[22],pch="|", cex=2)


##Phy III - Darwin's scenario
par(mar=c(5, 2, 0, 0))

plotRegimes(trI$tree, pal=function(x) rev(cols), show.tip.label=FALSE, direction="rightwards", y.lim=c(-1.5, 22.5), lwd=2)
tiplabels(pch=21, bg=cols[bayou:::.tipregime(pars, trI$tree)], cex=1.5,  col=cols[bayou:::.tipregime(pars, trI$tree)])
text(0.1, 21.5, labels="III.", cex=2)
text(1, 21.5, labels="Trait Y (D)")
text(0.5, 0, labels="Pagel", cex=2)
.Lpp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
points(.Lpp$xx[22]/2, .Lpp$yy[22],pch="|", cex=2)

par(mar=c(5, 0, 0, 2))
plotRegimes(trI$tree, pal=function(x) rev(cols2), show.tip.label=FALSE, direction="leftwards", y.lim=c(-1.5, 22.5), lwd=2)
tiplabels(pch=21, bg=cols2[bayou:::.tipregime(pars, trI$tree)], cex=1.5,  col=cols2[bayou:::.tipregime(pars, trI$tree)])
text(0.5, 21.5, labels="Trait X (D)")
.Lpp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
points(.Lpp$xx[22]+.Lpp$xx[22]/4, .Lpp$yy[22],pch="|", cex=2)

##Phy IV - BiSSE
par(mar=c(5, 1, 0, 0))
parsIV <- list(alpha=1, sig2=1, k=1, ntheta=2, theta=c(0, 4), sb=118, loc=0.25, t2=2)
trIV <- pars2simmap(parsIV, phyIV)
plotRegimes(trIV$tree, pal=function(x) rev(cols), show.tip.label=FALSE, direction="rightwards", y.lim=c(-7, 67), lwd=2)
text(0.1, 65, labels="IV.", cex=2)
text(1, 65, labels="Diversification Y")
text(0.5, -2, labels="BiSSE", cex=2)
.Lpp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
points(.Lpp$xx[71]/2, .Lpp$yy[71],pch="|", cex=2)

#tiplabels(pch=21, bg=cols[bayou:::.tipregime(pars, trI$tree)], col=cols[bayou:::.tipregime(pars, trI$tree)])
par(mar=c(5, 0, 0, 1))
plotRegimes(trIV$tree, pal=function(x) rev(cols2), show.tip.label=FALSE, direction="leftwards", y.lim=c(-7, 67), lwd=2)
tiplabels(pch=21, bg=cols2[bayou:::.tipregime(parsIV, trIV$tree)], col=cols2[bayou:::.tipregime(parsIV, trIV$tree)])
text(0.5, 65, labels="Trait X (D)")
points(1.1+.Lpp$xx[71]/2, .Lpp$yy[71],pch="|", cex=2)

dev.off()




phy0 <- felsTree(20, lambda=0)
phy1 <- felsTree(20, lambda=1)

P0 <- P1 <- NULL
phy0$edge.length <- phy0$edge.length/(max(branching.times(phy0)))
phy1$edge.length <- phy1$edge.length/(max(branching.times(phy1)))

for(S in 10^seq(-2, 3, length.out=10)){
  for(j in 1:200){
  dat0 <- sim.char(phy0, par=matrix(c(1, 0, 0, 1), ncol=2), model="BM", root=0)[,,1]
  dat1 <- sim.char(phy1, par=matrix(c(1, 0, 0, 1), ncol=2), model="BM", root=0)[,,1]
  
  dat0.B <- dat0
  dat1.B <- dat1
  #S <- 1
  
  dev <- MASS::mvrnorm(1, mu=c(0,0), Sigma=matrix(c(S, 0*S, 0*S, S), ncol=2))
  dat0.B[1:20, ] <- dat0.B[1:20,] + cbind(rep(dev[1], 20), rep(dev[2], 20))
  dat1.B[1:20, ] <- dat1.B[1:20,] + cbind(rep(dev[1], 20), rep(dev[2], 20))
  #plot(dat2)

  pic0.1 <- pic(dat0.B[,1], phy0)
  pic0.2 <- pic(dat0.B[,2], phy0)
  
  pic1.1 <- pic(dat1.B[,1], phy1)
  pic1.2 <- pic(dat1.B[,2], phy1)

  #.pic1 <- pic(dat[,1], phy)
  #.pic2 <- pic(dat[,2], phy)
  #summary(lm(.pic2 ~ .pic1 -1))
  tmp0 <- summary(lm(pic0.2 ~ pic0.1 -1))
  tmp1 <- summary(lm(pic1.2 ~ pic1.1 -1))
  P0 <- rbind(P0, data.frame(burstSize=S, pValue=tmp0$coefficients[4]))
  P1 <- rbind(P1, data.frame(burstSize=S, pValue=tmp1$coefficients[4]))
  }
}

#j <- sample(1:100000, 1)
set.seed(40160)
phy <- phy0
S=10^1.5
dat <- sim.char(phy, par=matrix(c(1, 0, 0, 1), ncol=2), model="BM", root=0)[,,1]
dat2 <- dat
#S <- 1
dev <- MASS::mvrnorm(1, mu=c(0,0), Sigma=matrix(c(S, 0*S, 0*S, S), ncol=2))
dat2[1:20, ] <- dat2[1:20,] + cbind(rep(dev[1], 20), rep(dev[2], 20))
#plot(dat2)

pic1 <- pic(dat2[,1], phy)
pic2 <- pic(dat2[,2], phy)

#.pic1 <- pic(dat[,1], phy)
#.pic2 <- pic(dat[,2], phy)
#summary(lm(.pic2 ~ .pic1 -1))
tmp <- lm(pic2 ~ pic1 -1)
#P <- rbind(P, data.frame(burstSize=S, pValue=tmp$coefficients[4]))

pdf("~/repos/pnh-ms/R/output/felsenstein/FelsensteinPlot.pdf", height=6, width=8)
par(mfrow=c(2,2))
par(mar=c(0,0,0,0))
plot(phy, show.tip.label=FALSE)
L <- get("last_plot.phylo", envir = .PlotPhyloEnv)
points(L$xx[61]/2, L$yy[61], pch="|", cex=2)
par(mar=c(5, 5, 0.25,2))
cat <- c("white", "black")[as.numeric(1:40 > 20)+1]
tiplabels(pch=21, bg=cat)
#plot(dat, xlab="Trait X", ylab="Trait Y", pch=21, bg=cat)
plot(dat2, xlab="Trait X", ylab="Trait Y", pch=21, bg=cat)#
plot(pic1, pic2, xlab="PIC X", ylab="PIC Y", pch=21, bg=c("white", "black")[c(2, rep(1, length(pic1)))])
abline(tmp)
plot(log(unique(P0[,1]),10), tapply(P0[,2], P0[,1], function(x) sum(x < 0.05)/200),
     ylab="Proportion significant at p < 0.05", xlab="log10(Shift Variance) - log10(BM Variance)")
lines(log(unique(P1[,1]),10), tapply(P1[,2], P1[,1], function(x) sum(x < 0.05)/200), lty=2)
lines(log(unique(P0[,1]),10), tapply(P0[,2], P0[,1], function(x) sum(x < 0.05)/200), lty=1)
points(log(unique(P1[,1]),10), tapply(P1[,2], P1[,1], function(x) sum(x < 0.05)/200), pch=24)
abline(h=0.05, lty=3)
legend(-2, 0.85, legend=c("Clade polytomies", "Fully bifurcating"), lty=c(1,2), pch=c(1, 24))
dev.off()


#plot(pic1, pic2)

#par(mfrow=c(1,3))
#plot(phy)
#plot(dat, xlab="Trait X", ylab="Trait Y")
#plot(dat2, xlab="Trait X", ylab="Trait Y")#


P0 <- P1 <- NULL
phy0$edge.length <- phy0$edge.length/(max(branching.times(phy0)))
phy1$edge.length <- phy1$edge.length/(max(branching.times(phy1)))

for(S in 10^seq(-2, 3, length.out=10)){
  for(j in 1:200){
    #S <- 1
    dev <- MASS::mvrnorm(1, mu=c(0,0), Sigma=matrix(c(S, 0*S, 0*S, S), ncol=2))
    
    dat0 <- cbind(c(rnorm(length(phy0$tip.label)/2, 0, 1),rnorm(length(phy0$tip.label)/2, dev[1], 1)), 
                  c(rnorm(length(phy0$tip.label)/2, 0, 1),rnorm(length(phy0$tip.label)/2, dev[2], 1)))#sim.char(phy0, par=matrix(c(1, 0, 0, 1), ncol=2), model="BM", root=0)[,,1]
    dat1 <- cbind(c(rnorm(length(phy1$tip.label)/2, 0, 1),rnorm(length(phy1$tip.label)/2, dev[1], 1)), 
                  c(rnorm(length(phy1$tip.label)/2, 0, 1),rnorm(length(phy1$tip.label)/2, dev[2], 1)))#sim.char(phy1, par=matrix(c(1, 0, 0, 1), ncol=2), model="BM", root=0)[,,1]
    
    dat0.B <- dat0
    dat1.B <- dat1
    
    #dat0.B[1:20, ] <- dat0.B[1:20,] + cbind(rep(dev[1], 20), rep(dev[2], 20))
    #dat1.B[1:20, ] <- dat1.B[1:20,] + cbind(rep(dev[1], 20), rep(dev[2], 20))
    #plot(dat2)
    
    pic0.1 <- pic(dat0.B[,1], phy0)
    pic0.2 <- pic(dat0.B[,2], phy0)
    
    pic1.1 <- pic(dat1.B[,1], phy1)
    pic1.2 <- pic(dat1.B[,2], phy1)
    
    #.pic1 <- pic(dat[,1], phy)
    #.pic2 <- pic(dat[,2], phy)
    #summary(lm(.pic2 ~ .pic1 -1))
    tmp0 <- summary(lm(dat0.B[,2] ~ dat0.B[,1]))
    tmp1 <- summary(lm(dat1.B[,2] ~ dat1.B[,1]))
    P0 <- rbind(P0, data.frame(burstSize=S, pValue=tmp0$coefficients[4]))
    P1 <- rbind(P1, data.frame(burstSize=S, pValue=tmp1$coefficients[4]))
  }
}





