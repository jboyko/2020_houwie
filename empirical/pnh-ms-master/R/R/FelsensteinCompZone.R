## The Felsenstein Problem
require(bayou)
require(treeplyr)
require(geiger)

felsTree <- function(n, stemL=1, lambda=1){
  phy1 <- reorder(sim.bdtree(b=1, d=0, stop="taxa", n=n), "postorder")
  phy2 <- reorder(sim.bdtree(b=1, d=0, stop="taxa", n=n), "postorder")
  phy1 <- rescale(phy1, model="lambda", lambda)
  phy2 <- rescale(phy1, model="lambda", lambda)
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
text(0.2, 39, labels="A.", cex=2)
L <- get("last_plot.phylo", envir = .PlotPhyloEnv)
points(L$xx[61]/2, L$yy[61], pch="|", cex=2)
par(mar=c(5, 4.5, 0.25,1))
cat <- c("white", "black")[as.numeric(1:40 > 20)+1]
tiplabels(pch=21, bg=cat)
#plot(dat, xlab="Trait X", ylab="Trait Y", pch=21, bg=cat)
plot(dat2, xlab="Trait X", ylab="Trait Y", pch=21, bg=cat, cex.lab=1.25)#
text(-7.85, 0.6, 2, labels="B.", cex=2)
plot(pic1, pic2, xlab="PIC X", ylab="PIC Y", pch=21, bg=c("white", "black")[c(2, rep(1, length(pic1)))], cex.lab=1.25)
text(-0.7, 2.7, 2, labels="C.", cex=2)
abline(tmp)
plot(log(unique(P0[,1]),10), tapply(P0[,2], P0[,1], function(x) sum(x < 0.05)/200),
     ylab="Prop. significant at p < 0.05", xlab="log10(Shift Variance) - log10(BM Variance)", cex.lab=1.25)
text(-1, 0.835, 2, labels="D.", cex=2)
lines(log(unique(P1[,1]),10), tapply(P1[,2], P1[,1], function(x) sum(x < 0.05)/200), lty=2)
lines(log(unique(P0[,1]),10), tapply(P0[,2], P0[,1], function(x) sum(x < 0.05)/200), lty=1)
points(log(unique(P1[,1]),10), tapply(P1[,2], P1[,1], function(x) sum(x < 0.05)/200), pch=24)
abline(h=0.05, lty=3)
legend(-2, 0.77, legend=c("Clade polytomies", "Fully bifurcating"), lty=c(1,2), pch=c(1, 24))


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





