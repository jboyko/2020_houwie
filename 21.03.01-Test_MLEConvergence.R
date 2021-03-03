# simulation tests for optimizing the optimization of hOUwie 
simulateData <- function(phy, q, alpha, sigma.sq, quiet=FALSE){
  phy$edge.length <- phy$edge.length/max(branching.times(phy)) # following Cresler et al. (2015) we set T to be 1
  Q = matrix(c(-q,q,q,-q), 2, 2)
  
  root.p = c(0,0) # we will sample the root with equal probability of either state
  root.p[sample(c(1,2), 1, prob =c(0.5, 0.5))] <- 1
  
  theta = c(4, 5) # following Cresler et al. (2015) we set delta theta to be 1 (tempor at 2)
  theta0 = rnorm(1, theta[which(root.p == 1)], sqrt(sigma.sq[which(root.p == 1)]/2*alpha[which(root.p == 1)])) # following Beaulieu et al (2012) we sample the root theta from the stationary distribution matchiing the root state
  full.data <- hOUwie.sim(phy, Q, root.p, alpha, sigma.sq, theta0, theta)
  obs.no.trans <- sum(unlist(lapply(full.data$simmap[[1]]$maps, function(x) length(x) - 1)))
  if(!quiet){
    cat("The observed number of transitions was found to be", obs.no.trans, "\n")
  }
  return(full.data)
}

nLikEval <- function(p.subset, nEval, nSimmaps, tree, trait){
  p <- c(mk = NA, alpha = NA, sig2 = NA, opt1 = 4, opt2 = 5)
  p[1] <- p.subset[1]
  p[2] <- p.subset[2]
  p[3] <- p.subset[3]
  lliks <- sapply(1:nEval, function(x) hOUwie(phy = tree, data = trait, rate.cat = 1, model.cor = "ER", model.ou = "OUM", p = p, nSim = nSimmaps, quiet = FALSE)$loglik)
  return(lliks)
}

# get a phylogeny rescaled to a height of 1
source("~/2020_hOUwie/hOUwie.R")
#source("hOUwie.R")
require(OUwie)
require(corHMM)
require(parallel)
require(lhs)

# simulate data
nTip <- 250
q <- 1 # the expected number of markov transitions in terms of br proportion
alpha <- c(0.5,0.5) 
sigma.sq <- c(0.5,0.5) 
phy <- sim.bdtree(b = 1, d = 0.5, stop = "taxa", n = nTip) # this will be scaled to H=1
phy <- drop.extinct(phy)
phy$edge.length <- phy$edge.length/max(branching.times(phy))
data.houwie <- simulateData(phy, q, alpha, sigma.sq)[[1]]

# set up a range of params to explore based on latin hypercube
n.par.combos <- 250 # the resolution of the heat map
RNGkind("L'Ecuyer-CMRG")
set.seed(1985)
X <- randomLHS(n = n.par.combos, k = 3)
X <- X * 1
X[,1] <- X[,1] * 2
colnames(X) <- c("p.mk", "alpha", "sig2")
pairs(X)
par.set <- vector("list", dim(X)[1])
for(i in 1:dim(X)[1]){
  par.set[[i]] <- X[i,]
}

# evaluate 
n.evaluations <- 50 # the accuracy of our mean and sd estimates
for(n.simmaps in c(50, 100, 500, 1000)){
  # set.seed(1985)
  out <- mclapply(par.set, function(x) nLikEval(x, n.evaluations, n.simmaps, phy, data.houwie), mc.cores = 50)
  save(out, par.set, file = paste0("out-", n.simmaps, "simmaps.Rsave"))
}
# n.simmaps <- 100 # how good will our hOUwie model be


out <- mclapply(par.set, function(x) nLikEval(x, nEval = 5, nSimmaps = 10, phy, data.houwie), mc.cores = 1)

# graph

# load("2020_hOUwie/out-10simmaps.Rsave")


out.sd <- unlist(lapply(out, function(x) sd( x[!x == -Inf])))
out.mu <- unlist(lapply(out, function(x) mean( x[!x == -Inf])))
true.par <- c("p.mk" = q, "alpha" = alpha[1], "sig2" = sigma.sq[1])

res <- do.call(rbind, par.set)
res <- cbind(res, out.mu, out.sd)
colnames(res)[c(4,5)] <- c("mean", "sd")
res <- res[apply(res, 1, function(x) !any(is.na(x))),]
res <- as.data.frame(res)

require(akima)
resolution <- 0.1
par(mfrow=c(3,3))
for(i in 1:3){
  for(j in 1:3){
    # organize what we are plotting
    i.name <- colnames(res)[i]
    j.name <- colnames(res)[j]
    if(i == j){
      # the diagonal
      plot(1, type = "n",xlab = "", ylab = "", xlim = c(0, 2), ylim = c(0, 2), axes = FALSE)
      text(x = 1, y = 1, i.name, cex = 3)
    }
    if(i > j){
      # the lower triangle
      data <- data.frame(x=res[,j],
                         y=res[,i],
                         distance=-res$sd)
      a <- interp(x=data$x, y=data$y, z=data$distance, 
                  xo=seq(min(data$x),max(data$x),by=resolution), 
                  yo=seq(min(data$y),max(data$y),by=resolution), duplicate="mean")
      image(a, xlab = j.name, ylab = i.name, col = hcl.colors(100, "viridis"), main = "-sd(llik)")
      contour(a, add = TRUE, col = "black")
    }
    if(i < j){
      # the upper triangle
      data <- data.frame(x=res[,i],
                         y=res[,j],
                         distance=res$mean)
      a <- interp(x=data$x, y=data$y, z=data$distance, 
                  xo=seq(min(data$x),max(data$x),by=resolution), 
                  yo=seq(min(data$y),max(data$y),by=resolution), duplicate="mean")
      image(a, xlab = i.name, ylab = j.name, col = hcl.colors(100, "viridis"), main = "mean(llik)")
      contour(a, add = TRUE, col = "black")
      true.pts <- true.par[match(c(i.name, j.name), names(true.par))]
      points(x = true.pts[1], y = true.pts[2], cex = 3, col = "red", pch = 13)
    }
    else{
    }
  }
}


