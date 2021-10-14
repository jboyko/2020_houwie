# utility functions

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Model simulating and fitting functions
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

ModelSetRun <- function(phy, pars, model.name, sim.index.cor, sim.index.ou, sim.rate.cat, fit.index.cor, fit.index.ou, fit.rate.cat, nmaps, iter, ub.ou=NULL, nstarts=1){
  nModels <- length(fit.index.ou)
  nState <- dim(sim.index.cor)[1]
  cat("Begining single run for", nModels, "models...\n")
  data.houwie <- generateData(phy, sim.index.cor, sim.index.ou, pars)
  while(length(unique(data.houwie$data[,2])) != nState){
    data.houwie <- generateData(phy, index.cor, index.ou, pars)
  }
  if(sim.rate.cat == 1){
    data <- data.houwie$data
  }
  if(sim.rate.cat == 2){
    data <- data.houwie$data
    data[data[,2] == 3, 2] <- 1
    data[data[,2] == 4, 2] <- 2
    data[,2] <- as.factor(as.numeric(data[,2]))
  }
  if(sim.rate.cat == 3){
    data <- data.houwie$data
    data[data[,2] == 3, 2] <- 1
    data[data[,2] == 5, 2] <- 1
    data[data[,2] == 4, 2] <- 2
    data[data[,2] == 6, 2] <- 2
    data[,2] <- as.factor(as.numeric(data[,2]))
  }
  obj <- fit.hOUwie.set(phy=phy, data=data, rate.cat=fit.rate.cat,
                        index.cor=fit.index.cor, root.p="yang", lb.cor=1e-3, 
                        ub.cor=10, index.ou=fit.index.ou, root.station=FALSE, 
                        get.root.theta=FALSE, mserr="none", ub.ou = ub.ou, 
                        nSim=nmaps, nCores=1, prune=FALSE)
  out <- list(simulating.data = data.houwie,
              model.fits = obj)
  file.name <- paste0("/space_2/jamesboyko/2020_hOUwie/ModelTesting/", model.name, "/Sim=", model.name, "-", "Pars=", paste(pars, collapse = "_"), "-nTip=", length(phy$tip.label), "-nMap=", nmaps,"-Iter=", iter, ".Rsave")
  save(out, file = file.name)
  cat("Done.\n")
  return(out)
}

getTotalErrorRate <- function(Rsaves, Type="Error"){
  TipPars <- WeightedParsTip <- WeightedParsBr <- WeightedParsRt <-AICTables <- ParTables <- list()
  for(i in 1:(length(Rsaves)/2)){
    if(length(grep(paste0("Iter=", i, ".Rsave"), Rsaves)) == 0){
      next
    }
    ToLoad <- Rsaves[grep(paste0("Iter=", i, ".Rsave"), Rsaves)]
    load(ToLoad[1])
    load(ToLoad[2])
    Simulators <- out$simulating.data
    WeightedParsTip[[i]] <- getMapWeightedPars(Simulators, "Tip")
    WeightedParsBr[[i]] <- getMapWeightedPars(Simulators, "Branch")
    WeightedParsRt[[i]] <- getMapWeightedPars(Simulators, "Root")
    AICTables[[i]] <- SumSet$model.table
    ParTables[[i]] <- SumSet$model.pars
  }
  
  if(Type == "Error"){
    AvgSignErrorTip <- matrix(NA, length(ParTables), 3, dimnames = list(rep("Tip", length(ParTables)), c("alpha", "sigma.sq", "theta"))) 
    AvgSignErrorBr <- matrix(NA, length(ParTables), 3, dimnames = list(rep("Branch", length(ParTables)), c("alpha", "sigma.sq", "theta")))
    AvgSignErrorRt <- matrix(NA, length(ParTables), 3, dimnames = list(rep("Total", length(ParTables)), c("alpha", "sigma.sq", "theta")))
    for(i in 1:length(ParTables)){
      if(is.null(WeightedParsTip[[i]])){
        next
      }
      AvgSignErrorTip[i,] <- colMeans(WeightedParsTip[[i]][,c(2,3,4)] - ParTables[[i]][,c(1,2,3)])
      AvgSignErrorBr[i,] <- colMeans(WeightedParsBr[[i]][,c(2,3,4)] - ParTables[[i]][,c(1,2,3)])
      AvgSignErrorRt[i,] <- colMeans(WeightedParsRt[[i]][,c(2,3,4)] - ParTables[[i]][,c(1,2,3)])
    }
    AvgSignErrorTip <- AvgSignErrorTip[!apply(AvgSignErrorTip, 1, function(x) any(is.na(x))),]
    MeltedAvgSignErrorTip <- melt(AvgSignErrorTip)
    AvgSignErrorBr <- AvgSignErrorBr[!apply(AvgSignErrorBr, 1, function(x) any(is.na(x))),]
    MeltedAvgSignErrorBr <- melt(AvgSignErrorBr)
    AvgSignErrorRt <- AvgSignErrorRt[!apply(AvgSignErrorRt, 1, function(x) any(is.na(x))),]
    MeltedAvgSignErrorRt <- melt(AvgSignErrorRt)
    TotalErrorRt <- rbind(MeltedAvgSignErrorTip, MeltedAvgSignErrorBr, MeltedAvgSignErrorRt)
    return(TotalErrorRt)
  }
  if(Type == "Raw"){
    AvgSignErrorTip <- list()
    AvgSignErrorBr <-list()
    AvgSignErrorRt <- list()
    for(i in 1:length(ParTables)){
      if(is.null(WeightedParsTip[[i]])){
        next
      }
      AvgSignErrorTip[[i]] <- cbind("Tip", cbind(WeightedParsTip[[i]][,c(2,3,4)], ParTables[[i]][,c(1,2,3)]))
      AvgSignErrorBr[[i]] <- cbind("Br", cbind(WeightedParsBr[[i]][,c(2,3,4)], ParTables[[i]][,c(1,2,3)]))
      AvgSignErrorRt[[i]] <- cbind("Tot", cbind(WeightedParsRt[[i]][,c(2,3,4)], ParTables[[i]][,c(1,2,3)]))
      colnames(AvgSignErrorTip[[i]]) <- colnames(AvgSignErrorBr[[i]]) <- colnames(AvgSignErrorRt[[i]]) <- c("Type","E_alpha", "E_sigma.sq", "E_optim", "P_alpha", "P_sigma.sq", "P_optim")
    }
    
    AvgValueRt <- do.call(rbind, AvgSignErrorRt)
    AvgValueBr <- do.call(rbind, AvgSignErrorBr)
    AvgValueTip <- do.call(rbind, AvgSignErrorTip)
    return(rbind(AvgValueRt, AvgValueBr, AvgValueTip))
  }
}

SingleModelTestRun <- function(phy, index.cor, index.ou, model.cor, model.ou, rate.cat, pars, nmaps, iter, ub.ou=NULL, nstarts=1){
  cat("Begining single run...\n")
  nState <- dim(index.cor)[1]
  nObs <- nState/rate.cat
  data.houwie <- generateData(phy, index.cor, index.ou, pars)
  while(length(unique(data.houwie$data[,2])) != nState){
    data.houwie <- generateData(phy, index.cor, index.ou, pars)
  }
  if(rate.cat == 2){
    data <- data.houwie$data
    data[data[,2] == 3, 2] <- 1
    data[data[,2] == 4, 2] <- 2
    data[,2] <- as.factor(as.numeric(data[,2]))
  }else{
    data <- data.houwie$data
  }
  starting.vals <- generateStartingValues(phy, index.cor, index.ou, data.houwie$data, nstarts)
  cat("\nBegining two step...\n")
  TwoStepFit <- list()
  TwoStepFit$loglik <- -Inf
  for(i in 1:dim(starting.vals)[1]){
    tmp <- fitTwoStep(phy=phy, data=data, rate.cat=rate.cat, 
                      model.cor=model.cor, root.p="yang", lb.cor=1e-3, ub.cor=21,
                      model.ou=model.ou, root.station=FALSE, get.root.theta=FALSE, mserr="none", ub.ou = ub.ou, ip = starting.vals[i,])
    if(tmp$loglik > TwoStepFit$loglik){
      TwoStepFit <- tmp
    }
  }
  cat("\nBegining non cens...\n")
  NonCensFit <- list()
  NonCensFit$loglik <- -Inf
  for(i in 1:dim(starting.vals)[1]){
    tmp <- fitNonCensored(phy=phy, data=data, rate.cat=rate.cat, 
                          model.cor=model.cor, root.p="yang", lb.cor=1e-3, ub.cor=21,
                          model.ou=model.ou, root.station=FALSE, get.root.theta=FALSE, mserr="none", ub.ou = ub.ou, ip = starting.vals[i,],
                          nSim=nmaps, nCores=1)
    if(tmp$loglik > NonCensFit$loglik){
      NonCensFit <- tmp
    }
  }
  cat("\nBegining houwie...\n")
  hOUwieFit <- list()
  hOUwieFit$loglik <- -Inf
  for(i in 1:dim(starting.vals)[1]){
    start.vals_i <- generateStartingValues.houwie(starting.vals[i,], data.houwie, index.ou, index.cor)
    tmp <- hOUwie(phy=phy, data=data, rate.cat=rate.cat,
                  index.cor=index.cor, root.p="yang", lb.cor=1e-3, ub.cor=21,
                  index.ou=index.ou, root.station=FALSE, get.root.theta=FALSE, mserr="none", ub.ou = ub.ou, ip = start.vals_i,
                  nSim=nmaps, nCores=1)
    if(tmp$loglik > hOUwieFit$loglik){
      hOUwieFit <- tmp
    }
  }
  cat("\nBegining truth\n")
  TruMapFit <- list()
  TruMapFit$loglik <- -Inf
  for(i in 1:dim(starting.vals)[1]){
    starting.vals_i <- starting.vals[i,]
    tmp <- OUwie(phy=data.houwie$simmap[[1]], data=data.houwie$data,
                 model=model.ou, root.station=FALSE, get.root.theta=FALSE, mserr="none", ub = ub.ou, starting.vals = starting.vals_i[2:length(starting.vals_i)],
                 simmap.tree=TRUE, scaleHeight=FALSE, algorithm="three.point", diagn=FALSE)
    if(tmp$loglik > TruMapFit$loglik){
      TruMapFit <- tmp
    }
  }
  
  obj <- list(data.houwie = data.houwie,
              TwoStepFit = TwoStepFit,
              NonCensFit  = NonCensFit,
              hOUwieFit = hOUwieFit,
              TruMapFit = TruMapFit)
  cat("Done.\n")
  
  file.name <- paste0("Fit=", model.cor, "_", model.ou, "-", "Pars=", paste(pars, collapse = "_"), "-Iter=", iter, ".Rsave")
  print(file.name)
  save(obj, file = file.name)
  return(obj)
}

# for evaluating parsiimony
SingleModelTestRunParsimony <- function(phy, index.cor, index.ou, model.cor, model.ou, rate.cat, pars, nmaps, iter, ub.ou=NULL, nstarts=1){
  cat("Begining single run...\n")
  nState <- dim(index.cor)[1]
  nObs <- nState/rate.cat
  data.houwie <- generateData(phy, index.cor, index.ou, pars)
  while(length(unique(data.houwie$data[,2])) != nState){
    data.houwie <- generateData(phy, index.cor, index.ou, pars)
  }
  if(rate.cat == 2){
    data <- data.houwie$data
    data[data[,2] == 3, 2] <- 1
    data[data[,2] == 4, 2] <- 2
    data[,2] <- as.factor(as.numeric(data[,2]))
  }else{
    data <- data.houwie$data
  }
  starting.vals <- generateStartingValues(phy, index.cor, index.ou, data.houwie$data, nstarts)

  cat("\nBegining houwie normal...\n")
  hOUwieFit <- list()
  hOUwieFit$loglik <- -Inf
  for(i in 1:dim(starting.vals)[1]){
    start.vals_i <- generateStartingValues.houwie(starting.vals[i,], data.houwie, index.ou, index.cor)
    tmp <- hOUwie(phy=phy, data=data, rate.cat=rate.cat,
                  index.cor=index.cor, root.p="yang", lb.cor=1e-3, ub.cor=21,
                  index.ou=index.ou, root.station=FALSE, get.root.theta=FALSE, mserr="none", ub.ou = ub.ou, ip = start.vals_i,
                  nSim=nmaps, nCores=1)
    if(tmp$loglik > hOUwieFit$loglik){
      hOUwieFit <- tmp
    }
  }
  
  cat("\nBegining houwie parsimony \n")
  hOUwieFitPars <- list()
  hOUwieFitPars$loglik <- -Inf
  for(i in 1:dim(starting.vals)[1]){
    start.vals_i <- generateStartingValues.houwie(starting.vals[i,], data.houwie, index.ou, index.cor)
    tmp <- hOUwie(phy=phy, data=data, rate.cat=rate.cat,
                  index.cor=index.cor, root.p="yang", lb.cor=1e-3, ub.cor=21,
                  index.ou=index.ou, root.station=FALSE, get.root.theta=FALSE, mserr="none", ub.ou = ub.ou, ip = start.vals_i,
                  nSim=nmaps, nCores=1, parsimony = TRUE)
    if(tmp$loglik > hOUwieFitPars$loglik){
      hOUwieFitPars <- tmp
    }
  }
  
  
  cat("\nBegining truth\n")
  TruMapFit <- list()
  TruMapFit$loglik <- -Inf
  for(i in 1:dim(starting.vals)[1]){
    starting.vals_i <- starting.vals[i,]
    tmp <- OUwie(phy=data.houwie$simmap[[1]], data=data.houwie$data,
                 model=model.ou, root.station=FALSE, get.root.theta=FALSE, mserr="none", ub = ub.ou, starting.vals = starting.vals_i[2:length(starting.vals_i)],
                 simmap.tree=TRUE, scaleHeight=FALSE, algorithm="three.point", diagn=FALSE)
    if(tmp$loglik > TruMapFit$loglik){
      TruMapFit <- tmp
    }
  }
  
  obj <- list(data.houwie = data.houwie,
              hOUwieFitPars  = hOUwieFitPars,
              hOUwieFit = hOUwieFit,
              TruMapFit = TruMapFit)
  cat("Done.\n")
  
  file.name <- paste0("Fit=", model.cor, "_", model.ou, "-", "Pars=", paste(pars, collapse = "_"), "-Iter=", iter, ".Rsave")
  print(file.name)
  save(obj, file = file.name)
  return(obj)
}


# a function that given a model structure and params, generates a hOUwie dataset.
generateData <- function(phy, index.cor, index.ou, pars, type ="CD", quiet=FALSE){
  phy$edge.length <- phy$edge.length/max(branching.times(phy)) # ensure tree height is 1
  k.cor <- max(index.cor, na.rm = TRUE) # number of corhmm params
  k.ou <- max(index.ou, na.rm = TRUE) # number of ouwie params
  
  if((k.cor + k.ou) != length(pars)){
    stop("Length of pars does not match index matrices.", call. = FALSE)
  }
  
  p.mk <- pars[1:k.cor]
  p.ou <- pars[(k.cor+1):length(pars)]
  
  # organize corhmm params
  index.cor[index.cor == 0] <- NA
  index.cor[is.na(index.cor)] <- max(index.cor, na.rm = TRUE) + 1
  Q <- matrix(0, dim(index.cor)[1], dim(index.cor)[1])
  Q[] <- c(p.mk, 0)[index.cor]
  diag(Q) <- -rowSums(Q)
  root.p = rep(0, dim(Q)[1]) # we will sample the root with equal probability of either state
  root.p[sample(1:dim(Q)[1], 1)] <- 1
  
  # organize ou params
  Rate.mat <- matrix(1, 3, dim(index.cor)[2])
  index.ou[is.na(index.ou)] <- max(index.ou, na.rm = TRUE) + 1
  Rate.mat[] <- c(p.ou, 1e-10)[index.ou]
  alpha = Rate.mat[1,]
  sigma.sq = Rate.mat[2,]
  theta = Rate.mat[3,]
  # theta0 = rnorm(1, theta[which(root.p == 1)], sqrt(sigma.sq[which(root.p == 1)]/2*alpha[which(root.p == 1)])) # following Beaulieu et al (2012) we sample the root theta from the stationary distribution matchiing the root state
  theta0 = theta[which(root.p == 1)]
  full.data <- hOUwie.sim(phy, Q, root.p, alpha, sigma.sq, theta0, theta)
  # if(type == "CD"){
  #   full.data <- hOUwie.sim(phy, Q, root.p, alpha, sigma.sq, theta0, theta)
  # }else{
  #   full.data <- hOUwie.sim(phy, Q, root.p, alpha, sigma.sq, theta0, theta)
  #   dat.cor <- rTraitDisc(phy, Q, states = 1:dim(Q)[1], root.value = sample(1:dim(Q)[1], 1, prob = root.p))
  #   while(!all(levels(dat.cor) %in% dat.cor)){
  #     dat.cor <- rTraitDisc(phy, Q, states = 1:dim(Q)[1], root.value = sample(1:dim(Q)[1], 1, prob = root.p))
  #   }
  #   full.data$data <- cbind(full.data$data, cid.reg = dat.cor)
  # }
  obs.no.trans <- sum(unlist(lapply(full.data$simmap[[1]]$maps, function(x) length(x) - 1)))
  if(!quiet){
    cat("The observed number of transitions was found to be", obs.no.trans, "\n")
  }
  simulators <- list(index.cor = index.cor, index.ou = index.ou, pars = pars)
  return(c(full.data, simulators))
}

generateStartingValues <- function(phy, index.cor, index.ou, data, nstarts){
  x <- data[,3]
  names(x) <- data[,1]
  k.cor <- max(index.cor, na.rm = TRUE)
  k.ou <- ifelse(any(is.na(index.ou)), 1, 2)
  ip.cor <- rep(10/sum(phy$edge.length), k.cor)
  ip.ou <- c(mean(pic(x, phy)^2), log(2)/max(branching.times(phy)))[k.ou:1]
  ip <- c()
  if(nstarts > 1){
    ip <- abs(sapply(c(ip.cor, ip.ou), function(y) rnorm(nstarts-1, y, 2*y)))
  }
  ip <- rbind(c(ip.cor, ip.ou), ip)
  return(ip)
}

generateStartingValues.houwie <- function(starting.vals, data.houwie, index.ou, index.cor){
  means.by.regime <- with(data.houwie$data, tapply(data.houwie$data[,3], data.houwie$data[,2], mean))
  k.cor <- max(index.cor, na.rm = TRUE)
  k.alpha <- length(na.omit(unique(index.ou[1,])))
  k.sigma <- length(na.omit(unique(index.ou[2,])))
  k.theta <- length(na.omit(unique(index.ou[3,])))
  ip.cor <- rep(starting.vals[1], k.cor)
  if(length(means.by.regime) != k.theta){
    ip.theta <- rep(mean(means.by.regime), k.theta)
  }else{
    ip.theta <- means.by.regime
  }
  if(k.alpha == 0){
    ip.alpha <- c()
    ip.sigma <- rep(starting.vals[2], k.sigma)
  }else{
    ip.alpha <- rep(starting.vals[2], k.alpha)
    ip.sigma <- rep(starting.vals[3], k.sigma)
  }
  start.vals <- c(ip.cor, ip.alpha, ip.sigma, ip.theta)
  return(start.vals)
}

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Model evaluation functions
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

getBigObj <- function(sub.folders){
  out <- list()
  count <- 1
  for(sub.folder in sub.folders){
    Rsaves <- dir(sub.folder, full.names = TRUE)
    big.obj <- list()
    for(j in 1:length(Rsaves)){
      load(Rsaves[j])
      big.obj[[j]] <- obj
    }
    out[[count]] <- big.obj
    count <- count + 1
  }
  return(out)
}

getBigObj2 <- function(sub.folders){
  big.list <- list()
  count <- 1
  for(sub.folder in sub.folders){
    Rsaves <- dir(sub.folder, full.names = TRUE)
    big.obj <- list()
    for(j in 1:length(Rsaves)){
      load(Rsaves[j])
      big.obj[[j]] <- out
    }
    big.list[[count]] <- big.obj
    count <- count + 1
  }
  return(big.list)
}

getTrueParsVector <- function(truepars, index.ou){
  if(any(index.ou[1,] > index.ou[2,])){
    index.ou[1,][index.ou[1,] > index.ou[2,]] <- 0
  }
  pars.count <- c(alpha=0, sigma=0, theta=0)
  out.vec <- c()
  for(i in 1:max(index.ou, na.rm = TRUE)){
    index_i <- index.ou == i
    name.par <- rownames(index_i)[apply(index_i, 1, any)]
    pars.count[apply(index_i, 1, any)] <- pars.count[apply(index_i, 1, any)] + 1
    name.par <- paste(name.par, pars.count[apply(index_i, 1, any)], sep = "_")
    vec_i <- truepars$pars.ou[index_i][1]
    names(vec_i) <- name.par
    out.vec <- c(out.vec, vec_i)
  }
  return(out.vec)
}

getParTable <- function(index.ou, big.obj){
  # determine how many of each parameter we have
  if(any(index.ou[1,] > index.ou[2,])){
    index.ou[1,][index.ou[1,] > index.ou[2,]] <- 0
  }
  pars.count <- c(alpha=0, sigma=0, theta=0)
  table <- c()
  for(i in 1:max(index.ou, na.rm = TRUE)){
    index_i <- index.ou == i
    name.par <- rownames(index_i)[apply(index_i, 1, any)]
    pars.count[apply(index_i, 1, any)] <- pars.count[apply(index_i, 1, any)] + 1
    name.par <- paste(name.par, pars.count[apply(index_i, 1, any)], sep = "_")
    name.model <- names(big.obj[[1]])[-1]
    table_i <- c()
    for(j in 1:length(name.model)){
      model.fits_j <- lapply(big.obj, "[[", name.model[j])
      if(length(grep("solution", names(model.fits_j[[1]]))) == 1){
        table_j <- data.frame(par = name.par, model = name.model[j], 
                              value = unlist(lapply(model.fits_j, function(x) x$solution[index_i][1])))
      }else{
        table_j <- data.frame(par = name.par, model = name.model[j], 
                              value = unlist(lapply(model.fits_j, function(x) x$solution.ou[index_i][1])))
      }
      table_i <- rbind(table_i, table_j)
    }
    table <- rbind(table, table_i)
  }
  return(table)
}

# from a simmap create an ancestral state reconstruciton at nodes
getReconTableFromSimmap <- function(simmap){
  nNodes <- max(simmap$edge)
  nTips <- length(simmap$tip.label)
  nStates <- dim(simmap$mapped.edge)[2]
  recon_table <- matrix(0, nNodes, nStates, dimnames = list(1:nNodes, 1:nStates))
  for(i in 1:nTips){
    edge_i <- simmap$maps[[which(simmap$edge[,2] == i)]]
    recon_table[i, as.numeric(names(edge_i[length(edge_i)]))] <- 1
  }
  for(i in (nTips+1):nNodes){
    edge_i <- simmap$maps[[which(simmap$edge[,1] == i)[1]]]
    recon_table[i, as.numeric(names(edge_i[1]))] <- 1
  }
  return(recon_table)
}


organizeGeneratingData <- function(simulators){
  k.cor <- max(simulators$index.cor, na.rm = TRUE) - 1 # number of corhmm params
  k.ou <- max(simulators$index.ou, na.rm = TRUE) - 1 # number of ouwie params
  p.mk <- simulators$pars[1:k.cor]
  p.ou <- simulators$pars[(k.cor+1):length(simulators$pars)]
  simulators$index.cor[simulators$index.cor == 0] <- NA
  simulators$index.cor[is.na(simulators$index.cor)] <- max(simulators$index.cor, na.rm = TRUE) + 1
  Q <- matrix(0, dim(simulators$index.cor)[1], dim(simulators$index.cor)[1])
  Q[] <- c(p.mk, 0)[simulators$index.cor]
  diag(Q) <- -rowSums(Q)
  Rate.mat <- matrix(1, 3, dim(simulators$index.cor)[2])
  simulators$index.ou[is.na(simulators$index.ou)] <- max(simulators$index.ou, na.rm = TRUE) + 1
  Rate.mat[] <- c(p.ou, 1e-10)[simulators$index.ou]
  rownames(Rate.mat) <- c("alpha", "sigma.sq", "theta")
  recon <- getReconTableFromSimmap(simulators$simmap[[1]])
  obj <- list(solution.disc = Q,
              solution.cont = Rate.mat,
              recon = recon)
  class(obj) <- "houwie"
  return(obj)
}

# weight a set of simulating parameters by their time spent in a particular state
getMapWeightedPars <- function(simulators, Weight="Tip"){
  pars <- organizeSimulators(simulators)
  map <- simulators$simmap[[1]]
  nTip <- length(map$tip.label)
  TipPaths <- lapply(1:nTip, function(x) OUwie:::getPathToRoot(simulators$simmap[[1]], x))
  Pars <- list()
  pars.mk <- -diag(pars$pars.cor)
  pars.ou <- pars$pars.ou
  if(Weight == "Tip"){
    for(i in 1:nTip){
      tmp_vec_i <- vector("numeric", length(pars.mk))
      EdgeIndex_i <- TipPaths[[i]][length(TipPaths[[i]])]
      TipState_i <- names(map$maps[[EdgeIndex_i]])[length(names(map$maps[[EdgeIndex_i]]))]
      tmp_vec_i[as.numeric(TipState_i)] <- 1
      MkRate_Tip_i <- sum(tmp_vec_i * pars.mk)
      OUPars_Tip_i <- colSums(t(pars.ou) * tmp_vec_i)
      Pars[[i]] <- c(rate = MkRate_Tip_i, OUPars_Tip_i)
    }
  }
  if(Weight == "Branch"){
    for(i in 1:nTip){
      EdgeIndex_i <- TipPaths[[i]][length(TipPaths[[i]])]
      MappedEdge_i <- map$mapped.edge[EdgeIndex_i,]/sum(map$mapped.edge[EdgeIndex_i,])
      MkRate_Tip_i <- sum(MappedEdge_i * pars.mk)
      OUPars_Tip_i <- colSums(t(pars.ou) * MappedEdge_i)
      Pars[[i]] <- c(rate = MkRate_Tip_i, OUPars_Tip_i)
    }
  }
  if(Weight == "Root"){
    for(i in 1:nTip){
      MappedEdge_i <- colSums(map$mapped.edge[TipPaths[[i]],])/sum(colSums(map$mapped.edge[TipPaths[[i]],]))
      MkRate_Tip_i <- sum(MappedEdge_i * pars.mk)
      OUPars_Tip_i <- colSums(t(pars.ou) * MappedEdge_i)
      Pars[[i]] <- c(rate = MkRate_Tip_i, OUPars_Tip_i)
    }
  }
  Pars <- do.call(rbind, Pars)
  return(Pars)
}

getParDifferences <- function(simulators, fit){
  k.cor <- max(simulators$index.cor, na.rm = TRUE) - 1 # number of corhmm params
  k.ou <- max(simulators$index.ou, na.rm = TRUE) - 1 # number of ouwie params
  p.mk <- simulators$pars[1:k.cor]
  p.ou <- simulators$pars[(k.cor+1):length(simulators$pars)]
  simulators$index.cor[simulators$index.cor == 0] <- NA
  simulators$index.cor[is.na(simulators$index.cor)] <- max(simulators$index.cor, na.rm = TRUE) + 1
  Q <- matrix(0, dim(simulators$index.cor)[1], dim(simulators$index.cor)[1])
  Q[] <- c(p.mk, 0)[simulators$index.cor]
  diag(Q) <- -rowSums(Q)
  Rate.mat <- matrix(1, 3, dim(simulators$index.cor)[2])
  simulators$index.ou[is.na(simulators$index.ou)] <- max(simulators$index.ou, na.rm = TRUE) + 1
  Rate.mat[] <- c(p.ou, 1e-10)[simulators$index.ou]
  rownames(Rate.mat) <- c("alpha", "sigma.sq", "theta")
}

# get the AIC table from a particular RSave for hOUwie fits
getAICParamTablefromRsave<- function(Rsave){
  load(Rsave)
  models <- c("CDBMS", "CDOUM", "CIDBMS", "CIDOUM", "HYBBMS", "HYBOUM")
  lapply(out$model.fits, "[[", "index.ou")
  lapply(out$model.fits, "[[", "solution.ou")
  lapply(out$model.fits, "[[", "index.cor")
  lapply(out$model.fits, "[[", "solution.cor")
  nPars <- unlist(lapply(out$model.fits, "[[", "param.count"))
  AICc <- unlist(lapply(out$model.fits, "[[", "AICc"))
  LogLik <- unlist(lapply(out$model.fits, "[[", "loglik"))
  names(nPars) <- names(LogLik) <- names(AICc) <- models
  dAICc <- AICc - min(AICc)
  AICcWt <- exp(-0.5 * dAICc)/sum(exp(-0.5 * dAICc))
  return(data.frame(LogLik = LogLik, nPars = nPars, AICc = AICc, dAICc = dAICc, AICcWt = AICcWt))
}

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Plotting functions
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

plotRSME <- function(in.obj, n.obj, type = "RMSE"){
  g <- list()
  count <- 1
  for(i in 1:n.obj){
    if(n.obj == 1){
      big.obj <- in.obj
    }else{
      big.obj <- in.obj[[i]]
    }
    truepars <- organizeSimulators(big.obj[[1]]$data.houwie)
    index.ou <- big.obj[[1]]$data.houwie$index.ou
    table.resids <- table.pars <- getParTable(index.ou, big.obj)
    sim.pars <- getTrueParsVector(truepars, index.ou)
    
    for(j in 1:length(sim.pars)){
      table.resids[names(sim.pars)[j] == table.resids[,1], 3] <- table.resids[names(sim.pars)[j] == table.resids[,1], 3] - sim.pars[j]
    }
    table.resids$model <- factor(table.resids$model)
    RMSE <- aggregate(table.resids$value, by=list(table.resids$model, table.resids$par), function(x) sqrt(mean(x^2)))
    SE <- aggregate(table.resids$value, by=list(table.resids$model, table.resids$par), function(x) sd(x)/length(x))
    RMSE.table <- data.frame(Model = RMSE[,1], Param = RMSE[,2], Value = RMSE[,3], SE = SE[,3])
    
    Model <- paste(big.obj[[1]]$hOUwieFit$model.cor, big.obj[[1]]$hOUwieFit$model.ou, sep="+")
    Model <- paste("Model:", Model)
    Pars <- paste(sim.pars, collapse = ", ")
    Pars <- paste("Generating Pars:", Pars)
    Rate <- truepars$pars.cor[2,1]
    Rate <- paste("Mk Rate:", Rate)
    table.resids[,3] <- table.resids[,3]/sd(table.resids[,3])
    
    if(type == "table"){
      g[[count]] <- table.resids
    }
    
    if(type == "RMSE"){
      g[[count]] <- 
        ggplot(RMSE.table, aes(x=Param, y=Value, fill=Model)) +
        theme(legend.position = "bottom") +
        geom_errorbar(aes(ymin=Value-SE, ymax=Value+SE), width=.2) + 
        labs(x = "Parameter", y = "RMSE", title = Model, subtitle = Pars, caption = Rate) + 
        scale_fill_brewer() +
        ylim(0, NA) + 
        theme_linedraw() +
        geom_point(shape=21, size = 4) 
    }
    if(type == "resid"){
      g[[count]] <-
        ggplot(table.resids, aes(x=par, y=value, fill=model)) +
        theme_linedraw() +
        labs(x = "Parameter", y = "Standardized Residuals", title = Model, subtitle = Pars, caption = Rate) +
        theme(legend.position="bottom") +
        scale_fill_brewer() +
        geom_boxplot()
    }
    count <- count + 1
  }
  return(g)
}

organizeDat <- function(dat, simmap){
  mapping <- unlist(lapply(simmap$maps, function(x) names(x[length(x)])))
  nTip <- length(simmap$tip.label)
  TipStates <- mapping[match(match(dat[,1], simmap$tip.label), simmap$edge[,2])]
  dat[,2] <- TipStates
  return(dat)
}




ModelSetRunOld <- function(phy, pars, model.name, sim.index.cor, sim.index.ou, sim.rate.cat, fit.index.cor, fit.index.ou, fit.rate.cat, nmaps, iter, ub.ou=NULL, nstarts=1){
  nModels <- length(fit.index.ou)
  nState <- dim(sim.index.cor)[1]
  cat("Begining single run for", nModels, "models...\n")
  data.houwie <- generateData(phy, sim.index.cor, sim.index.ou, pars)
  while(length(unique(data.houwie$data[,2])) != nState){
    data.houwie <- generateData(phy, index.cor, index.ou, pars)
  }
  if(sim.rate.cat == 2){
    data <- data.houwie$data
    data[data[,2] == 3, 2] <- 1
    data[data[,2] == 4, 2] <- 2
    data[,2] <- as.factor(as.numeric(data[,2]))
  }
  if(sim.rate.cat == 1){
    data <- data.houwie$data
  }
  if(sim.rate.cat == 3){
    data <- data.houwie$data
    data[data[,2] == 3, 2] <- 1
    data[data[,2] == 5, 2] <- 1
    data[data[,2] == 4, 2] <- 2
    data[data[,2] == 6, 2] <- 2
    data[,2] <- as.factor(as.numeric(data[,2]))
  }
  
  obj <- vector("list", nModels)
  for(j in 1:nModels){
    hOUwieFit <- list()
    hOUwieFit$loglik <- -Inf
    starting.vals <- generateStartingValues(phy, fit.index.cor[[j]], fit.index.ou[[j]], data.houwie$data, nstarts)
    for(i in 1:dim(starting.vals)[1]){
      start.vals_i <- generateStartingValues.houwie(starting.vals[i,], data.houwie, fit.index.ou[[j]], fit.index.cor[[j]])
      tmp <- hOUwie(phy=phy, data=data, rate.cat=fit.rate.cat[j],
                    index.cor=fit.index.cor[[j]], root.p="yang", lb.cor=1e-5, ub.cor=21,
                    index.ou=fit.index.ou[[j]], root.station=FALSE, get.root.theta=FALSE, mserr="none", ub.ou = ub.ou, ip = start.vals_i,
                    nSim=nmaps, nCores=1)
      if(tmp$loglik > hOUwieFit$loglik){
        hOUwieFit <- tmp
      }
    }
    obj[[j]] <- hOUwieFit
  }
  out <- list(simulating.data = data.houwie,
              model.fits = obj)
  file.name <- paste0("Sim=", model.name, "-", "Pars=", paste(pars, collapse = "_"), "-nTip=", length(phy$tip.label), "-nMap=", nmaps,"-Iter=", iter, ".Rsave")
  save(out, file = file.name)
  cat("Done.\n")
  return(out)
}

getRsaves <- function(Model, nTip=NULL, nMap=NULL){
  Rsaves <- dir(paste0("~/2020_hOUwie/ModelTesting/", Model), full.names = TRUE)
  if(is.null(nTip)){
    nTip = "*"
  }
  if(is.null(nMap)){
    nMap =  "*"
  }
  Rsaves <- Rsaves[grep(paste0("nTip=", nTip), Rsaves)]
  Rsaves <- Rsaves[grep(paste0("nMap=", nMap), Rsaves)]
  return(Rsaves)
}

plotDataSet <- function(data.houwie){
  xadd = 0.2
  plotSimmap(data.houwie$simmap[[1]], fsize = 0.01, colors = cols, xlim = c(0,1.2))
  dat.prop <- (data.houwie$data[,3] - min(data.houwie$data[,3]))/max(data.houwie$data[,3] - min(data.houwie$data[,3]))
  jitter <- 0.1 * xadd
  for(i in 1:length(dat.prop)){
    lines(list(x = c(1.01, 1.01 + (dat.prop[i] * xadd)) ,y = c(i,i)))
  }
}

getPVecFromModel <- function(hOUwie.model){
  phy <- hOUwie.model$phy
  data <- hOUwie.model$data
  rate.cat <- hOUwie.model$rate.cat
  nBins <- hOUwie.model$nBins
  collapse <- hOUwie.model$collapse
  model.cor <- hOUwie.model$model.cor
  index.cor <- hOUwie.model$index.cor
  root.p <- hOUwie.model$root.p
  model.ou <- hOUwie.model$model.ou
  index.ou <- hOUwie.model$index.ou
  root.station <- hOUwie.model$root.station
  get.root.theta <- hOUwie.model$get.root.theta
  mserr <- hOUwie.model$mserr
  Tmax <- max(branching.times(phy))
  hOUwie.dat <- organizeHOUwieDat(data, mserr)
  organizedData <- getHOUwieCombosAndData(data, rate.cat, collapse, nBins)
  nObs <- length(hOUwie.dat$ObservedTraits)
  #reorder phy
  phy <- reorder(phy, "pruningwise")
  null.cor <- FALSE
  if(is.null(model.cor)){
    model.cor <- "ER"
    null.cor <- TRUE
  }
  model.set.final <- corHMM:::rate.cat.set.corHMM.JDB(phy=phy,data=hOUwie.dat$data.cor,rate.cat=rate.cat, ntraits = nObs, model = model.cor)
  # get the appropriate OU model structure
  if(is.null(index.ou)){
    index.ou <- getOUParamStructure(model.ou, "three.point", root.station, get.root.theta, dim(model.set.final$Q)[1])
  }
  
  # this allows for custom rate matricies!
  if(!is.null(index.cor)){
    index.cor[index.cor == 0] <- NA
    rate <- index.cor
    model.set.final$np <- max(rate, na.rm=TRUE)
    rate[is.na(rate)]=max(rate, na.rm=TRUE)+1
    model.set.final$rate <- rate
    model.set.final$index.matrix <- index.cor
    model.set.final$Q <- matrix(0, dim(index.cor)[1], dim(index.cor)[2])
    ## for precursor type models ##
    col.sums <- which(colSums(index.cor, na.rm=TRUE) == 0)
    row.sums <- which(rowSums(index.cor, na.rm=TRUE) == 0)
    drop.states <- col.sums[which(col.sums == row.sums)]
    if(length(drop.states > 0)){
      model.set.final$liks[,drop.states] <- 0
    }
    ## need to do anything to the ouwie matrix?
    ###############################
  }
  index.cor <- model.set.final$rate
  index.cor[index.cor == max(index.cor)] <- 0  
  p.mk <- vector("numeric", max(index.cor))
  for(i in 1:max(index.cor)){
    p.mk[i] <- hOUwie.model$solution.cor[index.cor == i][1]
  }
  p.ou <- vector("numeric", max(index.ou, na.rm = TRUE))
  for(i in 1:max(index.ou, na.rm = TRUE)){
    p.ou[i] <- na.omit(hOUwie.model$solution.ou[index.ou == i])[1]
  }
  return(c(p.mk, p.ou))
}

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Utility functions for testing power of the models 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

load_model <- function(model, model_files){
  files_to_load <- model_files[grep(model, model_files)]
  model_res <- list()
  for(i in 1:length(files_to_load)){
    load(files_to_load[i])
    model_res[[i]] <- out
  }
  names(model_res) <- files_to_load
  return(model_res)
}

get_parameter_names_discrete <- function(discrete_index, continuous_index){
  all_names <- c("rate", "alpha", "sigma", "theta")
  n_p_disc <- max(discrete_index, na.rm = TRUE)
  n_p_alpha <- length(unique(na.omit(continuous_index[1,])))
  n_p_sigma <- length(unique(na.omit(continuous_index[2,])))
  n_p_theta <- length(unique(na.omit(continuous_index[3,])))
  tmp_names <- c(rep("rate", n_p_disc), rep("alpha", n_p_alpha), rep("sigma", n_p_sigma), rep("theta", n_p_theta))
  out_names <- tmp_names
  for(i in 1:4){
    out_names[tmp_names %in% all_names[i]] <- paste0(all_names[i], "_", seq_len(length(which(tmp_names %in% all_names[i]))))
  }
  return(out_names)
}

getSignError <- function(model_list_iter){
  true_pars <- model_list_iter$simulated_data$pars
  if(all(unlist(lapply(model_list_iter$model_res, function(x) x$loglik)) > 1e5)){
    return(c(NA, NA, NA))
  }
  best_fit <- model_list_iter$model_res[[which.max(unlist(lapply(model_list_iter$model_res, function(x) x$loglik))[unlist(lapply(model_list_iter$model_res, function(x) x$loglik)) < 1e5])]]
  est_pars <- best_fit$p
  par_names <- get_parameter_names_discrete(best_fit$index.disc, best_fit$index.cont)
  alpha_index <- grep("alpha", par_names)
  sigma_index <- grep("sigma", par_names)
  theta_index <- grep("theta", par_names)
  if(length(alpha_index) > 1){
    alpha <- diff(est_pars[alpha_index]) > 0
  }else{
    alpha <- NA
  }
  if(length(sigma_index) > 1){
    sigma <- diff(est_pars[sigma_index]) > 0
  }else{
    sigma <- NA
  }
  if(length(theta_index) > 1){
    theta <- diff(est_pars[theta_index]) > 0
  }else{
    theta <- NA
  }
  error_vector <- c(alpha_correct=alpha, sigma_correct=sigma, theta_correct=theta)
  return(error_vector)
}

getSignErrorList <- function(model_list){
  tmp_table <- do.call(rbind, lapply(model_list, function(x) getSignError(x)))
  new_tmp_table <- tmp_table[!apply(tmp_table, 1, function(x) all(is.na(x))),]
  if(dim(new_tmp_table)[1] == 0){
    return(c(alpha_correct=NA, sigma_correct=NA, theta_correct=NA))
  }else{
    return(colSums(new_tmp_table)/dim(new_tmp_table)[1])
  }
}

getRMSE <- function(model_list_iter){
  true_pars <- model_list_iter$simulated_data$pars
  if(all(unlist(lapply(model_list_iter$model_res, function(x) x$loglik)) > 1e5)){
    return(c(NA, NA, NA))
  }
  best_fit <- model_list_iter$model_res[[which.max(unlist(lapply(model_list_iter$model_res, function(x) x$loglik))[unlist(lapply(model_list_iter$model_res, function(x) x$loglik)) < 1e5])]]
  is_HMM <- dim(model_list_iter$simulated_data$index.cor)[2]==4
  est_pars <- best_fit$p
  par_names <- get_parameter_names_discrete(best_fit$index.disc, best_fit$index.cont)
  alpha_index <- grep("alpha", par_names)
  sigma_index <- grep("sigma", par_names)
  theta_index <- grep("theta", par_names)
  if(length(alpha_index) > 0){
    if(is_HMM){
      alpha <- sum((sort(true_pars[alpha_index]) - sort(est_pars[alpha_index]))^2)
    }else{
      alpha <- sum((true_pars[alpha_index] - est_pars[alpha_index])^2)
    }
  }else{
    alpha <- NA
  }
  if(length(sigma_index) > 0){
    if(is_HMM){
      sigma <- sum((sort(true_pars[sigma_index]) - sort(est_pars[sigma_index]))^2)
    }else{
      sigma <- sum((true_pars[sigma_index] - est_pars[sigma_index])^2)
    }
  }else{
    sigma <- NA
  }
  if(length(theta_index) > 0){
    if(is_HMM){
      theta <- sum((sort(true_pars[theta_index]) - sort(est_pars[theta_index]))^2)
    }else{
      theta <- sum((true_pars[theta_index] - est_pars[theta_index])^2)
    }
  }else{
    theta <- NA
  }
  error_vector <- c(alpha_se=alpha, sigma_se=sigma, theta_se=theta)
  return(error_vector)
}

getRMSEList <- function(model_list){
  tmp_table <- do.call(rbind, lapply(model_list, function(x) getRMSE(x)))
  new_tmp_table <- tmp_table[!apply(tmp_table, 1, function(x) all(is.na(x))),]
  if(dim(new_tmp_table)[1] == 0){
    return(c(alpha_correct=NA, sigma_correct=NA, theta_correct=NA))
  }else{
    return(sqrt(colSums(new_tmp_table)/dim(new_tmp_table)[1]))
  }
}

getPowerMeasures <- function(model_iter, focal_model, cd_vector){
  model_names <- paste0("M", 1:22)
  power_table <- getModelTable(model_iter$model_res)
  rownames(power_table) <- model_names
  if(any(power_table$lnLik > 1e5)){
    new_model_names <- model_names[!power_table$lnLik > 1e5]
    new_cd_vector <- cd_vector[!power_table$lnLik > 1e5]
    new_model_res <- model_iter$model_res[!power_table$lnLik > 1e5]
    new_power_table <- getModelTable(new_model_res)
    rownames(new_power_table) <- new_model_names
  }else{
    new_model_res <- model_iter$model_res
    new_power_table <- power_table
    rownames(new_power_table) <- model_names
    new_cd_vector <- cd_vector
  }
  focal_row <- new_power_table[rownames(new_power_table) == focal_model,]
  if(dim(focal_row)[1] == 0){
    return(c(NA, NA, NA))
  }else{
    obj <- data.frame(best = focal_row$dAIC == 0, dAIC = focal_row$AICwt, best_type = new_cd_vector[which.min(new_power_table$dAIC)])
    return(obj)
  }
}

getRootState <- function(simmap){
  nTip <- length(simmap$tip.label)
  root_state <- as.numeric(names(simmap$maps[[which(simmap$edge[,1] == nTip+1)[1]]]))[1]
  return(root_state)
}

getBestModel <- function(model_set){
  best_fit <- model_set[[which.max(unlist(lapply(model_set, function(x) x$loglik))[unlist(lapply(model_set, function(x) x$loglik)) < 1e5])]]
  return(best_fit)
}

#### #### #### #### #### #### #### #### #### #### #### #### 
# functions for the empirical dataset
#### #### #### #### #### #### #### #### #### #### #### #### 
getCSVs <- function(wd){
  CSVs <- dir(paste0(wd, "/trait_data/"))
  CSVs <- CSVs[grep("niche", CSVs)]
  CSVs <- paste0(wd, "/trait_data/", CSVs)
  return(CSVs)
}

getData <- function(csv){
  dat <- read.csv(csv)
  dat.temp.se <- data.frame(sp = dat$species, reg = dat$Fruit_type, temp = dat$mean_temp, se_temp = dat$within_sp_var_temp)
  dat.temp.se <- dat.temp.se[which(apply(dat.temp.se, 1, function(x) !any(is.na(x)))),]
  dat.temp <- data.frame(sp = dat$species, reg = dat$Fruit_type, temp = dat$mean_temp)
  dat.temp <- dat.temp[which(apply(dat.temp, 1, function(x) !any(is.na(x)))),]
  
  dat.prec.se <- data.frame(sp = dat$species, reg = dat$Fruit_type, prec = dat$mean_prec, se_prec = dat$within_sp_var_prec)
  dat.prec.se <- dat.prec.se[which(apply(dat.prec.se, 1, function(x) !any(is.na(x)))),]
  dat.prec <- data.frame(sp = dat$species, reg = dat$Fruit_type, prec = dat$mean_prec)
  dat.prec <- dat.prec[which(apply(dat.prec, 1, function(x) !any(is.na(x)))),]
  
  dat.pet.se <- data.frame(sp = dat$species, reg = dat$Fruit_type, pet = dat$mean_pet, se_pet = dat$within_sp_var_pet)
  dat.pet.se <- dat.pet.se[which(apply(dat.pet.se, 1, function(x) !any(is.na(x)))),]
  dat.pet <- data.frame(sp = dat$species, reg = dat$Fruit_type, pet = dat$mean_pet)
  dat.pet <- dat.pet[which(apply(dat.pet, 1, function(x) !any(is.na(x)))),]
  
  dat.arid.se <- data.frame(sp = dat$species, reg = dat$Fruit_type, arid = dat$mean_aridity, se_temp = dat$within_sp_var_aridity)
  dat.arid.se <- dat.arid.se[which(apply(dat.arid.se, 1, function(x) !any(is.na(x)))),]
  dat.arid <- data.frame(sp = dat$species, reg = dat$Fruit_type, arid = dat$mean_aridity)
  dat.arid <- dat.arid[which(apply(dat.arid, 1, function(x) !any(is.na(x)))),]
  return(list(dat.temp = dat.temp,
              dat.temp.se = dat.temp.se,
              dat.prec = dat.prec,
              dat.prec.se = dat.prec.se,
              dat.pet = dat.pet,
              dat.pet.se = dat.pet.se,
              dat.arid = dat.arid,
              dat.arid.se = dat.arid.se
  ))
}

quickRun <- function(data, phy, continuous_model, discrete_model_cid, discrete_model_cd, nSim=50, time_slice=20){
  if(dim(continuous_model)[2] == 2){
    discrete_model <- discrete_model_cd
    rate.cat <- 1
  }else{
    discrete_model <- discrete_model_cid
    rate.cat <- 2
  }
  if(dim(data)[2] == 4){
    mserr <- "known"
  }else{
    mserr <- "none"
  }
  fit <- hOUwie(phy = phy, data = data, rate.cat = rate.cat, nSim = nSim, time_slice = time_slice, discrete_model = discrete_model, continuous_model = continuous_model, recon = FALSE, mserr = mserr)
  return(fit)
}
