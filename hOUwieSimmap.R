# set of functions for the hidden rates OU model

# exported function with all the bells and whistles
#'@author James Boyko
#'@param phy a phylogeny of class phylo
#'@param data data.frame where [,1] = species name,  [,2:j] are discrete traits, [,j:m] are the continuous trait (and possible mserr). if mserr = "none" j=m-1, if mserr = "known j=m-2. m is the number of total columns in the dataset.
#'@param rate.cat the number rate categories in the corHMM model
#'@param model.cor a single corHMM model or a vector of corHMM models (ARD, ER, or SYM) 
#'@param index.cor custom rate/index matrix
#'@param lb.cor lower bound for trans rate. default is 1e-5
#'@param ub.cor upper bound for trans rate. default is 21
#'@param root.p designates the root prior (yang, maddfitz, flat, numeric vector)
#'@param model.ou designates the default OU model (BM1, BMS, OU1, OUM, OUMV, OUMA, OUMVA)
#'@param index.ou designates which ou params are being estimated
#'@param root.station logical indicating whether to assume a random starting point (TRUE) or a fixed starting point (FALSE) 
#'@param get.root.theta logical indicating whether the starting state, theta_0, should be estimated
#'@param mserr designates whether mserr is included (known) or not (none)
#'@param lb.ou designates the lower bounds for alpha, sigma.sq, optima. default is given in code
#'@param ub.ou designates the upper bounds for alpha, sigma.sq, optima. default is given in code
#'@param p a vector of params to calculate a fixed likelihood
#'@param ip a vector of initial params 
#'@param nSim the number of simmaps a single set of params is evaluated over
#'@param opts optioins for nloptr
#'@param nCores number of parallel cores
#'@param quiet a logical indicating whether to output user messages
hOUwie <- function(phy, data, rate.cat, discrete_model, continuous_model, nSim=1000, root.p="yang", dual = FALSE, collapse = TRUE, lb.cor=NULL, ub.cor=NULL, root.station=FALSE, get.root.theta=FALSE, mserr = "none", lb.ou=NULL, ub.ou=NULL, p=NULL, ip=NULL, opts=NULL, nCores=1, quiet=FALSE, parsimony=FALSE, sample.tips=TRUE){
  # if the data has negative values, shift it right - we will shift it back later
  if(mserr == "none"){
    if(any(data[,dim(data)[2]] < 0)){
      cat("Negative values detected... adding 50 to the trait mean\n")
      data[,dim(data)[2]] <- data[,dim(data)[2]] + 50
    }
  }else{
    if(any(data[,dim(data)[2]-1] < 0)){
      cat("Negative values detected... adding 50 to the trait mean\n")
      data[,dim(data)[2]-1] <- data[,dim(data)[2]-1] + 50
    }
  }
  # check that tips and data match
  # check for invariance of tip states and not that non-invariance isn't just ambiguity
  if(!is.null(phy$node.label)){
    if(!quiet){
      cat("Your phylogeny had node labels, these have been removed.\n")
    }
    phy$node.label <- NULL
  }
  
  # organize the data
  hOUwie.dat <- organizeHOUwieDat(data, mserr, collapse)
  nStates <- as.numeric(max(hOUwie.dat$data.cor[,2]))
  nCol <- dim(data)[2] - ifelse(mserr == "none", 2, 3)
  tip.paths <- lapply(1:length(phy$tip.label), function(x) OUwie:::getPathToRoot(phy, x))
  Tmax <- max(branching.times(phy))
  nObs <- length(hOUwie.dat$ObservedTraits)
  phy <- reorder.phylo(phy, "pruningwise")
  algorithm <- "three.point"
  
  if(class(discrete_model)[1] == "character"){
    index.disc <- getDiscreteModel(hOUwie.dat$data.cor, discrete_model, rate.cat, dual, collapse)
  }else{
    index.disc <- discrete_model
  }
  if(class(continuous_model)[1] == "character"){
    index.cont <- getOUParamStructure(continuous_model, "three.point", root.station, get.root.theta, nStates * rate.cat)
  }else{
    index.cont <- continuous_model
  }
  if(dim(index.disc)[2] > dim(index.cont)[2]){
    stop("Not all of your discrete states have OU parameters associated with them. Please check that your discrete index matrix matches your continuous index matrix.")
  }
  if(dim(index.cont)[2] > dim(index.disc)[2]){
    stop("You have specified more OU parameters than there are states in the discrete process. Please check that your discrete index matrix matches your continuous index matrix.")
  }
  
  if(is.null(lb.ou)){
    # the lower limit of alpha is defined as a halflife of 10000% of the max tree height
    # the lower limit of sigma is defined 10 times less than alpha
    # the lower limit of optim is defined 10 times lower than the minimum observation
    lb.alpha = log(2)/(100 * Tmax)
    lb.sigma = lb.alpha/10
    lb.optim = min(data[, 1+nCol+1])/10 
    lb.ou=c(lb.alpha,lb.sigma,lb.optim)
  }
  if(is.null(ub.ou)){
    # the upper limit of alpha is defined as a halflife of 1% of the max tree height
    # the upper limit of sigma is defined 10 times more than alpha
    # the upper limit of optim is defined 10 times more than the maximum observation
    ub.alpha = log(2)/(0.01 * Tmax)
    ub.sigma = 10 * ub.alpha
    ub.optim = max(data[, 1+nCol+1])*10 
    ub.ou=c(ub.alpha,ub.sigma,ub.optim)
  }
  if(is.null(lb.cor)){
    # the minimum dwell time is defined as 100 times the max tree height
    lb.cor= 1/(Tmax*100)
  }
  if(is.null(ub.cor)){
    # the maximum dwell time is defined as 0.01 times the max tree height
    ub.cor=1/(Tmax*0.01)
  }
  
  #Ensures that weird root state probabilities that do not sum to 1 are input:
  if(!is.null(root.p)){
    if(!is.character(root.p)){
      root.p <- root.p/sum(root.p)
    }
  }
  
  # a temporary corhmm model to set the rates up
  model.set.final <- corHMM:::rate.cat.set.corHMM.JDB(phy=phy,data=hOUwie.dat$data.cor,rate.cat=rate.cat, ntraits = nObs, model = "ER", rate.mat = index.disc, collapse = collapse)

  # this allows for custom rate matricies!
  order.test <- TRUE
  if(!is.null(index.disc)){
    order.test <- FALSE
    index.disc[index.disc == 0] <- NA
    rate <- index.disc
    model.set.final$np <- max(rate, na.rm=TRUE)
    rate[is.na(rate)]=max(rate, na.rm=TRUE)+1
    model.set.final$rate <- rate
    model.set.final$index.matrix <- index.disc
    model.set.final$Q <- matrix(0, dim(index.disc)[1], dim(index.disc)[2])
    ## for precursor type models ##
    col.sums <- which(colSums(index.disc, na.rm=TRUE) == 0)
    row.sums <- which(rowSums(index.disc, na.rm=TRUE) == 0)
    drop.states <- col.sums[which(col.sums == row.sums)]
    if(length(drop.states > 0)){
      model.set.final$liks[,drop.states] <- 0
    }
    ## need to do anything to the ouwie matrix?
    ###############################
  }
  # the number of parameters for each process
  n_p_alpha <- length(unique(na.omit(index.cont[1,])))
  n_p_sigma <- length(unique(na.omit(index.cont[2,])))
  n_p_theta <- length(unique(na.omit(index.cont[3,])))
  
  # default MLE search options
  if(is.null(opts)){
    opts <- list("algorithm"="NLOPT_LN_SBPLX", "maxeval"="1000000", "ftol_rel"=.Machine$double.eps^0.5)
  }
  # p is organized into 2 groups with the first set being corHMM and the second set being OUwie
  # organized as c(trans.rt, alpha, sigma.sq, theta)
  # evaluate likelihood
  if(!is.null(p)){
    if(!quiet){
      cat("Calculating likelihood from a set of fixed parameters", "\n")
      print(p)
    }
    if(max(index.cont, na.rm = TRUE) + max(model.set.final$index.matrix, na.rm = TRUE) != length(p)){
      message <- paste0("The number of parameters does not match the number required by the model structure. You have supplied ", length(p), ", but the model structure requires ", max(index.cont, na.rm = TRUE) + max(model.set.final$index.matrix, na.rm = TRUE), ".")
      stop(message, call. = FALSE)
      
    }
    out<-NULL
    est.pars<-log(p)
    out$solution <- log(p)
    out$objective <- hOUwie.dev(est.pars, phy=phy, rate.cat=rate.cat,data.cor=hOUwie.dat$data.cor, liks=model.set.final$liks, Q=model.set.final$Q, rate=model.set.final$rate, root.p=root.p, data.ou=hOUwie.dat$data.ou, index.ou=index.cont, algorithm=algorithm, mserr=mserr,nSim=nSim, nCores=nCores, tip.paths=tip.paths, order.test=order.test, fix.node=NULL, fix.state=NULL, parsimony = parsimony, sample.tips=sample.tips)
  }else{
    if(!quiet){
      cat("Starting a search of parameters with", nSim, "simmaps...\n")
    }
    out<-NULL
    # check for user input initial parameters 
    if(is.null(ip)){
      # means.by.regime <- with(hOUwie.dat$data.ou, tapply(hOUwie.dat$data.ou[,3], hOUwie.dat$data.ou[,2], mean))
      # if(length(unique(na.omit(index.ou[3,]))) == length(means.by.regime)){
      #   start.theta <- rep(means.by.regime, length(unique(index.ou[3,]))/length(means.by.regime))
      # }else{
      #   start.theta <- rep(mean(hOUwie.dat$data.ou[,3]), length(unique(na.omit(index.ou[3,]))))
      # }
      starts.alpha <- rep(log(2)/Tmax, n_p_alpha)
      starts.sigma <- rep(var(hOUwie.dat$data.ou[,3]), n_p_sigma)
      if(rate.cat > 1){
        bin_index <- cut(hOUwie.dat$data.ou[,3], rate.cat, labels = FALSE)
        combos <- expand.grid(1:max(hOUwie.dat$data.cor[,2]), 1:rate.cat)
        disc_tips <- vector("numeric", length(phy$tip.label))
        for(i in 1:dim(combos)[1]){
          disc_tips[hOUwie.dat$data.cor[,2] == combos[i,1] & bin_index == combos[i,2]] <- i
        }
      }else{
        disc_tips <- hOUwie.dat$data.cor[,2]
      }
      start.theta <- getIP.theta(hOUwie.dat$data.ou[,3], disc_tips, index.cont[3,])
      start.ou <- c(starts.alpha, starts.sigma, start.theta)
      start.cor <- rep(10/sum(phy$edge.length), model.set.final$np)
      # start.ou <- c(rep(log(2)/Tmax, length(unique(na.omit(index.ou[1,])))), 
      #               rep(var(hOUwie.dat$data.ou[,3]), length(unique(na.omit(index.ou[2,])))), 
      #               start.theta)
      starts = c(start.cor, start.ou)
    }else{
      starts <- ip
    }
    lower = log(c(rep(lb.cor, model.set.final$np), 
                  rep(lb.ou[1], length(unique(na.omit(index.cont[1,])))), 
                  rep(lb.ou[2], length(unique(na.omit(index.cont[2,])))), 
                  rep(lb.ou[3], length(unique(na.omit(index.cont[3,]))))))
    upper = log(c(rep(ub.cor, model.set.final$np), 
                  rep(ub.ou[1], length(unique(na.omit(index.cont[1,])))), 
                  rep(ub.ou[2], length(unique(na.omit(index.cont[2,])))), 
                  rep(ub.ou[3], length(unique(na.omit(index.cont[3,]))))))
    cat(c("TotalLnLik", "DiscLnLik", "ContLnLik"), "\n")
    out = nloptr(x0=log(starts), eval_f=hOUwie.dev, lb=lower, ub=upper, opts=opts, phy=phy, rate.cat=rate.cat,data.cor=hOUwie.dat$data.cor, liks=model.set.final$liks, Q=model.set.final$Q, rate=model.set.final$rate, root.p=root.p, data.ou=hOUwie.dat$data.ou, index.ou=index.cont, algorithm=algorithm, mserr=mserr, nSim=nSim, nCores=nCores, tip.paths=tip.paths, order.test=order.test, fix.node=NULL, fix.state=NULL, parsimony = parsimony, sample.tips=sample.tips, split.liks=FALSE)
    cat("\n")
    if(!quiet){
      cat("Finished.\n")
    }
  }
  # preparing output
  FinalLiks <- hOUwie.dev(out$solution, phy=phy, rate.cat=rate.cat,data.cor=hOUwie.dat$data.cor, liks=model.set.final$liks, Q=model.set.final$Q, rate=model.set.final$rate, root.p=root.p, data.ou=hOUwie.dat$data.ou, index.ou=index.cont, algorithm=algorithm, mserr=mserr,nSim=nSim, nCores=nCores, tip.paths=tip.paths, order.test=order.test, fix.node=NULL, fix.state=NULL, parsimony = parsimony, sample.tips=sample.tips, split.liks = TRUE)
  cat("\n")
  param.count <- max(model.set.final$index.matrix, na.rm = TRUE) + max(index.cont, na.rm = TRUE)
  nb.tip <- length(phy$tip.label)
  solution <- organizeHOUwiePars(out=out, rate=model.set.final$rate, Q=model.set.final$Q, index.ou=index.cont)
  if(rate.cat > 1){
    StateNames <- paste("(", rep(1:nObs, rate.cat), rep(LETTERS[1:rate.cat], each = nObs), ")", sep = "")
  }else{
    StateNames <- paste("(", rep(1:nObs, rate.cat), ")", sep = "")
  }
  rownames(solution$solution.cor) <- colnames(solution$solution.cor) <- StateNames
  colnames(solution$solution.ou) <- StateNames
  names(hOUwie.dat$ObservedTraits) <- 1:length(hOUwie.dat$ObservedTraits)
  obj <- list(
    loglik = -FinalLiks$TotalLik,
    DiscLik = FinalLiks[1],
    ContLik = FinalLiks[2],
    AIC = -2*FinalLiks$TotalLik + 2*param.count,
    AICc = -2*FinalLiks$TotalLik+ 2*param.count*(param.count/(nb.tip-param.count-1)),
    BIC = -2*FinalLiks$TotalLik + log(nb.tip) * param.count,
    param.count = param.count,
    solution.disc = solution$solution.cor,
    solution.cont = solution$solution.ou,
    index.disc = index.disc,
    index.cont = index.cont,
    RegimeMap = FinalLiks$BestMap,
    phy = phy,
    legend = hOUwie.dat$ObservedTraits,
    data = data, 
    hOUwie.dat = hOUwie.dat,
    rate.cat = rate.cat, 
    discrete_model=discrete_model, 
    continuous_model=continuous_model, 
    root.p=root.p, 
    lb.cor=lb.cor, 
    ub.cor=ub.cor,
    root.station=root.station, 
    get.root.theta=get.root.theta, 
    mserr = mserr, 
    lb.ou=lb.ou, 
    ub.ou=ub.ou,
    p=exp(out$solution), 
    ip=ip, 
    nSim=nSim, 
    opts=opts, 
    nCores=nCores, 
    quiet=quiet
    )
  class(obj) <- "houwie"
  return(obj)
}

# for a single set of parameters, evaluate the hOUwie likelihood
hOUwie.dev <- function(p, phy, rate.cat, 
                       data.cor, liks, Q, rate, root.p,
                       data.ou, index.ou, algorithm, mserr,
                       nSim, nCores, tip.paths=NULL, order.test=TRUE, parsimony=FALSE,
                       fix.node=NULL, fix.state=NULL, sample.tips=TRUE, split.liks=FALSE){
  # params are given in log form
  p <- exp(p)
  #print(p)
  # define which params are for the HMM
  k <- max(rate)-1
  p.mk <- p[1:k]
  # set the OU params
  p.ou <- p[(k+1):length(p)] 
  Rate.mat <- matrix(1, 3, dim(rate)[2])
  alpha.na <- is.na(index.ou[1,])
  index.ou[is.na(index.ou)] <- max(index.ou, na.rm = TRUE) + 1
  Rate.mat[] <- c(p.ou, 1e-10)[index.ou]
  alpha = Rate.mat[1,]
  sigma.sq = Rate.mat[2,]
  theta = Rate.mat[3,]
  if(algorithm == "invert"){
    theta = NULL
  }
  # fit the corHMM model. if rate.cat > 1, then ensure the order of the state mats is fast to slow.
  # returns negative loglikelihood (hence sign change)
  Mk.loglik <- -corHMM:::dev.corhmm(log(p.mk), phy, liks, Q, rate, root.p, rate.cat, order.test, FALSE)
  # set up the rate matrix
  Q[] <- c(p.mk, 0)[rate]
  diag(Q) <- -rowSums(Q)
  # fit the ancestral state reconstruction
  # simulate a set of simmaps
  # sample a set of tips based on the probability of a stationary distribution generating the value
  if(sample.tips & rate.cat > 1){
    data.cor.sample <- data.cor
    sample.index <- liks[1:length(phy$tip.label),]
    means <- theta
    tmp.denom <- rep(0.5, length(alpha))
    tmp.denom[!alpha.na] <- alpha[!alpha.na]
    vars <- sigma.sq/2*tmp.denom
    normal.params <- rbind(means, vars)
    sample.tip.probs <- apply(normal.params, 2, function(x) dnorm(data.ou[,3], x[1], sqrt(x[2])))
    # sample.tip.probs <- t(apply(sample.tip.probs, 1, function(x) exp(x - max(x))))
    sample.tip.probs <- sample.tip.probs * sample.index
    Check1 <- any(apply(sample.tip.probs, 1, function(x) all(x == 0)))
    # if our sample doesn't include all tip states we revert to the unsampled simmap
    
    if(Check1){
      simmap <- corHMM:::makeSimmap(phy, data.cor, Q, rate.cat, root.p = root.p, nSim = nSim, fix.node = fix.node, fix.state = fix.state, parsimony=parsimony)
    }else{
      data.cor.sample[,2] <- apply(sample.tip.probs, 1, function(x) sample(1:dim(Q)[1], 1, prob = x))
      Check2 <- length(unique(data.cor.sample[,2])) != dim(Q)[1]
      if(Check2){
        simmap <- corHMM:::makeSimmap(phy, data.cor, Q, rate.cat, root.p = root.p, nSim = nSim, fix.node = fix.node, fix.state = fix.state, parsimony=parsimony)
      }else{
        simmap <- corHMM:::makeSimmap(phy, data.cor.sample, Q, 1, root.p = root.p, nSim = nSim, fix.node = fix.node, fix.state = fix.state, parsimony=parsimony)
      }
    }
  }else{
    simmap <- corHMM:::makeSimmap(phy, data.cor, Q, rate.cat, root.p = root.p, nSim = nSim, fix.node = fix.node, fix.state = fix.state, parsimony=parsimony)
  }
  # fit the OU models to the simmaps
  # returns log likelihood
  OU.loglik <- lapply(simmap, function(x) OUwie.basic(x, data.ou, simmap.tree=TRUE, scaleHeight=FALSE, alpha=alpha, sigma.sq=sigma.sq, theta=theta, algorithm=algorithm, tip.paths=tip.paths, mserr=mserr))
  # to show that OUwie.basic is identical to OUwie.fixed
  # OUwie.basic(simmap[[1]], data.ou, simmap.tree=TRUE, scaleHeight=FALSE, alpha=alpha, sigma.sq=sigma.sq, theta=theta, algorithm=algorithm, tip.paths=tip.paths, mserr=mserr)
  # OUwie.fixed(simmap[[1]], data.ou, model = "OUM", simmap.tree=TRUE, scaleHeight=FALSE, alpha=alpha, sigma.sq=sigma.sq, theta=theta, algorithm="invert", tip.paths=tip.paths, mserr=mserr)
  # OUwie.fixed(simmap[[1]], data.ou, model = "OUM", simmap.tree=TRUE, scaleHeight=FALSE, alpha=alpha, sigma.sq=sigma.sq, theta=theta, algorithm="three.point", tip.paths=tip.paths)

  # get the likelihoods of the simmaps
  # OU.loglik <- max(unlist(OU.loglik))
  OU.loglik <- unlist(OU.loglik)
  # OU.loglik <- max(OU.loglik) + log(sum(exp(OU.loglik - max(OU.loglik)))) - log(nSim)
  Mk.loglik <- unlist(lapply(simmap, function(x) getMapProbability(x$maps, Q)))
  # Mk.loglik <- max(Mk.loglik) + log(sum(exp(Mk.loglik - max(Mk.loglik)))) - log(nSim)
  Total.loglik <- max(OU.loglik + Mk.loglik)
  Mk.loglik.tmp <- Mk.loglik[which.max(OU.loglik + Mk.loglik)]
  OU.loglik.tmp <- OU.loglik[which.max(OU.loglik + Mk.loglik)]
  simmap.tmp <- simmap[which.max(OU.loglik + Mk.loglik)]
  # OU.loglik <- log(mean(exp(unlist(OU.loglik)-comp)))+comp
  cat("\r", c(Total.loglik, Mk.loglik.tmp, OU.loglik.tmp), "     ")
  if(split.liks){
    return(c(TotalLik = Total.loglik, DiscLik = Mk.loglik.tmp, ContLik = OU.loglik.tmp, BestMap = simmap.tmp))
  }
  return(-(Total.loglik))
}
  # }else{
  #   # algebraic trick to prevent overflow and underflow while preserving as many accurate leading digits in the result as possible. The leading digits are preserved by pulling the maximum outside. The arithmetic is robust because subtracting the maximum on the inside makes sure that only negative numbers or zero are ever exponentiated, so there can be no overflow on those calculations. If there is underflow, we know the leading digits have already been returned as part of the max term on the outside. subtract log nsim for avg.
  #   OU.loglik <- unlist(OU.loglik)
  #   OU.loglik <- max(OU.loglik) + log(sum(exp(OU.loglik - max(OU.loglik)))) - log(nSim)
  #   # example
  #   # tmp <- c(-100, -200, -300, -100, -250)
  #   # max(tmp) + log(sum(exp(tmp - max(tmp)))) - log(length(tmp))
  #   # log(mean(exp(tmp)))
  #   # 
  #   print(c(Total = OU.loglik + Mk.loglik, DiscLik = Mk.loglik, ContLik = OU.loglik))
  #   if(split.liks){
  #     return(c(TotalLik = Total.loglik, DiscLik = Mk.loglik, ContLik = OU.loglik))
  #   }
  #   return(-(OU.loglik + Mk.loglik))
  # }
#}

# hOUwie's input data will be similar to OUwie, with the exception that there can be more than a single trait being evaluated thus it's defined as column 1 is species name, the last column is the continuous trait and anything in between are discrete characters (can also account for mserr if second last columns are continuous)
organizeHOUwieDat <- function(data, mserr, collapse = TRUE){
  # return a list of corHMM data and OU data
  if(mserr=="known"){
    data.cor <- data[, 1:(dim(data)[2]-2)]
    data.cor <- corHMM:::corProcessData(data.cor, collapse = collapse)
    data.ou <- data.frame(sp = data[,1], 
                          reg = data.cor$corData[,2], 
                          x = data[, dim(data)[2]-1],
                          err = data[, dim(data)[2]])
  }
  if(mserr=="none"){
    data.cor <- data[, 1:(dim(data)[2]-1)]
    data.cor <- corHMM:::corProcessData(data.cor)
    data.ou <- data.frame(sp = data[,1], 
                          reg = data.cor$corData[,2], 
                          x = data[, dim(data)[2]])
  }
  return(list(StateMats = data.cor$StateMats, 
              PossibleTraits = data.cor$PossibleTraits,
              ObservedTraits = data.cor$ObservedTraits,
              data.cor = data.cor$corData,
              data.ou = data.ou))
}

getIP.theta <- function(x, states, index){
  ip.theta <- vector("numeric", length(unique(index)))
  for(i in 1:length(unique(index))){
    state_i <- which(unique(index)[i] == index)
    ip.theta[i] <- mean(x[states %in% state_i])
  }
  ip.theta[is.nan(ip.theta)] <- mean(x)
  return(ip.theta)
}

getDiscreteModel <- function(data, model, rate.cat, dual, collapse){
  rate <- getStateMat4Dat(data, model, dual, collapse)$rate.mat
  if (rate.cat > 1) {
    StateMats <- vector("list", rate.cat)
    for (i in 1:rate.cat) {
      StateMats[[i]] <- rate
    }
    rate <- getFullMat(StateMats)
  }
  return(rate)
}

# gets the probability of a particular painting. limited to a single change at a shift point of 0.5
getMapProbability <- function(maps, Q){
  BranchProbs <- lapply(maps, function(x) log(probPath(x, Q)))
  LnLik_map <- sum(unlist(BranchProbs))
  return(LnLik_map)
}

# the probability of a particular path
probPath <- function(path, Q){
  nTrans <- length(path)
  P <- vector("numeric", length(path))
  for(i in sequence(nTrans-1)){
    state_i <- as.numeric(names(path)[1])
    state_j <- as.numeric(names(path)[2])
    time_i <- as.numeric(path[1])
    rate_i <- abs(Q[state_i,state_j])
    P[i] <- dexp(time_i, rate_i)
    path <- path[-1]
  }
  state_j <- as.numeric(names(path))
  time_j <- as.numeric(path)
  rate_j <- abs(Q[state_j,state_j])
  P[nTrans] <- 1 - pexp(rate_j, time_j)
  P <- prod(P)
  if(P == 0){
    P <- 1e-99
  }
  return(P)
}

# probPath <- function(path, Q){
#   nTrans <- length(path)
#   P <- vector("numeric", length(path))
#   for(i in sequence(nTrans-1)){
#     state_i <- as.numeric(names(path)[1])
#     state_j <- as.numeric(names(path)[2])
#     time_i <- as.numeric(path[1])
#     time_j <- as.numeric(path[2])
#     # time_j <- as.numeric(sum(path[-1]))
#     rate_i <- abs(Q[state_i,state_j])
#     rate_j <- abs(Q[state_j,state_j])
#     P[i] <- dexp(time_i, rate_i) * (1 - pexp(rate_j, time_j))
#     path <- path[-1]
#   }
#   state_j <- as.numeric(names(path))
#   time_j <- as.numeric(path)
#   rate_j <- abs(Q[state_j,state_j])
#   P[nTrans] <- 1 - pexp(rate_j, time_j)
#   P <- prod(P)
#   if(P == 0){
#     P <- 1e-99
#   }
#   return(P)
# }

# different OU models have different parameter structures. This will evaluate the appropriate one.
getOUParamStructure <- function(model, algorithm, root.station, get.root.theta, k){
  index.mat <- matrix(0, 2,k)
  if(algorithm == "three.point"){
    Rate.mat <- matrix(1, 3, k)
  }else{
    Rate.mat <- matrix(1, 2, k)
  }
  if (model == "BM1"){
    np <- 1
    index.mat[1,1:k] <- NA
    index.mat[2,1:k] <- 1
    if(algorithm == "three.point"){
      index.mat <- rbind(index.mat, rep(2,k))
    }
    param.count <- np+1
  }
  #The group mean model of Thomas et al, is trivial: set root.station to be TRUE:
  if (model == "BMS"){
    np <- k
    index.mat[1,1:k] <- NA
    index.mat[2,1:k] <- 1:np
    if(root.station==TRUE){
      if(algorithm == "three.point"){
        max.par.so.far <- max(index.mat, na.rm = TRUE)
        index.mat <- rbind(index.mat, (max.par.so.far + 1):(max.par.so.far + k))
      }
      param.count <- np+k
    }
    if(root.station==FALSE){
      if(algorithm == "three.point"){
        max.par.so.far <- max(index.mat, na.rm = TRUE)
        index.mat <- rbind(index.mat, max.par.so.far+1)
      }
      param.count <- np+k
    }
  }
  if (model == "OU1"){
    np <- 2
    index.mat[1,1:k] <- 1
    index.mat[2,1:k] <- 2
    if(algorithm == "three.point"){
      max.par.so.far <- max(index.mat)
      index.mat <- rbind(index.mat, rep(3,k))
    }
    param.count <- np + 1
    if(get.root.theta == TRUE){
      param.count <- param.count + 1
    }
  }
  if (model == "OUM"){
    np <- 2
    index.mat[1,1:k] <- 1
    index.mat[2,1:k] <- 2
    if(algorithm == "three.point"){
      max.par.so.far <- max(index.mat)
      index.mat <- rbind(index.mat, (max.par.so.far + 1):(max.par.so.far + k))
    }
    param.count <- np + k
    if(get.root.theta == TRUE){
      param.count <- param.count + 1
    }
  }
  if (model == "OUMV") {
    np <- k+1
    index.mat[1,1:k] <- 1
    index.mat[2,1:k] <- 2:(k+1)
    if(algorithm == "three.point"){
      max.par.so.far <- max(index.mat)
      index.mat <- rbind(index.mat, (max.par.so.far + 1):(max.par.so.far + k))
    }
    param.count <- np + k
    if(get.root.theta == TRUE){
      param.count <- param.count + 1
    }
  }
  if (model == "OUMA") {
    np <- k+1
    index.mat[1,1:k] <- 1:k
    index.mat[2,1:k] <- k+1
    if(algorithm == "three.point"){
      max.par.so.far <- max(index.mat)
      index.mat <- rbind(index.mat, (max.par.so.far + 1):(max.par.so.far + k))
    }
    param.count <- np+k
    if(get.root.theta == TRUE){
      param.count <- param.count + 1
    }
  }
  if (model == "OUMVA") {
    np <- k*2
    index.mat[1,1:k] <- 1:k
    index.mat[2,1:k] <- (k+1):(k*2)
    if(algorithm == "three.point"){
      max.par.so.far <- max(index.mat)
      index.mat <- rbind(index.mat, (max.par.so.far + 1):(max.par.so.far + k))
    }
    param.count <- np+k
    if(get.root.theta == TRUE){
      param.count <- param.count + 1
    }
  }
  if (model == "TrendyM"){
    index.mat <- matrix(0,3,k)
    np <- k+2 #We have k regimes, but 1 sigma, plus the root value
    index.mat[1,1:k] <- np+1
    index.mat[2,1:k] <- 1
    index.mat[3,1:k] <- 2:(np-1)
    param.count<-np
    root.station=TRUE
    trendy=1
  }
  if (model == "TrendyMS"){
    index.mat <- matrix(0,3,k)
    np <- (k*2)+1 #We have k regimes and assume k trends and k sigmas, plus the root value
    index.mat[1,1:k] <- np+1
    index.mat[2,1:k] <- 1:k
    index.mat[3,1:k] <- (k+1):(np-1)
    param.count <- np
    root.station=FALSE
    trendy=1
  }
  if(algorithm=="three.point"){
    rownames(index.mat) <- c("alpha", "sigma2", "theta")
  }
  return(index.mat)
}

# a basic version of OUwie that will evaluate params based on an index matrix
OUwie.basic<-function(phy, data, simmap.tree=TRUE, root.age=NULL, scaleHeight=FALSE, root.station=FALSE, get.root.theta=FALSE, shift.point=0.5, alpha, sigma.sq, theta, mserr="none", algorithm="three.point", tip.paths=NULL){
  
  # organize tip states based on what the simmap suggests
  mapping <- unlist(lapply(phy$maps, function(x) names(x[length(x)])))
  nTip <- length(phy$tip.label)
  TipStates <- mapping[match(match(data[,1], phy$tip.label), phy$edge[,2])]
  data[,2] <- TipStates
  
  #Makes sure the data is in the same order as the tip labels
  if(mserr=="none"){
    data <- data.frame(data[,2], data[,3], row.names=data[,1])
    data <- data[phy$tip.label,]
  }
  if(mserr=="known"){
    # algorithm = "invert"
    if(!dim(data)[2]==4){
      stop("You specified measurement error should be incorporated, but this information is missing.", call. = FALSE)
    }
    else{
      if(is.factor(data[,4]) == TRUE){
        stop("Check the format of the measurement error column. It's reading as a factor.",  call. = FALSE)
      }else{
        data <- data.frame(data[,2], data[,3], data[,4], row.names=data[,1])
        data <- data[phy$tip.label,]
      }
    }
  }
  
  #Values to be used throughout
  n <- max(phy$edge[,1])
  ntips <- length(phy$tip.label)

  # setup values when simmap (always simmap for hOUwie)
  k <- length(colnames(phy$mapped.edge))
  tot.states <- factor(colnames(phy$mapped.edge))
  tip.states <- factor(data[,1])
  data[,1] <- as.numeric(tip.states)
  #Obtains the state at the root
  root.edge.index <- which(phy$edge[,1] == ntips+1)
  root.state <- which(colnames(phy$mapped.edge)==names(phy$maps[[root.edge.index[2]]][1]))
  ##Begins the construction of the edges matrix -- similar to the ouch format##
  edges <- cbind(c(1:(n-1)),phy$edge,OUwie:::MakeAgeTable(phy, root.age=root.age))
  if(scaleHeight == TRUE){
    Tmax <- max(OUwie:::MakeAgeTable(phy, root.age=root.age))
    edges[,4:5]<-edges[,4:5]/Tmax
    root.age <-  1
    phy$maps <- lapply(phy$maps, function(x) x/Tmax)
  }
  edges <- edges[sort.list(edges[,3]),]
  #Resort the edge matrix so that it looks like the original matrix order
  edges <- edges[sort.list(edges[,1]),]

  if(algorithm == "three.point"){
    x <- data[,2]
    names(x) <- rownames(data)
  }else{
    x <- as.matrix(data[,2])
  }
  
  if(scaleHeight==TRUE){
    phy$edge.length <- phy$edge.length/Tmax
    Tmax <- 1
  }
  
  if(algorithm == "three.point"){
    if(simmap.tree == FALSE){
      map <- getMapFromNode(phy, tip.states, node.states, shift.point)
    }else{
      map <- phy$maps
    }
  }
  
  Rate.mat <- rbind(alpha, sigma.sq, theta)
  pars <- matrix(c(theta, sigma.sq, alpha), length(theta), 3, dimnames = list(1:length(sigma.sq), c("opt", "sig", "alp")))
  # if the simmap did not simulate every possible state in a given hmm
  if(dim(phy$mapped.edge)[2] != dim(Rate.mat)[2]){
    Rate.mat <- Rate.mat[,as.numeric(colnames(phy$mapped.edge))]
    pars <- pars[as.numeric(colnames(phy$mapped.edge)), ]
  }
  
  if(get.root.theta == TRUE){
    W <- OUwie:::weight.mat(phy, edges, Rate.mat, root.state=root.state, simmap.tree=simmap.tree, root.age=root.age, scaleHeight=scaleHeight, assume.station=FALSE, shift.point=shift.point)
  }else{
    W <- OUwie:::weight.mat(phy, edges, Rate.mat, root.state=root.state, simmap.tree=simmap.tree, root.age=root.age, scaleHeight=scaleHeight, assume.station=TRUE, shift.point=shift.point)
  }
  
  #Likelihood function for estimating model parameters
  if(get.root.theta == TRUE){
    expected.vals <- colSums(t(W) * c(theta0, pars[,1]))
    names(expected.vals) <- phy$tip.label
  }else{
    expected.vals <- colSums(t(W) * pars[,1])
    names(expected.vals) <- phy$tip.label
  }
  transformed.tree <- OUwie:::transformPhy(phy, map, pars, tip.paths)
  # generate a map from node based reconstructions
  if(mserr=="known"){
    TIPS <- transformed.tree$tree$edge[,2] <= length(transformed.tree$tree$tip.label)
    transformed.tree$tree$edge.length[TIPS] <- transformed.tree$tree$edge.length[TIPS] + (data[,3]^2/transformed.tree$diag/transformed.tree$diag)
  }
  comp <- NA
  try(comp <- phylolm::three.point.compute(transformed.tree$tree, x, expected.vals, transformed.tree$diag, check.precision = FALSE), silent=TRUE)
  if(is.na(comp[1])){
    return(-1e10)
  }else{
    nTips <- length(phy$tip.label)
    logl <- -as.numeric(nTips * log(2 * pi) + comp$logd + (comp$PP - 2 * comp$QP + comp$QQ))/2
    return(logl)
  }
}

# organize output
organizeHOUwiePars <- function(out, rate, Q, index.ou){
  pars <- exp(out$solution)
  k <- max(rate)-1
  p.mk <- pars[1:k]
  p.ou <- pars[(k+1):length(pars)] 
  # ouwie pars
  Rate.mat <- matrix(1, 3, dim(rate)[2])
  index.ou[is.na(index.ou)] <- max(index.ou, na.rm = TRUE) + 1
  Rate.mat[] <- c(p.ou, NA)[index.ou]
  rownames(Rate.mat) <- rownames(index.ou)
  Q[] <- c(p.mk, NA)[rate]
  diag(Q) <- -rowSums(Q)
  if(!is.null(colnames(rate))){
    colnames(Rate.mat) <- colnames(rate)
    colnames(Q) <- rownames(Q) <- colnames(rate)
  }
  # corhmm pars 
  return(list(solution.ou = Rate.mat,
              solution.cor = Q))
}

simCharacterHistory <- function(phy, Q, root.freqs, Q2 = NA, NoI = NA){
  
  #Randomly choose starting state at root using the root.values as the probability:
  root.value <- sample.int(dim(Q)[2], 1, FALSE, prob=root.freqs/sum(root.freqs))
  #Reorder the phy:
  phy <- reorder.phylo(phy, "postorder")
  ntips <- length(phy$tip.label)
  N <- dim(phy$edge)[1]
  ROOT <- ntips + 1 #perhaps use an accessor to get the root node id
  
  #Generate vector that contains the simulated states:
  CharacterHistory <- integer(ntips + phy$Nnode)
  CharacterHistory[ROOT] <- as.integer(root.value)
  anc <- phy$edge[, 1]
  des <- phy$edge[, 2]
  edge.length <- phy$edge.length
  diag(Q) = 0
  diag(Q) = -rowSums(Q)
  
  # setting up the alternative Q matrix at the node of interest
  if(!any(is.na(Q2))){
    diag(Q2) = 0
    diag(Q2) = -rowSums(Q2)
  }
  if(!is.na(NoI)){
    NewQDesc <- getDescendants(phy, NoI)
  }
  
  #standard simulation protocol
  if(any(is.na(Q2)) | is.na(NoI)){
    for (i in N:1) {
      p <- expm(Q * edge.length[i], method="Ward77")[CharacterHistory[anc[i]], ]
      CharacterHistory[des[i]] <- sample.int(dim(Q)[2], size = 1, FALSE, prob = p)
    }
  }
  
  # simulating a clade under a different (Q2) evolutionary model
  if(!any(is.na(Q2)) & !is.na(NoI)){
    for (i in N:1) {
      if(anc[i] %in% NewQDesc){
        p <- expm(Q2 * edge.length[i], method="Ward77")[CharacterHistory[anc[i]], ]
        CharacterHistory[des[i]] <- sample.int(dim(Q2)[2], size = 1, FALSE, prob = p)
      }else{
        p <- expm(Q * edge.length[i], method="Ward77")[CharacterHistory[anc[i]], ]
        CharacterHistory[des[i]] <- sample.int(dim(Q)[2], size = 1, FALSE, prob = p)
      }
    }
  }
  
  TipStates <-  CharacterHistory[1:ntips]
  names(TipStates) <- phy$tip.label
  NodeStates <- CharacterHistory[ROOT:(N+1)]
  names(NodeStates) <- ROOT:(N+1)
  
  res <- list(TipStates = TipStates, NodeStates = NodeStates)
  return(res)
  #return(CharacterHistory)
}

# simulate a hOUwie model
hOUwie.sim <- function(phy, Q, root.freqs, alpha, sig2, theta0, theta, nMap=1){
  # simulate an Mk dataset
  dat.cor <- simCharacterHistory(phy, Q, root.freqs)$TipStates
  while(!all(levels(dat.cor) %in% dat.cor)){
    dat.cor <- simCharacterHistory(phy, Q, root.freqs)$TipStates
  }
  dat.cor <- data.frame(sp=names(dat.cor), d=dat.cor)
  # simulate a stochastic map with true Q
  simmap <- corHMM:::makeSimmap(phy, dat.cor, Q, 1, nSim=nMap)
  # lik <- corHMM:::getSimmapLik(simmap, Q)
  # simulate the ou dataset
  dat.ou <- lapply(simmap, function(x) OUwie.sim(x, simmap.tree = TRUE, alpha = alpha, sigma.sq = sig2, theta0 = theta0, theta = theta)[,2])
  dat.ou <- colMeans(do.call(rbind, dat.ou))
  # dat.ou <- OUwie.sim(simmap, simmap.tree = TRUE, alpha = alpha, sigma.sq = sig2, theta0 = theta0, theta = theta)
  # return true params and data
  data <- data.frame(sp = phy$tip.label, reg = dat.cor[,2], x = dat.ou)
  return(list(data = data, simmap = simmap))
}

# print a houwie object
print.houwie <- function(x, ...){
  ntips <- Ntip(x$phy)
  output <- data.frame(x$loglik,x$AIC,x$AICc,x$BIC, ntips, x$param.count, row.names="")
  names(output) <- c("lnL","AIC","AICc","BIC","nTaxa","nPars")
  cat("\nFit\n")
  print(output)
  cat("\nLegend\n")
  print(x$legend)
  cat("\nRegime Rate matrix\n")
  print(x$solution.disc)
  cat("\nOU Estimates\n")
  print(x$solution.cont)
  cat("\n")
  if(!any(is.na(x$solution.cont[1,]))){
    cat("\nHalf-life (another way of reporting alpha)\n")
    print(log(2)/x$solution.cont[1,])
  }
  cat("\nDon't forget: your parameters have units!\n")
}

# fits a model based on the best corhmm model, simmaps, ouwie
fitNonCensored <- function(phy, data, rate.cat, 
                   model.cor, root.p="yang", lb.cor=1e-5, ub.cor=10,
                   model.ou, root.station=FALSE, get.root.theta=FALSE, mserr = "none", ub.ou=NULL, 
                   nSim=1000, nCores=1, quiet=FALSE, ip=NULL){
  # a single attempt to fit an OU model to a simmap where the data is organized to have simmap tips
  singleRun <- function(dat, simmap, model, mserr, ub.ou, ip){
    data <- organizeDat(dat, simmap)
    out <- try(OUwie(simmap, data, model, simmap.tree = TRUE, algorithm = "three.point", scaleHeight = FALSE, mserr = mserr, quiet = TRUE, ub = ub.ou, starting.vals = ip))
    return(out)
  }
  #Ensures that weird root state probabilities that do not sum to 1 are input:
  algorithm="three.point"
  # organize the data into the corHMM data and the OUwie data
  hOUwie.dat <- organizeHOUwieDat(data, mserr)
  # fit the corhmm model
  cor.fit <- corHMM(phy = phy, data = hOUwie.dat$data.cor, rate.cat = rate.cat, model = model.cor, node.states = "none", root.p = root.p, lower.bound = lb.cor, upper.bound = ub.cor, ip = ip[1])
  # generate maps based on MLE corhmm
  maps <- makeSimmap(tree = phy, data = hOUwie.dat$data.cor, model = cor.fit$solution, rate.cat = rate.cat, root.p = root.p, nSim = nSim)
  # fit ouwie to all maps
  ou.fits <- mclapply(maps, function(x) singleRun(dat = hOUwie.dat$data.ou, simmap = x, model = model.ou, mserr = mserr, ub.ou = ub.ou, ip = ip[2:length(ip)]), mc.cores = nCores)
  # the average solution is our estimate
  mean.ou.solution <- Reduce("+", lapply(ou.fits, function(x) x$solution))/nSim
  # probbaly can't do this, but take average likelihood and add to corhmm's
  loglik <- Reduce("+", lapply(ou.fits, function(x) x$loglik))/nSim + cor.fit$loglik
  param.count <- max(ou.fits[[1]]$index.mat, na.rm = TRUE) + max(cor.fit$index.mat, na.rm = TRUE) - 1
  nb.tip <- length(phy$tip.label)
  
  obj <- list(
    loglik = loglik,
    AIC = -2*loglik + 2*param.count,
    AICc = -2*loglik+ 2*param.count*(param.count/(nb.tip-param.count-1)),
    BIC = -2*loglik + log(nb.tip) * param.count,
    param.count = param.count,
    solution.cor = cor.fit$solution,
    solution.ou = mean.ou.solution,
    phy = phy, 
    data = data, 
    hOUwie.dat = hOUwie.dat,
    rate.cat = rate.cat, 
    model.cor=model.cor, 
    root.p=root.p, 
    lb.cor=lb.cor, 
    ub.cor=ub.cor,
    model.ou=model.ou, 
    root.station=root.station, 
    get.root.theta=get.root.theta, 
    mserr = mserr, 
    nSim=nSim, 
    nCores=nCores)
  class(obj) <- "houwie"
  return(obj)
}

# fit a corhmm model then fit an ou model, bingo
fitTwoStep <- function(phy, data, rate.cat, 
                       model.cor, root.p="yang", lb.cor=1e-5, ub.cor=10,
                       model.ou, root.station=FALSE, get.root.theta=FALSE, mserr = "none", ub.ou=NULL, ip=NULL){
  #Ensures that weird root state probabilities that do not sum to 1 are input:
  algorithm="three.point"
  # organize the data into the corHMM data and the OUwie data
  hOUwie.dat <- organizeHOUwieDat(data, mserr)
  nObs <- length(hOUwie.dat$ObservedTraits)
  # a way to speed up the three.point function

  cor.fit <- corHMM(phy = phy, data = hOUwie.dat$data.cor, rate.cat = rate.cat, model = model.cor, node.states = "none", root.p = root.p, lower.bound = lb.cor, upper.bound = ub.cor, ip = ip[1])
  
  node.states <- apply(ace(hOUwie.dat$data.ou[,2], phy, type = "discrete")$lik.anc, 1, function(x) which(round(x) == 1))
  phy$node.label <- node.states
  ou.fit <- OUwie(phy = phy, data = hOUwie.dat$data.ou, model = model.ou, simmap.tree = FALSE, root.station = root.station, get.root.theta = get.root.theta, algorithm = algorithm, mserr = mserr, ub = ub.ou, starting.vals = ip[2:length(ip)])
  
  loglik <- cor.fit$loglik + ou.fit$loglik
  param.count <- ou.fit$param.count + max(cor.fit$index.mat, na.rm = TRUE)
  nb.tip <- length(phy$tip.label)
  
  obj <- list(
    loglik = loglik,
    AIC = -2*loglik + 2*param.count,
    AICc = -2*loglik+ 2*param.count*(param.count/(nb.tip-param.count-1)),
    BIC = -2*loglik + log(nb.tip) * param.count,
    param.count = param.count,
    solution.cor = cor.fit$solution,
    solution.ou = ou.fit$solution,
    phy = phy, 
    data = data, 
    hOUwie.dat = hOUwie.dat,
    rate.cat = rate.cat, 
    model.cor=model.cor, 
    root.p=root.p, 
    lb.cor=lb.cor, 
    ub.cor=ub.cor,
    model.ou=model.ou, 
    root.station=root.station, 
    get.root.theta=get.root.theta, 
    mserr = mserr)
  class(obj) <- "houwie"
  return(obj)
}


organizeDat <- function(dat, simmap){
  mapping <- unlist(lapply(simmap$maps, function(x) names(x[length(x)])))
  nTip <- length(simmap$tip.label)
  TipStates <- mapping[match(match(dat$sp, simmap$tip.label), simmap$edge[,2])]
  dat$reg <- TipStates
  return(dat)
}

getModelTable <- function(model.list, type="AIC"){
  # checks
  if(class(model.list) != "list"){
    stop("Input object must be of class list with each element as a separet fit model to the same dataset.", call. = FALSE)
  }
  if(length(model.list) == 1){
    stop("Two or models are needed to conduct model averaging.", call. = FALSE)
  }
  if(!all(unlist(lapply(model.list, function(x) class(x))) == "houwie")){
    stop("Not all models are of class houwie.", call.=FALSE)
  }
  if(var(unlist(lapply(model.list, function(x) dim(x$data)[1]))) != 0){
    stop("The number of rows in your data are not the same for all models. Models cannot be averaged if they are not evaluating the same dataset", call.=FALSE)
  }
  
  ParCount <- unlist(lapply(model.list, function(x) x$param.count))
  nTip <- length(model.list[[1]]$phy$tip.label)
  AIC <- simplify2array(lapply(model.list, "[[", type))
  dAIC <- AIC - min(AIC)
  AICwt <- exp(-0.5 * dAIC)/sum(exp(-0.5 * dAIC))
  LogLik <- simplify2array(lapply(model.list, "[[", "loglik"))
  out <- data.frame(np = ParCount, lnLik = LogLik, AIC = AIC, dAIC = dAIC, AICwt = AICwt)
  colnames(out) <- gsub("AIC", type, colnames(out))
  return(out)
}

# prunes redundant models. modiified from https://github.com/thej022214/hisse/blob/master/R/PruneRedundantModels.R
PruneRedundantModels <- function(..., precision=1e-5) {
  models <- list(...)
  ## Check if we have a single list with all the models or the arguments are a list.
  if( !inherits(models[[1]], what = c("houwie"))){
    ## This is a list of models. Need to extract.
    models <- models[[1]]
  }
  ## Check if elements of the list are houwie
  mod.class <- sapply(models, function(x) inherits(x, what = "houwie"))

  if(!all(mod.class)){
    ## Strange! Break.
    stop( "list of models need to be only hOUwie fits.", .call=FALSE)
  }

  mod.nparameters <- simplify2array(lapply(models, "[[", "param.count"))
  index <- order(mod.nparameters, decreasing=FALSE)
  models <- models[index]
  mod.loglik <- as.numeric(simplify2array(lapply(models, "[[", "loglik")))
  models_to_delete <- c()
  isTrueAllEqual <- function(...) {
    return(isTRUE(all.equal(...)))
  }
  if(length(models)>1) {
    for(i in 2:(length(models))) {
      if(any(sapply(mod.loglik[1:(i-1)], isTrueAllEqual, mod.loglik[i], tolerance=precision))) {
        models_to_delete <- c(models_to_delete, i)
      }
    }
  }
  if(length(models_to_delete)>0) {
    cat(length(models_to_delete), "of", length(models), "models have been removed due to model redundancy.\n")
    models <- models[-models_to_delete]
  }
  return(models)
}

## wrapper function that will evaluate some number of nodes via simmaps or marginal 
# hOUwieRecon <- function(nodes="all", type="marginal", nMap=50, nCores=1, houwie.obj){
#   phy <- houwie.obj$phy
#   data <- houwie.obj$data
#   mserr <- houwie.obj$mserr
#   rate.cat <- houwie.obj$rate.cat
#   model.cor <- houwie.obj$model.cor
#   index.cor <- houwie.obj$index.cor
#   model.ou <- houwie.obj$model.ou
#   index.ou <- houwie.obj$index.ou
#   solution.cor <- houwie.obj$solution.cor
#   solution.ou <- houwie.obj$solution.ou
#   root.p <- houwie.obj$root.p
#   root.station <- houwie.obj$root.station
#   get.root.theta <- houwie.obj$get.root.theta
#   weighted <- houwie.obj$weighted
#   
#   if(is.character(nodes)){
#     if(nodes == "all"){
#       nodes.to.eval <- 1:max(phy$edge)
#     }
#     if(nodes == "external"){
#       nodes.to.eval <- 1:length(phy$tip.label)
#     }
#     if(nodes == "internal"){
#       nodes.to.eval <- (length(phy$tip.label)+1):max(phy$edge)
#     }
#   }
#   if(is.numeric(nodes)){
#     nodes.to.eval <- nodes
#   }
#   
#   algorithm="three.point"
#   # organize the data into the corHMM data and the OUwie data
#   # TO DO: add a way to shift negative continuous variables to positive then shift back
#   hOUwie.dat <- organizeHOUwieDat(data, mserr)
#   nObs <- length(hOUwie.dat$ObservedTraits)
#   #reorder phy
#   phy <- reorder.phylo(phy, "pruningwise")
#   # a way to speed up the three.point function
#   tip.paths <- lapply(1:length(phy$tip.label), function(x) OUwie:::getPathToRoot(phy, x))
#   # get the appropriate OU model structure
#   if(is.null(index.ou)){
#     index.ou <- getOUParamStructure(model.ou, "three.point", root.station, get.root.theta, dim(model.set.final$Q)[1])
#   }
#   # get the appropriate corhmm model structure
#   if(is.null(model.cor)){
#     model.cor <- "ER"
#   }
#   model.set.final <- corHMM:::rate.cat.set.corHMM.JDB(phy=phy,data=hOUwie.dat$data.cor,rate.cat=rate.cat, ntraits = nObs, model = model.cor)
#   # this allows for custom rate matricies!
#   if(!is.null(index.cor)){
#     order.test <- FALSE
#     index.cor[index.cor == 0] <- NA
#     rate <- index.cor
#     model.set.final$np <- max(rate, na.rm=TRUE)
#     rate[is.na(rate)]=max(rate, na.rm=TRUE)+1
#     model.set.final$rate <- rate
#     model.set.final$index.matrix <- index.cor
#     model.set.final$Q <- matrix(0, dim(index.cor)[1], dim(index.cor)[2])
#     ## for precursor type models ##
#     col.sums <- which(colSums(index.cor, na.rm=TRUE) == 0)
#     row.sums <- which(rowSums(index.cor, na.rm=TRUE) == 0)
#     drop.states <- col.sums[which(col.sums == row.sums)]
#     if(length(drop.states > 0)){
#       model.set.final$liks[,drop.states] <- 0
#     }
#     ## need to do anything to the ouwie matrix?
#   }
# 
#   p.mk <- sapply(1:max(index.cor, na.rm = TRUE), function(x) c(na.omit(c(solution.cor)))[match(x, c(na.omit(c(index.cor))))])
#   p.ou <- sapply(1:max(index.ou, na.rm = TRUE), function(x) c(na.omit(c(solution.ou)))[match(x, c(na.omit(c(index.ou))))])
#   p <- c(p.mk, p.ou)
#   ReconTable <- matrix(0, length(nodes.to.eval), dim(model.set.final$liks)[2], dimnames = list(nodes.to.eval))
#   for(i in 1:length(nodes.to.eval)){
#     fix.node <- nodes.to.eval[i]
#     ReconTable[i,] <- NodeSpecific.hOuwieRecon(pars=p, phy=phy, rate.cat=rate.cat, data.cor=hOUwie.dat$data.cor, liks=model.set.final$liks, Q=model.set.final$Q, rate=model.set.final$rate, root.p=root.p, data.ou=hOUwie.dat$data.ou, index.ou=index.ou, algorithm=algorithm, mserr=mserr, nSim=nMap, nCores=nCores, tip.paths=tip.paths, weighted=weighted, fix.node)
#   }
#   
#   if(rate.cat > 1){
#     StateNames <- paste("(", rep(1:nObs, rate.cat), rep(LETTERS[1:rate.cat], each = nObs), ")", sep = "")
#   }else{
#     StateNames <- paste("(", rep(1:nObs, rate.cat), ")", sep = "")
#   }
#   colnames(ReconTable) <- StateNames
#   
#   return(ReconTable)
# }

## houwie node specific regime reconstruction

## evaluates the regime specific state and node 
# NodeSpecific.hOuwieRecon <- function(pars, phy, rate.cat, data.cor, liks, Q, rate, root.p, data.ou, index.ou, algorithm, mserr, nSim, nCores, tip.paths, weighted, fix.node){
#   est.pars<-log(pars)
#   ReconVec <- liks[fix.node,]
#   # are we reconstructing a tip?
#   if(any(ReconVec == 1)){
#     for(i in which(ReconVec == 1)){
#       ReconVec[i] <- -hOUwie.dev(est.pars, phy=phy, rate.cat=rate.cat, data.cor=data.cor, liks=liks, Q=Q, rate=rate, root.p=root.p, data.ou=data.ou, index.ou=index.ou, algorithm=algorithm, mserr=mserr, nSim=nSim, nCores=nCores, tip.paths=tip.paths, weighted=weighted, order.test=FALSE, fix.node = fix.node, fix.state = i)
#     }
#   }else{
#     for(i in 1:length(ReconVec)){
#       ReconVec[i] <- -hOUwie.dev(est.pars, phy=phy, rate.cat=rate.cat, data.cor=data.cor, liks=liks, Q=Q, rate=rate, root.p=root.p, data.ou=data.ou, index.ou=index.ou, algorithm=algorithm, mserr=mserr, nSim=nSim, nCores=nCores, tip.paths=tip.paths, weighted=weighted, order.test=FALSE, fix.node = fix.node, fix.state = i)
#     }
#   }
#   ReconVec[ReconVec == 0] <- NA
#   ReconVec <- exp(ReconVec)/sum(exp(ReconVec), na.rm = TRUE)
#   ReconVec[is.na(ReconVec)] <- 0
#   return(ReconVec)
# }

## wrapper function used mainly for model averaging
getTipRecons <- function(models, nSim){
  model.recon <- vector("list", length(models))
  for(i in 1:length(models)){
    if(models[[i]]$rate.cat == 1){
      model.recon[[i]] <- corHMM:::rate.cat.set.corHMM.JDB(phy=models[[i]]$phy,data=models[[i]]$hOUwie.dat$data.cor,rate.cat=models[[i]]$rate.cat, ntraits = length(models[[i]]$hOUwie.dat$ObservedTraits), model = "ER")$liks[1:length(models[[i]]$phy$tip.label), ]
    }else{
      model.recon[[i]] <- hOUwieRecon("external", "marginal", nSim, 1, models[[i]])
    }
  }
  return(model.recon)
}

# model average some params of significance 
getModelAvgParams <- function(models, model.recon, AICwt){
  # ou params and corhmm params
  ou.pars <- vector("list", length(models))
  cor.pars <- vector("list", length(models))
  for(i in 1:length(models)){
    models[[i]]$solution.ou[is.na(models[[i]]$solution.ou)] <- 1e-10
    ou.pars[[i]] <- apply(models[[i]]$solution.ou, 1, function(x) colSums(x * t(model.recon[[i]]))) * AICwt[i]
    cor.pars[[i]] <- colSums(rowSums(models[[i]]$solution.cor, na.rm = TRUE) * t(model.recon[[i]]), na.rm = TRUE) * AICwt[i]
  }
  
  cor.pars <- 1/Reduce("+", cor.pars)
  ou.pars <- Reduce("+", ou.pars)
  
  ObsValue <- models[[1]]$hOUwie.dat$data.ou[,3]
  names(ObsValue) <- models[[1]]$hOUwie.dat$data.ou[,1]
  ObsValue[models[[1]]$phy$tip.label]
  
  pars <- data.frame(alpha = ou.pars[,1], 
                     sigma = ou.pars[,2],
                     optim = ou.pars[,3],
                     half.life = log(2)/ou.pars[,1], 
                     station.var = ou.pars[,2]/2*ou.pars[,1],
                     dist.2.opt = ou.pars[,3] - ObsValue,
                     wait.time = cor.pars, 
                     row.names = models[[1]]$phy$tip.label)
  return(pars)
}

## evaluates a set of hOUwie models 
fit.hOUwie.set <- function(phy, data, rate.cats, index.cor.list, index.ou.list, root.p="yang", lb.cor=NULL, ub.cor=NULL, root.station=FALSE, get.root.theta=FALSE, mserr = "none", lb.ou=NULL, ub.ou=NULL, p=NULL, ip=NULL, nSim=1000, opts=NULL, nCores=1, weighted=FALSE, quiet=TRUE, prune=TRUE){
  # checks
  if(var(c(length(rate.cats), length(index.cor.list), length(index.ou.list))) != 0){
    stop("The number of rate categories, discrete character models, and continuous character models don't match.", call. = FALSE)
  }
  nModels <- length(rate.cats)
  model.fits <- vector("list", nModels)
  # set the model names
  if(is.null(names(rate.cats))){
    if(is.null(names(index.cor.list))){
      if(is.null(names(index.ou.list))){
        names(model.fits) <- paste0("Model", 1:nModels)
      }else{
        names(model.fits) <- names(index.ou.list)
      }
    }else{
      names(model.fits) <- names(index.cor.list)
    }
  }else{
    names(model.fits) <- names(rate.cats)
  }
  # fit the models
  cat("Begining optimization of", nModels, "hOUwie models with", nSim, "simmaps...\n")
  for(i in 1:nModels){
    model.fits[[i]] <- hOUwie(phy, data, rate.cats[i], index.cor=index.cor.list[[i]], index.ou=index.ou.list[[i]], root.p = root.p, lb.cor=lb.cor, ub.cor=ub.cor, root.station=root.station, get.root.theta=get.root.theta, mserr=mserr, lb.ou=lb.ou, ub.ou=ub.ou, p=p, ip=ip, nSim=nSim, opts=opts, nCores=nCores, weighted=weighted, quiet=quiet)
    cat("\r", paste0(round(i/nModels*100), "%"), "complete...           ")
  }
  cat("\n Model fitting complete. Now assessing models and averaging parameters.\n")
  
  # check for any redundant models
  if(prune){
    models.pruned <- PruneRedundantModels(model.fits)
  }else{
    models.pruned <- model.fits
  }
  
  # evaluate the AIC and summarize the models
  model.table <- getModelTable(models.pruned, "AICc")
  
  # reconstruct the ancestors
  model.recon <- getTipRecons(models.pruned, nSim)
  
  # average the paramaters of the models - slow because requires tip reconstruction
  model.pars <- getModelAvgParams(models.pruned, model.recon, model.table[,5])
  
  cat("Done.\n")
  obj = list(model.fits=model.fits, 
             models.pruned = models.pruned,
             model.table = model.table,
             model.recon = model.recon,
             model.pars = model.pars)
  class(obj) <- "houwie.set"
  return(obj)
}

summarize.hOUwie.set <- function(model.fits, prune=TRUE, nSim = NULL){
  if(is.null(nSim)){
    nSim <- round(mean(unlist(lapply(model.fits, "[[", "nSim"))))
  }
  
  # check for any redundant models
  if(prune){
    models.pruned <- PruneRedundantModels(model.fits)
  }else{
    models.pruned <- model.fits
  }
  
  # evaluate the AIC and summarize the models
  model.table <- getModelTable(models.pruned, "AICc")
  
  # reconstruct the ancestors
  model.recon <- getTipRecons(models.pruned, nSim)
  
  # average the paramaters of the models - slow because requires tip reconstruction
  model.pars <- getModelAvgParams(models.pruned, model.recon, model.table[,5])
  
  cat("Done.\n")
  obj = list(model.fits = model.fits, 
             models.pruned = models.pruned,
             model.table = model.table,
             model.recon = model.recon,
             model.pars = model.pars)
  class(obj) <- "houwie.set"
  return(obj)
}

print.houwie.set <- function(x, ...){
  x$model.table[,2] <- round(x$model.table[,2], 3)
  x$model.table[,3] <- round(x$model.table[,3], 3)
  x$model.table[,4] <- round(x$model.table[,4], 3)
  x$model.table[,5] <- round(x$model.table[,5], 3)
  cat("Model Fit Table\n")
  print(x$model.table)
  cat("\n")
}


# forward simulation of a branch. kind of a cheap way to do discrete time stuff by just shrinking time, but the markov property should ensure that this works just fine
forwardSimBranch <- function(theta0, state0, Q, alpha, sigma, theta, time.unit, time.steps){
  nStates <- dim(Q)[1]
  state.vec <- c(state0, rep(NA, time.steps-1))
  conti.vec <- c(theta0, rep(NA, time.steps-1))
  MkPMat <- matrix(NA, nStates, nStates)
  for(i in 1:nStates){
    tmp_vec_i <- rep(0, nStates)
    tmp_vec_i[i] <- 1
    MkPMat[i,] <- c(expm(Q * time.unit) %*% tmp_vec_i)
  }
  
  for(i in 2:time.steps){
    state_i <- state.vec[i-1]
    conti_i <- conti.vec[i-1]
    conti.vec[i] <- rOU(1, conti_i, time.unit, alpha[state_i], theta[state_i], sigma[state_i])
    state.vec[i] <- sample.int(nStates, size = 1, prob = MkPMat[state_i,])
  }
  out <- data.frame(D = state.vec, C = conti.vec)
  return(out)
}