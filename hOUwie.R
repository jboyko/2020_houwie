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
#'@param weighted a logicial indicating whether to take the average OU likelihood over all simmaps (FALSE) or just take the maximum likelihood
#'@param quiet a logical indicating whether to output user messages
hOUwie <- function(phy, data, rate.cat, nBins=10, collapse=TRUE,type="joint",
                   model.cor=NULL, index.cor=NULL, root.p="yang", lb.cor=NULL, ub.cor=NULL,
                   model.ou=NULL, index.ou=NULL, root.station=FALSE, get.root.theta=FALSE, mserr = "none", lb.ou=NULL, ub.ou=NULL, order.test=TRUE, nSim=1000, 
                   p=NULL, ip=NULL, opts=NULL, nCores=1, weighted=FALSE, quiet=FALSE, parsimony=FALSE, split.LnLiks=TRUE, all.roots=FALSE, sample.recons=FALSE, shift.point=0.5){
  # check that tips and data match
  # check for invariance of tip states and not that non-invariance isn't just ambiguity
  if(is.null(model.cor) & is.null(index.cor)){
    stop("No model for discrete evolution has been specified. Please select either a default model with the model.cor argument or specify an index matrix with index.cor.", call. = FALSE)
  }
  
  if(is.null(model.ou) & is.null(index.ou)){
    stop("No model for continuous evolution has been specified. Please select either a default model with the model.ou argument or specify an index matrix with index.ou.", call. = FALSE)
  }
  
  if(!is.null(index.cor) & !is.null(index.ou)){
    if(dim(index.cor)[2] > dim(index.ou)[2]){
      stop("Not all of your discrete states have OU parameters associated with them. Please check that your discrete index matrix (index.cor) matches your continuous index matrix (index.ou).")
    }
    if(dim(index.ou)[2] > dim(index.cor)[2]){
      stop("You have specified more OU parameters than there are states in the discrete process. Please check that your discrete index matrix (index.cor) matches your continuous index matrix (index.ou).")
    }
  }
  if(!is.null(phy$node.label)){
    if(!quiet){
      cat("Your phylogeny had node labels, these have been removed.\n")
    }
    phy$node.label <- NULL
  }
  
  
  Tmax <- max(branching.times(phy))
  if(mserr == "none"){
    nDiscrete <- dim(data)[2] - 2
  }else{
    nDiscrete <- dim(data)[2] - 3
  }
  if(is.null(lb.ou)){
    # the lower limit of alpha is defined as a halflife of 10000% of the max tree height
    # the lower limit of sigma is defined 10 times less than alpha
    # the lower limit of optim is defined 10 times lower than the minimum observation
    lb.alpha = log(2)/(100 * Tmax)
    lb.sigma = lb.alpha/10
    lb.optim = min(data[, 1+nDiscrete+1])/10 
    lb.ou=c(lb.alpha,lb.sigma,lb.optim)
  }
  if(is.null(ub.ou)){
    # the upper limit of alpha is defined as a halflife of 1% of the max tree height
    # the upper limit of sigma is defined 10 times more than alpha
    # the upper limit of optim is defined 10 times more than the maximum observation
    ub.alpha = log(2)/(0.01 * Tmax)
    ub.sigma = 10 * ub.alpha
    ub.optim = max(data[, 1+nDiscrete+1])*10 
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
  
  
  algorithm="three.point"
  # organize the data into the corHMM data and the OUwie data
  # TO DO: add a way to shift negative continuous variables to positive then shift back
  hOUwie.dat <- organizeHOUwieDat(data, mserr)
  organizedData <- getHOUwieCombosAndData(data, rate.cat, collapse, nBins)
  nObs <- length(hOUwie.dat$ObservedTraits)
  #reorder phy
  phy <- reorder(phy, "pruningwise")
  # a way to speed up the three.point function
  tip.paths <- lapply(1:length(phy$tip.label), function(x) OUwie:::getPathToRoot(phy, x))
  
  #scale the tree to a root height of 1
  # phy$edge.length <- phy$edge.length/max(branching.times(phy))
  # a temporary corhmm model to set the rates up
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
    if(max(index.ou, na.rm = TRUE) + max(model.set.final$index.matrix, na.rm = TRUE) != length(p)){
      message <- paste0("The number of parameters does not match the number required by the model structure. You have supplied ", length(p), ", but the model structure requires ", max(index.ou, na.rm = TRUE) + max(model.set.final$index.matrix, na.rm = TRUE), ".")
      stop(message, call. = FALSE)
      
    }
    out<-NULL
    est.pars<-log(p)
    out$solution <- log(p)
    out$objective <- hOUwie.dev(est.pars, phy=phy, rate.cat=rate.cat,organizedData=organizedData, data.cor=hOUwie.dat$data.cor, liks=model.set.final$liks, Q=model.set.final$Q, rate=model.set.final$rate, root.p=root.p, data.ou=hOUwie.dat$data.ou, index.ou=index.ou, algorithm=algorithm, mserr=mserr,nSim=nSim, nCores=nCores, tip.paths=tip.paths, weighted=weighted, order.test=order.test, fix.node=NULL, fix.state=NULL, parsimony = parsimony,type=type, split.LnLiks=FALSE, all.roots=all.roots, sample.recons=sample.recons, shift.point=shift.point)
  }else{
    if(!quiet){
      if(type == "simmap"){
        cat("Starting a search of parameters with", nSim, "simmaps...\n")
      }
      if(type == "joint"){
        cat("Starting a search of parameters with joint character reconstructions...\n")
      }
    }
    out<-NULL
    # check for user input initial parameters 
    if(is.null(ip)){
      means.by.regime <- with(hOUwie.dat$data.ou, tapply(hOUwie.dat$data.ou[,3], hOUwie.dat$data.ou[,2], mean))
      if(length(unique(na.omit(index.ou[3,]))) == length(means.by.regime)){
        start.theta <- rep(means.by.regime, length(unique(index.ou[3,]))/length(means.by.regime))
      }else{
        start.theta <- rep(mean(hOUwie.dat$data.ou[,3]), length(unique(na.omit(index.ou[3,]))))
      }
      start.cor <- rep(10/sum(phy$edge.length), model.set.final$np)
      start.ou <- c(rep(log(2)/Tmax, length(unique(na.omit(index.ou[1,])))), 
                    rep(var(hOUwie.dat$data.ou[,3]), length(unique(na.omit(index.ou[2,])))), 
                    start.theta)
      starts = c(start.cor, start.ou)
    }else{
      starts <- ip
    }
    lower = log(c(rep(lb.cor, model.set.final$np), 
                  rep(lb.ou[1], length(unique(na.omit(index.ou[1,])))), 
                  rep(lb.ou[2], length(unique(na.omit(index.ou[2,])))), 
                  rep(lb.ou[3], length(unique(na.omit(index.ou[3,]))))))
    upper = log(c(rep(ub.cor, model.set.final$np), 
                  rep(ub.ou[1], length(unique(na.omit(index.ou[1,])))), 
                  rep(ub.ou[2], length(unique(na.omit(index.ou[2,])))), 
                  rep(ub.ou[3], length(unique(na.omit(index.ou[3,]))))))
    cat("Loglik:\n")
    out = nloptr(x0=log(starts), eval_f=hOUwie.dev, lb=lower, ub=upper, opts=opts, phy=phy, rate.cat=rate.cat, organizedData=organizedData, data.cor=hOUwie.dat$data.cor, liks=model.set.final$liks, Q=model.set.final$Q, rate=model.set.final$rate, root.p=root.p, data.ou=hOUwie.dat$data.ou, index.ou=index.ou, algorithm=algorithm, mserr=mserr, nSim=nSim, nCores=nCores, tip.paths=tip.paths, weighted=weighted, order.test=order.test, fix.node=NULL, fix.state=NULL, parsimony = parsimony, type=type, split.LnLiks=FALSE, all.roots=all.roots, sample.recons=sample.recons, shift.point=shift.point)
  }
  # preparing output
  # params are independent corhmm rates, alpha, sigma, theta, and 1 intercept
  if(split.LnLiks == TRUE){
    p <- out$solution
    SplitLiks <- hOUwie.dev(p, phy=phy, rate.cat=rate.cat,organizedData=organizedData, data.cor=hOUwie.dat$data.cor, liks=model.set.final$liks, Q=model.set.final$Q, rate=model.set.final$rate, root.p=root.p, data.ou=hOUwie.dat$data.ou, index.ou=index.ou, algorithm=algorithm, mserr=mserr,nSim=nSim, nCores=nCores, tip.paths=tip.paths, weighted=weighted, order.test=order.test, fix.node=NULL, fix.state=NULL, parsimony = parsimony,type=type, split.LnLiks=TRUE, all.roots=all.roots, sample.recons=sample.recons, shift.point=shift.point)
  }else{
    SplitLiks <- NULL
  }
  if(!quiet){
    cat("\nFinished.\n")
  }
  if(null.cor){
    model.cor <- NULL
  }
  param.count <- max(model.set.final$index.matrix, na.rm = TRUE) + max(index.ou, na.rm = TRUE)
  nb.tip <- length(phy$tip.label)
  solution <- organizeHOUwiePars(out=out, rate=model.set.final$rate, Q=model.set.final$Q, index.ou=index.ou)
  if(rate.cat > 1){
    StateNames <- paste("(", rep(1:nObs, rate.cat), rep(LETTERS[1:rate.cat], each = nObs), ")", sep = "")
  }else{
    StateNames <- paste("(", rep(1:nObs, rate.cat), ")", sep = "")
  }
  rownames(solution$solution.cor) <- colnames(solution$solution.cor) <- StateNames
  colnames(solution$solution.ou) <- StateNames
  names(hOUwie.dat$ObservedTraits) <- 1:length(hOUwie.dat$ObservedTraits)
  obj <- list(
    loglik = -out$objective,
    OULnLik = SplitLiks[1], 
    MkLnLik = SplitLiks[2],
    AIC = 2*out$objective + 2*param.count,
    AICc = 2*out$objective+ 2*param.count*(param.count/(nb.tip-param.count-1)),
    BIC = 2*out$objective + log(nb.tip) * param.count,
    param.count = param.count,
    solution.cor = solution$solution.cor,
    solution.ou = solution$solution.ou,
    phy = phy,
    legend = hOUwie.dat$ObservedTraits,
    data = data, 
    hOUwie.dat = hOUwie.dat,
    rate.cat = rate.cat, 
    model.cor=model.cor, 
    index.cor=index.cor, 
    root.p=root.p, 
    lb.cor=lb.cor, 
    ub.cor=ub.cor,
    model.ou=model.ou, 
    index.ou=index.ou, 
    root.station=root.station, 
    get.root.theta=get.root.theta, 
    mserr = mserr, 
    lb.ou=lb.ou, 
    ub.ou=ub.ou,
    p = exp(out$solution),
    ip=ip, 
    nSim=nSim, 
    nBins=nBins,
    collapse=collapse,
    opts=opts, 
    nCores=nCores, 
    weighted=weighted, 
    quiet=quiet
  )
  class(obj) <- "houwie"
  return(obj)
}

# for a single set of parameters, evaluate the hOUwie likelihood
hOUwie.dev <- function(p, phy, rate.cat, organizedData,type="joint",
                       data.cor, liks, Q, rate, root.p,
                       data.ou, index.ou, algorithm, mserr,
                       nSim, nCores, tip.paths=NULL, weighted = FALSE, order.test=TRUE,parsimony=FALSE,
                       split.LnLiks=FALSE,fix.node=NULL, fix.state=NULL, all.roots=FALSE, sample.recons=FALSE, shift.point=0.5){
  # params are given in log form
  p <- exp(p)
  # define which params are for the HMM
  k <- max(rate)-1
  p.mk <- p[1:k]
  # set the OU params
  p.ou <- p[(k+1):length(p)] 
  Rate.mat <- matrix(1, 3, dim(rate)[2])
  index.ou[is.na(index.ou)] <- max(index.ou, na.rm = TRUE) + 1
  Rate.mat[] <- c(p.ou, 1e-10)[index.ou]
  alpha = Rate.mat[1,]
  sigma.sq = Rate.mat[2,]
  theta = Rate.mat[3,]
  if(algorithm == "invert"){
    theta = NULL
  }
  # if there are more than one rate categories we test the order to help consistency
  if(rate.cat > 1 & order.test){
    # for each state and for each parameter ensure that parameter values for A > B > C etc..
    order.check.vec <- vector("logical", organizedData$nStates/rate.cat)
    for(i in 1:(organizedData$nStates/rate.cat)){
      Rate.mat_i <- Rate.mat[,grep(i, names(organizedData$namedStates))]
      order.check.vec[i] <- check.order(Rate.mat_i)
    }
    if(!all(order.check.vec)){
      return(1e10)
    }
  }
  # fit the corHMM model. if rate.cat > 1, then ensure the order of the state mats is fast to slow.
  # returns negative loglikelihood (hence sign change)
  Mk.loglik <- -corHMM:::dev.corhmm(log(p.mk), phy, liks, Q, rate, root.p, rate.cat, order.test, FALSE)
  # set up the rate matrix
  Q[] <- c(p.mk, 0)[rate]
  diag(Q) <- -rowSums(Q)
  # fit the ancestral state reconstruction
  # simulate a set of simmaps
  if(type == "pupko"){
    index.cor <- rate
    index.cor[index.cor == max(index.cor)] <- 0
    nTip <- length(phy$tip.label)
    JointRecon <- hOUwieRecon.dev(log(p), phy=phy, organizedData=organizedData, rate.cat=rate.cat, index.cor=index.cor, index.ou=index.ou, all.roots = all.roots, sample.recons = sample.recons, shift.point=shift.point, return.lik = TRUE)
    cat("LnLik:", JointRecon, "pars:", p, "\n")
    return(-JointRecon)
  }
  if(type == "joint"){
    index.cor <- rate
    index.cor[index.cor == max(index.cor)] <- 0
    nTip <- length(phy$tip.label)
    if(sample.recons){
      JointRecon <- hOUwieRecon.dev(log(p), phy=phy, organizedData=organizedData, rate.cat=rate.cat, index.cor=index.cor, index.ou=index.ou, all.roots = all.roots, sample.recons = sample.recons, shift.point=shift.point)
      JointLiks <- JointRecon[[length(JointRecon)]]
      JointReconOnly <- JointRecon[-length(JointRecon)]
      map <- vector("list", length(JointReconOnly))
      for(ReconIndex in 1:length(JointReconOnly)){
        tip.states <- organizedData$AllCombos[JointReconOnly[[ReconIndex]][1:nTip], 1] 
        node.states <- organizedData$AllCombos[JointReconOnly[[ReconIndex]][(nTip+1):max(phy$edge[,1])], 1]
        map[[ReconIndex]] <- getJointMapFromNode(phy, tip.states, node.states, 0)
      }
      OU.loglik <- lapply(map, function(x) OUwie.basic(x, data.ou, simmap.tree=TRUE, scaleHeight=FALSE, alpha=alpha, sigma.sq=sigma.sq, theta=theta, algorithm=algorithm, tip.paths=tip.paths, mserr=mserr))
    }else{
      if(!all.roots){
        JointRecon <- hOUwieRecon.dev(log(p), phy=phy, organizedData=organizedData, rate.cat=rate.cat, index.cor=index.cor, index.ou=index.ou, all.roots = all.roots, sample.recons = sample.recons, shift.point=shift.point)
        tip.states <- organizedData$AllCombos[JointRecon$Recon[1:nTip], 1] 
        node.states <- organizedData$AllCombos[JointRecon$Recon[(nTip+1):max(phy$edge[,1])], 1]
        # if not all hidden states are present there is no need for a hidden states model
        if(!all(organizedData$namedStates %in% tip.states)){
          return(1e10)
        }
        #cat("Root state:", node.states[1], "Root index:", JointRecon$Recon[1], "\n")
        map <- getJointMapFromNode(phy, tip.states, node.states, 0)
        OU.loglik <- OUwie.basic(map, data.ou, simmap.tree=TRUE, scaleHeight=FALSE, alpha=alpha, sigma.sq=sigma.sq, theta=theta, algorithm=algorithm, tip.paths=tip.paths, mserr=mserr)
      }else{
        JointRecon <- hOUwieRecon.dev(log(p), phy=phy, organizedData=organizedData, rate.cat=rate.cat, index.cor=index.cor, index.ou=index.ou, all.roots = all.roots, sample.recons = sample.recons, shift.point=shift.point)
        JointLiks <- JointRecon[[length(JointRecon)]]
        JointReconOnly <- JointRecon[-length(JointRecon)]
        map <- vector("list", length(JointReconOnly))
        for(ReconIndex in 1:length(JointReconOnly)){
          tip.states <- organizedData$AllCombos[JointReconOnly[[ReconIndex]][1:nTip], 1] 
          node.states <- organizedData$AllCombos[JointReconOnly[[ReconIndex]][(nTip+1):max(phy$edge[,1])], 1]
          map[[ReconIndex]] <- getJointMapFromNode(phy, tip.states, node.states, 0)
        }
        OU.loglik <- lapply(map, function(x) OUwie.basic(x, data.ou, simmap.tree=TRUE, scaleHeight=FALSE, alpha=alpha, sigma.sq=sigma.sq, theta=theta, algorithm=algorithm, tip.paths=tip.paths, mserr=mserr))
        # # addig in a root prior
        # ScaledJointLiks <- (JointLiks - max(JointLiks))
        # ScaledJointLiks <- exp(ScaledJointLiks)/sum(exp(ScaledJointLiks))
        # OU.loglik <- unlist(OU.loglik) + log(ScaledJointLiks)
      }
    }
  }
  if(type == "simmap"){
    simmap <- corHMM:::makeSimmap(phy, data.cor, Q, rate.cat, root.p = root.p, nSim = nSim, fix.node = fix.node, fix.state = fix.state, parsimony=parsimony)
    # fit the OU models to the simmaps
    # returns log likelihood
    OU.loglik <- mclapply(simmap, function(x) OUwie.basic(x, data.ou, simmap.tree=TRUE, scaleHeight=FALSE, alpha=alpha, sigma.sq=sigma.sq, theta=theta, algorithm=algorithm, tip.paths=tip.paths, mserr=mserr), mc.cores = nCores)
  }
  
  # to show that OUwie.basic is identical to OUwie.fixed
  # OUwie.basic(simmap[[1]], data.ou, simmap.tree=TRUE, scaleHeight=FALSE, alpha=alpha, sigma.sq=sigma.sq, theta=theta, algorithm=algorithm, tip.paths=tip.paths, mserr=mserr)
  # OUwie.fixed(simmap[[1]], data.ou, model = "OUM", simmap.tree=TRUE, scaleHeight=FALSE, alpha=alpha, sigma.sq=sigma.sq, theta=theta, algorithm="invert", tip.paths=tip.paths, mserr=mserr)
  # OUwie.fixed(simmap[[1]], data.ou, model = "OUM", simmap.tree=TRUE, scaleHeight=FALSE, alpha=alpha, sigma.sq=sigma.sq, theta=theta, algorithm="three.point", tip.paths=tip.paths)
  
  if(weighted == TRUE){
    # get the likelihoods of the simmaps
    OU.loglik <- max(unlist(OU.loglik))
    # OU.loglik <- log(mean(exp(unlist(OU.loglik)-comp)))+comp
    return(-(OU.loglik + Mk.loglik))
  }else{
    # algebraic trick to prevent overflow and underflow while preserving as many accurate leading digits in the result as possible. The leading digits are preserved by pulling the maximum outside. The arithmetic is robust because subtracting the maximum on the inside makes sure that only negative numbers or zero are ever exponentiated, so there can be no overflow on those calculations. If there is underflow, we know the leading digits have already been returned as part of the max term on the outside. subtract log nsim for avg.
    OU.loglik <- unlist(OU.loglik) 
    OU.loglik <- max(OU.loglik) + log(sum(exp(OU.loglik - max(OU.loglik)))) - log(length(OU.loglik))
    # example
    # tmp <- c(-100, -200, -300, -100, -250)
    # max(tmp) + log(sum(exp(tmp - max(tmp)))) - log(length(tmp))
    # log(mean(exp(tmp)))
    # 
    llik <- OU.loglik + Mk.loglik
    # cat("Mk:", Mk.loglik, "OU:", OU.loglik, "\n")
    cat("LnLik:", llik, "Mk:", Mk.loglik, "OU:", OU.loglik, "pars:", p, "\n")
    # cat("\r", llik)
    if(split.LnLiks == FALSE){
      return(-llik)
    }else{
      return(c(OULnLik = OU.loglik, MkLnLik = Mk.loglik))
    }
  }
}

# checks the order of matrix left to right
check.order <- function(mat, decreasing = FALSE){
  order.vec <- vector("logical", dim(mat)[1])
  for(i in 1:dim(mat)[1]){
    order.vec[i] <- identical(mat[i,], sort(mat[i,], decreasing = decreasing))
  }
  return(all(order.vec))
}


# organize data for the joint reconstruction
getHOUwieCombosAndData <- function(data, rate.cat, collapse, nBins){
  # organize the data
  nCol <- dim(data)[2]
  # process the corHMM data to combine 
  discrete.data <- corHMM:::corProcessData(data[c(1:(nCol-1))], collapse = collapse)
  hOUwie.dat <- cbind(discrete.data$corData, x = data[dim(data)[2]])
  # the number of discrete states
  # the number of discrete bins of the continuous trait findInterval is another possibility
  bin_index <- cut(hOUwie.dat[,3], nBins, labels = FALSE)
  bin <- cut(hOUwie.dat[,3], nBins)
  namesBin <- levels(bin)
  valueBin <- unlist(lapply(strsplit(namesBin, ","), function(x) mean(as.numeric(gsub("\\(|\\]", "", x)))))
  boundBin <- do.call(rbind, lapply(strsplit(namesBin, ","), function(x) c(as.numeric(gsub("\\(|\\]", "", x)))))
  hOUwie.dat <- cbind(hOUwie.dat, bin = bin, valueBin = valueBin[bin_index], lowerBin = boundBin[bin_index,1], upperBin = boundBin[bin_index,2])
  # find all possible state combinations once the bins and states have been defined
  nStates <- as.numeric(max(discrete.data$corData[,2])) * rate.cat
  AllCombos <- expand.grid(list(State = 1:nStates, Bin = valueBin))
  # nCombos <- dim(AllCombos)[1]
  nObs <- nStates/rate.cat
  ObsStateMatrix <- matrix(1:nStates, nObs, rate.cat, dimnames = list(1:nObs, LETTERS[1:rate.cat]))
  namedStates <- 1:nStates
  names(namedStates) <- paste("(", rep(1:nObs, rate.cat), rep(LETTERS[1:rate.cat], each = nObs), ")", sep = "")
  FullDat <- list(nStates = nStates,
                  nObs = nObs, 
                  PossibleTraits = discrete.data$PossibleTraits,
                  ObservedTraits = discrete.data$ObservedTraits,
                  namedStates = namedStates,
                  ObsStateMatrix = ObsStateMatrix,
                  discreteData = discrete.data$corData,
                  hOUwieData = hOUwie.dat,
                  AllCombos = AllCombos)
  return(FullDat)
}

## takes a node based reconstruction and returns a map (identical to a map from simmap)
getJointMapFromNode <- function(phy, tipstates, nodestates, shift.point){
  Map <- vector("list", dim(phy$edge)[1])
  Data <- c(tipstates, nodestates)
  NodeStates <- cbind(Data[phy$edge[,1]], Data[phy$edge[,2]])
  for(i in 1:dim(phy$edge)[1]){
    from <- as.character(NodeStates[i,1])
    to <- as.character(NodeStates[i,2])
    if(from == to){
      tmp <- phy$edge.length[i]
      names(tmp) <- from
      Map[[i]] <- tmp
    }else{
      shift.time <- shift.point * phy$edge.length[i]
      tmp <- c(phy$edge.length[i] - shift.time, shift.time)
      names(tmp) <- c(from, to)
      Map[[i]] <- tmp
    }
  }
  mapped.edge <- corHMM:::convertSubHistoryToEdge(phy, Map)
  phy$maps <- Map
  phy$mapped.edge <- mapped.edge
  attr(phy, "map.order") <- "right-to-left"
  if (!inherits(phy, "simmap")) 
    class(phy) <- c("simmap", setdiff(class(phy), "simmap"))
  return(phy)
}

# joint recon warapper funciton 
hOUwieRecon <- function(hOUwie.model=NULL, phy=NULL, data=NULL, rate.cat=NULL, nSim=1000, nBins=10, collapse=TRUE,type="joint",model.cor=NULL, index.cor=NULL, root.p="yang", lb.cor=NULL, ub.cor=NULL, model.ou=NULL, index.ou=NULL, root.station=FALSE, get.root.theta=FALSE, mserr = "none", lb.ou=NULL, ub.ou=NULL, order.test=TRUE, p=NULL, ip=NULL, opts=NULL, nCores=1, weighted=FALSE, quiet=FALSE, parsimony=FALSE, split.LnLiks=TRUE){
  # check that tips and data match
  # check for invariance of tip states and not that non-invariance isn't just ambiguity
  if(!is.null(hOUwie.model)){
    if(class(hOUwie.model) != "houwie"){
      stop("hOUwie model has been given, but is not of class hOUwie. Either give a model that has been fit which is of class hOUwie or the necessary prerequisites for a state reconstruction.", call. = FALSE)
    }
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
    p <- hOUwie.model$p
  }else{
    if(is.null(model.cor) & is.null(index.cor)){
      stop("No model for discrete evolution has been specified. Please select either a default model with the model.cor argument or specify an index matrix with index.cor.", call. = FALSE)
    }
    
    if(is.null(model.ou) & is.null(index.ou)){
      stop("No model for continuous evolution has been specified. Please select either a default model with the model.ou argument or specify an index matrix with index.ou.", call. = FALSE)
    }
    
    if(!is.null(index.cor) & !is.null(index.ou)){
      if(dim(index.cor)[2] > dim(index.ou)[2]){
        stop("Not all of your discrete states have OU parameters associated with them. Please check that your discrete index matrix (index.cor) matches your continuous index matrix (index.ou).")
      }
      if(dim(index.ou)[2] > dim(index.cor)[2]){
        stop("You have specified more OU parameters than there are states in the discrete process. Please check that your discrete index matrix (index.cor) matches your continuous index matrix (index.ou).")
      }
    }
    if(!is.null(phy$node.label)){
      if(!quiet){
        cat("Your phylogeny had node labels, these have been removed.\n")
      }
      phy$node.label <- NULL
    }
  }
  
  Tmax <- max(branching.times(phy))
  if(mserr == "none"){
    nDiscrete <- dim(data)[2] - 2
  }else{
    nDiscrete <- dim(data)[2] - 3
  }
  
  #Ensures that weird root state probabilities that do not sum to 1 are input:
  if(!is.null(root.p)){
    if(!is.character(root.p)){
      root.p <- root.p/sum(root.p)
    }
  }
  
  algorithm="three.point"
  # organize the data into the corHMM data and the OUwie data
  # TO DO: add a way to shift negative continuous variables to positive then shift back
  hOUwie.dat <- organizeHOUwieDat(data, mserr)
  organizedData <- getHOUwieCombosAndData(data, rate.cat, collapse, nBins)
  nObs <- length(hOUwie.dat$ObservedTraits)
  #reorder phy
  phy <- reorder(phy, "pruningwise")
  # a way to speed up the three.point function
  tip.paths <- lapply(1:length(phy$tip.label), function(x) OUwie:::getPathToRoot(phy, x))
  
  #scale the tree to a root height of 1
  # phy$edge.length <- phy$edge.length/max(branching.times(phy))
  # a temporary corhmm model to set the rates up
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
  # default MLE search options
  # p is organized into 2 groups with the first set being corHMM and the second set being OUwie
  # organized as c(trans.rt, alpha, sigma.sq, theta)
  # evaluate likelihood
  index.cor <- model.set.final$rate
  index.cor[index.cor == max(index.cor)] <- 0
  JointRecon <- hOUwieRecon.dev(log(p), phy=phy, organizedData=organizedData, rate.cat=rate.cat, index.cor=index.cor, index.ou=index.ou)
  nTip <- length(phy$tip.label)
  tip.states <- organizedData$AllCombos[JointRecon$Recon[1:nTip], ] 
  rownames(tip.states) <- 1:nTip
  node.states <- organizedData$AllCombos[JointRecon$Recon[(nTip+1):max(phy$edge[,1])], ]
  rownames(node.states) <- (nTip+1):max(phy$edge[,1])
  # preparing output
  # params are independent corhmm rates, alpha, sigma, theta, and 1 intercept
  if(rate.cat > 1){
    StateNames <- paste("(", rep(1:nObs, rate.cat), rep(LETTERS[1:rate.cat], each = nObs), ")", sep = "")
  }else{
    StateNames <- paste("(", rep(1:nObs, rate.cat), ")", sep = "")
  }
  obj <- list(TipStates = data.frame(TipState=tip.states[,1], TipStateName=StateNames[tip.states[,1]], Bin=tip.states[,2]),
              NodeStates = data.frame(NodeState=node.states[,1], NodeStateName=StateNames[node.states[,1]], Bin=node.states[,2]))
  return(obj)
}

# the joint reconstruction for a continuous and discrete variable
hOUwieRecon.dev <- function(p, phy, organizedData, rate.cat, index.cor, index.ou, return.lik=FALSE, all.roots=FALSE, sample.recons=FALSE, nSamples=100, shift.point=0){
  # organizing the parameters
  p <- exp(p)
  # define which params are for the HMM
  k <- max(index.cor)
  p.mk <- p[1:k]
  Q <- index.cor
  Q[Q!=0] <- p.mk[index.cor]
  diag(Q) <- -rowSums(Q)
  # set the OU params
  p.ou <- p[(k+1):length(p)] 
  Rate.mat <- matrix(1, 3, organizedData$nStates)
  index.ou[is.na(index.ou)] <- max(index.ou, na.rm = TRUE) + 1
  Rate.mat[] <- c(p.ou, 1e-10)[index.ou]
  alpha = Rate.mat[1,]
  sigma.sq = Rate.mat[2,]
  theta = Rate.mat[3,]
  
  nTip <- length(phy$tip.label)
  # the pupko et al algorithms the calculation of two quantities L and C
  # L can be stored in an array since its dimension will be the same for every node 
  # each entry in the L array contains the joint likelihood of a particular parental transition
  # there will be nState by nBin potential combinations of states, therefore each entry in L will be nBin by nState
  nCombos <- dim(organizedData$AllCombos)[1]
  nStates <- organizedData$nStates
  AllCombos <- organizedData$AllCombos
  ObsStateMatrix <- organizedData$ObsStateMatrix
  dat <- organizedData$hOUwieData
  L_z <- array(data=-Inf, dim = c(nCombos, nCombos, dim(phy$edge)[1]))
  # each evaluation of a particular L_z(i) will have a best joint likelihood and that will be defined in C_z
  C_z <- array(data=0, dim = c(nCombos, nCombos, dim(phy$edge)[1]))
  
  # step 1 of pupko: initialize the tips
  for(tip_index in 1:nTip){
    focal_edge_index <- match(tip_index, phy$edge[,2])
    RootAsAnc <- phy$edge[focal_edge_index,1] == nTip+1
    focal_data <- dat[match(phy$tip.label[tip_index], organizedData$hOUwieData[,1]),]
    time_edge <- phy$edge.length[focal_edge_index]
    edge_Mk <- matrix(0, length(unique(AllCombos[,1])), length(unique(AllCombos[,1])), dimnames = list(unique(AllCombos[,1]), unique(AllCombos[,1])))
    # for(Mk_index in 1:nStates){
    #   state <- unique(AllCombos[,1])[Mk_index]
    #   Mk_vec <- rep(0, nStates)
    #   Mk_vec[state] <- 1
    #   edge_Mk[Mk_index,] <- log(c(expm(Q * time_edge) %*% Mk_vec))
    # }
    # P_mat <- log((expm(Q * time_edge)))
    for(i in 1:nCombos){
      for(j in 1:rate.cat){
        obsState <- ObsStateMatrix[dat[tip_index,2],j]
        value_i <- which(AllCombos[,1] == obsState & AllCombos[,2] == focal_data$valueBin) # the observed all combo
        # calculate the Pij(ty) for each possible starting state i
        # start with the Mk likelihood
        state_i <- AllCombos[i,1]
        if(RootAsAnc){
          # init_i <- theta[state_i]
          init_i <- AllCombos[i,2]
        }else{
          init_i <- AllCombos[i,2]
        }
        if(state_i == obsState){
          Mk_llik <- log(1 - (1 - exp(time_edge * Q[obsState,obsState])))
        }else{
          Mk_llik <- log(Q[state_i,obsState]) + log(1 - (1 - exp(time_edge * Q[obsState,obsState])))
        } # P_ij for the Mk process. the probability that we transition at time 0 from i to j then remain in branch j for the entire branch length. the probability of a transtion of time 0 is dexp(0, rate), which is just rate.
        # next calculate the OU likelihood
        OU_llik <- dOU(dat[tip_index,3], init_i, time_edge, alpha[obsState], theta[obsState], sqrt(sigma.sq[obsState]), TRUE)
        
        L_z[i,value_i,focal_edge_index] <- Mk_llik + OU_llik
        # L_z[i,value_i,focal_edge_index] <- OU_llik
        # L_z[,,focal_edge_index]
      }
      C_z[i,which.max(L_z[i,,focal_edge_index]),focal_edge_index] <- 1
    }
  }
  # step 2 of pupko: calculations for the internal nodes, excluding the root
  pruningwise.index <- unique(phy$edge[reorder(phy, "pruningwise", index.only = TRUE), 1])
  for(node_index in pruningwise.index[-length(pruningwise.index)]){
    focal_anc_index <- which(phy$edge[,2] %in% node_index)
    focal_dec_index <- which(phy$edge[,1] %in% node_index)
    RootAsAnc <- phy$edge[focal_anc_index,1] == nTip+1
    # time is based on the parent to focal length. dec values calculated already
    time_edge <- phy$edge.length[focal_anc_index]
    # for(Mk_index in 1:nStates){
    #   state <- unique(AllCombos[,1])[Mk_index]
    #   Mk_vec <- rep(0, nStates)
    #   Mk_vec[state] <- 1
    #   edge_Mk[Mk_index,] <- log(c(expm(Q * time_edge) %*% Mk_vec))
    # }
    # P_mat <- log((expm(Q * time_edge)))
    # i indexes the parental state, j indexes the descendants 
    for(i in 1:nCombos){
      # Mk prereqs for the parent j
      state_i <- AllCombos[i,1]
      if(RootAsAnc){
        # init_i <- theta[state_i]
        init_i <- AllCombos[i,2]
      }else{
        init_i <- AllCombos[i,2]
      }
      # Mk_vec_i <- edge_Mk[state_i,] 
      # the probability that we end up in state j at node z with the parent as state i. 
      for(j in 1:nCombos){
        # given the parental node defines the distribution, we evaluate over the possible k ending values
        state_j <- AllCombos[j,1]
        # Mk_llik <- Mk_vec_i[state_j] # P_ij for the Mk process. the probability that we transition from i to j over time_edge is just the probability of being state_j
        if(state_i == state_j){
          Mk_llik <- log(1 - (1 - exp(time_edge * Q[state_j,state_j])))
        }else{
          Mk_llik <- log(Q[state_i,state_j]) + log(1 - (1 - exp(time_edge * Q[state_j,state_j])))
        } # P_ij for the Mk process. the probability that we transition at time 0 from i to j then remain in branch j for the entire branch length. the probability of a transtion of time 0 is dexp(0, rate), which is just rate.
        # next calculate the OU likelihood
        OU_llik <- dOU(AllCombos[j,2], init_i, time_edge, alpha[state_j], theta[state_j], sigma.sq[state_j], TRUE) # P_ij for the OU process. we observe the particular j bin, we start in the i bin, and the OU parameters are defined by the i bin's state
        # L_z[i,j,focal_anc_index] <- OU_llik + sum(apply(L_z[j,,focal_dec_index], 2, max))
        L_z[i,j,focal_anc_index] <- Mk_llik + OU_llik + sum(apply(L_z[j,,focal_dec_index], 2, max)) # L_z(i) is the maximum of row j. in pupko he doesn't hold all of the values only keeping the maximum. this would be achieved by collapsing the collumns of L_z and only keeping the max
      }
      # the C_z(i) is the maximum of the ith row
      C_z[i,which.max(L_z[i,,focal_anc_index]),focal_anc_index] <- 1
    }
  }
  # step 4 of pupko: evaluate the root
  focal_dec_index <- which(phy$edge[,1] %in% (nTip+1))
  C_root <- L_root <- vector("numeric", nCombos)
  for(k in 1:nCombos){
    state_k <- AllCombos[k,1]
    #P_k <- root.p[state_k]
    #L_root[k] <- log(P_k) + sum(apply(L_z[k,,focal_dec_index], 2, max))
    L_root[k] <- sum(apply(L_z[k,,focal_dec_index], 2, max))
  }
  if(return.lik == TRUE){
    return(max(L_root))
  }
  # for possible ways to reconstruct...
  if(sample.recons){
    AllRootsRecons <- vector("list", nSamples + 1)
    AllRootsRecons[[nSamples + 1]] <- L_root
    for(i in 1:nSamples){
      C_root <- vector("numeric", nCombos)
      PRoot_i <- getProbFromLliks(L_root)
      C_root[sample(1:length(PRoot_i), 1, prob = PRoot_i)] <- 1
      JointComboRecon <- vector("numeric", max(phy$edge))
      names(JointComboRecon) <- 1:max(phy$edge)
      JointComboRecon[nTip+1] <- which(C_root == 1)
      # step 5 of pupko: assign joint recon to nodes
      preorder.index <- reorder(phy, "pruningwise")$edge[,2]
      preorder.index <- preorder.index[length(preorder.index):1]
      for(node_index in preorder.index){
        focal_anc_index <- which(phy$edge[,2] %in% node_index)
        focal_parent <- phy$edge[focal_anc_index,1]
        # C_z_i <- C_z[JointComboRecon[focal_parent],,focal_anc_index] # denote by i the reconstructed ancestor at y
        L_z_i <- L_z[JointComboRecon[focal_parent],,focal_anc_index]
        P_z_i <- getProbFromLliks(L_z_i)
        JointComboRecon[node_index] <- sample(1:length(P_z_i), 1, prob = P_z_i)
      }
      AllRootsRecons[[i]] <- JointComboRecon
    }
    return(AllRootsRecons)
  }else{
    if(all.roots){
      AllRootsRecons <- vector("list", length(L_root) + 1)
      AllRootsRecons[[length(L_root) + 1]] <- L_root
      for(i in 1:length(L_root)){
        C_root <- vector("numeric", nCombos)
        C_root[i] <- 1
        JointComboRecon <- vector("numeric", max(phy$edge))
        names(JointComboRecon) <- 1:max(phy$edge)
        JointComboRecon[nTip+1] <- which(C_root == 1)
        # step 5 of pupko: assign joint recon to nodes
        preorder.index <- reorder(phy, "pruningwise")$edge[,2]
        preorder.index <- preorder.index[length(preorder.index):1]
        for(node_index in preorder.index){
          focal_anc_index <- which(phy$edge[,2] %in% node_index)
          focal_parent <- phy$edge[focal_anc_index,1]
          C_z_i <- C_z[JointComboRecon[focal_parent],,focal_anc_index] # denote by i the reconstructed ancestor at y
          JointComboRecon[node_index] <- which(C_z_i == 1)
        }
        AllRootsRecons[[i]] <- JointComboRecon
      }
      return(AllRootsRecons)
    }else{
      C_root[which.max(L_root)] <- 1
      JointComboRecon <- vector("numeric", max(phy$edge))
      names(JointComboRecon) <- 1:max(phy$edge)
      JointComboRecon[nTip+1] <- which(C_root == 1)
      
      # step 5 of pupko: assign joint recon to nodes
      preorder.index <- reorder(phy, "pruningwise")$edge[,2]
      preorder.index <- preorder.index[length(preorder.index):1]
      for(node_index in preorder.index){
        focal_anc_index <- which(phy$edge[,2] %in% node_index)
        focal_parent <- phy$edge[focal_anc_index,1]
        C_z_i <- C_z[JointComboRecon[focal_parent],,focal_anc_index] # denote by i the reconstructed ancestor at y
        JointComboRecon[node_index] <- which(C_z_i == 1)
      }
      llik <- max(L_root)
      return(list(LnLik = llik, Recon = JointComboRecon))
    }
  }
}

getProbFromLliks <- function(lliks){
  P <- exp(lliks - max(lliks))/sum(exp(lliks - max(lliks)))
  return(P)
}

# hOUwie's input data will be similar to OUwie, with the exception that there can be more than a single trait being evaluated thus it's defined as column 1 is species name, the last column is the continuous trait and anything in between are discrete characters (can also account for mserr if second last columns are continuous)
organizeHOUwieDat <- function(data, mserr){
  # return a list of corHMM data and OU data
  if(mserr=="known"){
    data.cor <- data[, 1:(dim(data)[2]-2)]
    data.cor <- corHMM:::corProcessData(data.cor)
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
    return(10000000)
  }else{
    logl <- -as.numeric(Ntip(phy) * log(2 * pi) + comp$logd + (comp$PP - 2 * comp$QP + comp$QQ))/2
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

# simulate a hOUwie model
hOUwie.sim <- function(phy, Q, root.freqs, alpha, sig2, theta0, theta, nMap=1){
  # simulate an Mk dataset
  dat.cor <- rTraitDisc(phy, Q, states = 1:dim(Q)[1], root.value = sample(1:dim(Q)[1], 1, prob = root.freqs))
  while(!all(levels(dat.cor) %in% dat.cor)){
    dat.cor <- rTraitDisc(phy, Q, states = 1:dim(Q)[1], root.value = sample(1:dim(Q)[1], 1, prob = root.freqs))
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
  cat("Fit\n")
  print(output)
  cat("\n")
  cat("Legend\n")
  print(x$legend)
  cat("\n")
  cat("\nRegime Rate matrix\n")
  print(x$solution.cor)
  cat("\n")
  cat("\nOU Estimates\n")
  print(x$solution.ou)
  cat("\n")
  if(!any(is.na(x$solution.ou[1,]))){
    cat("\nHalf-life (another way of reporting alpha)\n")
    print(log(2)/x$solution.ou[1,])
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
#   phy <- reorder(phy, "pruningwise")
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
NodeSpecific.hOuwieRecon <- function(pars, phy, rate.cat, data.cor, liks, Q, rate, root.p, data.ou, index.ou, algorithm, mserr, nSim, nCores, tip.paths, weighted, fix.node){
  est.pars<-log(pars)
  ReconVec <- liks[fix.node,]
  # are we reconstructing a tip?
  if(any(ReconVec == 1)){
    for(i in which(ReconVec == 1)){
      ReconVec[i] <- -hOUwie.dev(est.pars, phy=phy, rate.cat=rate.cat, data.cor=data.cor, liks=liks, Q=Q, rate=rate, root.p=root.p, data.ou=data.ou, index.ou=index.ou, algorithm=algorithm, mserr=mserr, nSim=nSim, nCores=nCores, tip.paths=tip.paths, weighted=weighted, order.test=FALSE, fix.node = fix.node, fix.state = i)
    }
  }else{
    for(i in 1:length(ReconVec)){
      ReconVec[i] <- -hOUwie.dev(est.pars, phy=phy, rate.cat=rate.cat, data.cor=data.cor, liks=liks, Q=Q, rate=rate, root.p=root.p, data.ou=data.ou, index.ou=index.ou, algorithm=algorithm, mserr=mserr, nSim=nSim, nCores=nCores, tip.paths=tip.paths, weighted=weighted, order.test=FALSE, fix.node = fix.node, fix.state = i)
    }
  }
  ReconVec[ReconVec == 0] <- NA
  ReconVec <- exp(ReconVec)/sum(exp(ReconVec), na.rm = TRUE)
  ReconVec[is.na(ReconVec)] <- 0
  return(ReconVec)
}

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
