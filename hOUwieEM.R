hOUwie.EM <- function(phy, data, discrete_model, continuous_model, rate.cat, shift.point=0.5, mserr = "none", dual = FALSE, collapse = TRUE, root.station = FALSE, get.root.theta = FALSE, lb.disc = NULL, ub.disc = NULL, lb.cont = NULL, ub.cont = NULL, opts = NULL, niter = 1000, ip = NULL, nodeproposals = 2, tipproposals = 1){
  
  # source("~/2020_hOUwie/hOUwieEM.R")
  # require(OUwie)
  # require(corHMM)
  # data(tworegime)
  # phy <- tree
  # data <- trait
  # discrete_model <- "ER"
  # continuous_model <- "OUM"
  # rate.cat <- 2
  # shift.point <- 0.5
  # mserr <- "none"
  # dual = FALSE
  # collapse <- TRUE
  # root.station <- FALSE
  # get.root.theta <- FALSE
  # lb.disc = NULL
  # ub.disc = NULL
  # lb.cont = NULL
  # ub.cont = NULL
  # opts = NULL
  # niter = 50
  # ip = NULL
  # nodeproposals = 3
  # tipproposals = 1

  # if the data has negative values, shift it right
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
  
  # organize the data
  hOUwie.dat <- organizeHOUwieDat(data, mserr, collapse)
  nStates <- as.numeric(max(hOUwie.dat$data.cor[,2]))
  tip.paths <- lapply(1:length(phy$tip.label), function(x) OUwie:::getPathToRoot(phy, x))
  
  # establish the models being evaluated 
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
  
  # establish the upper and lower bounds
  Tmax <- max(branching.times(phy))
  if(mserr == "none"){
    nDiscrete <- dim(data)[2] - 2
  }else{
    nDiscrete <- dim(data)[2] - 3
  }
  if(is.null(lb.cont)){
    # the lower limit of alpha is defined as a halflife of 10000% of the max tree height
    # the lower limit of sigma is defined 10 times less than alpha
    # the lower limit of optim is defined 10 times lower than the minimum observation
    lb.alpha = 1e-10
    lb.sigma = 1e-10
    lb.optim = min(data[, 1+nDiscrete+1])/10
    lb.cont=c(lb.alpha,lb.sigma,lb.optim)
  }
  if(is.null(ub.cont)){
    # the upper limit of alpha is defined as a halflife of 1% of the max tree height
    # the upper limit of sigma is defined 10 times more than alpha
    # the upper limit of optim is defined 10 times more than the maximum observation
    ub.alpha = log(2)/(0.01 * Tmax)
    ub.sigma = 10 * ub.alpha
    ub.optim = max(data[, 1+nDiscrete+1])*10 
    ub.cont=c(ub.alpha,ub.sigma,ub.optim)
  }
  if(is.null(lb.disc)){
    # the minimum dwell time is defined as 100 times the max tree height
    lb.disc = 1/(Tmax*1000)
  }
  if(is.null(ub.disc)){
    # the maximum dwell time is defined as 0.01 times the max tree height
    ub.disc = 1/(Tmax*0.001)
  }
  
  # default ML search options
  if(is.null(opts)){
    opts <- list("algorithm"="NLOPT_LN_SBPLX", "maxeval"="1000000", "ftol_rel"=.Machine$double.eps^0.1)
  }

  # the number of parameters for each process
  n_p_discr <- max(index.disc)
  n_p_alpha <- length(unique(na.omit(index.cont[1,])))
  n_p_sigma <- length(unique(na.omit(index.cont[2,])))
  n_p_theta <- length(unique(na.omit(index.cont[3,])))
  
  # combine parameter number with the upper and lower bound values
  upper.cont <- c(rep(ub.cont[1], n_p_alpha), rep(ub.cont[2], n_p_sigma), rep(ub.cont[3], n_p_theta))
  lower.cont <- c(rep(lb.cont[1], n_p_alpha), rep(lb.cont[2], n_p_sigma), rep(lb.cont[3], n_p_theta))
  upper.disc <- rep(ub.disc, n_p_discr)
  lower.disc <- rep(lb.disc, n_p_discr)
  
  # initial generation of a node reconstruction
  # get marginal recon of discrete states | discrete_model
  cat("Initializing with a corHMM fit...\n")
  corHMM_init <- quiet(corHMM(phy = phy, data = hOUwie.dat$data.cor, rate.cat = rate.cat, rate.mat = index.disc, lower.bound = lb.disc, upper.bound = ub.disc, get.tip.states = rate.cat > 1))
  marginal_nodes <- corHMM_init$states
  marginal_tips <- corHMM_init$tip.states
  # possible hidden states based on bins
  # all possible combinations of observed and "hidden" states
  if(rate.cat > 1){
    bin_index <- cut(hOUwie.dat$data.ou[,3], rate.cat, labels = FALSE)
    combos <- expand.grid(1:max(hOUwie.dat$data.cor[,2]), 1:rate.cat)
    disc_tips <- vector("numeric", dim(marginal_tips)[1])
    for(i in 1:dim(combos)[1]){
      disc_tips[hOUwie.dat$data.cor[,2] == combos[i,1] & bin_index == combos[i,2]] <- i
    }
  }
  stochasticmap <- NULL
  count <- 1
  while(is.null(stochasticmap)){
    disc_nodes <- apply(marginal_nodes, 1, function(x) sample(1:length(x), 1, replace = FALSE, x))
    disc_tips <- apply(marginal_tips, 1, function(x) sample(1:length(x), 1, replace = FALSE, x))
    # if(rate.cat == 1){
    #   disc_tips <- apply(marginal_tips, 1, function(x) sample(1:length(x), 1, replace = FALSE, x))
    # }
    stochasticmap <- getStochmapFromNode(phy, disc_tips, disc_nodes, shift.point, index.disc)
    print(count)
    count <- count + 1
  }
  
  # initial paramaters for discrete
  ip.disc <- rep(mean(corHMM_init$solution[!is.na(corHMM_init$solution)]), n_p_discr)
  
  # get OU param estimates | continuous model
  # OUwie_init <- quiet(OUwie(phy = stochasticmap, data = hOUwie.dat$data.ou, model = "OU1", algorithm = "three.point", simmap.tree = TRUE, root.station = root.station, get.root.theta = get.root.theta))
  
  # initial parameters for the continuous
  starts.alpha <- rep(log(2)/Tmax, n_p_alpha)
  starts.sigma <- rep(var(hOUwie.dat$data.ou[,3]), n_p_sigma)
  start.theta <- getIP.theta(hOUwie.dat$data.ou[,3], disc_tips, index.cont[3,])
  ip.cont <- c(starts.alpha, starts.sigma, start.theta)
  # means.by.regime <- with(hOUwie.dat$data.ou, tapply(hOUwie.dat$data.ou[,3], hOUwie.dat$data.ou[,2], mean))
  # if(n_p_theta == length(means.by.regime)){
  #   start.theta <- means.by.regime
  # }else{
  #   start.theta <- rep(mean(hOUwie.dat$data.ou[,3]), n_p_theta)
  # }

  # some vectors
  best_model <- vector("numeric", 3 + n_p_discr + n_p_alpha + n_p_sigma + n_p_theta)
  names(best_model) <- c("total_likelihood", "discretelikelihood", "continuouslikelihood")
  best_model[1] <- -Inf
  current_best_id <- 0
  all_nodes <- c()
  state_independ_cont <- all(apply(index.cont, 1, function(x) length(unique(x)) == 1))
  
  # may change from sequence to a tolerance of total likelihood
  for(i in sequence(niter)){
    # add the current reconstruction to the vector of nodes already examined
    all_nodes[i] <- paste0(c(disc_nodes,disc_tips), collapse = "")

    # discrete_result <- optimize(discrete model given discrete_internal, start at best_model$discrete_params
    # discrete_result <- nloptr(x0=log(ip.disc), eval_f=dev.discretelikelihood, lb=log(lower.disc), ub=log(upper.disc), opts=opts, maps = stochasticmap$maps, index.disc = index.disc)
    discrete_result <- lapply(stochasticmap, function(x) nloptr(x0=log(ip.disc), eval_f=dev.discretelikelihood, lb=log(lower.disc), ub=log(upper.disc), opts=opts, maps = x$maps, index.disc = index.disc))
    
    
    # continuous_result <- optimize(continuous model given discrete_internal, start at best_model$continuous_params)
    if(state_independ_cont){
      # the continuous result won't change based on discrete regime paintings if the model is OU1 or BM1
      if(i == 1){
        # continuous_result <- nloptr(x0=log(ip.cont), eval_f=dev.continuouslikelihood, lb=log(lower.cont), ub=log(upper.cont), opts=opts, stochasticmap = stochasticmap, data.ou = hOUwie.dat$data.ou, index.cont = index.cont, tip.paths = tip.paths, mserr = mserr)
        continuous_result <- lapply(stochasticmap, function(x) nloptr(x0=log(ip.cont), eval_f=dev.continuouslikelihood, lb=log(lower.cont), ub=log(upper.cont), opts=opts, stochasticmap = x, data.ou = hOUwie.dat$data.ou, index.cont = index.cont, tip.paths = tip.paths, mserr = mserr))
      }
    }else{
      # continuous_result <- nloptr(x0=log(ip.cont), eval_f=dev.continuouslikelihood, lb=log(lower.cont), ub=log(upper.cont), opts=opts, stochasticmap = stochasticmap, data.ou = hOUwie.dat$data.ou, index.cont = index.cont, tip.paths = tip.paths, mserr = mserr)
      continuous_result <- lapply(stochasticmap, function(x) nloptr(x0=log(ip.cont), eval_f=dev.continuouslikelihood, lb=log(lower.cont), ub=log(upper.cont), opts=opts, stochasticmap = x, data.ou = hOUwie.dat$data.ou, index.cont = index.cont, tip.paths = tip.paths, mserr = mserr))
    }
    # OUwie(phy = stochasticmap, data = hOUwie.dat$data.ou, model = "OUM", algorithm = "three.point", simmap.tree = TRUE, root.station = root.station, get.root.theta = get.root.theta)
    # OUwie.fixed(stochasticmap, data = hOUwie.dat$data.ou, model = continuous_model, algorithm = "three.point", simmap.tree = TRUE, root.station = root.station, get.root.theta = get.root.theta, alpha = exp(continuous_result$solution)[1], sigma.sq = exp(continuous_result$solution)[2:3], theta = exp(continuous_result$solution)[4:5])

    # total_likelihood <- discrete_result$likelihood + continuous_result$likelihood
    # total_likelihood <- (-discrete_result$objective) + (-continuous_result$objective)
    total_likelihood <- -unlist(lapply(discrete_result, function(x) x$objective)) -unlist(lapply(continuous_result, function(x) x$objective))
    best_map <- which.max(total_likelihood)
    
    if(best_model[1] < total_likelihood[best_map]){
      # track the current best model
      best_model <- c(total_likelihood[best_map], -discrete_result[[best_map]]$objective, -continuous_result[[best_map]]$objective, exp(discrete_result[[best_map]]$solution), exp(continuous_result[[best_map]]$solution))
      # the current ip come from the best model 
      ip.disc <- exp(discrete_result[[best_map]]$solution)
      ip.cont <- exp(continuous_result[[best_map]]$solution)
      # the current best discrete sample
      best_disc.nodes <- disc_nodes
      best_disc.tips <- disc_tips
    }
    
    print(c(LnLik = best_model[1], LnLikDisc = best_model[2], LnLikCont = best_model[3]))
    # a new sample is generated and we start all over
    nsamples.nodes = round(runif(1) * nodeproposals) + 1
    nsamples.tips = round(runif(1) * tipproposals)
    type = c("unif", "entropy")[round(runif(1) + 1)]
    # the new sample is checked that if has no duplicates
    stochasticmap <- NULL
    while(is.null(stochasticmap)){
      duplicate.node.states <- TRUE
      count <- 1
      while(duplicate.node.states){
        disc_nodes <- sample.discretestates(best_disc.nodes, marginal_nodes, nsamples.nodes, type)
        if(rate.cat > 1){
          disc_tips <- sample.discretestates(best_disc.tips, marginal_tips, nsamples.tips, type)
        }
        duplicate.node.states <- paste0(c(disc_nodes, disc_tips), collapse="") %in% all_nodes
        count <- count + 1
        if(count == 20){
          type = "unif"
        }
        if(count == niter){
          return(list(best_model, best_disc.nodes, best_disc.tips))
        }
      }
      stochasticmap <- getStochmapFromNode(phy, disc_tips, disc_nodes, shift.point, index.disc)
    }
  }
  return(list(best_model, best_disc.nodes, best_disc.tips))
}

# given a current reconstruction sample a new set of states that is unique
sample.discretestates <- function(states, probability, nsamples, type="unif"){
  if(type == "entropy"){
    probability[probability == 0] <- 1e-10
    entropy <- apply(probability, 1, function(x) -sum(log(x) * x))
    p_sample <- entropy/sum(entropy)
  }
  if(type == "unif"){
    p_sample <- rep(1/length(states), length(states))
  }
  nodestochange <- sample(sequence(length(states)), size = nsamples, prob = p_sample)
  p_nodestochange <- matrix(0, nsamples, dim(probability)[2])
  p_nodestochange[] <- probability[nodestochange,]
  state_nodestochanges <- states[nodestochange]
  for(i in sequence(nsamples)){
    p_nodestochange[i,state_nodestochanges[i]] <- 0
    state_nodestochanges[i] <- sample(sequence(length(p_nodestochange[i,])), 1, prob = p_nodestochange[i,])
  }
  states[nodestochange] <- state_nodestochanges
  return(states)
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

# a wrapper function for the continuous likelihood
dev.continuouslikelihood <- function(pars, stochasticmap, data.ou, index.cont, tip.paths, mserr){
  p <- exp(pars)
  index.cont[is.na(index.cont)] <- max(index.cont, na.rm = TRUE) + 1
  Rate.mat <- index.cont
  Rate.mat[] <- c(p, 1e-10)[index.cont]
  alpha = Rate.mat[1,]
  sigma.sq = Rate.mat[2,]
  theta = Rate.mat[3,]
  continuouslikelihood <- OUwie.basic(stochasticmap, data.ou, simmap.tree=TRUE, scaleHeight=FALSE, alpha=alpha, sigma.sq=sigma.sq, theta=theta, algorithm="three.point", tip.paths=tip.paths, mserr=mserr)
  return(-continuouslikelihood)
}

# a wrapper function for the discrete likelihood
dev.discretelikelihood <- function(pars, maps, index.disc){
  p <- exp(pars)
  Q <- index.disc
  for(i in 1:max(index.disc)){
    Q[Q == i] <- p[i]
  }
  diag(Q) <- -rowSums(Q)
  discretelikelihood <- getMapProbability(maps, Q)
  return(-discretelikelihood)
}

# what to do about the cases where we've sampled an impossible ancestral state and descendant state in a single step?
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
    P <- 1e-100
  }
  return(P)
}

probPathOld <- function(path, Q){
  state_i <- as.numeric(names(path)[1])
  state_j <- as.numeric(names(path)[length(path)])
  if(state_i != state_j){
    time_i <- as.numeric(path[1])
    time_j <- as.numeric(path[length(path)])
    rate_i <- abs(Q[state_i,state_j])
    rate_j <- abs(Q[state_j,state_j])
    P <- dexp(time_i, rate_i) * (1 - pexp(rate_j, time_j))
  }else{
    time_i <- as.numeric(path[1])
    rate_i <- abs(Q[state_i,state_j])
    P <- (1 - pexp(rate_i, time_i))
  }
  return(P)
}

# generates a discrete model index if not supplied by the user
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

# organize data to be consistent with OUwie and corHMM
organizeHOUwieDat <- function(data, mserr, collapse){
  # return a list of corHMM data and OU data
  if(mserr=="known"){
    data.cor <- data[, 1:(dim(data)[2]-2)]
    data.cor <- corHMM:::corProcessData(data.cor, collapse)
    data.ou <- data.frame(sp = data[,1], 
                          reg = data.cor$corData[,2], 
                          x = data[, dim(data)[2]-1],
                          err = data[, dim(data)[2]])
  }
  if(mserr=="none"){
    data.cor <- data[, 1:(dim(data)[2]-1)]
    data.cor <- corHMM:::corProcessData(data.cor, collapse)
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

# check the order of the parameters, useful for hidden states
check.order <- function(mat, decreasing = FALSE){
  order.vec <- vector("logical", dim(mat)[1])
  for(i in 1:dim(mat)[1]){
    order.vec[i] <- identical(mat[i,], sort(mat[i,], decreasing = decreasing))
  }
  return(all(order.vec))
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

#shhhhhh
quiet <- function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}

## takes a node based reconstruction and returns a map (identical to a map from simmap)
getStochmapFromNode <- function(phy, tipstates, nodestates, shift.point, index.disc, mid.trans=TRUE, nSamples=10){
  Map <- vector("list", dim(phy$edge)[1])
  Data <- c(tipstates, nodestates)
  NodeStates <- cbind(Data[phy$edge[,1]], Data[phy$edge[,2]])
  mid.trans.index <- c()
  for(i in 1:dim(phy$edge)[1]){
    from <- as.numeric(NodeStates[i,1])
    to <- as.numeric(NodeStates[i,2])
    if(from == to){
      tmp <- phy$edge.length[i]
      names(tmp) <- from
      Map[[i]] <- tmp
    }else{
      if(index.disc[from,to] == 0){
        if(!mid.trans){
          # if we only want one transition max per branch return NULL if not
          return(NULL)
        }
        # if the transition is disallowed we identify where we can go in two steps
        # if a two step transition is allowed, return a list of all possible two steps
        mid <- which(index.disc[from,] != 0)
        ManyMap <- vector("list", length(mid))
        mid.trans.index <- c(mid.trans.index, i)
        for(j in 1:length(mid)){
          tmp <- rep(phy$edge.length[i]/3, 3)
          names(tmp) <- c(from, mid[j], to)
          ManyMap[[j]] <- tmp
        }
        Map[[i]] <- ManyMap
      }else{
        shift.time <- shift.point * phy$edge.length[i]
        tmp <- c(phy$edge.length[i] - shift.time, shift.time)
        names(tmp) <- c(from, to)
        Map[[i]] <- tmp
      }
    }
  }
  # because we can have multiple mappings we will take each of those and produce a map
  if(length(mid.trans.index) > 0){
    # apply(expand.grid(mid.trans.index, c("1", "2")), 1, function(x) paste(x, collapse = "_"))
    possible.combos <- lapply(mid.trans.index, function(x) paste(x, c("1", "2"), sep = "_"))
    full.set.combos <- t(sapply(1:10, function(x) unlist(lapply(possible.combos, function(x) sample(x, 1)))))
    # full.set.combos <- expand.grid(possible.combos)
    # we create as many maps as potential combinations then we make a map specific to each combo
    Maps <- vector("list", dim(full.set.combos)[1])
    for(i in 1:dim(full.set.combos)[1]){
      Maps[[i]] <- Map
      combo_i <- full.set.combos[i,]
      for(j in 1:length(combo_i)){
        # in essence what we're doing: Maps[[k]][[2]] <- Maps[[k]][[2]][[1]]
        # i is the particular mapping, j is the edge of the map we're changing to combo_i
        edge_i_combo_j <- as.numeric(strsplit(as.vector(combo_i[j]), "_")[[1]])
        Maps[[i]][[edge_i_combo_j[1]]] <- Maps[[i]][[edge_i_combo_j[1]]][[edge_i_combo_j[2]]]
      }
    }
  }else{
    Maps <- list(Map)
  }
  Phys <- vector("list", length(Maps))
  for(i in 1:length(Maps)){
    Phys[[i]] <- phy
    mapped.edge <- corHMM:::convertSubHistoryToEdge(phy, Maps[[i]])
    Phys[[i]]$maps <- Maps[[i]]
    Phys[[i]]$mapped.edge <- mapped.edge
    attr(Phys[[i]], "map.order") <- "right-to-left"
    if (!inherits(Phys[[i]], "simmap")) 
      class(Phys[[i]]) <- c("simmap", setdiff(class(Phys[[i]]), "simmap"))
    
  }
  # mapped.edge <- corHMM:::convertSubHistoryToEdge(phy, Map)
  # phy$maps <- Map
  # phy$mapped.edge <- mapped.edge
  # attr(phy, "map.order") <- "right-to-left"
  # if (!inherits(phy, "simmap")) 
  #   class(phy) <- c("simmap", setdiff(class(phy), "simmap"))
  return(Phys)
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

# siple simulation function
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
