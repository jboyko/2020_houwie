# set of functions for the hidden rates OU model

# exported function with all the bells and whistles
#'@author James Boyko
#'@param phy a phylogeny of class phylo
#'@param data data.frame where [,1] = species name,  [,2:j] are discrete traits, [,j:m] are the continuous trait (and possible mserr). if mserr = "none" j=m-1, if mserr = "known j=m-2. m is the number of total columns in the dataset.
#'@param rate.cat the number rate categories for the discrete model (rate.cat > 1 means hidden markov model)
#'@param discrete_model designates a default discrete models (ARD, ER, or SYM) or a custom index matrix
#'@param continuous_model designates a default continuous model (BM1, BMS, OU1, OUM, OUMV, OUMA, OUMVA) or a custom index matrix
#'@param root.p designates the root prior (yang, maddfitz, flat, or numeric vector)
#'@param collapse a boolean indicating whether multiple traits should be collapsed to exlude non-observed trait combinations (e.g. if 00, 11, and 10 are all present, if collapse = FALSE the possible combination of 01 would be part of the model)
#'@param dual a boolean indicating whether transitions which require two changes are allowed (e.g. by default 00 to 11 is not allowed, but can be if dual = TRUE)
#'@param get.root.theta a boolean indicating whether the starting state, theta_0, should be estimated
#'#'@param root.station a booleean indicating whether to assume a random starting point (TRUE) or a fixed starting point (FALSE) 
#'@param mserr designates whether mserr is included (known) or not (none)
#'@param lb_discrete_model lower bound for trans rate. default is 1e-5
#'@param ub_discrete_model upper bound for trans rate. default is 21
#'#'@param lb_continuous_model designates the lower bounds for alpha, sigma.sq, optima. default is given in code
#'@param ub_continuous_model designates the upper bounds for alpha, sigma.sq, optima. default is given in code
#'@param optimizer designates the optimizer. one of "nlopt" for nloptr, "sann" for simulated annealing
#'@param opts options for nloptr
#'@param ip "fast" for analytic, "good" for optimized values, or a vector of initial params 
#'@param p a vector of params to calculate a fixed likelihood
#'@param nSim the number of simmaps a single set of params is evaluated over
#'@param quiet a boolean indicating whether to output user messages
#'@param recon a boolean indicating whether marginal node state reconstruction should be preformed
#'@param nodes which nodes should be reconstructed ("all", "internal", "external", or a numeric vector)
#'@param sample_tips a boolean which indicates whether the continuous values should inform the discrete tip values when initializing the stochastic mapping (default is TRUE)
#'@param n_starts a numeric value specifying how many optimizations are to occur
#'@param ncores a nunmeric value indicating how many cores are to be used, only useful if n_starts is greater than 2.
hOUwie <- function(phy, data, rate.cat, discrete_model, continuous_model, time_slice=NULL, nSim=1000, root.p="yang", dual = FALSE, collapse = TRUE, root.station=FALSE, get.root.theta=FALSE, mserr = "none", lb_discrete_model=NULL, ub_discrete_model=NULL, lb_continuous_model=NULL, ub_continuous_model=NULL, recon=FALSE, nodes="internal", p=NULL, ip=NULL, optimizer="nlopt", opts=NULL, quiet=FALSE, sample_tips=TRUE, n_starts = 1, ncores = 1){
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
  
  if(ncores > n_starts){
    cat("You have specified more cores are to be used than the number of starts. Setting ncores to be equal to the number of optimizations.\n")
    ncores <- n_starts
  }

  # organize the data
  phy <- reorder.phylo(phy, "pruningwise")
  hOUwie.dat <- organizeHOUwieDat(data, mserr, collapse)
  nStates <- as.numeric(max(hOUwie.dat$data.cor[,2]))
  nCol <- dim(data)[2] - ifelse(mserr == "none", 2, 3)
  Tmax <- max(branching.times(phy))
  tip.paths <- lapply(1:length(phy$tip.label), function(x) OUwie:::getPathToRoot(phy, x))

  if(class(discrete_model)[1] == "character"){
    index.disc <- getDiscreteModel(hOUwie.dat$data.cor, discrete_model, rate.cat, dual, collapse)
    index.disc[index.disc == 0] <- NA
  }else{
    index.disc <- discrete_model
    index.disc[index.disc == 0] <- NA
  }
  if(class(continuous_model)[1] == "character"){
    index.cont <- getOUParamStructure(continuous_model, "three.point", root.station, get.root.theta, nStates * rate.cat)
  }else{
    continuous_model[continuous_model == 0] <- NA
    index.cont <- continuous_model
  }
  if(dim(index.disc)[2] > dim(index.cont)[2]){
    stop("Not all of your discrete states have OU parameters associated with them. Please check that your discrete index matrix matches your continuous index matrix.")
  }
  if(dim(index.cont)[2] > dim(index.disc)[2]){
    stop("You have specified more OU parameters than there are states in the discrete process. Please check that your discrete index matrix matches your continuous index matrix.")
  }
  if(class(root.p[1]) != "character"){
    if(dim(index.disc)[2] != length(root.p)){
      stop("You have entered a custom root prior whose length does not equal the number of states in your discrete model.")
    }
  }
  
  if(is.null(lb_continuous_model)){
    # the lower limit of alpha is defined as a halflife of 10000% of the max tree height
    # the lower limit of sigma is defined 10 times less than alpha
    # the lower limit of optim is defined 10 times lower than the minimum observation
    lb.alpha = 1e-10
    lb.sigma = lb.alpha/10
    lb.optim = min(data[, 1+nCol+1])/10 
    lb_continuous_model=c(lb.alpha,lb.sigma,lb.optim)
  }
  if(is.null(ub_continuous_model)){
    # the upper limit of alpha is defined as a halflife of 1% of the max tree height
    # the upper limit of sigma is defined 10 times more than alpha
    # the upper limit of optim is defined 10 times more than the maximum observation
    ub.alpha = log(2)/(0.01 * Tmax)
    ub.sigma = 10 * ub.alpha
    ub.optim = max(data[, 1+nCol+1])*10 
    ub_continuous_model=c(ub.alpha,ub.sigma,ub.optim)
  }
  if(is.null(lb_discrete_model)){
    # the minimum dwell time is defined as 100 times the max tree height
    lb_discrete_model = 1/(Tmax*100)
  }
  if(is.null(ub_discrete_model)){
    ub_discrete_model = 1/(Tmax*0.01)
  }
  #Ensures that weird root state probabilities that do not sum to 1 are input:
  if(!is.null(root.p)){
    if(!is.character(root.p)){
      root.p <- root.p/sum(root.p)
    }
  }
  if(is.null(time_slice)){
    time_slice <- Tmax/10
  }
  
  # the number of parameters for each process
  n_p_trans <- max(index.disc, na.rm = TRUE)
  n_p_alpha <- length(unique(na.omit(index.cont[1,])))
  n_p_sigma <- length(unique(na.omit(index.cont[2,])))
  n_p_theta <- length(unique(na.omit(index.cont[3,])))
  
  # an internal data structure (internodes liks matrix) for the dev function
  edge_liks_list <- getEdgeLiks(phy, hOUwie.dat$data.cor, nStates, rate.cat, time_slice)
  
  # default MLE search options
  if(is.null(opts)){
    if(optimizer == "nlopt" | optimizer == "twostep"){
      opts <- list("algorithm"="NLOPT_LN_SBPLX", "maxeval"="1000000", "ftol_rel"=.Machine$double.eps^0.5)
    }
    if(optimizer == "sann"){
      opts <- list(max.call=10000, smooth=FALSE)
    }
  }
  # p is organized into 2 groups with the first set being corHMM and the second set being OUwie
  # organized as c(trans.rt, alpha, sigma.sq, theta)
  # evaluate likelihood
  if(!is.null(p)){
    if(!quiet){
      cat("Calculating likelihood from a set of fixed parameters.\n")
      print(p)
    }
    if(max(index.cont, na.rm = TRUE) + max(index.disc, na.rm = TRUE) != length(p)){
      message <- paste0("The number of parameters does not match the number required by the model structure. You have supplied ", length(p), ", but the model structure requires ", max(index.cont, na.rm = TRUE) + max(index.disc, na.rm = TRUE), ".")
      stop(message, call. = FALSE)
    }
    out<-NULL
    pars <- out$solution <- log(p)
    out$objective <- hOUwie.dev(p = log(p), phy=phy, data=hOUwie.dat$data.ou, 
                                rate.cat=rate.cat, mserr=mserr, 
                                index.disc=index.disc, index.cont=index.cont, root.p=root.p,
                                edge_liks_list=edge_liks_list, nSim=nSim, tip.paths=tip.paths, 
                                sample_tips=sample_tips, split.liks=FALSE)
  }else{
    out<-NULL
    lower = log(c(rep(lb_discrete_model, n_p_trans), 
                  rep(lb_continuous_model[1], length(unique(na.omit(index.cont[1,])))), 
                  rep(lb_continuous_model[2], length(unique(na.omit(index.cont[2,])))), 
                  rep(lb_continuous_model[3], length(unique(na.omit(index.cont[3,]))))))
    upper = log(c(rep(ub_discrete_model, n_p_trans), 
                  rep(ub_continuous_model[1], length(unique(na.omit(index.cont[1,])))), 
                  rep(ub_continuous_model[2], length(unique(na.omit(index.cont[2,])))), 
                  rep(ub_continuous_model[3], length(unique(na.omit(index.cont[3,]))))))
    # cat(c("TotalLnLik", "DiscLnLik", "ContLnLik"), "\n")
    # check for user input initial parameters 
    if(is.null(ip)){
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
      starts.alpha <- rep(log(2)/Tmax, n_p_alpha)
      starts.sigma <- rep(var(hOUwie.dat$data.ou[,3]), n_p_sigma)
      start.theta <- getIP.theta(hOUwie.dat$data.ou[,3], disc_tips, index.cont[3,])
      start.cor <- rep(10/sum(phy$edge.length), n_p_trans)
      starts.basic = c(start.cor, starts.alpha, starts.sigma, start.theta)
      if(!quiet){
        cat("Generating an initial parameter estimation point based on independent models.\n This may take some time, but is generally worth it.\n")
      }
      spec_opts <- list("algorithm"="NLOPT_LN_SBPLX", "maxeval"="1000", "ftol_rel"=.Machine$double.eps^0.5)
      discrete_init <- silence(corHMM(phy=phy, data=hOUwie.dat$data.cor, rate.cat=rate.cat, rate.mat=index.disc, node.states="joint", opts = spec_opts))
      map <- OUwie:::getMapFromNode(phy, hOUwie.dat$data.cor[,2], discrete_init$states, 0.5)
      phy_ouwie <- getMapFromSubstHistory(list(map), phy)
      tmp_continuous <- nloptr(x0 = log(c(starts.alpha, starts.sigma, start.theta)), eval_f = OUwie.basic.dev, lb=lower[-seq(n_p_trans)], ub=upper[-seq(n_p_trans)], opts=spec_opts, phy = phy_ouwie[[1]], data = hOUwie.dat$data.ou, mserr = mserr, index.cont = index.cont, tip.paths = tip.paths)
      # continuous_init <- OUwie(phy=phy_ouwie[[1]], data=hOUwie.dat$data.ou, model="OUM", simmap.tree=TRUE, algorithm="three.point")
      starts.ou <- exp(tmp_continuous$solution)
      start.cor <- sapply(seq(max(discrete_init$index.mat, na.rm = TRUE)), function(x) na.omit(discrete_init$solution[discrete_init$index.mat == x])[1])
      starts <- c(start.cor, starts.ou)
      starts[starts < exp(lower)] <- starts.basic[starts < exp(lower)]
      starts[starts > exp(upper)] <- starts.basic[starts > exp(upper)]
    }else{
      starts <- ip
    }
    if(!quiet){
      cat("Starting a thorough search of parameters with", nSim, "simmaps using the", optimizer, "optimization protocol...\n")
    }
    multiple_starts <- generateMultiStarting(starts, index.disc, index.cont, n_starts, exp(lower), exp(upper))
    if(optimizer == "nlopt"){
      # out = nloptr(x0=log(starts), eval_f=hOUwie.dev, lb=lower, ub=upper, opts=opts,
      #              phy=phy, data=hOUwie.dat$data.ou,
      #              rate.cat=rate.cat, mserr=mserr,
      #              index.disc=index.disc, index.cont=index.cont, root.p=root.p,
      #              edge_liks_list=edge_liks_list, nSim=nSim, tip.paths=tip.paths,
      #              sample_tips=sample_tips, split.liks=FALSE)
      multi_out <- mclapply(multiple_starts, function(x) nloptr(x0=log(x), eval_f=hOUwie.dev, lb=lower, ub=upper, opts=opts, phy=phy, data=hOUwie.dat$data.ou, rate.cat=rate.cat, mserr=mserr,index.disc=index.disc, index.cont=index.cont, root.p=root.p,edge_liks_list=edge_liks_list, nSim=nSim, tip.paths=tip.paths, sample_tips=sample_tips, split.liks=FALSE), mc.cores = ncores)
      multi_logliks <- unlist(lapply(multi_out, function(x) x$objective))
      search_summary <- c(best_loglik = -min(multi_logliks), mean_loglik = -mean(multi_logliks), sd_logliks = sd(multi_logliks))
      if(!quiet){
        cat("\nOptimization complete. Optimization summary:\n")
        print(search_summary)
      }
      out <- multi_out[[which.min(multi_logliks)]]
      pars <- out$solution
    }
    if(optimizer == "sann"){
      # out = GenSA(par=log(starts), fn=hOUwie.dev, lower=lower, upper=upper, control=opts, 
      #              phy=phy, data=hOUwie.dat$data.ou, 
      #              rate.cat=rate.cat, mserr=mserr, 
      #              index.disc=index.disc, index.cont=index.cont, root.p=root.p,
      #              edge_liks_list=edge_liks_list, nSim=nSim, tip.paths=tip.paths, 
      #              sample_tips=sample_tips, split.liks=FALSE)
      multi_out <- mclapply(multiple_starts, function(x) GenSA(par=log(x), fn=hOUwie.dev, lower=lower, upper=upper, control=opts, phy=phy, data=hOUwie.dat$data.ou, rate.cat=rate.cat, mserr=mserr, index.disc=index.disc, index.cont=index.cont, root.p=root.p, edge_liks_list=edge_liks_list, nSim=nSim, tip.paths=tip.paths, sample_tips=sample_tips, split.liks=FALSE), mc.cores = ncores)
      multi_logliks <- unlist(lapply(multi_out, function(x) x$value))
      search_summary <- c(best_loglik = -min(multi_logliks), mean_loglik = -mean(multi_logliks), sd_logliks = sd(multi_logliks))
      if(!quiet){
        cat("Optimization complete. Optimization summary:")
        print(search_summary)
      }
      out <- multi_out[[which.min(multi_logliks)]]
      pars <- out$par
    }
  }
  # preparing output
  liks_houwie <- hOUwie.dev(p = pars, phy=phy, data=hOUwie.dat$data.ou, rate.cat=rate.cat, mserr=mserr, index.disc=index.disc, index.cont=index.cont, root.p=root.p, edge_liks_list=edge_liks_list, nSim=nSim, tip.paths=tip.paths, sample_tips=sample_tips, split.liks=TRUE)
  houwie_obj <- getHouwieObj(liks_houwie, pars=exp(pars), phy=phy, data=data, hOUwie.dat=hOUwie.dat, rate.cat=rate.cat, mserr=mserr, index.disc=index.disc, index.cont=index.cont, root.p=root.p, nSim=nSim, sample_tips=sample_tips, nStates=nStates, discrete_model=discrete_model, continuous_model=continuous_model, time_slice=time_slice, root.station=root.station, get.root.theta=get.root.theta,lb_discrete_model,ub_discrete_model,lb_continuous_model,ub_continuous_model, ip=ip, opts=opts, quiet=quiet)
  # adding independent model if included
  if(is.null(p)){
    liks_indep <- hOUwie.dev(p = log(starts), phy=phy, data=hOUwie.dat$data.ou, rate.cat=rate.cat, mserr=mserr, index.disc=index.disc, index.cont=index.cont, root.p=root.p, edge_liks_list=edge_liks_list, nSim=nSim, tip.paths=tip.paths, sample_tips=sample_tips, split.liks=TRUE)
    houwie_obj$init_model <- getHouwieObj(liks_indep, pars=starts, phy=phy, data=data, hOUwie.dat=hOUwie.dat, rate.cat=rate.cat, mserr=mserr, index.disc=index.disc, index.cont=index.cont, root.p=root.p, nSim=nSim, sample_tips=sample_tips, nStates=nStates, discrete_model=discrete_model, continuous_model=continuous_model, time_slice=time_slice, root.station=root.station, get.root.theta=get.root.theta,lb_discrete_model,ub_discrete_model,lb_continuous_model,ub_continuous_model, ip=ip, opts=opts, quiet=quiet)
  }
  # conducting ancestra state resconstruction
  if(recon){
    houwie_recon <- hOUwieRecon(houwie_obj, nodes)
    houwie_obj$recon <- houwie_recon
  }
  return(houwie_obj)
}

hOUwieRecon <- function(houwie_obj, nodes="all"){
  # if the class is houwie_obj
  phy <- houwie_obj$phy
  hOUwie.dat <- houwie_obj$hOUwie.dat
  root.p <- houwie_obj$root.p
  mserr <- houwie_obj$mserr
  rate.cat <- houwie_obj$rate.cat
  index.cont <- houwie_obj$index.cont
  index.disc <- houwie_obj$index.disc
  p <- houwie_obj$p
  time_slice <- houwie_obj$time_slice
  state_names <- colnames(houwie_obj$solution.disc)
  sample_tips <- houwie_obj$sample_tips
  nSim <- houwie_obj$nSim
  # organize the data
  phy <- reorder.phylo(phy, "pruningwise")
  nTip <- length(phy$tip.label)
  Tmax <- max(branching.times(phy))
  nStates <- as.numeric(max(hOUwie.dat$data.cor[,2]))
  tip.paths <- lapply(1:length(phy$tip.label), function(x) OUwie:::getPathToRoot(phy, x))
  # an internal data structure (internodes liks matrix) for the dev function
  edge_liks_list <- getEdgeLiks(phy, hOUwie.dat$data.cor, nStates, rate.cat, time_slice)
  if(is.character(nodes[1])){
    if(nodes == "internal"){
      nodes_to_fix <- min(phy$edge[,1]):max(phy$edge[,1])
    }
    if(nodes == "external"){
      nodes_to_fix <- 1:nTip
    }
    if(nodes  == "all"){
      nodes_to_fix <- 1:max(phy$edge[,1])
    }
  }else{
    nodes_to_fix <- nodes
  }
  recon_matrix <- matrix(-Inf, length(nodes_to_fix), nStates * rate.cat, dimnames = list(nodes_to_fix, state_names))
  for(i in 1:length(nodes_to_fix)){
    node_i <- nodes_to_fix[i]
    anc_edges_to_fix <- which(phy$edge[,1] == node_i)
    dec_edges_to_fix <- which(phy$edge[,2] == node_i)
    if(node_i <= nTip){
      possible_states <- which(edge_liks_list[[dec_edges_to_fix]][1,] == 1)
    }else{
      possible_states <- 1:(nStates*rate.cat)
    }
    check_unreasonable_recon <- TRUE
    while(check_unreasonable_recon){
      for(state_j in 1:(nStates*rate.cat)){
        if(!state_j %in% possible_states){
          next
        }
        edge_liks_list_i <- edge_liks_list
        fix_vector <- numeric(nStates * rate.cat)
        fix_vector[state_j] <- 1
        for(k in dec_edges_to_fix){
          edge_liks_list_i[[k]][1,] <- fix_vector
        }
        for(k in anc_edges_to_fix){
          last_row <- dim(edge_liks_list_i[[k]])[1]
          edge_liks_list_i[[k]][last_row,] <- fix_vector
        }
        fixed_loglik <- -hOUwie.dev(p = log(p), phy=phy, data=hOUwie.dat$data.ou, rate.cat=rate.cat, mserr=mserr, index.disc=index.disc, index.cont=index.cont, root.p=root.p, edge_liks_list=edge_liks_list_i, nSim=nSim, tip.paths=tip.paths, sample_tips=sample_tips, split.liks=FALSE)
        recon_matrix[i, state_j] <- fixed_loglik
      }
      recon_loglik <- max(recon_matrix[i, ]) + log(sum(exp(recon_matrix[i, ] - max(recon_matrix[i, ]))))
      check_unreasonable_recon <- recon_loglik < houwie_obj$loglik * 1.05
    }
  }
  recon_matrix <- t(apply(recon_matrix, 1, function(x) x - max(x)))
  recon_matrix <- round(exp(recon_matrix)/rowSums(exp(recon_matrix)), 10)
  return(recon_matrix)
}

hOUwie.dev <- function(p, phy, data, rate.cat, mserr,
                       index.disc, index.cont, root.p, 
                       edge_liks_list, nSim, tip.paths=NULL, 
                       sample_tips=TRUE, split.liks=FALSE){
  p <- exp(p)
  k <- max(index.disc, na.rm = TRUE)
  p.mk <- p[1:k]
  p.ou <- p[(k+1):length(p)] 
  Rate.mat <- matrix(1, 3, dim(index.disc)[2])
  alpha.na <- is.na(index.cont[1,])
  index.cont[is.na(index.cont)] <- max(index.cont, na.rm = TRUE) + 1
  Rate.mat[] <- c(p.ou, 1e-10)[index.cont]
  alpha = Rate.mat[1,]
  sigma.sq = Rate.mat[2,]
  theta = Rate.mat[3,]
  rate <- index.disc
  rate[is.na(rate)] <- k + 1
  Q <- matrix(0, dim(rate)[1], dim(rate)[2])
  Q[] <- c(p.mk, 0)[rate]
  diag(Q) <- -rowSums(Q)
  if(rate.cat > 1){
    if(sample_tips){
      means <- theta
      vars <- sigma.sq
      normal.params <- rbind(means, vars)
      sample.tip.probs <- apply(normal.params, 2, function(x) dnorm(data[,3], x[1], sqrt(x[2])))
      for(i in 1:length(phy$tip.label)){
        if(all(sample.tip.probs[i,] == 0)){
          sample.tip.probs[i,] <- rep(1, length(sample.tip.probs[i,]))
        }
        sample_tip_i <- edge_liks_list[[i]][1,] * sample.tip.probs[i,]
        edge_liks_list[[i]][1,] <- sample_tip_i/sum(sample_tip_i)
      }
    }
  }
  conditional_probs <- getConditionalInternodeLik(phy, Q, edge_liks_list)
  if(any(is.na(conditional_probs$root_state))){
    return(-1e6)
  }
  if(class(root.p)[1] == "character"){
    if(root.p == "yang"){
      root_liks <- c(MASS:::Null(Q))
      root_liks <- root_liks/sum(root_liks)
    }
    if(root.p == "flat"){
      root_liks <- rep(1/dim(Q)[1], dim(Q)[1])
    }
    if(root.p == "maddfitz"){
      root_liks <- conditional_probs$root_state/sum(conditional_probs$root_state)
    }
  }else{
    root_liks <- root.p/sum(root.p)
  }
  # initial sample
  internode_maps_and_discrete_probs <- getInternodeMap(phy, Q, conditional_probs$edge_liks_list, conditional_probs$root_state, root_liks, nSim, check_vector = NULL)
  internode_maps <- internode_maps_and_discrete_probs$maps
  internode_samples <- internode_maps_and_discrete_probs$state_samples
  # if additional samples are needed (didn't reach nSim), they are taken ~randomly
  if(length(internode_samples) < nSim){
    additional_sims <- nSim - length(internode_samples)
    check_vector <- unlist(lapply(internode_samples, function(x) paste0(unlist(x), collapse="")))
    random_internode_maps_and_discrete_probs <- getInternodeMap(phy, Q * 100, edge_liks_list, conditional_probs$root_state, root_liks, additional_sims, check_vector = check_vector)
    internode_maps <- c(internode_maps, random_internode_maps_and_discrete_probs$maps)
    internode_samples <- c(internode_samples, random_internode_maps_and_discrete_probs$state_samples)
  }
  # calculte the discrete probabilities based on the given Q matrix (Pij already calculated)
  discrete_probs <- lapply(internode_samples, function(x) getStateSampleProb(state_sample = x, Pij = internode_maps_and_discrete_probs$Pij, root_liks = root_liks, root_edges = internode_maps_and_discrete_probs$root_edges))
  llik_discrete <- unlist(discrete_probs)
  # generate maps
  simmaps <- getMapFromSubstHistory(internode_maps, phy)
  # if there is no character dependence the map has no influence on continuous likleihood
  character_dependence_check <- all(apply(index.cont, 1, function(x) length(unique(x)) == 1))
  if(character_dependence_check){
    llik_continuous <- OUwie.basic(simmaps[[1]], data, simmap.tree=TRUE, scaleHeight=FALSE, alpha=alpha, sigma.sq=sigma.sq, theta=theta, algorithm="three.point", tip.paths=tip.paths, mserr=mserr)
  }else{
    llik_continuous <- unlist(lapply(simmaps, function(x) OUwie.basic(x, data, simmap.tree=TRUE, scaleHeight=FALSE, alpha=alpha, sigma.sq=sigma.sq, theta=theta, algorithm="three.point", tip.paths=tip.paths, mserr=mserr)))
  }
  # combine probabilities being careful to avoid underflow
  llik_houwies <- llik_discrete + llik_continuous
  llik_houwie <- max(llik_houwies) + log(sum(exp(llik_houwies - max(llik_houwies))))
  llik_discrete <- max(llik_discrete) + log(sum(exp(llik_discrete - max(llik_discrete))))
  llik_continuous <- max(llik_continuous) + log(sum(exp(llik_continuous - max(llik_continuous))))
  if(split.liks){
    expected_vals <- lapply(simmaps, function(x) OUwie.basic(x, data, simmap.tree=TRUE, scaleHeight=FALSE, alpha=alpha, sigma.sq=sigma.sq, theta=theta, algorithm="three.point", tip.paths=tip.paths, mserr=mserr,return.expected.vals=TRUE))
    expected_vals <- colSums(do.call(rbind, expected_vals) * exp(llik_houwies - max(llik_houwies))/sum(exp(llik_houwies - max(llik_houwies))))
      return(list(TotalLik = llik_houwie, DiscLik = llik_discrete, ContLik = llik_continuous, expected_vals = expected_vals))
  }
  print(c(llik_houwie, llik_discrete, llik_continuous))
  print(p)
  return(-llik_houwie)
}

getModelParams <- function(houwie_obj){
  if(class(houwie_obj) != "houwie"){
    stop("Object must be of class houwie.")
  }
  if(is.null(houwie_obj$recon)){
    stop("hOUwie object must include a state reconstruction for model averaged parameters.")
  }
  parameter_matrix <- houwie_obj$solution.cont
  parameter_matrix[is.na(parameter_matrix)] <- 1e-10
  diag(houwie_obj$solution.disc) <- NA
  parameter_matrix <- rbind(parameter_matrix, wait.times = 1/rowSums(houwie_obj$solution.disc, na.rm = TRUE))
  parameters_by_node <- t(apply(houwie_obj$recon, 1, function(x) colSums(x * t(parameter_matrix))))
  return(parameters_by_node)
}

getEdgeLiks <- function(phy, data, n.traits, rate.cat, time_slice){
  edge_liks_list <- vector("list", dim(phy$edge)[1])
  nTip <- length(phy$tip.label)
  for(edge_i in 1:dim(phy$edge)[1]){
    # +2 because we slice the middle of the branch and need 2 terminal nodes (ancestor and descendent)
    n_slice <- (phy$edge.length[edge_i] %/% time_slice) + 2
    edge_liks_list[[edge_i]] <- matrix(1, n_slice, n.traits * rate.cat)
    if(phy$edge[edge_i,2] <= nTip){
      tmp <- numeric(n.traits)
      species_i <- phy$tip.label[phy$edge[edge_i,2]]
      state_i <- data[data[,1] == species_i, 2]
      state_i_index<- as.numeric(unlist(strsplit(as.character(state_i), "&")))
      tmp[state_i_index] <- 1
      edge_liks_list[[edge_i]][1,] <- rep(tmp, rate.cat)
    }
  }
  return(edge_liks_list)
}

getConditionalInternodeLik <- function(phy, Q, edge_liks_list){
  nTip <- length(phy$tip.label)
  external_index <- which(phy$edge[,2] <= nTip)
  # external edges
  for(edge_i in external_index){
  # move rootward along all tips to their root
    n_edges <- dim(edge_liks_list[[edge_i]])[1] - 1
    time_edge <- phy$edge.length[edge_i]
    p_mat_i <- expm(Q * (time_edge/n_edges), method=c("Ward77"))
    for(inter_edge_i in 2:(n_edges+1)){
      dec_states <- edge_liks_list[[edge_i]][inter_edge_i-1,]
      v <- edge_liks_list[[edge_i]][inter_edge_i,] * c(p_mat_i %*% dec_states )
      edge_liks_list[[edge_i]][inter_edge_i,] <- v/sum(v)
    }
  }
  # internal edges
  anc <- unique(phy$edge[,1])
  # remove the root
  root <- anc[length(anc)]
  anc <- anc[-length(anc)]
  for(anc_i in anc){
    # for the start of an internal node, combine the decs
    edge_i <- which(phy$edge[,2] == anc_i)
    dec_combo_index_i <- which(phy$edge[,1] == anc_i)
    v <- 1
    for(j in dec_combo_index_i){
      liks_j <- edge_liks_list[[j]]
      v <- v * liks_j[dim(liks_j)[1],]
    }
    v <- edge_liks_list[[edge_i]][1,] * v
    edge_liks_list[[edge_i]][1,] <- v/sum(v)
    n_edges <- dim(edge_liks_list[[edge_i]])[1] - 1
    time_edge <- phy$edge.length[edge_i]
    p_mat_i <- expm(Q * (time_edge/n_edges), method=c("Ward77"))
    for(inter_edge_i in 2:(n_edges+1)){
      dec_states <- edge_liks_list[[edge_i]][inter_edge_i-1,]
      v <- edge_liks_list[[edge_i]][inter_edge_i,] * c(p_mat_i %*% dec_states )
      edge_liks_list[[edge_i]][inter_edge_i,] <- v/sum(v)
    }
  }
  # do the root
  dec_combo_index_i <- which(phy$edge[,1] == root)
  v <- 1
  for(j in dec_combo_index_i){
    liks_j <- edge_liks_list[[j]]
    v <- v * edge_liks_list[[j]][dim(liks_j)[1],]
  }
  root_state <- v/sum(v)
  return(list(root_state = root_state,
              edge_liks_list = edge_liks_list))
}

getInternodeMap <- function(phy, Q, edge_liks_list, root_state, root_liks, nSim, check_vector=NULL){
  # set-up
  nStates <- dim(Q)[1]
  nTip <- length(phy$tip.label)
  # a potential speedup is to calculate all Pij (bollback eq.3) for all branches first
  Pij <- array(0, c(dim(Q)[1], dim(Q)[2], length(phy$edge.length)))
  # reduced edge.lengths since we are including internodes
  number_of_nodes_per_edge <- unlist(lapply(edge_liks_list, function(x) dim(x)[1]))
  number_of_edges_per_edge <- number_of_nodes_per_edge - 1
  reduced_edge_length <- phy$edge.length/number_of_edges_per_edge
  for(i in 1:length(phy$edge.length)){
    Pij[,,i] <- expm(Q * reduced_edge_length[i])
  }
  # the probability of a descendent being in state j given starting in the row of the Pj matrix
  Pj <- vector("list", length(phy$edge.length))
  for(i in 1:length(Pj)){
    Pj[[i]] <- array(0, c(dim(Q)[1], dim(Q)[2], number_of_nodes_per_edge[i]))
    for(j in 1:number_of_nodes_per_edge[i]){
      Pj[[i]][,,j] <- sweep(Pij[,,i], MARGIN = 2, edge_liks_list[[i]][j,], '*') 
    }
  }
  # simulate nSim substitution histories
  rev.pruning.order <- rev(reorder.phylo(phy, "pruningwise", index.only = TRUE))
  sub_histories <- vector("list", nSim)
  root_edges <- which(phy$edge[,1] == nTip + 1)
  edge_index <- phy$edge
  Map_i <- mapply(function(x, y) rep(x, y), x=reduced_edge_length/2, y=number_of_edges_per_edge*2, SIMPLIFY = FALSE)
  if(!is.null(check_vector)){
    state_samples <- vector("list", nSim)
    sim_counter <- 0
    while(sim_counter < nSim){
      state_sample <- getInternodeStateSample(Pj, root_state, root_edges, rev.pruning.order, edge_index, nStates, number_of_nodes_per_edge)
      current_mapping_id <- paste0(unlist(state_sample), collapse="")
      if(!current_mapping_id %in% check_vector){
        sim_counter <- sim_counter + 1
        state_samples[[sim_counter]] <- state_sample
        check_vector <- c(check_vector, current_mapping_id)
      }
    }
  }else{
    state_samples <- lapply(1:nSim, function(x) getInternodeStateSample(Pj, root_state, root_edges, rev.pruning.order, edge_index, nStates, number_of_nodes_per_edge))
    mapping_ids <- unlist(lapply(state_samples, function(x) paste0(unlist(x), collapse="")))
    state_samples <- state_samples[!duplicated(mapping_ids, nmax = 1)]
  }
  maps <- lapply(state_samples, function(x) getMapFromStateSample(Map_i, x))
  return(list(state_samples=state_samples, maps = maps, root_edges=root_edges,
              Pij = Pij, Pj = Pj))
}

getMapFromStateSample <- function(map, state_sample){
  for(edge_i in 1:length(map)){
    state_transitions <- rep(state_sample[[edge_i]][-c(1, length(state_sample[[edge_i]]))], each = 2)
    state_samples_i <- c(state_sample[[edge_i]][1],state_transitions,state_sample[[edge_i]][length(state_sample[[edge_i]])])
    names(map[[edge_i]]) <- state_samples_i
  }
  return(map)
}

getInternodeStateSample <- function(Pj, root_state, root_edge, rev.pruning.order, edge_index, nStates, number_of_nodes_per_edge){
  # each map will have edges split into equal time portions
  state_samples <- lapply(number_of_nodes_per_edge, function(x) numeric(x))
  root_sample <- sample(1:nStates, 1, prob = root_state)
  for(i in root_edge){
    state_samples[[i]][1] <- root_sample
  }
  # sample the nodes along a branch the last dec node goes into the next map
  for(edge_i in rev.pruning.order){
    from <- state_samples[[edge_i]][1]
    count <- 2
    n_inter_nodes <- length(state_samples[[edge_i]])
    for(inter_edge_i in (n_inter_nodes-1):1){
      to <- sample(1:nStates, 1, prob = Pj[[edge_i]][from,,inter_edge_i])
      from <- state_samples[[edge_i]][count] <- to
      count <- count + 1
    }
    ancestor_to_add <- edge_index[edge_i,2]
    anc_edge <- which(ancestor_to_add == edge_index[,1])
    for(i in anc_edge){
      state_samples[[i]][1] <- to
    }
  }
  return(state_samples)
}

# get path probability internal
getPathStateProb <- function(path_states, p_mat){
  P <- vector("numeric", length(path_states)-1)
  for(i in 1:(length(path_states)-1)){
    P[i] <- p_mat[path_states[1],path_states[2]]
    path_states <- path_states[-1]
  }
  return(sum(log(P)))
}

getStateSampleProb <- function(state_sample, Pij, root_liks, root_edges){
  path_probs <- numeric(length(state_sample))
  root_sample <- state_sample[[root_edges[1]]][1] # the root sample
  for(i in 1:length(state_sample)){
    path_probs[i] <- getPathStateProb(state_sample[[i]], Pij[,,i])
  }
  llik <- sum(path_probs) + log(root_liks[root_sample])
  return(llik)
}

# probability of a particular stochastic map
getMapProb <- function(simmap, Q=NULL, root_prior, p_mats=NULL){
  map <- simmap$maps
  root_state <- as.numeric(names(map[[which.min(simmap$edge[,1])]][1]))
  if(!is.null(p_mats)){
    pathway_liks <- unlist(mapply(function(x,y) getPathProb(path = x, Q=Q, p_mat=y), x = map, y = p_mats))
  }else{
    pathway_liks <- unlist(lapply(map, function(x) getPathProb(path = x, Q=Q, p_mat=NULL)))
  }
  llik <- sum(c(pathway_liks, log(root_prior)[root_state]))
  # llik <- sum(pathway_liks)
  return(llik)
}

# take substition histories and make them simmaps
getMapFromSubstHistory <- function(maps, phy){
  mapped.edge <- lapply(maps, function(x) corHMM:::convertSubHistoryToEdge(phy, x))
  obj <- vector("list", length(maps))
  for (i in 1:length(maps)){
    tree.simmap <- phy
    tree.simmap$maps <- maps[[i]]
    tree.simmap$mapped.edge <- mapped.edge[[i]]
    attr(tree.simmap, "map.order") <- "right-to-left"
    if (!inherits(tree.simmap, "simmap")) 
      class(tree.simmap) <- c("simmap", setdiff(class(tree.simmap), 
                                                "simmap"))
    obj[[i]] <- tree.simmap
  }
  if (length(maps) > 1) {
    class(obj) <- c("multiSimmap", "multiPhylo")
  }
  return(obj)
  
}
# a basic optimization for OUwie basic
OUwie.basic.dev <- function(p, phy, data, mserr, index.cont, tip.paths=NULL){
  p <- exp(p)
  index.cont[is.na(index.cont)] <- max(index.cont, na.rm = TRUE) + 1
  Rate.mat <- matrix(1, 3, dim(index.cont)[2])
  Rate.mat[] <- c(p, 1e-10)[index.cont]
  alpha = Rate.mat[1,]
  sigma.sq = Rate.mat[2,]
  theta = Rate.mat[3,]
  llik_continuous <- OUwie.basic(phy, data, simmap.tree=TRUE, scaleHeight=FALSE, alpha=alpha, sigma.sq=sigma.sq, theta=theta, algorithm="three.point", tip.paths=tip.paths, mserr=mserr)
  return(-llik_continuous)
}


# probability of the continuous parameter
OUwie.basic<-function(phy, data, simmap.tree=TRUE, root.age=NULL, scaleHeight=FALSE, root.station=FALSE, get.root.theta=FALSE, shift.point=0.5, alpha, sigma.sq, theta, mserr="none", algorithm="three.point", tip.paths=NULL, return.expected.vals=FALSE){
  
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
  if(return.expected.vals){
    return(expected.vals)
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


getAllContinuousModelStructures <- function(k, type = "OU"){
  # index.mat <- matrix(0, 3, k, dimnames = list(c("alpha", "sigma.sq", "theta"), c(1:k)))
  # we want all unique combinations of a parameter. then we can add a single all same
  # how many combinations are there of 1:k numbers?
  potential_combos <- apply(partitions:::setparts(k), 2, function(x) paste(x, collapse="_"))
  # this technically isn't all the possible alpha combinations, but for sim purposes we're fine.
  if(type == "BM"){
    alpha.combos <- paste(rep(0, k), collapse="_")
    theta.combos <- paste(rep(1, k), collapse="_")
  }
  if(type == "OU"){
    alpha.combos <- potential_combos
    theta.combos <- potential_combos
  }
  if(type == "BMOU"){
    if(k > 2){
      stop("BMOU must be manually created if k > 1 atm. Sorry.")
    }
    # needed_numerals <- 1:((2^k)-2)
    needed_numerals <- 1
    alpha.combos <- apply(sapply(needed_numerals, function(x) as.numeric(intToBits(x)[1:k])), 2, function(x) paste(x, collapse="_")) # currently doesn't allow for BM mixed with OUA
    # theta.combos <- potential_combos
    theta.combos <- paste(rep(1, k), collapse="_")
  }
  sigma.sq.combos <- potential_combos
  all_combos <- expand.grid(list(alpha.combos, sigma.sq.combos, theta.combos))
  index_mats <- array(NA, c(3, k, dim(all_combos)[1]), dimnames = list(c("alpha", "sigma.sq", "theta"), c(1:k)))
  for(i in 1:dim(all_combos)[1]){
    alpha_i <- as.numeric(unlist(strsplit(as.character(all_combos[i,1]), "_")))
    alpha_i[alpha_i == 0] <- NA
    sigma_i <- max(c(0, alpha_i), na.rm = TRUE) + as.numeric(unlist(strsplit(as.character(all_combos[i,2]), "_")))
    theta_i <- max(sigma_i) + as.numeric(unlist(strsplit(as.character(all_combos[i,3]), "_")))
    index_mats[,,i] <- rbind(alpha_i, sigma_i, theta_i)
  }
  return(index_mats)
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

organizeHOUwiePars <- function(pars, index.disc, index.cont){
  k <- max(index.disc, na.rm = TRUE)
  p.mk <- pars[1:k]
  p.ou <- pars[(k+1):length(pars)] 
  # ouwie pars
  Rate.mat <- matrix(1, 3, dim(index.disc)[2])
  index.cont[is.na(index.cont)] <- max(index.cont, na.rm = TRUE) + 1
  Rate.mat[] <- c(p.ou, NA)[index.cont]
  rownames(Rate.mat) <- rownames(index.cont)
  rate <- index.disc
  rate[is.na(rate)] <- k + 1
  Q <- matrix(0, dim(rate)[1], dim(rate)[2])
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
  while(!all(1:dim(Q)[1] %in% dat.cor)){
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

# script for generating all the possible underlying mappings and looking at joint probablity. using this we can look at the bias produced by looking only at the discrete mappings.
fixEdgeLiksLiks <- function(edge_liks_list, combo, phy, n_tips, n_nodes, n_internodes, nStates, rate.cat){
  # fix the externals
  for(j in 1:n_tips){
    tip_index <- n_nodes + n_internodes + which(phy$edge[phy$edge[,2] <= n_tips,2] == j)
    tip_state <- combo[tip_index]
    dec_edges_to_fix <- which(phy$edge[,2] == j)
    fix_vector <- numeric(nStates * rate.cat)
    fix_vector[tip_state] <- 1
    for(k in dec_edges_to_fix){
      edge_liks_list[[k]][1,] <- fix_vector
    }
  }
  # fix the internals
  for(j in 1:n_nodes){
    node_index <- unique(phy$edge[,1])[j]
    node_state <- combo[j]
    anc_edges_to_fix <- which(phy$edge[,1] == node_index)
    dec_edges_to_fix <- which(phy$edge[,2] == node_index)
    fix_vector <- numeric(nStates * rate.cat)
    fix_vector[node_state] <- 1
    for(k in dec_edges_to_fix){
      edge_liks_list[[k]][1,] <- fix_vector
    }
    for(k in anc_edges_to_fix){
      last_row <- dim(edge_liks_list[[k]])[1]
      edge_liks_list[[k]][last_row,] <- fix_vector
    }
  }
  # fix the inter nodes
  if(n_internodes > 0){
    for(j in 1:n_internodes){
      internode_index_list <- which(unlist(lapply(edge_liks_list, function(x) dim(x)[1] - 2)) >= 1)[j]
      internode_state <- combo[j + n_nodes]
      internode_index_edge <- which(apply(edge_liks_list[[internode_index_list]], 1, function(x) sum(x) > 1))[1]
      fix_vector <- numeric(nStates * rate.cat)
      fix_vector[internode_state] <- 1
      edge_liks_list[[internode_index_list]][internode_index_edge,] <- fix_vector
    }
  }
  return(edge_liks_list)
}

getAllJointProbs<- function(phy, data, rate.cat, time_slice, Q, alpha, sigma.sq, theta, quiet=TRUE){
  # prerequisites
  hOUwie.dat <- organizeHOUwieDat(data, "none", TRUE)
  nStates <- as.numeric(max(hOUwie.dat$data.cor[,2]))
  tip.paths <- lapply(1:length(phy$tip.label), function(x) OUwie:::getPathToRoot(phy, x))
  # generate the edge_liks_list
  edge_liks_list <- getEdgeLiks(phy, hOUwie.dat$data.cor, nStates, rate.cat, time_slice)
  # determine all the possible ways to fix the node states
  # how many nodes and internodes are we fixing?
  n_tips <- length(phy$tip.label)
  n_nodes <- n_tips - 1
  n_internodes <- sum(unlist(lapply(edge_liks_list, function(x) dim(x)[1] - 2)))
  # what are the possible internal combinations?
  internal_possibilities <- rep(list(1:(nStates*rate.cat)), n_nodes + n_internodes)
  external_possibilities <- lapply(edge_liks_list[phy$edge[,2] <= n_tips], function(x) which(x[1,] == 1))
  all_combinations <- expand.grid(c(internal_possibilities, external_possibilities))
  # the joint probability table
  joint_probability_table <- matrix(NA, dim(all_combinations)[1], 3, dimnames = list(1:dim(all_combinations)[1], c("disc", "cont", "total")))
  simmap_list <- vector("list", dim(all_combinations)[1])
  if(!quiet){
    cat("Begining to calcualte all possible map combinations...\n")
  }
  # for each possibility, generate an edge_liks_list to match
  for(i in 1:dim(all_combinations)[1]){
    if(!quiet){
      cat("\r", i, "of", dim(all_combinations)[1], "complete.         ")
    }
    combo_i <- as.numeric(all_combinations[i,])
    root_state <- numeric(nStates * rate.cat)
    root_state[combo_i[which(unique(phy$edge[,1]) == n_tips + 1)[1]]] <- 1
    edge_liks_list_i <- fixEdgeLiksLiks(edge_liks_list, combo_i, phy, n_tips, n_nodes, n_internodes, nStates, rate.cat)
    root_liks <- rep(1, nStates * rate.cat)/(nStates * rate.cat)
    # calculate the discrete probability of the edge_liks_list
    tmp <- getInternodeMap(phy, Q, edge_liks_list_i, root_state, root_liks, 1)
    internode_maps <- tmp$maps
    internode_samples <- tmp$state_samples
    llik_discrete <- unlist(lapply(internode_samples, function(x) getStateSampleProb(state_sample = x, Pij = tmp$Pij, root_liks = root_liks, root_edges = tmp$root_edges)))
    # generate a stochstic map
    simmaps <- getMapFromSubstHistory(internode_maps, phy)
    llik_continuous <- unlist(lapply(simmaps, function(x) OUwie.basic(x, data, simmap.tree=TRUE, scaleHeight=FALSE, alpha=alpha, sigma.sq=sigma.sq, theta=theta, algorithm="three.point", tip.paths=tip.paths, mserr="none")))
    simmap_list[[i]] <- simmaps[[1]]
    joint_probability_table[i,] <- c(llik_discrete, llik_continuous, llik_discrete + llik_continuous)
  }
  if(!quiet){
    cat("\n")
  }
  class(simmap_list) <- c("multiSimmap", "multiPhylo")
  return(list(joint_probability_table=joint_probability_table, simmap_list=simmap_list, all_combinations=all_combinations))
}

# print a houwie object
print.houwie <- function(x, ...){
  ntips <- Ntip(x$phy)
  output <- data.frame(x$loglik,x$DiscLik, x$ContLik, x$AIC,x$AICc,x$BIC, ntips, x$param.count, row.names="")
  names(output) <- c("lnLTot","lnLDisc", "lnLCont", "AIC","AICc","BIC","nTaxa","nPars")
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

silence <- function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}

# generate random starting values for multiple starts
generateMultiStarting <- function(starts, index.disc, index.cont, n_starts, lower, upper){
  n_p_trans <- max(index.disc, na.rm = TRUE)
  n_p_alpha <- length(unique(na.omit(index.cont[1,])))
  n_p_sigma <- length(unique(na.omit(index.cont[2,])))
  n_p_theta <- length(unique(na.omit(index.cont[3,])))
  multiple_starts <- vector("list", n_starts)
  multiple_starts[[1]] <- starts
  if(n_starts > 1){
    for(i in 2:n_starts){
      multiple_starts[[i]] <- numeric(length(starts))
      for(j in 1:length(starts)){
        if(starts[j]/10 <= lower[j] | starts[j]*10 >= upper[j]){
          lb <- (lower[j]+lower[j]*0.01)*10
          ub <- (upper[j]-upper[j]*0.01)/10
        }else{
          lb <- starts[j]
          ub <- starts[j]
        }
        multiple_starts[[i]][j] <- exp(runif(1, log(lb/10), log(ub*10)))
      }
    }
  }
  return(multiple_starts)
}

# organize the houwie output
getHouwieObj <- function(liks_houwie, pars, phy, data, hOUwie.dat, rate.cat, mserr, index.disc, index.cont, root.p, nSim, sample_tips, nStates, discrete_model, continuous_model, time_slice, root.station, get.root.theta,lb_discrete_model,ub_discrete_model,lb_continuous_model,ub_continuous_model,ip, opts, quiet){
  param.count <- max(index.disc, na.rm = TRUE) + max(index.cont, na.rm = TRUE)
  nb.tip <- length(phy$tip.label)
  solution <- organizeHOUwiePars(pars=pars, index.disc=index.disc, index.cont=index.cont)
  if(rate.cat > 1){
    StateNames <- paste("(", rep(1:nStates, rate.cat), rep(LETTERS[1:rate.cat], each = nStates), ")", sep = "")
  }else{
    StateNames <- paste("(", rep(1:nStates, rate.cat), ")", sep = "")
  }
  rownames(solution$solution.cor) <- colnames(solution$solution.cor) <- StateNames
  colnames(solution$solution.ou) <- StateNames
  names(hOUwie.dat$ObservedTraits) <- 1:length(hOUwie.dat$ObservedTraits)
  obj <- list(
    loglik = liks_houwie$TotalLik,
    DiscLik = liks_houwie$DiscLik,
    ContLik = liks_houwie$ContLik,
    AIC = -2*liks_houwie$TotalLik + 2*param.count,
    AICc = -2*liks_houwie$TotalLik+ 2*param.count*(param.count/(nb.tip-param.count-1)),
    BIC = -2*liks_houwie$TotalLik + log(nb.tip) * param.count,
    param.count = param.count,
    solution.disc = solution$solution.cor,
    solution.cont = solution$solution.ou,
    recon = NULL,
    index.disc = index.disc,
    index.cont = index.cont,
    phy = phy,
    legend = hOUwie.dat$ObservedTraits,
    expected_vals = liks_houwie$expected_vals,
    data = data, 
    hOUwie.dat = hOUwie.dat,
    rate.cat = rate.cat, 
    discrete_model=discrete_model, 
    continuous_model=continuous_model, 
    root.p=root.p, 
    time_slice=time_slice,
    root.station=root.station, 
    get.root.theta=get.root.theta, 
    mserr = mserr, 
    sample_tips = sample_tips,
    lb.cor=lb_discrete_model, 
    ub.cor=ub_discrete_model,
    lb.ou=lb_continuous_model, 
    ub.ou=ub_continuous_model,
    p=pars, 
    ip=ip,
    init_model = NULL,
    nSim=nSim, 
    opts=opts,
    quiet=quiet
  )
  class(obj) <- "houwie"
  return(obj)
}

