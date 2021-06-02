hOUwie <- function(phy, data, rate.cat, nBins,
                   model.cor=NULL, index.cor=NULL, root.p="yang", lb.cor=NULL, ub.cor=NULL, collapse=TRUE, dual=FALSE,
                   model.ou=NULL, index.ou=NULL, lb.ou=NULL, ub.ou=NULL,
                   p=NULL, ip=NULL, opts=NULL, quiet=FALSE){
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
  
  # preliminaries
  nTip <- length(phy$tip.label)
  Tmax <- max(branching.times(phy))
  organizedData <- getHOUwieCombosAndData(data, rate.cat, collapse, nBins)
  
  # setting the upper and lower bounds
  if(is.null(lb.ou)){
    # the lower limit of alpha is defined as a halflife of 10000% of the max tree height
    # the lower limit of sigma is defined 10 times less than alpha
    # the lower limit of optim is defined 10 times lower than the minimum observation
    lb.alpha = 1e-10
    lb.sigma = lb.alpha/10
    lb.optim = min(data[, dim(data)[2]])/10 
    lb.ou=c(lb.alpha,lb.sigma,lb.optim)
  }
  if(is.null(ub.ou)){
    # the upper limit of alpha is defined as a halflife of 1% of the max tree height
    # the upper limit of sigma is defined 10 times more than alpha
    # the upper limit of optim is defined 10 times more than the maximum observation
    ub.alpha = log(2)/(0.01 * Tmax)
    ub.sigma = 10 * ub.alpha
    ub.optim = max(data[, dim(data)[2]])*10 
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
  
  # if the index matrices are null, that means we need to create based on default model inputs
  if(is.null(index.cor)){
    index.cor <- corHMM:::getStateMat4Dat(organizedData$discreteData, model.cor, dual)$rate.mat
    if(rate.cat > 1){
      StateMats <- vector("list", rate.cat)
      for(i in 1:rate.cat){
        StateMats[[i]] <- index.cor
      }
      index.cor <- corHMM:::getFullMat(StateMats)
    }
  }
  
  # TO DO: add an argument for automatic CID, CD, or HYB model
  if(is.null(index.ou)){
    index.ou <- getOUParamStructure(model.ou, "three.point", FALSE, FALSE, organizedData$nStates)
  }
  
  # if user options for nloptr are not given
  if(is.null(opts)){
    opts <- list("algorithm"="NLOPT_LN_SBPLX", "maxeval"="1000000", "ftol_rel"=.Machine$double.eps^0.5)
  }
  
  if(!is.null(p)){
    if(!quiet){
      cat("Calculating likelihood from a set of fixed parameters", "\n")
      print(p)
    }
    if(max(index.ou, na.rm = TRUE) + max(index.cor, na.rm = TRUE) != length(p)){
      message <- paste0("The number of parameters does not match the number required by the model structure. You have supplied ", length(p), ", but the model structure requires ", max(index.ou, na.rm = TRUE) + max(index.cor, na.rm = TRUE), ".")
      stop(message, call. = FALSE)
    }
    out<-NULL
    est.pars<-log(p)
    out$solution <- log(p)
    out$objective <- hOUwie.dev(est.pars, phy=phy, organizedData=organizedData, rate.cat=rate.cat, index.cor=index.cor, index.ou=index.ou)
  }else{
    out<-NULL
    # check for user input initial parameters 
    if(is.null(ip)){
      means.by.regime <- with(organizedData$hOUwieData, tapply(organizedData$hOUwieData[,3], organizedData$hOUwieData[,2], mean))
      if(length(unique(na.omit(index.ou[3,]))) == length(means.by.regime)){
        start.theta <- rep(means.by.regime, length(unique(index.ou[3,]))/length(means.by.regime))
      }else{
        start.theta <- rep(mean(organizedData$hOUwieData[,3]), length(unique(na.omit(index.ou[3,]))))
      }
      start.cor <- rep(10/sum(phy$edge.length), max(index.cor, na.rm = TRUE))
      start.ou <- c(rep(log(2)/Tmax, length(unique(na.omit(index.ou[1,])))), 
                    rep(var(organizedData$hOUwieData[,3]), length(unique(na.omit(index.ou[2,])))), 
                    start.theta)
      starts = c(start.cor, start.ou)
    }else{
      starts <- ip
    }
    lower = log(c(rep(lb.cor, max(index.cor, na.rm = TRUE)), 
                  rep(lb.ou[1], length(unique(na.omit(index.ou[1,])))), 
                  rep(lb.ou[2], length(unique(na.omit(index.ou[2,])))), 
                  rep(lb.ou[3], length(unique(na.omit(index.ou[3,]))))))
    upper = log(c(rep(ub.cor, max(index.cor, na.rm = TRUE)), 
                  rep(ub.ou[1], length(unique(na.omit(index.ou[1,])))), 
                  rep(ub.ou[2], length(unique(na.omit(index.ou[2,])))), 
                  rep(ub.ou[3], length(unique(na.omit(index.ou[3,]))))))
    cat("LogLik:\n")
    out = nloptr(x0=log(starts), eval_f=hOUwie.dev, lb=lower, ub=upper, opts=opts, phy=phy, organizedData=organizedData, rate.cat=rate.cat, index.cor=index.cor, index.ou=index.ou)
    cat("\n")
  }
  return(out)
}

avgLogs <- function(x){
  x <- x[!is.na(x)]
  out <- max(x) + log(sum(exp(x - max(x))))
  return(out)
}

hOUwie.dev <- function(p, phy, organizedData, rate.cat, index.cor, index.ou){
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

  # the pupko et al algorithms the calculation of two quantities L and C
  # L can be stored in an array since its dimension will be the same for every node
  # each entry in the L array contains the joint likelihood of a particular parental transition
  # there will be nState by nBin potential combinations of states, therefore each entry in L will be nBin by nState
  nCombos <- dim(organizedData$AllCombos)[1]
  nStates <- organizedData$nStates
  AllCombos <- organizedData$AllCombos
  ObsStateMatrix <- organizedData$ObsStateMatrix
  dat <- organizedData$hOUwieData
  L_z <- array(data=NA, dim = c(nCombos, nCombos, dim(phy$edge)[1]))
  # step 1 of pupko: initialize the tips
  for(tip_index in 1:nTip){
    focal_edge_index <- match(tip_index, phy$edge[,2])
    RootAsAnc <- phy$edge[focal_edge_index,1] == nTip+1
    focal_data <- dat[match(phy$tip.label[tip_index], organizedData$hOUwieData[,1]),]
    time_edge <- phy$edge.length[focal_edge_index]
    edge_Mk <- matrix(0, length(unique(AllCombos[,1])), length(unique(AllCombos[,1])), dimnames = list(unique(AllCombos[,1]), unique(AllCombos[,1])))
    for(Mk_index in 1:nStates){
      state <- unique(AllCombos[,1])[Mk_index]
      Mk_vec <- rep(0, nStates)
      Mk_vec[state] <- 1
      edge_Mk[Mk_index,] <- log(c(expm(Q * time_edge) %*% Mk_vec))
    }
    for(i in 1:nCombos){
      # establish the possible tip combos
      for(j in 1:rate.cat){
        obsState <- ObsStateMatrix[dat[tip_index,2],j]
        value_i <- which(AllCombos[,1] == obsState & AllCombos[,2] == focal_data$valueBin) # the observed bin
        # calculate the Pij(ty) for each possible starting state i
        # start with the Mk likelihood
        state_i <- AllCombos[i,1]
        # if we are at the root, the value is just theta for the particular state
        if(RootAsAnc){
          init_i <- theta[state_i]
        }else{
          init_i <- AllCombos[i,2]
        }
        Mk_llik <- edge_Mk[state_i,obsState] # P_ij for the Mk process. the probability that we transition from i to j over time_edge is just the probability of being state_j
        # next calculate the OU likelihood
        OU_llik <- dOU(dat[tip_index,3], init_i, time_edge, alpha[state_i], theta[state_i], sqrt(sigma.sq[state_i]), TRUE)
        # OU_llik <- dOU(focal_data$valueBin, AllCombos[i,2], time_edge, alpha[state_i], theta[state_i], sigma.sq[state_i], TRUE)
        L_z[i,value_i,focal_edge_index] <- Mk_llik + OU_llik
      }
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
    edge_Mk <- matrix(0, length(unique(AllCombos[,1])), length(unique(AllCombos[,1])), dimnames = list(unique(AllCombos[,1]), unique(AllCombos[,1])))
    for(Mk_index in 1:nStates){
      state <- unique(AllCombos[,1])[Mk_index]
      Mk_vec <- rep(0, nStates)
      Mk_vec[state] <- 1
      edge_Mk[Mk_index,] <- log(c(expm(Q * time_edge) %*% Mk_vec))
    }
    # i indexes the parental state, j indexes the descendants
    for(i in 1:nCombos){
      # Mk prereqs for the parent j
      state_i <- AllCombos[i,1]
      if(RootAsAnc){
        init_i <- theta[state_i]
      }else{
        init_i <- AllCombos[i,2]
      }
      Mk_vec_i <- edge_Mk[state_i,]
      # the probability that we end up in state j at node z with the parent as state i.
      for(j in 1:nCombos){
        # given the parental node defines the distribution, we evaluate over the possible k ending values
        state_j <- AllCombos[j,1]
        Mk_llik <- Mk_vec_i[state_j] # P_ij for the Mk process. the probability that we transition from i to j over time_edge is just the probability of being state_j
        OU_llik <- dOU(AllCombos[j,2], init_i, time_edge, alpha[state_i], theta[state_i], sqrt(sigma.sq[state_i]), TRUE) # P_ij for the OU process. we observe the particular j bin, we start in the i bin, and the OU parameters are defined by the i bin's state
        # L_z[i,j,focal_anc_index] <- Mk_llik + OU_llik + sum(apply(L_z[j,,focal_dec_index], 2, max))
        L_z[i,j,focal_anc_index] <- Mk_llik + OU_llik + sum(apply(L_z[j,,focal_dec_index], 2, avgLogs)) # L_z(i) is the maximum of row j. in pupko he doesn't hold all of the values only keeping the maximum. this would be achieved by collapsing the collumns of L_z and only keeping the max
      }
    }
  }
  # step 4 of pupko: evaluate the root
  focal_dec_index <- which(phy$edge[,1] %in% (nTip+1))
  C_root <- L_root <- vector("numeric", nCombos)
  for(k in 1:nCombos){
    state_k <- AllCombos[k,1]
    #P_k <- root.p[state_k]
    #L_root[k] <- log(P_k) + sum(apply(L_z[k,,focal_dec_index], 2, max))
    # L_root[k] <- sum(apply(L_z[k,,focal_dec_index], 2, max))
    L_root[k] <- sum(apply(L_z[k,,focal_dec_index], 2, avgLogs))
    
  }
  # print(L_root)
  # llik <- max(L_root)
  llik <- max(L_root) + log(sum(exp(L_root - max(L_root))))
  cat("\r", llik)
  return(-llik)
}

# hOUwie.dev.DEVELOP <- function(p, phy, organizedData, rate.cat, index.cor, index.ou, recon=FALSE){
#   # organizing the parameters
#   p <- exp(p)
#   # define which params are for the HMM
#   k <- max(index.cor)
#   p.mk <- p[1:k]
#   Q <- index.cor
#   Q[Q!=0] <- p.mk[index.cor]
#   diag(Q) <- -rowSums(Q)
#   # set the OU params
#   p.ou <- p[(k+1):length(p)]
#   Rate.mat <- matrix(1, 3, organizedData$nStates)
#   index.ou[is.na(index.ou)] <- max(index.ou, na.rm = TRUE) + 1
#   Rate.mat[] <- c(p.ou, 1e-10)[index.ou]
#   alpha = Rate.mat[1,]
#   sigma.sq = Rate.mat[2,]
#   theta = Rate.mat[3,]
# 
#   # the pupko et al algorithms the calculation of two quantities L and C
#   # L can be stored in an array since its dimension will be the same for every node
#   # each entry in the L array contains the joint likelihood of a particular parental transition
#   # there will be nState by nBin potential combinations of states, therefore each entry in L will be nBin by nState
#   nCombos <- max(organizedData$ComboMatrix)
#   nStates <- organizedData$nStates
#   AllCombos <- organizedData$AllCombos
#   namedStates <- organizedData$namedStates
#   ObsStateMatrix <- organizedData$ObsStateMatrix
#   dat <- organizedData$hOUwieData
#   L_z <- array(data=-Inf, dim = c(nCombos, nCombos, dim(phy$edge)[1]))
#   # step 1 of pupko: initialize the tips
#   for(tip_index in 1:nTip){
#     focal_edge_index <- match(tip_index, phy$edge[,2])
#     RootAsAnc <- phy$edge[focal_edge_index,1] == nTip+1
#     time_edge <- phy$edge.length[focal_edge_index]
#     edge_Mk <- matrix(0, length(unique(AllCombos[,1])), length(unique(AllCombos[,1])), dimnames = list(unique(AllCombos[,1]), unique(AllCombos[,1])))
#     # establish the possible tip combos
#     focal_data <- dat[match(phy$tip.label[tip_index], organizedData$hOUwieData[,1]),]
#     TipStates <- as.numeric(strsplit(dat[tip_index,2], "&")[[1]])
#     TipStates <- paste("(", rep(TipStates, rate.cat), rep(LETTERS[1:rate.cat], each = length(TipStates)), ")", sep = "")
#     TipStates <- namedStates[match(TipStates, names(namedStates))]
#     nTipStates <- length(TipStates)
#     if(dat[tip_index,3] == "?"){
#       TipCont <- unique(AllCombos$Bin)
#     }else{
#       TipCont <- dat[tip_index,3]
#     }
#     TipContIndex <- focal_data$valueBin
#     nTipCont <- length(TipCont)
#     nTipCombos <- nTipCont * nTipStates
#     tipCombos <- expand.grid(TipStates, TipCont)
# 
#     # calculate all possible Mk outcomes beforehand (time deosnt change for a particular tip)
#     for(Mk_index in 1:nStates){
#       state <- unique(AllCombos[,1])[Mk_index]
#       Mk_vec <- rep(0, nStates)
#       Mk_vec[state] <- 1
#       edge_Mk[Mk_index,] <- log(c(expm(Q * time_edge) %*% Mk_vec))
#     }
#     for(i in 1:nCombos){
#       for(j in 1:nTipCombos){
#         obsState <- as.numeric(tipCombos[j,1])
#         binValue <- as.numeric(tipCombos[j,2])
#         if(TipContIndex == "?"){
#           TipContIndex = binValue
#         }
#         value_j <- which(AllCombos[,1] == obsState & AllCombos[,2] == TipContIndex) # the observed bin
#         # calculate the Pij(ty) for each possible starting state i
#         # start with the Mk likelihood
#         state_i <- AllCombos[i,1]
#         # if we are at the root, the value is just theta for the particular state
#         if(RootAsAnc){
#           init_i <- theta[state_i]
#         }else{
#           init_i <- AllCombos[i,2]
#         }
#         Mk_llik <- edge_Mk[state_i,obsState] # P_ij for the Mk process. the probability that we transition from i to j over time_edge is just the probability of being state_j
#         # next calculate the OU likelihood
#         OU_llik <- dOU(binValue, init_i, time_edge, alpha[state_i], theta[state_i], sigma.sq[state_i], TRUE)
#         # OU_llik <- dOU(focal_data$valueBin, AllCombos[i,2], time_edge, alpha[state_i], theta[state_i], sigma.sq[state_i], TRUE)
# 
#         L_z[i,value_j,focal_edge_index] <- Mk_llik + OU_llik
#       }
#     }
#   }
#   # step 2 of pupko: calculations for the internal nodes, excluding the root
#   pruningwise.index <- unique(phy$edge[reorder(phy, "pruningwise", index.only = TRUE), 1])
#   for(node_index in pruningwise.index[-length(pruningwise.index)]){
#     focal_anc_index <- which(phy$edge[,2] %in% node_index)
#     focal_dec_index <- which(phy$edge[,1] %in% node_index)
#     RootAsAnc <- phy$edge[focal_anc_index,1] == nTip+1
#     # time is based on the parent to focal length. dec values calculated already
#     time_edge <- phy$edge.length[focal_anc_index]
#     for(Mk_index in 1:nStates){
#       state <- unique(AllCombos[,1])[Mk_index]
#       Mk_vec <- rep(0, nStates)
#       Mk_vec[state] <- 1
#       edge_Mk[Mk_index,] <- log(c(expm(Q * time_edge) %*% Mk_vec))
#     }
#     # i indexes the parental state, j indexes the descendants
#     for(i in 1:nCombos){
#       # Mk prereqs for the parent j
#       state_i <- AllCombos[i,1]
#       if(RootAsAnc){
#         init_i <- theta[state_i]
#       }else{
#         init_i <- AllCombos[i,2]
#       }
#       Mk_vec_i <- edge_Mk[state_i,]
#       # the probability that we end up in state j at node z with the parent as state i.
#       for(j in 1:nCombos){
#         # given the parental node defines the distribution, we evaluate over the possible k ending values
#         state_j <- AllCombos[j,1]
#         Mk_llik <- Mk_vec_i[state_j] # P_ij for the Mk process. the probability that we transition from i to j over time_edge is just the probability of being state_j
#         OU_llik <- dOU(AllCombos[j,2], init_i, time_edge, alpha[state_i], theta[state_i], sigma.sq[state_i], TRUE) # P_ij for the OU process. we observe the particular j bin, we start in the i bin, and the OU parameters are defined by the i bin's state
#         L_z[i,j,focal_anc_index] <- Mk_llik + OU_llik + sum(apply(L_z[j,,focal_dec_index], 2, max)) # L_z(i) is the maximum of row j. in pupko he doesn't hold all of the values only keeping the maximum. this would be achieved by collapsing the collumns of L_z and only keeping the max
#       }
#     }
#   }
#   # step 4 of pupko: evaluate the root
#   focal_dec_index <- which(phy$edge[,1] %in% (nTip+1))
#   C_root <- L_root <- vector("numeric", nCombos)
#   for(k in 1:nCombos){
#     state_k <- AllCombos[k,1]
#     #P_k <- root.p[state_k]
#     #L_root[k] <- log(P_k) + sum(apply(L_z[k,,focal_dec_index], 2, max))
#     L_root[k] <- sum(apply(L_z[k,,focal_dec_index], 2, max))
#   }
#   # print(L_root)
#   llik <- max(L_root)
#   cat("\r", llik)
#   return(-llik)
# }


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

hOUwieRecon <- function(phy, data, rate.cat, nBins, p,
                        model.cor=NULL, index.cor=NULL, root.p="yang", collapse=TRUE, dual=FALSE,
                        model.ou=NULL, index.ou=NULL, lb.ou=NULL, ub.ou=NULL){
  organizedData <- getHOUwieCombosAndData(data, rate.cat, collapse, nBins)
  if(is.null(index.cor)){
    index.cor <- corHMM:::getStateMat4Dat(organizedData$discreteData, model.cor, dual)$rate.mat
    if(rate.cat > 1){
      StateMats <- vector("list", rate.cat)
      for(i in 1:rate.cat){
        StateMats[[i]] <- index.cor
      }
      index.cor <- corHMM:::getFullMat(StateMats)
    }
  }
  # TO DO: add an argument for automatic CID, CD, or HYB model
  if(is.null(index.ou)){
    index.ou <- getOUParamStructure(model.ou, "three.point", FALSE, FALSE, organizedData$nStates)
  }
  JointRecon <- hOUwieRecon.dev(log(p), phy=phy, organizedData=organizedData, rate.cat=rate.cat, index.cor=index.cor, index.ou=index.ou)
  ExteriorReconDF <- data.frame(State = organizedData$AllCombos[JointRecon,1], Value = organizedData$AllCombos[JointRecon,2], Index = JointRecon)
  InteriorReconDF <- ExteriorReconDF[(length(phy$tip.label)+1):dim(ExteriorReconDF)[1],]
  ExteriorReconDF <- ExteriorReconDF[1:length(phy$tip.label),]
  return(list(TipStates = ExteriorReconDF,
              NodeStates = InteriorReconDF,
              AllCombos = organizedData$AllCombos))
}

hOUwieRecon.dev <- function(p, phy, organizedData, rate.cat, index.cor, index.ou){
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
    for(Mk_index in 1:nStates){
      state <- unique(AllCombos[,1])[Mk_index]
      Mk_vec <- rep(0, nStates)
      Mk_vec[state] <- 1
      edge_Mk[Mk_index,] <- log(c(expm(Q * time_edge) %*% Mk_vec))
    }
    for(i in 1:nCombos){
      for(j in 1:rate.cat){
        obsState <- ObsStateMatrix[dat[tip_index,2],j]
        value_i <- which(AllCombos[,1] == obsState & AllCombos[,2] == focal_data$valueBin) # the observed bin
        # calculate the Pij(ty) for each possible starting state i
        # start with the Mk likelihood
        state_i <- AllCombos[i,1]
        if(RootAsAnc){
          init_i <- theta[state_i]
        }else{
          init_i <- AllCombos[i,2]
        }
        Mk_llik <- edge_Mk[state_i,obsState] # P_ij for the Mk process. the probability that we transition from i to j over time_edge is just the probability of being state_j
        # next calculate the OU likelihood
        OU_llik <- dOU(dat[tip_index,3], init_i, time_edge, alpha[state_i], theta[state_i], sigma.sq[state_i], TRUE)
        L_z[i,value_i,focal_edge_index] <- Mk_llik + OU_llik
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
    for(Mk_index in 1:nStates){
      state <- unique(AllCombos[,1])[Mk_index]
      Mk_vec <- rep(0, nStates)
      Mk_vec[state] <- 1
      edge_Mk[Mk_index,] <- log(c(expm(Q * time_edge) %*% Mk_vec))
    }
    # i indexes the parental state, j indexes the descendants 
    for(i in 1:nCombos){
      # Mk prereqs for the parent j
      state_i <- AllCombos[i,1]
      if(RootAsAnc){
        init_i <- theta[state_i]
      }else{
        init_i <- AllCombos[i,2]
      }
      Mk_vec_i <- edge_Mk[state_i,] 
      # the probability that we end up in state j at node z with the parent as state i. 
      for(j in 1:nCombos){
        # given the parental node defines the distribution, we evaluate over the possible k ending values
        state_j <- AllCombos[j,1]
        Mk_llik <- Mk_vec_i[state_j] # P_ij for the Mk process. the probability that we transition from i to j over time_edge is just the probability of being state_j
        OU_llik <- dOU(AllCombos[j,2], init_i, time_edge, alpha[state_i], theta[state_i], sigma.sq[state_i], TRUE) # P_ij for the OU process. we observe the particular j bin, we start in the i bin, and the OU parameters are defined by the i bin's state
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
  print(llik)
  return(JointComboRecon)
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
