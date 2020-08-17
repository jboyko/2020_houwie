
transf.branch.lengths <-
  function(phy, model = c("BM","OUrandomRoot","OUfixedRoot","lambda","kappa","delta","EB","trend"), parameters = NULL, check.pruningwise = TRUE, check.ultrametric=TRUE, D=NULL, check.names = TRUE){
    
    n <- length(phy$tip.label)
    N <- dim(phy$edge)[1]
    ROOT <- n + 1L
    anc <- phy$edge[, 1]
    des <- phy$edge[, 2]
    externalEdge = (des <= n)
    ## Default parameters
    parameters.default = c(0,1,1,1,0,0)
    names(parameters.default) = c("alpha", "lambda", "kappa", "delta", "rate", "sigma2_error")
    
    p = list(alpha = parameters[1],
             lambda = parameters[2],
             kappa = parameters[3],
             delta = parameters[4],
             rate = parameters[5],
             sigma2_error = parameters[6]) # note that sigma2_error = true_sigma2_error/sigma2
    
    root.edge = 0 # default, holds for most models. Assumes original tree has no root edge.
    diagWeight = rep(1,n)
    errEdge = rep(p$sigma2_error,n)
    
    ## BM model
    if (model %in% c("BM","trend")) {
      edge.length = phy$edge.length
    }	
    ## OU models
    OU = c("OUrandomRoot","OUfixedRoot")
    if (model %in% OU) {
      if (check.ultrametric){
        D = numeric(n) # adjustments to external branck lengths
        if (!is.ultrametric(phy)){
          flag = 1
          dis = pruningwise.distFromRoot(phy) # has all nodes
          D = max(dis[1:n]) - dis[1:n]
          D = D - mean(D)
          phy$edge.length[externalEdge] <- phy$edge.length[externalEdge] + D[des[externalEdge]]
        }
        ## phy is now ultrametric
      } 
      times <- pruningwise.branching.times(phy) # has internal nodes only
      Tmax <- max(times)
      alpha = p$alpha
      errEdge = errEdge*exp(-2*alpha*D[des[externalEdge]]) # adjust measurement errors for OU models
      ## OUrandomRoot model	
      if (model=="OUrandomRoot") {
        distFromRoot <-  exp(-2*alpha*times) # fixit: divide by 2 alpha??
        d1 = distFromRoot[anc-n] # distFromRoot has internal nodes only, not the n external nodes.
        d2 = numeric(N)
        d2[externalEdge]  = exp(-2*alpha*D[des[externalEdge]])
        d2[!externalEdge] = distFromRoot[des[!externalEdge]-n]
      }
      ## OUfixedRoot model
      if (model=="OUfixedRoot") {	
        distFromRoot <-  exp(-2*alpha*times)*(1 - exp(-2*alpha*(Tmax-times))) # fixit: divide by 2 alpha?
        d1 = distFromRoot[anc-n]
        d2 = numeric(N)
        d2[externalEdge] = exp(-2*alpha*D[des[externalEdge]]) * (1-exp(-2*alpha*(Tmax-D[des[externalEdge]])))
        d2[!externalEdge]= distFromRoot[des[!externalEdge]-n]
      }
      edge.length = d2 - d1
      root.edge = min(distFromRoot)
      diagWeight = exp(alpha*D)
    }
    ## lambda model
    if (model=="lambda") {
      lambda = p$lambda
      distFromRoot <- pruningwise.distFromRoot(phy)
      edge.length = phy$edge.length * lambda 
      edge.length[externalEdge] = edge.length[externalEdge] + (1-lambda)*distFromRoot[des[externalEdge]]
    }
    ## kappa model
    if (model=="kappa") {
      kappa = p$kappa
      edge.length = phy$edge.length^kappa
    }
    ## delta model
    if (model=="delta") {
      delta = p$delta
      distFromRoot <- pruningwise.distFromRoot(phy)
      depth = max(distFromRoot)
      edge.length = (distFromRoot[des]^delta - distFromRoot[anc]^delta)*depth^(1-delta)
    }
    ## early burst model
    if (model=="EB") {
      rate = p$rate
      if (rate==0) edge.length = phy$edge.length
      else {
        distFromRoot <- pruningwise.distFromRoot(phy)
        edge.length = (exp(rate*distFromRoot[des])-exp(rate*distFromRoot[anc]))/rate
      }			
    }
    
    edge.length[externalEdge] = edge.length[externalEdge] + errEdge # add measurement errors to the tree
    phy$edge.length = edge.length
    phy$root.edge = root.edge
    names(diagWeight) = phy$tip.label
    return(list(tree = phy, diagWeight = diagWeight))
  }
