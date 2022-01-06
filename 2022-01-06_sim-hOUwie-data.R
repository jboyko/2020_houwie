#### #### #### #### #### #### #### #### #### #### #### #### 
# imports
#### #### #### #### #### #### #### #### #### #### #### #### 

setwd("~/2020_hOUwie/")

source("hOUwieNode.R")
require(OUwie)
require(corHMM)
require(parallel)
require(phytools)
require(expm)
require(POUMM)
require(geiger)
require(partitions)

#### #### #### #### #### #### #### #### #### #### #### #### 
# Functions
#### #### #### #### #### #### #### #### #### #### #### #### 

# a function which generates a parameter set given a continuous model
generateParameters <- function(continuous_model, alpha, sigma2, theta, discrete_model, rate){
  k.alpha <- length(unique(na.omit(continuous_model[1,])))
  k.sigma <- length(unique(na.omit(continuous_model[2,])))
  k.theta <- length(unique(na.omit(continuous_model[3,])))
  k.discrete <- max(discrete_model, na.rm=TRUE)
  # par_alpha <- runif(k.alpha, minAlpha, maxAlpha)
  # par_sigma <- runif(k.sigma, minSigma2, maxSigma2)
  # par_theta <- runif(k.theta, minTheta, maxTheta)
  par_alpha <- alpha[1:k.alpha]
  par_sigma <- sigma2[1:k.sigma]
  par_theta <- theta[1:k.theta]
  par_rates <- rate[1:k.discrete]
  p <- c(par_rates, par_alpha, par_sigma, par_theta)
  return(p)
}

makeData <- function(nTip, continuous_model, discrete_model, alpha, sigma2, theta, rate){
  phy <- sim.bdtree(b = 1, d = 0, stop = "taxa", n = nTip) 
  phy <- drop.extinct(phy)
  phy$edge.length <- phy$edge.length/max(branching.times(phy))
  pars <- generateParameters(continuous_model, alpha, sigma2, theta, discrete_model, rate)
  full_data <- try(generateData(phy, discrete_model, continuous_model, pars))
  while(class(full_data) == "try-error"){
    full_data <- try(generateData(phy, discrete_model, continuous_model, pars))
  }
  return(full_data)
}

generateData <- function(phy, index.cor, index.ou, pars, type ="CD", root.p="yee", quiet=FALSE){
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
  if(root.p=="sample"){
    root.p = rep(0, dim(Q)[1]) # we will sample the root with equal probability of either state
    root.p[sample(1:dim(Q)[1], 1)] <- 1
  }else{
    root.p = rep(0, dim(Q)[1])
    root.p[1] <- 1
  }
  # organize ou params
  Rate.mat <- matrix(1, 3, dim(index.cor)[2])
  index.ou[is.na(index.ou)] <- max(index.ou, na.rm = TRUE) + 1
  Rate.mat[] <- c(p.ou, 1e-10)[index.ou]
  alpha = Rate.mat[1,]
  sigma.sq = Rate.mat[2,]
  theta = Rate.mat[3,]
  theta0 = theta[which(root.p == 1)]
  full.data <- hOUwie.sim(phy, Q, root.p, alpha, sigma.sq, theta0, theta)
  obs.no.trans <- sum(unlist(lapply(full.data$simmap[[1]]$maps, function(x) length(x) - 1)))
  if(!quiet){
    cat("The observed number of transitions was found to be", obs.no.trans, "\n")
  }
  simulators <- list(index.cor = index.cor, index.ou = index.ou, pars = pars)
  return(c(full.data, simulators))
}

generateDataset <- function(nTip, continuous_model_cd, continuous_model_cid, discrete_model_cd, discrete_model_cid, alpha, sigma2, theta, rate){
  # generate a dataset for CD
  dat_cd <- makeData(nTip, continuous_model_cd, discrete_model_cd, alpha, sigma2, theta, rate)
  # ensure sampling of both observed states
  tip_count <- length(which(dat_cd$data[,2] == 2))
  while(!(tip_count > nTip*0.25 & tip_count < nTip*.75)){
    dat_cd <- makeData(nTip, continuous_model_cd, discrete_model_cd, alpha, sigma2, theta, rate)
    tip_count <- length(which(dat_cd$data[,2] == 2))
  }
  # generate a dataset for CID
  dat_cid <- makeData(nTip, continuous_model_cid, discrete_model_cid, alpha, sigma2, theta, rate)
  # ensure sampling of both hidden states
  tip_count <- length(which(dat_cid$data[,2] == 3 | dat_cid$data[,2] == 4))
  while(tip_count > nTip*0.25 & tip_count < nTip*.75){
    dat_cid <- makeData(nTip, continuous_model_cid, discrete_model_cid, alpha, sigma2, theta, rate)
    tip_count <- length(which(dat_cid$data[,2] == 3 | dat_cid$data[,2] == 4))
  }
  out <- list(dat_cd=dat_cd, dat_cid=dat_cid)
  return(out)
}


#### #### #### #### #### #### #### #### #### #### #### #### 
# prerequisites
#### #### #### #### #### #### #### #### #### #### #### #### 

nTip <- 25
alpha <- c(10,10)
sigma2 <- alpha/5
theta <- c(12,24)
rate <- .1

#### #### #### #### #### #### #### #### #### #### #### #### 
# the 22 2-state models we want to fit - put into list form
#### #### #### #### #### #### #### #### #### #### #### #### 

continuous_models_cd_ou <- getAllContinuousModelStructures(2, "OU")
continuous_models_cd_bm <- getAllContinuousModelStructures(2, "BM")
continuous_models_cd_bmou <- getAllContinuousModelStructures(2, "BMOU")
continuous_models_cid_ou <- continuous_models_cd_ou[,c(1,1,2,2),]
continuous_models_cid_bm <- continuous_models_cd_bm[,c(1,1,2,2),]
continuous_models_cid_bmou <- continuous_models_cd_bmou[,c(1,1,2,2),]
continuous_models_cd_ou <- lapply(seq(dim(continuous_models_cd_ou)[3]), function(x) continuous_models_cd_ou[,,x])
continuous_models_cd_bm <- lapply(seq(dim(continuous_models_cd_bm)[3]), function(x) continuous_models_cd_bm[,,x])
continuous_models_cd_bmou <- lapply(seq(dim(continuous_models_cd_bmou)[3]), function(x) continuous_models_cd_bmou[,,x])
continuous_models_cid_ou <- lapply(seq(dim(continuous_models_cid_ou)[3]), function(x) continuous_models_cid_ou[,,x])
continuous_models_cid_bm <- lapply(seq(dim(continuous_models_cid_bm)[3]), function(x) continuous_models_cid_bm[,,x])
continuous_models_cid_bmou <- lapply(1:dim(continuous_models_cid_bmou)[3], function(x) continuous_models_cid_bmou[,,x])
all_model_structures <- c(continuous_models_cd_bm, continuous_models_cid_bm[-1], continuous_models_cd_ou, continuous_models_cid_ou[-1], continuous_models_cd_bmou, continuous_models_cid_bmou)
model_names <- c("CID_BM1", "CD_BMV", "CID+_BMV", "CID_OU1", "CD_OUA", "CD_OUV", "CD_OUVA", "CD_OUM",
                 "CD_OUMA", "CD_OUMV", "CD_OUMVA", "CID+_OUA", "CID+_OUV", "CID+_OUVA", "CID+_OUM", "CID+_OUMA",
                 "CID+_OUMV", "CID+_OUMVA", "CD_OUBM1", "CD_OUBMV", "CID+_OUBM1", "CID+_OUBMV")
names(all_model_structures) <- model_names

# the 2 discrete models being evaluated
discrete_model_cd <- equateStateMatPars(getRateCatMat(2), 1:2)
discrete_model_cid <- getFullMat(list(discrete_model_cd, discrete_model_cd), getRateCatMat(2))
discrete_model_cid <- equateStateMatPars(discrete_model_cid, c(1,2,3,4))

#### #### #### #### #### #### #### #### #### #### #### #### 
# run
#### #### #### #### #### #### #### #### #### #### #### #### 

# models 1 and 4 are standard CID and should always be included
# models 2 and 3 pair as BMV
# models 5 and 12 pair as OUA
# models 6 and 13 pair as OUV
# models 7 and 14 pair as OUVA
# models 8 and 15 pair as OUM
# models 9 and 16 pair as OUMA
# models 10 and 17 pair as OUMV
# models 11 and 18 pair as OUMVA
# models 19 and 21 pair as OUBM1
# models 20 and 22 pair as OUBMV

# for each model pair i want to generate a dataset consistent with the CD and one consistent with CID+
generateDataset(nTip, continuous_model_cd, continuous_model_cid, discrete_model_cd, discrete_model_cid, alpha, sigma2, theta, rate)



continuous_model_cd <- all_model_structures[[8]]
continuous_model_cid <- all_model_structures[[15]]
OU1_classic <- all_model_structures[[4]]

out <- list()
map_number <- c(100, 250, 500, 1000)
for(i in 1:length(map_number)){
  nmaps <- map_number[i]
  out[[i]] <- mclapply(1:20, function(x) run_hmm_help(nmaps), mc.cores = 5)
}

# nmap_100 <- mclapply(1:20, function(x) run_hmm_help(100), mc.cores = 5)
# nmap_200 <- mclapply(1:20, function(x) run_hmm_help(200), mc.cores = 5)
# nmap_1000 <- mclapply(1:20, function(x) run_hmm_help(1000), mc.cores = 5)
# nmap_10000 <- mclapply(1:20, function(x) run_hmm_help(10000), mc.cores = 20)
