setwd("~/2020_hOUwie/")

source("hOUwieNode.R")
source("Utils.R")

require(OUwie)
require(corHMM)
require(parallel)
require(expm)

#### #### #### #### #### #### #### #### #### #### #### #### 
# prerequisites
#### #### #### #### #### #### #### #### #### #### #### #### 

nTip = 100
nSim = 50
alpha = 1.5
sigma_sq = 0.75
theta_1 = 2
theta_2 = 4
theta_0 = 2
rate = 0.5

#### #### #### #### #### #### #### #### #### #### #### #### 
# the 40 2-state models we want to fit - put into list form
#### #### #### #### #### #### #### #### #### #### #### #### 

continuous_model <- getAllContinuousModelStructures(2, "OU")[,,5]

# the 2 discrete models being evaluated
discrete_model <- equateStateMatPars(getRateCatMat(2), 1:2)

# deciding what tree size we would like to use
phy <- sim.bdtree(b = 1, d = 0, stop = "taxa", n = nTip) 
phy <- drop.extinct(phy)
phy$edge.length <- phy$edge.length/max(branching.times(phy))

#### #### #### #### #### #### #### #### #### #### #### #### 
# Functions
#### #### #### #### #### #### #### #### #### #### #### #### 

SingleRun <- function(opt_combo, full_data, iter){
  n_starts <- n_cores <- as.numeric(opt_combo[1])
  ip <- as.character(opt_combo[1,2])
  optimizer <- as.character(opt_combo[1,3])
  if(optimizer == "nlopt_ln"){
    opts = list("algorithm"="NLOPT_LN_SBPLX", "maxeval"="1000", "ftol_rel"=.Machine$double.eps^0.5)
  }
  if(optimizer == "nlopt_gn"){
    opts = list("algorithm"="NLOPT_GN_DIRECT_L", "maxeval"="1000", "ftol_rel"=.Machine$double.eps^0.5)
  }
  if(optimizer == "sann"){
    opts=list(max.call=1000, smooth=FALSE)
  }
  fit <- hOUwie(phy = phy, data = data, rate.cat = 1, nSim = nSim, discrete_model = discrete_model, continuous_model = continuous_model, recon = FALSE, optimizer = gsub("_.*", "", optimizer), opts=opts, ip = ip, n_starts = n_starts, ncores = n_cores)
  obj <- list(sim = full_data, fit = fit)
  file_name <- paste0("optim_test/", optimizer, "-", n_starts, "-", ip, "-", iter,  ".Rsave")
  save(obj, file = file_name)
  return(obj)
}


#### #### #### #### #### #### #### #### #### #### #### #### 
# running
#### #### #### #### #### #### #### #### #### #### #### #### 

n_start_options = c(1,10)
ip_options = c("fast", "good")
optimizer_options = c("nlopt_ln", "nlopt_gn", "sann")
opimizer_combinations <- expand.grid(n_start_options, ip_options, optimizer_options)
opimizer_combinations <- lapply(seq_len(nrow(opimizer_combinations)), function(i) opimizer_combinations[i,])

yee <- vector("list", 10)
for(iter in 1:10){
  phy <- sim.bdtree(b = 1, d = 0, stop = "taxa", n = nTip) 
  phy <- drop.extinct(phy)
  phy$edge.length <- phy$edge.length/max(branching.times(phy))
  full_data <- try(generateData(phy, discrete_model, continuous_model, c(rate, alpha, sigma_sq, theta_1, theta_2)))
  while(class(full_data) == "try-error"){
    full_data <- try(generateData(phy, discrete_model, generating_model, pars))
  }
  yee[[iter]] <- mclapply(opimizer_combinations, function(x) SingleRun(x, full_data, iter), mc.cores = 12)
}




