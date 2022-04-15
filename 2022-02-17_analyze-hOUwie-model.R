#### #### #### #### #### #### #### #### #### #### #### #### 
# imports
#### #### #### #### #### #### #### #### #### #### #### #### 

setwd("~/2020_hOUwie/")

source("hOUwieNode.R")
source("Utils.R")

require(OUwie)
require(corHMM)
require(parallel)
require(phytools)
require(expm)
require(POUMM)
require(geiger)
require(partitions)
require(ggplot2)
require(gridExtra)
require(reshape)

#### #### #### #### #### #### #### #### #### #### #### #### 
# Functions
#### #### #### #### #### #### #### #### #### #### #### #### 
getFocalFiles <- function(model_type, ntip, nmap){
  focal_dir <- paste0("simulated_fit/", model_type, "/", ntip)
  fit_files <- dir(focal_dir, full.names = TRUE)
  fit_files <- fit_files[grep(paste0("nmap=", nmap, ".Rsave"), fit_files)]
  return(fit_files)
}

summarizeFile <- function(focal_file){
  load(focal_file)
  # summarize character dependent models
  cd_sim <- out$cd_out
  model_table_cd <- getModelTable(cd_sim)
  pars_cd <- c(cd_sim$cd_fit$p[1], t(cd_sim$cd_fit$solution.cont))
  index_mat <- cd_sim$cd_fit$index.cont
  names(pars_cd) <- c("rate", "alpha_1", "alpha_2", "sigma.sq_1", "sigma.sq_2", "theta_1", "theta_2")
  if(mean(cd_sim$bm1_fit$data[,3]) > 25){
    pars_cd[c(6,7)] <- pars_cd[c(6,7)] - 50
  }
  AICwts <- model_table_cd$AICwt
  names(AICwts) <- c("BM1", "OU1", "CD", "CID+")
  cd_out <- c(k=model_table_cd$np[3], loglik=model_table_cd$lnLik[3], AICwt=AICwts, best_model=model_table_cd$dAIC[3]==0, pars_cd)
  
  # summarize character independent models
  cid_sim <- out$cid_out
  model_table_cid <- getModelTable(cid_sim)
  pars_cid <- c(cid_sim$cid_fit$p[1], t(cid_sim$cid_fit$solution.cont[,c(1,3)]))
  names(pars_cid) <- c("rate", "alpha_1", "alpha_2", "sigma.sq_1", "sigma.sq_2", "theta_1", "theta_2")
  if(mean(cid_sim$bm1_fit$data[,3]) > 25){
    pars_cid[c(6,7)] <- pars_cid[c(6,7)] - 50
  }
  AICwts <- model_table_cid$AICwt
  names(AICwts) <- c("BM1", "OU1", "CD", "CID+")
  cid_out <- c(k=model_table_cid$np[4], loglik=model_table_cid$lnLik[4], AICwt=AICwts, best_model=model_table_cid$dAIC[4]==0, pars_cid)
  
  return(list(cd_out=cd_out, cid_out=cid_out))
}

organizeTable <- function(big_table){
  RMSE <- sqrt(colMeans(big_table[,c(4:10, 14:20)]^2))
  cd_results <- c(big_table[1,1], colMeans(big_table[,2:3]), RMSE[1:7])
  cid_results <- c(big_table[11,1], colMeans(big_table[,12:13]), RMSE[8:14])
  tmp_1 <- cbind(as.data.frame("CD"), t(as.data.frame(cd_results)))
  colnames(tmp_1)[1] <- "model_class"
  tmp_2 <- cbind(as.data.frame("CID"), t(as.data.frame(cid_results)))
  colnames(tmp_2)[1] <- "model_class"
  out <- cbind(tmp_1, tmp_2)
  return(out)
}

getError <- function(table_row, par_table){
  focal_pars <- par_table[,table_row[1] == colnames(par_table)]
  if(any(as.numeric(table_row[14:15]) > 25)){
    table_row[14:15] <- table_row[14:15] - 50
  }
  error <- as.numeric(table_row[9:15]) - focal_pars
  return(error)
}

reconstructCIDFits <- function(focal_file){
  load(focal_file)
  if(is.null(out$cd_out$cid_fit$recon) | is.null(out$cid_out$cid_fit$recon)){
    cid_fit_cd <- out$cd_out$cid_fit
    cid_fit_cid <- out$cid_out$cid_fit
    cid_fit_cid$sample_nodes <- cid_fit_cd$sample_nodes <- TRUE
    cid_fit_cid$adaptive_sampling <- cid_fit_cd$adaptive_sampling <- TRUE
    cid_recon_cd <- hOUwieRecon(cid_fit_cd, nodes = "external")
    cid_recon_cid <- hOUwieRecon(cid_fit_cid, nodes = "external")
    out$cd_out$cid_fit$recon <- cid_recon_cd
    out$cid_out$cid_fit$recon <- cid_recon_cid
    save(out, file = focal_file)
  }else{
    print("already reconstructed.")
    return(NULL)
  }
}

#### #### #### #### #### #### #### #### #### #### #### #### 
# run
#### #### #### #### #### #### #### #### #### #### #### #### 

### testing
# model_types <- c("BMV", "OUV", "OUA", "OUM", "OUVA", "OUMV", "OUMA", "OUMVA", "OUBM1", "OUBMV")
# for(model_i in model_types){
#   print(model_i)
#   focal_files <- getFocalFiles(model_i, 100, 25)
#   sapply(focal_files, reconstructCIDFits)
# }
# focal_files <- getFocalFiles("OUM", 25, 25)
# focal_file <- focal_files[10]
# # reconstructCIDFits(focal_files[10])
# load(focal_files[10])
# lapply(out$cd_out, "[[", "AIC")
# lapply(out$cid_out, "[[", "AIC")
# 
# phy <- out$cid_out$bm1_fit$phy
# dat <- out$cid_out$bm1_fit$dat
# discrete_model <- out$cid_out$cid_fit$discrete_model
# continuous_model <- out$cid_out$cid_fit$continuous_model
# 
# new_fit <- hOUwie(phy = phy, data = dat, rate.cat = 2, nSim = 100, time_slice = 1.1, discrete_model = discrete_model, continuous_model = continuous_model, recon = FALSE, sample_tips = FALSE, sample_nodes = TRUE, adaptive_sampling = TRUE, optimizer = "nlopt_ln")
# 
# houwie_res <- out$cd_out$cid_fit
# houwie_res$adaptive_sampling <- TRUE
# houwie_res$sample_nodes <- TRUE
# test <- hOUwieRecon(houwie_res, "external")
# 
# out$cd_out$cd_fit
# colSums(t(test) * houwie_res$solution.cont[2,])
# 
# 
# debug(hOUwieRecon)
# 
# 
# getModelAvgParams(houwie_res)

# for each model type i want to generate a dataset consistent with the CD and one consistent with CID+
model_types <- c("BMV", "OUV", "OUA", "OUM", "OUVA", "OUMV", "OUMA", "OUMVA", "OUBM1", "OUBMV")

ntips <- c(25, 100, 250)
nmaps <- c(25, 100, 250)

# model_type <- "OUA"
# ntip <- 25
# nmap <- 25
#

bigger_table <- list()
count <- 1
total_count <- length(model_types) * length(nmaps) * length(ntips)
for(k in 1:length(nmaps)){
  nmap <- nmaps[k]
  for(j in 1:length(ntips)){
    ntip <- ntips[j]
    for(i in 1:length(model_types)){
      model_type <- model_types[i]
      focal_files <- getFocalFiles(model_type, ntip, nmap)
      if(length(focal_files) == 0){
        next
      }
      list_of_errors <- lapply(focal_files, summarizeFile)
      cd_errors <- do.call(rbind, lapply(list_of_errors, "[[", "cd_out"))
      cd_errors <- cbind(model_type = model_type, model_class="CD", nTip=ntip, nMap=nmap, as.data.frame(cd_errors))
      cid_errors <- do.call(rbind, lapply(list_of_errors, "[[", "cid_out"))
      cid_errors <- cbind(model_type = model_type, model_class="CID", nTip=ntip, nMap=nmap, as.data.frame(cid_errors))
      big_table <- rbind(cd_errors, cid_errors)
      bigger_table[[count]] <- big_table
      cat("\r", round(count/total_count, 3) * 100, "% complete ...")
      count <- count + 1
    }
  }
}

biggest_table <- do.call(rbind, bigger_table)
rownames(biggest_table) <- NULL

# start here

load("biggest_table.Rsave")
model_types <- c("BMV", "OUV", "OUA", "OUM", "OUVA", "OUMV", "OUMA", "OUMVA", "OUBM1", "OUBMV")
biggest_table <- as.data.frame(biggest_table)
biggest_table$model_type <- factor(biggest_table$model_type, model_types)
biggest_table$model_class <- as.factor(biggest_table$model_class)
biggest_table$nTip <- as.factor(biggest_table$nTip)
biggest_table$nMap <- as.factor(biggest_table$nMap)


alpha <- c(3, 1.5)
sigma.sq <- c(0.35, 1)
theta <- c(2, 0.75)
rate <- .1
pars <- c(rate, alpha, sigma.sq, theta)
par_index_table <- matrix(data = 1, nrow = 7, ncol = 10, dimnames = list(c("rate", "alpha_1", "alpha_2", "sigma.sq_1", "sigma.sq_2", "theta_1", "theta_2"), model_types))
par_index_table[2,1] <- NA
par_index_table[3,c(1,2,4,6,9,10)] <- NA
par_index_table[5,c(3,4,7,9)] <- NA
par_index_table[7,c(1,2,3,5,9,10)] <- NA
# removing parameters that are not estimates (e.g., theta_1 = theta_2 for BMV)
for(i in 1:nrow(biggest_table)){
  cat("\r", round(i/nrow(biggest_table) * 100), "%")
  focal_par_index <- par_index_table[,biggest_table[i,1] == colnames(par_index_table)]
  start_col <- which(names(biggest_table[i,]) == "rate")
  end_col <- length(names(biggest_table[i,]))
  if(biggest_table[i, 2] == "CID"){
    biggest_table[i, 13:14] <- sort(biggest_table[i, 13:14], na.last = TRUE, decreasing = TRUE)
    biggest_table[i, 15:16] <- sort(biggest_table[i, 15:16])
    biggest_table[i, 17:18] <- sort(biggest_table[i, 17:18], decreasing = TRUE)
  }
  biggest_table[i, start_col:end_col] <- biggest_table[i, start_col:end_col] * focal_par_index
  if(biggest_table[i,6] > 1e10){
    biggest_table[i,] <- NA
  }
}
# remove failed optimizations
length(which(apply(biggest_table, 1, function(x) all(is.na(x)))))/dim(biggest_table)[1]
biggest_table <- biggest_table[!apply(biggest_table, 1, function(x) all(is.na(x))),]

# working with an error table
error_table <- biggest_table
# calculate RMSE
error_table[,12:18] <- t((t((biggest_table[,12:18]))) - (pars))
error_table <- error_table[,c(1:4,12:18)] # define what i want here (removing likelihood, AIC, etc.)

## plotting
plot_table <- melt(error_table, id = c("model_class","model_type","nTip","nMap"))
par_cols <- c("#6a3d9a", "#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c")
## by model type
# aggregate(error_table[-c(1:4)], by=list(error_table$model_type), FUN=median, na.rm=TRUE)
# a <- ggplot(plot_table, aes(x = model_type, y = value, fill = variable)) +
#   geom_boxplot(outlier.shape = NA) +
#   scale_fill_manual(values = par_cols) + 
#   ggtitle("a) Model Type") + 
#   ylab("Log10 diff from sim pars") + 
#   coord_cartesian(ylim=c(-2,2))
a <- ggplot(plot_table, aes(x = variable, y = value, fill = variable)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = par_cols) + 
  ggtitle("a) Model Type") + 
  ylab("Difference from generating parameters") + 
  facet_wrap(~model_class) +
  theme_classic() +
  geom_hline(yintercept = 0) +
  coord_cartesian(ylim=c(-3,3))

b <- ggplot(plot_table, aes(x = model_type, y = value, fill = variable)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = par_cols) +
  ggtitle("a) Model Type") + 
  ylab("") + 
  facet_wrap(~model_class) +
  theme_classic() +
  geom_hline(yintercept = 0) +
  coord_cartesian(ylim=c(-3,3))

c <- ggplot(plot_table, aes(x = nTip, y = value, fill = variable)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = par_cols) +
  ggtitle("b) Number of Taxa") + 
  ylab("Difference from generating parameters") + 
  facet_wrap(~model_class) +
  theme_classic() +
  geom_hline(yintercept = 0) +
  coord_cartesian(ylim=c(-3,3))

d <- ggplot(plot_table, aes(x = nMap, y = value, fill = variable)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = par_cols) +
  ggtitle("c) Number of maps per iteration") + 
  ylab("") + 
  facet_wrap(~model_class) +
  theme_classic() +
  geom_hline(yintercept = 0) +
  coord_cartesian(ylim=c(-3,3))

final_plot <- grid.arrange(b,c,d, nrow=3)
ggsave(filename = "figures/raw/model_error_rates.png", plot = final_plot, height = 8, width = 12, units = "in")

# ggplot(plot_table, aes(x = variable, y = value)) +
#   geom_boxplot() +
#   ylab("Raw diff from sim pars") + 
#   coord_cartesian(ylim=c(-5,5))

# PARAMATERS - RMSLE
error_table <- biggest_table
# calculate RMSE
error_table[,12:18] <- sqrt(t(((t(biggest_table[,12:18]))) - (pars))^2)
error_table <- error_table[,c(1:4,12:18)] # define what i want here (removing likelihood, AIC, etc.)

# by tip number
aggregate(error_table[-c(1:4)], by=list(error_table$nTip, error_table$model_class), FUN=mean, na.rm=TRUE)
# by map number
aggregate(error_table[-c(1:4)], by=list(error_table$nMap, error_table$model_class), FUN=mean, na.rm=TRUE)
# by model type
aggregate(error_table[-c(1:4)], by=list(error_table$model_type, error_table$model_class), FUN=mean, na.rm=TRUE)
# by model class
aggregate(error_table[-c(1:4)], by=list(error_table$model_class), FUN=mean, na.rm=TRUE)


rmse_table <- aggregate(error_table[-c(1:4)], by=list(error_table$model_type, error_table$model_clas), FUN=mean, na.rm=TRUE)
rmse_table[,-c(1,2)] <- round(rmse_table[,-c(1,2)], 2)
colnames(rmse_table)[1:2] <- c("model_type", "model_class")
write.csv(rmse_table, file = "tables/error-table.csv", row.names = FALSE)

# PARAMATERS - log errors
error_table <- biggest_table
# calculate RMSE
error_table[,12:18] <- t((t((biggest_table[,12:18]))) - (pars))
error_table <- error_table[,c(1:4,12:18)] # define what i want here (removing likelihood, AIC, etc.)
start_col <- which(colnames(error_table) == "rate")
end_col <- length(colnames(error_table))
# what percent is removed from setting a threshold of error (misestimattions)
# round(length(which(apply(error_table[,c(start_col:end_col)], 1, function(x) any(x[!is.na(x)] > 500))))/nrow(error_table) * 100, 2) # percent removed
# error_table <- error_table[!apply(error_table[,c(start_col:end_col)], 1, function(x) any(x[!is.na(x)] > 500)),]

# by tip number
aggregate(error_table[-c(1:4)], by=list(error_table$nTip, error_table$model_class), FUN=mean, na.rm=TRUE)
# by map number
aggregate(error_table[-c(1:4)], by=list(error_table$nMap, error_table$model_class), FUN=mean, na.rm=TRUE)
# by model type
aggregate(error_table[-c(1:4)], by=list(error_table$model_type, error_table$model_class), FUN=mean, na.rm=TRUE)
# by model class
aggregate(error_table[-c(1:4)], by=list(error_table$model_class), FUN=median, na.rm=TRUE)


# MODEL COMPARISON
error_table <- biggest_table
# calculate RMSE
error_table <- error_table[,c(1:4,7:11)] # define what i want here (removing likelihood, AIC, etc.)

# by tip number
aggregate(error_table[-c(1:4)], by=list(error_table$nTip, error_table$model_class), FUN=mean, na.rm=TRUE)
# by map number
aggregate(error_table[-c(1:4)], by=list(error_table$nMap, error_table$model_class), FUN=mean, na.rm=TRUE)
# by model type
aggregate(error_table[-c(1:4)], by=list(error_table$model_type, error_table$model_class), FUN=mean, na.rm=TRUE)
# by model class
aggregate(error_table[-c(1:4)], by=list(error_table$model_class), FUN=mean, na.rm=TRUE)

model_table <- aggregate(error_table[-c(1:4)], by=list(error_table$model_type, error_table$model_class), FUN=mean, na.rm=TRUE)
model_table[,-c(1,2)] <- round(model_table[,-c(1,2)], 2)
colnames(model_table)[1:2] <- c("model_type", "model_class")
write.csv(model_table, file = "tables/model-table.csv", row.names = FALSE)

tip_table <- aggregate(error_table[-c(1:4)], by=list(error_table$nTip, error_table$model_class), FUN=mean, na.rm=TRUE)
tip_table[,-c(1,2)] <- round(tip_table[,-c(1,2)], 2)
colnames(tip_table)[1:2] <- c("model_type", "model_class")
write.csv(tip_table, file = "tables/model-tip-table.csv", row.names = FALSE)


## plotting
error_table <- biggest_table
# calculate RMSE
error_table <- error_table[,c(1:4,7:10)] # define what i want here (removing likelihood, AIC, etc.)
plot_table <- melt(error_table, id = c("model_class","model_type","nTip","nMap"))
par_cols <- c("#6a3d9a", "#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c")
## by model type
# aggregate(error_table[-c(1:4)], by=list(error_table$model_type), FUN=median, na.rm=TRUE)
# a <- ggplot(plot_table, aes(x = model_type, y = value, fill = variable)) +
#   geom_boxplot(outlier.shape = NA) +
#   scale_fill_manual(values = par_cols) + 
#   ggtitle("a) Model Type") + 
#   ylab("Log10 diff from sim pars") + 
#   coord_cartesian(ylim=c(-2,2))

model_weight_plot <- ggplot(plot_table, aes(x = variable, y = value, fill = variable)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = par_cols) +
  ylab("AICwt") + 
  facet_wrap(~model_class + model_type) +
  theme_classic() +
  geom_hline(yintercept = 0) +
  coord_cartesian(ylim=c(0,1))

ggsave(filename = "figures/raw/model_weight.png", plot = model_weight_plot, height = 8, width = 12, units = "in")
