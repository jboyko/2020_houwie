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
  cd_out <- c(k=model_table_cd$np[3], loglik=model_table_cd$lnLik[3], AICwt=model_table_cd$AICwt[3],best_model=model_table_cd$dAIC[3]==0, pars_cd)
  
  # summarize character independent models
  cid_sim <- out$cid_out
  model_table_cid <- getModelTable(cid_sim)
  pars_cid <- c(cid_sim$cid_fit$p[1], t(cid_sim$cid_fit$solution.cont[,c(1,3)]))
  names(pars_cid) <- c("rate", "alpha_1", "alpha_2", "sigma.sq_1", "sigma.sq_2", "theta_1", "theta_2")
  if(mean(cid_sim$bm1_fit$data[,3]) > 25){
    pars_cid[c(6,7)] <- pars_cid[c(6,7)] - 50
  }
  cid_out <- c(k=model_table_cid$np[4], loglik=model_table_cid$lnLik[4], AICwt=model_table_cid$AICwt[4], best_model=model_table_cid$dAIC[4]==0, pars_cid)
  
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
  error <- focal_pars - as.numeric(table_row[9:15])
  return(error)
}


#### #### #### #### #### #### #### #### #### #### #### #### 
# run
#### #### #### #### #### #### #### #### #### #### #### #### 
# for each model type i want to generate a dataset consistent with the CD and one consistent with CID+
model_types <- c("BMV", "OUV", "OUA", "OUM", "OUVA", "OUMV", "OUMA", "OUMVA", "OUBM1", "OUBMV")

ntips <- c(25, 100, 250)
nmaps <- c(25, 100)

# model_type <- "OUA"
# ntip <- 25
# nmap <- 25
# 

bigger_table <- list()
count <- 1
for(k in 1:length(nmaps)){
  for(j in 1:length(ntips)){
    for(i in 1:length(model_types)){
      print(count)
      nmap <- nmaps[k]
      ntip <- ntips[j]
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
      count <- count + 1
    }
  }
}

biggest_table <- do.call(rbind, bigger_table)
rownames(biggest_table) <- NULL

load("biggest_table.Rsave")
biggest_table <- as.data.frame(biggest_table)
biggest_table$model_type <- as.factor(biggest_table$model_type)
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
  focal_par_index <- par_index_table[,biggest_table[i,1] == colnames(par_index_table)]
  start_col <- which(names(biggest_table[i,]) == "rate")
  end_col <- length(names(biggest_table[i,]))
  biggest_table[i, start_col:end_col] <- biggest_table[i, start_col:end_col] * focal_par_index
}

# working with an error table
error_table <- biggest_table
# calculate RMSE
error_table[,9:15] <- (biggest_table[,9:15] - pars)
error_table <- error_table[,c(1:4,9:15)] # define what i want here (removing likelihood, AIC, etc.)
start_col <- which(colnames(error_table) == "rate")
end_col <- length(colnames(error_table))
# what percent is removed from setting a threshold of error (misestimattions)
round(length(which(apply(error_table[,c(start_col:end_col)], 1, function(x) any(x[!is.na(x)] > 500))))/nrow(error_table) * 100, 2) # percent removed
error_table <- error_table[!apply(error_table[,c(start_col:end_col)], 1, function(x) any(x[!is.na(x)] > 500)),]
# aggregate the table into the broadest categories
rmse_table <-  aggregate(error_table[-c(1:4)], by=list(error_table$model_class, error_table$model_type, error_table$nTip, error_table$nMap), FUN=mean, na.rm=TRUE)

colnames(rmse_table)[1:4] <- c("model_class","model_type","nTip","nMap")
rmse_table[rmse_table$model_type == "OUMVA",]
rmse_table[rmse_table$model_type == "OUM",]


## plotting
plot_table <- melt(error_table, id = c("model_class","model_type","nTip","nMap"))
par_cols <- c("#6a3d9a", "#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c")
## by model type
# aggregate(error_table[-c(1:4)], by=list(error_table$model_type), FUN=median, na.rm=TRUE)
a <- ggplot(plot_table, aes(x = model_type, y = value, fill = variable)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = par_cols) + 
  ggtitle("a) Model Type") + 
  ylab("Raw diff from sim pars") + 
  coord_cartesian(ylim=c(-5,5))

b <- ggplot(plot_table, aes(x = model_class, y = value, fill = variable)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = par_cols) +
  ggtitle("b) Model Class") + 
  ylab("Raw diff from sim pars") + 
  coord_cartesian(ylim=c(-5,5))

c <- ggplot(plot_table, aes(x = nTip, y = value, fill = variable)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = par_cols) +
  ggtitle("c) Number of Taxa") + 
  ylab("Raw diff from sim pars") + 
  coord_cartesian(ylim=c(-5,5))

d <- ggplot(plot_table, aes(x = nMap, y = value, fill = variable)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = par_cols) +
  ggtitle("d) Number of maps per iteration") + 
  ylab("Raw diff from sim pars") + 
  coord_cartesian(ylim=c(-5,5))

grid.arrange(a,b,c,d, nrow=2)





