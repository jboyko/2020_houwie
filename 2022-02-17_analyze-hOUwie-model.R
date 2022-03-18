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
  error <- as.numeric(table_row[9:15]) - focal_pars
  return(error)
}


#### #### #### #### #### #### #### #### #### #### #### #### 
# run
#### #### #### #### #### #### #### #### #### #### #### #### 
# for each model type i want to generate a dataset consistent with the CD and one consistent with CID+
# model_types <- c("BMV", "OUV", "OUA", "OUM", "OUVA", "OUMV", "OUMA", "OUMVA", "OUBM1", "OUBMV")
# 
# ntips <- c(25, 100, 250)
# nmaps <- c(25, 100)
# 
# # model_type <- "OUA"
# # ntip <- 25
# # nmap <- 25
# # 
# 
# bigger_table <- list()
# count <- 1
# for(k in 1:length(nmaps)){
#   for(j in 1:length(ntips)){
#     for(i in 1:length(model_types)){
#       print(count)
#       nmap <- nmaps[k]
#       ntip <- ntips[j]
#       model_type <- model_types[i]
#       focal_files <- getFocalFiles(model_type, ntip, nmap)
#       if(length(focal_files) == 0){
#         next
#       }
#       list_of_errors <- lapply(focal_files, summarizeFile)
#       cd_errors <- do.call(rbind, lapply(list_of_errors, "[[", "cd_out"))
#       cd_errors <- cbind(model_type = model_type, model_class="CD", nTip=ntip, nMap=nmap, as.data.frame(cd_errors))
#       cid_errors <- do.call(rbind, lapply(list_of_errors, "[[", "cid_out"))
#       cid_errors <- cbind(model_type = model_type, model_class="CID", nTip=ntip, nMap=nmap, as.data.frame(cid_errors))
#       big_table <- rbind(cd_errors, cid_errors)
#       bigger_table[[count]] <- big_table
#       count <- count + 1
#     }
#   }
# }
# 
# biggest_table <- do.call(rbind, bigger_table)
# rownames(biggest_table) <- NULL

# start here

load("biggest_table.Rsave")
biggest_table <- as.data.frame(biggest_table)
biggest_table$model_type <- as.factor(biggest_table$model_type)
biggest_table$model_class <- as.factor(biggest_table$model_class)
biggest_table$nTip <- as.factor(biggest_table$nTip)
biggest_table$nMap <- as.factor(biggest_table$nMap)

model_types <- c("BMV", "OUV", "OUA", "OUM", "OUVA", "OUMV", "OUMA", "OUMVA", "OUBM1", "OUBMV")

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
  biggest_table[i, start_col:end_col] <- biggest_table[i, start_col:end_col] * focal_par_index
}

# working with an error table
error_table <- biggest_table
# calculate RMSE
error_table[,9:15] <- t(log10(t(biggest_table[,9:15])) - log10(pars))
error_table <- error_table[,c(1:4,9:15)] # define what i want here (removing likelihood, AIC, etc.)
start_col <- which(colnames(error_table) == "rate")
end_col <- length(colnames(error_table))
# what percent is removed from setting a threshold of error (misestimattions)
round(length(which(apply(error_table[,c(start_col:end_col)], 1, function(x) any(x[!is.na(x)] > 500))))/nrow(error_table) * 100, 2) # percent removed
error_table <- error_table[!apply(error_table[,c(start_col:end_col)], 1, function(x) any(x[!is.na(x)] > 500)),]
# aggregate the table into the broadest categories
# rmse_table <-  aggregate(error_table[-c(1:4)], by=list(error_table$model_class, error_table$model_type, error_table$nTip, error_table$nMap), FUN=mean, na.rm=TRUE)
# 
# colnames(rmse_table)[1:4] <- c("model_class","model_type","nTip","nMap")
# rmse_table[rmse_table$model_type == "OUMVA",]
# rmse_table[rmse_table$model_type == "OUM",]
# 
# rmse_table

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
  ylab("Log10 diff from sim pars") + 
  coord_cartesian(ylim=c(-2,2))

b <- ggplot(plot_table, aes(x = model_class, y = value, fill = variable)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = par_cols) +
  ggtitle("b) Model Class") + 
  ylab("Log10 diff from sim pars") + 
  coord_cartesian(ylim=c(-2,2))

c <- ggplot(plot_table, aes(x = nTip, y = value, fill = variable)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = par_cols) +
  ggtitle("c) Number of Taxa") + 
  ylab("Log10 diff from sim pars") + 
  coord_cartesian(ylim=c(-2,2))

d <- ggplot(plot_table, aes(x = nMap, y = value, fill = variable)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = par_cols) +
  ggtitle("d) Number of maps per iteration") + 
  ylab("Log10 diff from sim pars") + 
  coord_cartesian(ylim=c(-2,2))

grid.arrange(a,b,c,d, nrow=2)

# ggplot(plot_table, aes(x = variable, y = value)) +
#   geom_boxplot() +
#   ylab("Raw diff from sim pars") + 
#   coord_cartesian(ylim=c(-5,5))

# table by model type and class
rmse_table_type <-  aggregate(error_table[-c(1:4)], by=list(error_table$model_class, error_table$model_type), FUN=mean, na.rm=TRUE)
rmse_table_type[,-c(1,2)] <- round(rmse_table_type[,-c(1,2)], 4)

# sign errors
sign_error_table <- data.frame(model_type = biggest_table$model_type,
                               model_class = biggest_table$model_class,
                               nTip = biggest_table$nTip,
                               nMap = biggest_table$nMap,
                               alpha_error = abs(biggest_table$alpha_1 - biggest_table$alpha_2) - 1.5,
                               sigma.sq_error = abs(biggest_table$sigma.sq_2 - biggest_table$sigma.sq_1) - 0.65,
                               theta_error = abs(biggest_table$theta_1 - biggest_table$theta_2) - 1.25)

sign_plot_table <- melt(sign_error_table, id = c("model_class","model_type","nTip","nMap"))
# sign_error_table[,c(5,6,7)] <- abs(sign_error_table[,c(5,6,7)])
# sign_error_table[,c(5,6,7)] <- sign_error_table[,c(5,6,7)] + abs(min(sign_error_table[,c(5,6,7)], na.rm = TRUE)) * 1.1
# sign_plot_table <- sign_plot_table[!is.na(sign_plot_table$value),]
# sign_plot_table$value <- factor(sign_plot_table$value, levels = c(FALSE, TRUE))

# ggplot(data = sign_error_table, aes(x = alpha_error, y = sigma.sq_error))+
#   geom_point() +
#   theme_classic() +
#   scale_x_continuous(trans = 'log10') +
#   scale_y_continuous(trans = 'log10')
sign_error_table_type <-  aggregate(sign_error_table[-c(1:4)], by=list(sign_error_table$model_class, sign_error_table$model_type), FUN=mean, na.rm=TRUE)
sign_error_table_type[,-c(1,2)] <- round(sign_error_table_type[,-c(1,2)], 4)



a <- ggplot(sign_plot_table, aes(x = model_type, y = value, fill = variable)) +
  geom_boxplot(outlier.shape = NA) +
  ggtitle("a) Model Type") + 
  scale_fill_brewer() + 
  theme_classic() + 
  geom_hline(yintercept = 0) + 
  coord_cartesian(ylim=c(-5,5))

b <- ggplot(sign_plot_table, aes(x = model_class, y = value, fill = variable)) +
  geom_boxplot(outlier.shape = NA) +
  ggtitle("b) Model Class") + 
  scale_fill_brewer() + 
  theme_classic() + 
  geom_hline(yintercept = 0) + 
  coord_cartesian(ylim=c(-5,5))

c <- ggplot(sign_plot_table, aes(x = nTip, y = value, fill = variable)) +
  geom_boxplot(outlier.shape = NA) +
  ggtitle("c) Number of Taxa") +   
  scale_fill_brewer() + 
  theme_classic() + 
  geom_hline(yintercept = 0) + 
  coord_cartesian(ylim=c(-5,5))

d <- ggplot(sign_plot_table, aes(x = nMap, y = value, fill = variable)) +
  geom_boxplot(outlier.shape = NA) +
  ggtitle("d) Number of maps per iteration") + 
  scale_fill_brewer() + 
  theme_classic() + 
  geom_hline(yintercept = 0) + 
  coord_cartesian(ylim=c(-5,5))

grid.arrange(a,b,c,d, nrow=2)


load("biggest_table.Rsave")

model_types <- c("BMV", "OUV", "OUA", "OUM", "OUVA", "OUMV", "OUMA", "OUMVA", "OUBM1", "OUBMV")
s_vars <- c(NA, 
  abs(diff(sigma.sq/(2*alpha[1]))), 
  abs(diff(sigma.sq[1]/(2*alpha))),
  NA,
  abs(diff(sigma.sq/(2*alpha))), 
  abs(diff(sigma.sq/(2*alpha[1]))),
  abs(diff(sigma.sq[1]/(2*alpha))),
  abs(diff(sigma.sq/(2*alpha))),
  NA,
  NA)
names(s_vars) <- model_types

v_1 <- biggest_table$sigma.sq_1/(2*biggest_table$alpha_1)
v_2 <- biggest_table$sigma.sq_2/(2*biggest_table$alpha_2)

s_var_table <- cbind(biggest_table[,1:4], abs(v_1 - v_2))
colnames(s_var_table)[5] <- c("stat_var")
s_var_table[,5] <- abs(s_var_table[,5])

s_var_table <- s_var_table[!s_var_table$model_type =="BMV",]
s_var_table <- s_var_table[!s_var_table$model_type =="OUM",]
s_var_table <- s_var_table[!s_var_table$model_type =="OUBM1",]
s_var_table <- s_var_table[!s_var_table$model_type =="OUBMV",]

for(i in 1:length(model_types)){
  s_var_table[s_var_table$model_type ==model_types[i],5] <- s_var_table[s_var_table$model_type ==model_types[i],5] - s_vars[i]
}

s_var_table$model_type <- as.factor(s_var_table$model_type)
s_var_table$model_class <- as.factor(s_var_table$model_class)
s_var_table$nTip <- as.factor(s_var_table$nTip)
s_var_table$nMap <- as.factor(s_var_table$nMap)

a <- ggplot(s_var_table, aes(x = model_type, y = stat_var)) +
  geom_boxplot(outlier.shape = NA) +
  ggtitle("a) Model Type") + 
  scale_fill_brewer() + 
  theme_classic() + 
  geom_hline(yintercept = 0) + 
  coord_cartesian(ylim=c(-5,5))

b <- ggplot(s_var_table, aes(x = model_class, y = stat_var)) +
  geom_boxplot(outlier.shape = NA) +
  ggtitle("b) Model Class") + 
  scale_fill_brewer() + 
  theme_classic() + 
  geom_hline(yintercept = 0) + 
  coord_cartesian(ylim=c(-5,5))

c <- ggplot(s_var_table, aes(x = nTip, y = stat_var)) +
  geom_boxplot(outlier.shape = NA) +
  ggtitle("c) Number of Taxa") +   
  scale_fill_brewer() + 
  theme_classic() + 
  geom_hline(yintercept = 0) + 
  coord_cartesian(ylim=c(-5,5)) +
  facet_wrap(~model_class)

d <- ggplot(s_var_table, aes(x = nMap, y = stat_var)) +
  geom_boxplot(outlier.shape = NA) +
  ggtitle("d) Number of maps per iteration") + 
  scale_fill_brewer() + 
  theme_classic() + 
  geom_hline(yintercept = 0) + 
  coord_cartesian(ylim=c(-5,5)) +
  facet_wrap(~model_class)

grid.arrange(a,b,c,d, nrow=2)






