set.seed(1985)
setwd("~/2020_houwie/")
source("hOUwieNode.R")

require(ape) 
require(magick) 
require(expm) 
require(POUMM)
require(ggplot2)
require(reshape2)
require(corHMM)
require(gridExtra)

# functions

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

generateParameters <- function(continuous_model, alpha, sigma.sq, theta, discrete_model, rate){
  par_alpha <- alpha[continuous_model[1,]]
  par_alpha[is.na(par_alpha)] <- 1e-10
  par_sigma <- sigma.sq[continuous_model[2,] - min(continuous_model[2,]) + 1]
  par_sigma[is.na(par_sigma)] <- 1e-10
  par_theta <- theta[continuous_model[3,] - min(continuous_model[3,]) + 1]
  par_theta[is.na(par_theta)] <- 1e-10
  k.discrete <- max(discrete_model, na.rm=TRUE)
  par_rates <- rate[seq_len(k.discrete)]
  discrete_model[discrete_model == 0] <- NA
  discrete_model[is.na(discrete_model)] <- max(discrete_model, na.rm = TRUE) + 1
  Q <- matrix(0, dim(discrete_model)[1], dim(discrete_model)[1])
  Q[] <- c(par_rates, 0)[discrete_model]
  diag(Q) <- -rowSums(Q)
  p <- list(alpha = par_alpha, sigma.sq = par_sigma, theta = par_theta, Q = Q)
  return(p)
}


# continuous models
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
cd_model_structures <- all_model_structures[grep("CD_", model_names)]
# discrete models
discrete_model_cd <- equateStateMatPars(getRateCatMat(2), 1:2)
discrete_model_cid <- getFullMat(list(discrete_model_cd, discrete_model_cd), getRateCatMat(2))

# parameters
alpha <- c(3, 1.5)
sigma.sq <- c(0.35, 1)
theta <- c(2, 0.75)
root.p <- c(1,0,0,0)
theta0 <- theta[1]
rate <- .1
time.unit <- 0.1
time.steps <- 100

names(cd_model_structures)
cols <- c("#d53e4f", "#3288bd")
all_plots <- list()
odd_numbers <- seq(from = 1, to = 10, by = 2)
for(j in 1:10){
  pars <- generateParameters(cd_model_structures[[j]], alpha, sigma.sq, theta, discrete_model_cd, rate)
  model_type <- paste0(letters[j], ") ", gsub("CD_", "", names(cd_model_structures)[j]))
  if(j %in% odd_numbers){
    y_label <- "Continuous trait value"
  }else{
    y_label <- ""
  }
  if(j %in% c(9, 10)){
    x_label_a <- "Time"
    x_label_b <- "Count"
  }else{
    x_label_a <- ""
    x_label_b <- ""
  }
  # df <- cbind(X=1:time.steps, forwardSimBranch(pars$theta[1], 1, pars$Q, pars$alpha, pars$sigma.sq, pars$theta, time.unit, time.steps))
  df_list <- lapply(1:20, function(x) forwardSimBranch(pars$theta[1], 1, pars$Q, pars$alpha, pars$sigma.sq, pars$theta, time.unit, time.steps))
  for(i in 1:length(df_list)){
    tmp_line <- cbind(I=i, X=1:time.steps, df_list[[i]])
    tmp_line$Tr <- 0
    tmp_line$Tr[1] <- 0
    for(k in 2:length(tmp_line$D)){
      if(tmp_line$D[k] == tmp_line$D[k-1]){
        tmp_line$Tr[k] <- 0
      }else{
        tmp_line$Tr[k] <- ifelse(paste0(tmp_line$D[c(k-1, k)], collapse = "") == "12", 1, 2)
      }
    }
    df_list[[i]] <- tmp_line
  }
  plot_line <- do.call(rbind, df_list)
  plot_line$D <- as.factor(plot_line$D)
  plot_line$I <- as.factor(plot_line$I)
  plot_point <- plot_line[plot_line$Tr == "1" | plot_line$Tr == "2",]
  focal_line <- sample(which(as.vector(summary(plot_point$I)) >= 2), 1)
  plot_focal_line <- plot_line[plot_line$I == focal_line,]
  plot_focal_point <- plot_point[plot_point$I == focal_line,]
  
  a <- ggplot(data=plot_line, aes(x=X, y = C, group = I)) +
    geom_line(alpha = 0.25, color = cols[plot_line$D]) + 
    geom_point(data = plot_point,aes(x = X, y = C), color = "black", fill=cols[c(2,1)][plot_point$Tr], size=2, alpha=0.25, shape = c(21, 22)[plot_point$Tr]) +
    # geom_point(data = plot_point, aes(x = X, y = C), color = "black", fill=cols[c(2,1)][plot_point$Tr], size=1, alpha=0.25,  shape = 25) +
    geom_line(data = plot_focal_line, alpha = 1, color = cols[plot_focal_line$D]) +
    geom_point(data = plot_focal_point,aes(x = X, y = C), color = "black", fill=cols[c(2,1)][plot_focal_point$Tr], size=2, alpha=1, shape = c(21, 22)[plot_focal_point$Tr]) +
    ggtitle(model_type) +
    labs(y = y_label, x = x_label_a) +
    # geom_point(alpha = 0.05, size = 0.1, color = plot_line$D) +
    theme_bw()
  ylims <- layer_scales(a)$y$range$range
  
  plot_hist <- do.call(rbind, lapply(df_list, function(x) x[time.steps,]))
  plot_hist$D <- as.factor(plot_hist$D)
  b <- ggplot(plot_hist, aes(x=C, fill=D)) +
    geom_histogram(alpha=0.7, position="identity", aes(y = ..density..), color="black") +
    geom_density(alpha=0.7) +
    scale_fill_manual(values = cols) + 
    labs(y = x_label_b, x = "") +
    ggtitle("") +
    coord_flip() +
    xlim(ylims) +
    theme(legend.position = "none")
  
  all_plots[[j]] <- grid.arrange(a, b, nrow=1, widths=c(2,1))
}

final_plot <- grid.arrange(all_plots[[1]], all_plots[[2]], all_plots[[3]], all_plots[[4]], all_plots[[5]],
                           all_plots[[6]], all_plots[[7]], all_plots[[8]], all_plots[[9]], all_plots[[10]], nrow = 5)

ggsave(filename = "figures/raw/simulated_hOUwie_models.jpg", plot = final_plot, width = 15, height = 15, units = "in")





