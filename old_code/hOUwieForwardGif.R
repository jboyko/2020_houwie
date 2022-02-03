require(ape)
require(magick)
require(expm)
require(POUMM)

# forward simulation of a branch. kind of a cheap way to do discrete time stuff by just shrinking time, but the markov property should ensure that this works just fine
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

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### GIF 1 of presentation - the earliest version not used in the end
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
# creates a forward simulation gif for hOUwie
# intial conditions
theta0 = 0
state0 = 1
# model parameters
Q = rbind(c(-1,1), c(1, -1))
alpha = c(0, 0)
sigma = c(2, 0.5)
theta = c(10, -10)
# resolution parameters
time.unit <- 0.1
time.steps <- 100
# branching information
BranchingDF <- data.frame(time.step = c(0, 50), n = c(2,3))

df.ForwardBranch <- forwardSimBranch(theta0, state0, Q, alpha, sigma, theta, time.unit, time.steps)
dir_out <- "~/presentation/2021_Evolution/gif_1/"
cols <- c("orange", "purple")

for(j in 2:time.steps){
  # png(file = paste0(dir_out, "gif_1-", j-1, ".png"), width = 1500, height = 500)
  plot(1, type = "n", xlab = "", ylab = "", xlim = c(0, time.steps), ylim = c(-20, 20), axes=FALSE)
  for(i in 2:j){
    df_i <- df.ForwardBranch[i:(i-1),]
    # discrete state
    rect(rownames(df_i)[1], 15, rownames(df_i)[2], 10, col = cols[df_i$D[1]], border = NA)
    # # # discrete theta
    # lines_i <- list(x = mean(as.numeric(rownames(df_i))), y = (theta[df_i$D][1] * 2))
    # lines(lines_i, col = cols[df_i$D[1]], type = "h", lwd = 2)
    # continuous trait a
    lines_i <- list(x = rownames(df_i), y = df_i$C - 7.5)
    lines(lines_i, col = cols[df_i$D[1]], lwd = 4)
  }
  # dev.off()
  Sys.sleep(0.1)
}

## list file names and read in
imgs <- list.files(dir_out, full.names = TRUE)
imgs <- imgs[sapply(1:99, function(x) grep(paste0("-", x, ".png"), imgs))]
img_list <- lapply(imgs, image_read)

## join the images together
img_joined <- image_join(img_list)

## animate at x frames per second
img_animated <- image_animate(img_joined, fps = 10)

## save to disk
image_write(image = img_animated,
            path = "~/presentation/2021_Evolution/gif_1.gif")


#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### IMG 1 of presentation
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
theta0 = 0
state0 = 1
# model parameters
Q = rbind(c(-1,1), c(1, -1))
alpha = c(0, 0)
sigma = c(1, 1)
theta = c(10, -10)
# resolution parameters
time.unit <- 0.1
time.steps <- 100
# branching information
BranchingDF <- data.frame(time.step = c(0, 50), n = c(2,3))

df.ForwardBranch <- lapply(1:1000, function(x) forwardSimBranch(theta0, state0, Q, alpha, sigma, theta, time.unit, time.steps))

tmp <- df.ForwardBranch[[1]]
png("~/presentation/2021_Evolution/IMG1_A.png", width = 15.69, height = 5, units = "in", res = 96)
plot(x = rownames(tmp), y = tmp$C, xlab = "", ylab = "", xlim = c(0, time.steps), axes=FALSE, type = "l", ylim = c(-10, 10), col = rgb(0, 0, 0, 0.15), )
for(i in 2:length(df.ForwardBranch)){
  tmp <- df.ForwardBranch[[i]]
  lines_i <- list(x = rownames(tmp), y = tmp$C)
  lines(lines_i, col = rgb(0, 0, 0, 0.15))
}
dev.off()
lines_i <- list(x = rownames(tmp), y = colMeans(do.call(rbind, lapply(df.ForwardBranch, function(x) x$C))))
lines(lines_i, col = rgb(1, 0, 0), lwd = 2)

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### IMG 2 of presentation
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
theta0 = -8
state0 = 1
# model parameters
Q = rbind(c(-1,1), c(1, -1))
alpha = c(0.5, 0.5)
sigma = c(2, 2)
theta = c(5, 5)
# resolution parameters
time.unit <- 0.1
time.steps <- 100

df.ForwardBranch <- lapply(1:1000, function(x) forwardSimBranch(theta0, state0, Q, alpha, sigma, theta, time.unit, time.steps))

tmp <- df.ForwardBranch[[1]]
png("~/presentation/2021_Evolution/IMG2_A.png", width = 15.69, height = 5, units = "in", res = 96)
plot(x = rownames(tmp), y = tmp$C, xlab = "", ylab = "", xlim = c(0, time.steps), axes=FALSE, type = "l", ylim = c(-10, 12), col = rgb(0, 0, 0, 0.15), )
for(i in 2:length(df.ForwardBranch)){
  tmp <- df.ForwardBranch[[i]]
  lines_i <- list(x = rownames(tmp), y = tmp$C)
  lines(lines_i, col = rgb(0, 0, 0, 0.15))
}
dev.off()
lines_i <- list(x = rownames(tmp), y = colMeans(do.call(rbind, lapply(df.ForwardBranch, function(x) x$C))))
lines(lines_i, col = rgb(1, 0, 0), lwd = 2)

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### GIF 2 of presentation
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
# creates a forward simulation gif for hOUwie
# intial conditions
theta0 = 0
state0 = 1
# model parameters
Q = rbind(c(-1,1), c(1, -1))/2
alpha = c(2, 1)
sigma = c(2, 2)
theta = c(10, -10)
# resolution parameters
time.unit <- 0.1
time.steps <- 200

df.ForwardBranch <- forwardSimBranch(theta0, state0, Q, alpha, sigma, theta, time.unit, time.steps)
dir_out <- "~/presentation/2021_Evolution/gif_2/"
cols <- c("orange", "purple")

# uncomment the png and dev.off to produce images that can then be read in
for(j in 2:time.steps){
  # png(file = paste0(dir_out, "gif_1-", j-1, ".png"), width = 15.69, height = 5, units = "in", res = 96)
  plot(1, type = "n", xlab = "", ylab = "", xlim = c(0, time.steps), ylim = c(-25, 20), axes=FALSE)
  for(i in 2:j){
    df_i <- df.ForwardBranch[i:(i-1),]
    # discrete state
    rect(rownames(df_i)[1], 20, rownames(df_i)[2], 15, col = cols[df_i$D[1]], border = NA)
    # # # discrete theta
    lines_i <- list(x = rownames(df_i), y = theta[df_i$D] - 10)
    lines(lines_i, col = "light grey", lwd = 3)
    # continuous trait a
    lines_i <- list(x = rownames(df_i), y = df_i$C - 10)
    lines(lines_i, col = cols[df_i$D[1]], lwd = 4)
  }
  # dev.off()
  Sys.sleep(0.1)
}

## list file names and read in
imgs <- list.files(dir_out, full.names = TRUE)
imgs <- imgs[sapply(1:(time.steps-1), function(x) grep(paste0("-", x, ".png"), imgs))]
img_list <- lapply(imgs, image_read)

## join the images together
img_joined <- image_join(img_list)

## animate at x frames per second
img_animated <- image_animate(img_joined, fps = 10)

## save to disk
image_write(image = img_animated,
            path = "~/presentation/2021_Evolution/gif_2.gif")

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### IMG 3 of presentation
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
theta0 = -8
state0 = 1
# model parameters
Q = rbind(c(-1,1), c(1, -1))
alpha = c(2, 2)
sigma = c(2, 2)
theta = c(5, 5)
# resolution parameters
time.unit <- 0.1
time.steps <- 100

df.ForwardBranch <- lapply(1:1000, function(x) forwardSimBranch(theta0, state0, Q, alpha, sigma, theta, time.unit, time.steps))

tmp <- df.ForwardBranch[[1]]
png("~/presentation/2021_Evolution/IMG3.png", width = 15.69, height = 5, units = "in", res = 96)
plot(x = rownames(tmp), y = tmp$C, xlab = "", ylab = "", xlim = c(0, time.steps), axes=FALSE, type = "l", ylim = c(-10, 12), col = rgb(0, 0, 0, 0.15), )
for(i in 2:length(df.ForwardBranch)){
  tmp <- df.ForwardBranch[[i]]
  lines_i <- list(x = rownames(tmp), y = tmp$C)
  lines(lines_i, col = rgb(0, 0, 0, 0.15))
}
dev.off()
lines_i <- list(x = rownames(tmp), y = colMeans(do.call(rbind, lapply(df.ForwardBranch, function(x) x$C))))
lines(lines_i, col = rgb(1, 0, 0), lwd = 2)


#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### IMG 4 of presentation
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
theta0 = -8
state0 = 1
# model parameters
Q = rbind(c(-1,1), c(1, -1))
alpha = c(0.5, 0.5)
sigma = c(3.5, 3.5)
theta = c(5, 5)
# resolution parameters
time.unit <- 0.1
time.steps <- 100

df.ForwardBranch <- lapply(1:1000, function(x) forwardSimBranch(theta0, state0, Q, alpha, sigma, theta, time.unit, time.steps))

tmp <- df.ForwardBranch[[1]]
png("~/presentation/2021_Evolution/IMG4.png", width = 15.69, height = 5, units = "in", res = 96)
plot(x = rownames(tmp), y = tmp$C, xlab = "", ylab = "", xlim = c(0, time.steps), axes=FALSE, type = "l", ylim = c(-10, 12), col = rgb(0, 0, 0, 0.15), )
for(i in 2:length(df.ForwardBranch)){
  tmp <- df.ForwardBranch[[i]]
  lines_i <- list(x = rownames(tmp), y = tmp$C)
  lines(lines_i, col = rgb(0, 0, 0, 0.15))
}
dev.off()
lines_i <- list(x = rownames(tmp), y = colMeans(do.call(rbind, lapply(df.ForwardBranch, function(x) x$C))))
lines(lines_i, col = rgb(1, 0, 0), lwd = 2)


#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### IMG 5 of presentation
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
theta0 = 8
state0 = 1
# model parameters
Q = rbind(c(-1,1), c(1, -1))
alpha = c(0.75, 0.75)
sigma = c(2, 2)
theta = c(-5, -5)
# resolution parameters
time.unit <- 0.1
time.steps <- 100

df.ForwardBranch <- lapply(1:1000, function(x) forwardSimBranch(theta0, state0, Q, alpha, sigma, theta, time.unit, time.steps))

tmp <- df.ForwardBranch[[1]]
png("~/presentation/2021_Evolution/IMG5.png", width = 15.69, height = 5, units = "in", res = 96)
plot(x = rownames(tmp), y = tmp$C, xlab = "", ylab = "", xlim = c(0, time.steps), axes=FALSE, type = "l", ylim = c(-10, 12), col = rgb(0, 0, 0, 0.15), )
for(i in 2:length(df.ForwardBranch)){
  tmp <- df.ForwardBranch[[i]]
  lines_i <- list(x = rownames(tmp), y = tmp$C)
  lines(lines_i, col = rgb(0, 0, 0, 0.15))
}
dev.off()
lines_i <- list(x = rownames(tmp), y = colMeans(do.call(rbind, lapply(df.ForwardBranch, function(x) x$C))))
lines(lines_i, col = rgb(1, 0, 0), lwd = 2)




# 
# 
# #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
# #### GIF 2 of presentation
# #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
# 
# nInitSpecies <- BranchingDF[BranchingDF$time.step==0,2]
# df.ForwardBranch <- lapply(1:nInitSpecies, function(x) forwardSimBranch(theta0, state0, Q, alpha, sigma, theta, time.unit, time.steps))
# initSplitStates <- df.ForwardBranch[[2]][50,]
# tmp <- df.ForwardBranch[[2]]
# tmp[51:100,] <- forwardSimBranch(initSplitStates$C, initSplitStates$D, Q, alpha, sigma, theta, time.unit, 50)
# df.ForwardBranch <- c(df.ForwardBranch, list(tmp))
# #df.ForwardBranch <- forwardSimBranch(theta0, state0, Q, alpha, sigma, theta, time.unit, time.steps)
# 
# 
# dir_out <- "~/presentation/2021_Evolution/gif_2/"
# cols <- c("orange", "purple")
# 
# 
# for(j in 2:time.steps){
#   png(file = paste0(dir_out, "gif_2-", j-1, ".png"), width = 1450, height = 660)
#   plot(1, type = "n", xlab = "", ylab = "", xlim = c(0, time.steps), ylim = c(-10, 10), axes=FALSE)
#   for(i in 2:j){
#     # continuous trait a
#     df_i <- df.ForwardBranch[[1]][i:(i-1),]
#     lines_i <- list(x = rownames(df_i), y = df_i$C)
#     lines(lines_i, col = cols[df_i$D[1]], lwd = 3)
#     # continuous trait b
#     df_i <- df.ForwardBranch[[2]][i:(i-1),]
#     lines_i <- list(x = rownames(df_i), y = df_i$C)
#     lines(lines_i, col = cols[df_i$D[1]], lwd = 3)
#     # continuous trait c
#     if(j >= 52){
#       df_i <- df.ForwardBranch[[3]][i:(i-1),]
#       lines_i <- list(x = rownames(df_i), y = df_i$C)
#       lines(lines_i, col = cols[df_i$D[1]], lwd = 3)
#     }
#     
#     # # discrete theta
#     # lines_i <- list(x = rownames(df_i), y = theta[df_i$D] + 8)
#     # lines(lines_i, col = cols[df_i$D[1]], lwd = 3)
#     # # discrete state
#     # lines_i <- list(x = rownames(df_i), y = c(20, 20))
#     # lines(lines_i, col = cols[df_i$D[1]], lwd = 3)
#   }
#   dev.off()
#   # Sys.sleep(0.1)
# }
# 
# ## list file names and read in
# imgs <- list.files(dir_out, full.names = TRUE)
# imgs <- imgs[sapply(1:99, function(x) grep(paste0("-", x, ".png"), imgs))]
# img_list <- lapply(imgs, image_read)
# 
# ## join the images together
# img_joined <- image_join(img_list)
# 
# ## animate at 2 frames per second
# img_animated <- image_animate(img_joined, fps = 10)
# 
# ## save to disk
# image_write(image = img_animated,
#             path = "~/presentation/2021_Evolution/gif_2.gif")
# 
# #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
# #### GIF 4 of presentation
# #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
# 
# df.ForwardBranch <- forwardSimBranch(theta0, state0, Q, alpha, sigma, theta, time.unit, time.steps)
# 
# for(j in 2:time.steps){
#   # png(file = paste0(dir_out, "gif_1-", j-1, ".png"), width = 1450, height = 660)
#   plot(1, type = "n", xlab = "", ylab = "", xlim = c(0, time.steps), ylim = c(-20, 20), axes=FALSE)
#   for(i in 2:j){
#     # continuous trait a
#     df_i <- df.ForwardBranch[i:(i-1),]
#     lines_i <- list(x = rownames(df_i), y = df_i$C - 12)
#     lines(lines_i, col = cols[df_i$D[1]], lwd = 3)
#     # discrete theta
#     lines_i <- list(x = rownames(df_i), y = theta[df_i$D] + 8)
#     lines(lines_i, col = cols[df_i$D[1]], lwd = 3)
#     # discrete state
#     lines_i <- list(x = rownames(df_i), y = c(20, 20))
#     lines(lines_i, col = cols[df_i$D[1]], lwd = 3)
#   }
#   # dev.off()
#   Sys.sleep(0.1)
# }
# 
# ## list file names and read in
# imgs <- list.files(dir_out, full.names = TRUE)
# imgs <- imgs[sapply(1:99, function(x) grep(paste0("-", x, ".png"), imgs))]
# img_list <- lapply(imgs, image_read)
# 
# ## join the images together
# img_joined <- image_join(img_list)
# 
# ## animate at 2 frames per second
# img_animated <- image_animate(img_joined, fps = 10)
# 
# ## save to disk
# image_write(image = img_animated,
#             path = "~/presentation/2021_Evolution/gif_2.gif")
# 
# 
# 
# df.ForwardBranch <- forwardSimBranch(theta0, state0, Q, alpha, sigma, theta, time.unit, time.steps)
# 
# for(j in 2:time.steps){
#   # png(file = paste0(dir_out, "gif_1-", j-1, ".png"), width = 1000, height = 500)
#   plot(1, type = "n", xlab = "", ylab = "", xlim = c(0, time.steps), ylim = c(-50, 50), axes=FALSE)
#   for(i in 2:j){
#     # discrete state
#     rect(rownames(df_i)[1], 40, rownames(df_i)[2], 50, col = cols[df_i$D[1]], border = NA)
#     # # discrete theta
#     lines_i <- list(x = mean(as.numeric(rownames(df_i))), y = (theta[df_i$D][1] * 2))
#     lines(lines_i, col = cols[df_i$D[1]], type = "h", lwd = 2)
#     # continuous trait a
#     df_i <- df.ForwardBranch[i:(i-1),]
#     lines_i <- list(x = rownames(df_i), y = df_i$C - 40)
#     lines(lines_i, col = cols[df_i$D[1]], lwd = 5)
#   }
#   # dev.off()
#   Sys.sleep(0.1)
# }
# 
# ## list file names and read in
# imgs <- list.files(dir_out, full.names = TRUE)
# imgs <- imgs[sapply(1:99, function(x) grep(paste0("-", x, ".png"), imgs))]
# img_list <- lapply(imgs, image_read)
# 
# ## join the images together
# img_joined <- image_join(img_list)
# 
# ## animate at 2 frames per second
# img_animated <- image_animate(img_joined, fps = 10)
# 
# ## save to disk
# image_write(image = img_animated,
#             path = "~/presentation/2021_Evolution/gif_1.gif")
# 
# 
# #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
# #### GIF 1 of presentation
# #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
# # creates a forward simulation gif for hOUwie
# # intial conditions
# theta0 = 0
# state0 = 1
# # model parameters
# Q = rbind(c(-1,1), c(1, -1))
# alpha = c(0.5, 0.5)
# sigma = c(0, 0)
# theta = c(10, 10)
# # resolution parameters
# time.unit <- 0.1
# time.steps <- 100
# # branching information
# 
# df.ForwardBranch <- forwardSimBranch(theta0, state0, Q, alpha, sigma, theta, time.unit, time.steps)
# 
# plot(x = rownames(df.ForwardBranch), y = df.ForwardBranch$C, xlab = "", ylab = "", xlim = c(0, time.steps), axes=FALSE, type = "l")
# 
# 
# dir_out <- "~/presentation/2021_Evolution/gif_1/"
# cols <- c("orange", "purple")
# 
# for(j in 100:time.steps){
#   # png(file = paste0(dir_out, "gif_1-", j-1, ".png"), width = 1500, height = 500)
#   plot(1, type = "n", xlab = "", ylab = "", xlim = c(0, time.steps), ylim = c(-20, 20), axes=FALSE)
#   for(i in 2:j){
#     # discrete state
#     rect(rownames(df_i)[1], 15, rownames(df_i)[2], 10, col = cols[df_i$D[1]], border = NA)
#     # # # discrete theta
#     # lines_i <- list(x = mean(as.numeric(rownames(df_i))), y = (theta[df_i$D][1] * 2))
#     # lines(lines_i, col = cols[df_i$D[1]], type = "h", lwd = 2)
#     # continuous trait a
#     df_i <- df.ForwardBranch[i:(i-1),]
#     lines_i <- list(x = rownames(df_i), y = df_i$C - 7.5)
#     lines(lines_i, col = cols[df_i$D[1]], lwd = 4)
#   }
#   # dev.off()
#   # Sys.sleep(0.1)
# }
# 
# ## list file names and read in
# imgs <- list.files(dir_out, full.names = TRUE)
# imgs <- imgs[sapply(1:99, function(x) grep(paste0("-", x, ".png"), imgs))]
# img_list <- lapply(imgs, image_read)
# 
# ## join the images together
# img_joined <- image_join(img_list)
# 
# ## animate at x frames per second
# img_animated <- image_animate(img_joined, fps = 10)
# 
# ## save to disk
# image_write(image = img_animated,
#             path = "~/presentation/2021_Evolution/gif_1.gif")
# 
# 
# # df.ForwardBranch <- lapply(1:1000, function(x) forwardSimBranch(theta0, round(runif(1))+1, Q, alpha, sigma, theta, time.unit, time.steps))
# # 
# # par(mfrow=c(1,3))
# # plot(density(unlist(lapply(df.ForwardBranch, function(x) x[100,2]))), main = "high rates")
# 
