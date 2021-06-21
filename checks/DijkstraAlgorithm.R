require(igraph)
model <- CID.cor
init <- 1
final <- 4

# for a given, vertex i
# see where vertex i leads
# if vertex i leads to final, done
# if vertex i does not lead to final, visit nodes
# add vertex i to path and assign new nodes as vertex i


nStates <- dim(model)[1]
Dist <- matrix(Inf, nStates, nStates)
Next <- matrix(NA, nStates, nStates)
for(V_index in sequence(nStates)){
  Dist[V_index, V_index] <- 0
  Next[V_index, V_index] <- V_index
}
for(V_index in sequence(nStates)){
  To <- which(model[V_index, ] > 0)
  Dist[V_index, To] <- 1
  Next[V_index, To] <- To
}
for(k in sequence(nStates)){
  for(i in sequence(nStates)){
    for(j in sequence(nStates)){
      if(Dist[i,j] > Dist[i,k] + Dist[k,j]){
        Dist[i,j] = Dist[i,k] + Dist[k,j]
        Next[i,j] = Next[i,k]
      }
    }
  }
}
if(is.na(Next[init, final])){
  stop("Impossible transition on branch detected...")
}else{
  path = init
  u = init
  v = final
  while(u != v){
    u <- Next[u, v]
    path <- c(path, u)
  }
}


vertex.set <- c()

count <- 1
distance <- rep(Inf, dim(model)[1])
unvisited <- sequence(dim(model)[1])
distance[init] <- 0
current.node <- init

cbind(current.node, which(model[current.node,] > 0))

distance[which(model[current.node,] > 0)] <- count
distance[current.node]

visited <- init
current.path <- init
potential.visits <- which(model[init,] > 0)
not.visited <- potential.visits[!potential.visits %in% visited]
for(i in 1:length(not.visited)){
  which(model[not.visited[i],] > 0)
}

