# In computer science, the Floyd–Warshall algorithm (also known as Floyd's algorithm, the Roy–Warshall algorithm, the Roy–Floyd algorithm, or the WFI algorithm) is an algorithm for finding shortest paths in a directed weighted graph with positive or negative edge weights (but with no negative cycles).[1][2] A single execution of the algorithm will find the lengths (summed weights) of shortest paths between all pairs of vertices. Although it does not return details of the paths themselves, it is possible to reconstruct the paths with simple modifications to the algorithm. Versions of the algorithm can also be used for finding the transitive closure of a relation R {\displaystyle R} R, or (in connection with the Schulze voting system) widest paths between all pairs of vertices in a weighted graph. 
model <- CID.cor
init <- 1
final <- 4

FloydWalshAlg(model, init, final)

FloydWalshAlg <- function(model, init, final){
  nStates <- dim(model)[1]
  Dist <- matrix(Inf, nStates, nStates)
  Next <- matrix(NA, nStates, nStates)
  for(V_index in sequence(nStates)){
    Dist[V_index, V_index] <- 0
    Next[V_index, V_index] <- V_index
  }
  for(V_index in sequence(nStates)){
    To <- which(model[V_index, ] > 0)
    Dist[V_index, To] <- 1/model[V_index, To]
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
    return(NULL)
  }else{
    path = init
    u = init
    v = final
    while(u != v){
      u <- Next[u, v]
      path <- c(path, u)
    }
  }
  return(path)
}
