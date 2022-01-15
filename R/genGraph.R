
#' Function to generate bow free graph containing multi-edges
#' The distribution of the error terms can be changed in the function "graphSample".
#' @param p is the number of maximum vertices. The generated one has
#'        vertices <= p
#' @param d is the number of maximum edges in the graph. The generated one
#'        has edges<= the minimum of (d, p(p-1)/2)
#' @param n is the generated sample size. 
#' @return 
#' \itemize{
#' \item Y is the sample matrix with dim (n,k) where k<=p
#' \item directedGraph is a k by k matrix representing the directed part of the graph
#' \item multiedge is list of vector represent multi-treks
#' }
randomGraph<- function(p,d,n,g){
  d = min(d, p*(p-1)/2)
  vector = c()
  for(i in 1:(p-1)){
    for(j in (i+1):p){
      vector = cbind(vector, c(i,j))
    }
  }
  dim = dim(vector)
  edgeSelectedIndex = sort(sample(1:dim[2],d))
  edgeSelected = vector[, edgeSelectedIndex]
  vertexWithPar = unique(edgeSelected[2, ])
  vertexWithoutPar = setdiff(unique(edgeSelected[1,]), vertexWithPar)
  
  #k is the number of vertices randomly removed
  k = sample(1:length(vertexWithoutPar), 1)
  vertexRemoved = sample(vertexWithoutPar, k)
  edgeSelected = removeBow(edgeSelected, vertexRemoved)
  
  vertexRemained = setdiff(1:p, vertexRemoved)
  
  Y = graphSample(edgeSelected, vertexWithPar, vertexWithoutPar, vertexRemained, p, n)
  graph = graphPlot(edgeSelected, vertexRemoved, vertexRemained, p)
  directGraph = graph$directGraph
  biEdge = graph$biEdge
  multiEdge = graph$multiEdge
  
  return(list(Y=Y, directGraph=directGraph, biEdge = biEdge, multiEdge=multiEdge))
}

#' This function is a helper function
#' errors can be changed.
#' 
#' @param edge is 2xd matrix containing selected edges
#' @param u is vertices with parents
#' @param v is vertices in edge without parents
#' @param r is vertices remained after randomly removing vertices
#' @param p is the number of initial vertices
#' @param n is the sample size
#' @return a sample matrix n by r (r is the length of the parameter r)

graphSample <- function(edge, u, v, r, p, n){
  
  ########## here is setting the distribution of error
  err = c()
  for(i in 1:p){
#    err = cbind(err, runif(n,min=-10,max=10))
#    err = cbind(err, rt(n,df=10))
#    err = cbind(err, rgamma(n,shape=2,rate=4)-2/4)
    err = cbind(err, rchisq(n,df=2)-2)
  }
  ##########
  
  Y = c()
  for(i in 1:p){
    if (!(i %in% u)){
      Y = cbind(Y, err[,i])
    }else{
      y = err[,i]
      vertices = unique((edge[,edge[2,]==i,drop = FALSE])[1,])
      for(v in vertices){
        k = runif(1,0.6,1)
        y = y + (k)*Y[,v]
      }
      Y = cbind(Y, y)
    }
  }
  return(Y[,r])
}

#' This function is a helper function. 
#' 
#' @param edge is 2xd matrix containing selected edges
#' @param u is vertices removed
#' @param r is vertices remained after randomly removing vertices
#' @param p is the number of initial vertices
#' @return list(directedGraph, multitreks) directed graph is a r by r matrix, 
#'         multitreks is a list of vector
graphPlot <- function(edge, u, r, p){
  G = diag(p)
  vertexWithChild = unique(edge[1, ])
  for (i in vertexWithChild){
    entryIndex = (edge[,edge[1,]==i, drop = FALSE])[2,]
    G[i, entryIndex]=1
  }
  G = G[r,r]
  
  bG = diag(length(r))
  
  multiEdge = list()
  for (i in u){
    entryIndex = (edge[,edge[1,]==i, drop = FALSE])[2,]
    k=1
    for (j in entryIndex){
      entryIndex[k]=which(r==j)
      k = k+1
    }
    if (length(entryIndex)>1){
      multiEdge[[length(multiEdge)+1]]=entryIndex
    }
    for(t in entryIndex){
      temp = entryIndex[entryIndex>t]
      for (s in temp){
        bG[t,s]=1
        bG[s,t]=1
      }
    }
  }
  return(list(directGraph = G, biEdge = bG, multiEdge = multiEdge))
}

#' This is a helper function
#' 
#' @param edge is the set of selected edges
#' @param vertices is the set of vertices to be removed
#' @return vector of edges containing no bow
removeBow <- function(edge, vertices){
  
  #etr is edge to be returned
  etr = edge
  for(v in vertices){
    vector = (edge[,edge[1,]==v, drop = FALSE])[2,]
    
    vec = c()
    for (i in vector){
      temp = vector[vector > i]
      for (j in temp){
        etr = etr[, (etr[1,]!=i|etr[2,]!=j), drop = FALSE]
      }
    }
  }
  return(etr)
}

#' This compare if two list have the same item
compareList <- function(l1, l2){
  if(length(l1)>0){
    for(i in l1){
      if (!(Position(function(x) identical(x,i),l2,nomatch = 0)>0)){
        return(FALSE)
      }
    }
  }
  if(length(l1)==length(l2)){
    return(TRUE)
  }else{
    return(FALSE)
  }
}



