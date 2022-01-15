#'  Determines multidirected edges from a data set and the corresponding estimated graphical model.
#'  The estimated graphical model is the object returned by ngBap::bang.
#'  Assume the estimated graphical model is correct.
#'  Assume all the sample variables have zero expectation.

library(partitions)
library(ngBap)

#' @param Y p x n matrix of observations with variable in row, sample in column.
#' @param model Graphical model estimation from Y. This is the object returned by ngBap::bang.
#' @param tolerance A non-negative number such that cumulant!=0 is recognized by abs(cumulant)>tolerance.
#' @return multiEdges A list of bidirected edges and multidirected edges.
#'  Each element is a set of the sinks of the edge.
#'  Each sink is represented by its row index in Y.
MBANG=function(Y,model,tolerance=0.05)
{
  #Remove all the directed edges and then normalize the data to make standard deviation = 1.
  Ybiedge=scale(rmvdedge(Y,model),center = F,scale = T)
  #Find multidiected edges.
  P=1:nrow(Y)
  multiEdges=findmultiedge(Ybiedge,model,trek=list(),R=c(),P=P,Q=c(),tolerance=tolerance)

  #Assume the BANG result is correct, add any omitted bidirected edge back to the list.
  for(S in multiEdges)
  {
    #Remove the vertices in S from the bidirected edge matrix.
    model$bEdge=change_matrix(model,S)
  }
  biEdges=findbiedge(model$bEdge) #Assume all the remainging edges are bidirected.

  return(append(multiEdges,biEdges))
}



#' Removes all the directed edges from the sample matrix Y based on the corresponding graphical model.
#' @param Y p x n matrix of observations with variable in row, sample in column.
#' @param model The corresponding graphical model estimation from Y, using the function "bang".
#' @return Ybiedge The sample matrix with directed effects removed.
rmvdedge=function(Y,model)
{
  B=model$directEffect
  Ybiedge=Y-B%*%Y
  return(Ybiedge)
}



#' Finds multidirected edges in the model.
#' Assumes the model generating the matrix Ybiedge has no directed effect.
#' R, P, Q are disjoint subsets of the vertices of the graphical model.
#' @param Ybiedge p x n matrix of observations with variable in row, sample in column.
#' @param model Graphical model estimation returned by ngBap::bang.
#' @param trek The existing multidirected edges.
#' @param R A set such that each multiedge found by this function contains R as a subset.
#' @param P A set containing all the possible new vertices having a common hidden variable with the set R.
#' @param Q A set such that each multiedge found by this function contains none of the element in Q.
#' @param tolerance A non-negative number such that cumulant!=0 is recognized by abs(cumulant)>tolerance.
#' @return trek The list of multidirected edges.
findmultiedge=function(Ybiedge,model,trek=list(),R=c(),P,Q=c(),tolerance=0.05)
{
  if(length(P)==0 & length(Q)==0)
  {
    trek[[length(trek)+1]]<-R
    return (trek)
  }

  for(v in P)
  {
    edge=model$bEdge
    Nv=setdiff(which(edge[,v]!=0),c(v))
    #Only consider vertices with more than 1 neighbourhoods.
    if(length(Nv)>1 && maxcumulant(Ybiedge[union(R,v),],tolerance))
    {
      trek=findmultiedge(Ybiedge,model,trek,union(R,c(v)),intersect(P,Nv),intersect(Q,Nv),tolerance)
    }
    P=setdiff(P,c(v))
    Q=union(Q,c(v))
  }
  return(trek)
}


#Converts the bidirected edge matrix to a list of bidirected edges.
#Each element in the list is the two sinks of the bidirected edge.
#' @param matrix The bidirected edge matrix returned by the function "bang".
#' @return trek The list of bidirected edges.
findbiedge=function(matrix)
{
  trek=list()
  for(i in 1:(ncol(matrix)-1))
  {
    for(j in (i+1):ncol(matrix))
    {
      if(matrix[i,j]==1)
      {
        trek[[length(trek)+1]]<-c(i,j)
      }
    }
  }
  return(trek)
}



#----------------------------------------------------------------------------------------------------
#Below are the helper functions.


#Removes the vertices in the set R from the bidirected edge matrix.
change_matrix=function(model,R)
{
 edge=model$bEdge
 for(i in R)
 {
   rowtochange=edge[i,]
   rowtochange[R]=0
   #Update the results to the matrix.
   edge[i,]=rowtochange
 }
 return(edge)
}


#If the cumulant of data is greater than the tolerance, return True directly.
#Otherwise, it repeats each row of data and finds the corresponding cumulants.
#Then, it returns True if the maximum of the absolute values of cumulants is greater than the tolerance.
#Returns False otherwise.
maxcumulant=function(data,tolerance)
{
  maxval=abs(cumulant(data))
  if(maxval>tolerance)
  {
    return (TRUE)
  }
  else
  {
    if(is.vector(data))
    {
      data=as.matrix(t(data))
    }
    for(i in 1:nrow(data))
    {
      maxval=max(maxval,abs(cumulant(rbind(data,data[i,]))))
    }
    return (maxval>tolerance)
  }
}



#' Calculates the cumulant of the data.
#' @param data n x m matrix (or vector) of observations with variable in row, sample in column.
#' @return sumofcum The cumulant of the data.
cumulant=function(data)
{
  if(is.vector(data))
  {
    data=as.matrix(t(data))
  }

  x=1:nrow(data)
  parts = listParts(length(x))
  allpartition = rapply(parts, function(ii) x[ii], how="replace")
  sumofcum=0
  for(y in allpartition)
  {
    prodmom=c()
    for(element in y)
    {
      prodmom=c(prodmom,moment(data[element,]))
    }
    l=length(y)
    toadd=((-1)^(l-1))*factorial(l-1)*prod(prodmom)
    sumofcum=sumofcum+toadd
  }
  return (sumofcum)
}


#' Calculates the moment of the data.
#' @param data n x m matrix (or vector) of observations with variable in row, sample in column.
#' @return mom The moment of the data.
moment=function(data)
{
  if(is.vector(data))
  {
    data=as.matrix(t(data))
  }

  mom=mean(apply(data, 2, prod))
  return(mom)
}

