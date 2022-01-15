library(ngBap)
library(multiedge)
library(partitions)

testModel <- function(p,d,n, bangR, bothR){
  b = bangR
  bR = bothR
  model = randomGraph(p,d,n)
  Y = model$Y
  bResult = bang(Y, K = 4, level = .01, verbose = F)
  fResult = MBANG(t(Y), bResult,tolerance = 0.5)
  directEdge = t(bResult$dEdge+diag(dim(Y)[2]))
  if(all(directEdge == model$directGraph) & all(bResult$bEdge == model$biEdge)){
    b = b + 1
    if(compareList(fResult, model$multiEdge)){
      bR = bR+1
    }
  }

  corrmulti=length(which(fResult %in% model$multiEdge))
  totalmulti=length(model$multiEdge)
  if(totalmulti!=0){
    percent=corrmulti/totalmulti
    }
  if(totalmulti==0){
    if(length(fResult)==0){
      percent=1
    }
    else{
      percent=0
    }
  }
  return(list(bangR=b, bothR=bR,percent=percent))
}


#n = 10000
#n = 50000
n = 25000
bangR = 0
bothR = 0
totalNum = 0
percent=0
p=7
d=12
for(i in 1:100){
  temp = testModel(p,d,n, bangR, bothR)
  bangR = temp$bangR
  bothR = temp$bothR
  totalNum = totalNum + 1
  percent=percent+temp$percent
}
print(c(bangR,bothR, totalNum,percent))



#n = 10000
#n = 50000
n = 25000
bangR = 0
bothR = 0
totalNum = 0
percent=0
p=7
d=8
for(i in 1:100){
  temp = testModel(p,d,n, bangR, bothR)
  bangR = temp$bangR
  bothR = temp$bothR
  totalNum = totalNum + 1
  percent=percent+temp$percent
}
print(c(bangR,bothR, totalNum,percent))



#n = 10000
#n = 50000
n = 25000
bangR = 0
bothR = 0
totalNum = 0
percent=0
p=7
d=5
for(i in 1:100){
  temp = testModel(p,d,n, bangR, bothR)
  bangR = temp$bangR
  bothR = temp$bothR
  totalNum = totalNum + 1
  percent=percent+temp$percent
}
print(c(bangR,bothR, totalNum,percent))

