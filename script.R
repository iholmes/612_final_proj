n.host <- 100
hostMat <- matrix(NA, nrow=n.host, ncol=2)


hostMat[,1] <- runif(100, 0, 1)

hostMat[,2] <- 1-hostMat[,1]

n.mic.gen <- 10

# number of offspring per generation

nOff1 <- 2
nOff2 <- 20

# starting population size

micStart <- 1e3
micK <- 1e5

# interaction term
alpha12 <- 0.5
alpha21 <- 0.05

internalLV.show <- function(vector){
  
  tmp <- matrix(NA, nrow=10, ncol=2)
  tmp[1,1] <- vector[1]*micStart
  tmp[1,2] <- vector[2]*micStart
  
  
  for(i in 2:10){
    tmp[i,1] <-tmp[i-1,1]*(exp(nOff1*(1-(tmp[i-1,1]+alpha12*tmp[i-1,2])/micK)))
    
    tmp[i,2] <-tmp[i-1,2]*(exp(nOff1*(1-(tmp[i-1,2]+alpha21*tmp[i-1,1])/micK)))
  }
  return(tmp)
}

show <- internalLV.show(c(50, 50))

plot(x=1:10, y=seq(0, 130000, 130000/9), type="n", xlab="time step", ylab="population")
lines(show[,1], col="blue")
lines(show[,2], col="red")

nextSeed <- function(oneGenVec, tProb1, tProb2){
  oneTrans <- oneGenVec[1]*tProb1*micStart/((oneGenVec[1]*tProb1+ oneGenVec[2]*tProb2)+1)
  
  twoTrans <- oneGenVec[2]*tProb2*micStart/((oneGenVec[1]*tProb1+ oneGenVec[2]*tProb2)+1)
  
  return(c(oneTrans, twoTrans))
}

make.input <- function(alphamin, alphamax, tProbmin, tProbmax, nOffmin, nOffmax, step.size){
  alpha <- seq(alphamin, alphamax, step.size)
  tprob <- seq(tProbmin, tProbmax, step.size)
  nOff <- seq(nOffmin, nOffmax)
  alphmat <- cbind(rep(alpha, each=length(alpha)), rep(alpha, length(alpha)))
  
  probmat1 <- cbind(rep(alphmat[,1], length(tprob)), rep(alphmat[,2], length(tprob)), 
                    rep(tprob, each=nrow(alphmat)))
  probmat2 <- cbind(rep(probmat1[,1], length(tprob)), rep(probmat1[,2], length(tprob)), rep(probmat1[,3], length(tprob)),rep(tprob, each=nrow(probmat1)))
  
  off1mat <- cbind(rep(probmat2[,1], length(nOff)), rep(probmat2[,2], length(nOff)), rep(probmat2[,3], length(nOff)),rep(probmat2[,4], length(nOff)),rep(nOff, each=nrow(probmat2)))
  
  off2mat <- cbind(rep(off1mat[,1], length(nOff)), rep(off1mat[,2], length(nOff)), rep(off1mat[,3], length(nOff)),
                   rep(off1mat[,4], length(nOff)),rep(off1mat[,5], length(nOff)), rep(nOff, nrow(off1mat)))
  
  off2mat <- unique(off2mat)
  return(off2mat)
}

inMat <- make.input(0.1, 0.4, 0.1, 0.4, 1, 3, 0.1)

n.gen=10

paramTest <- function(inVec, hostMat){
  
  oneMat <- matrix(NA, nrow=n.host, ncol=n.gen)
  oneMat[,1] <- hostMat[,1]*micStart
  
  twoMat <- matrix(NA, nrow=n.host, ncol=n.gen)
  twoMat[,1] <- hostMat[,2]*micStart
  
  internalLV <- function(vec){
    
    tmp <- matrix(NA, nrow=10, ncol=2)
    tmp[1,1] <- vec[1]
    tmp[1,2] <- vec[2]
    
    
    for(i in 2:10){
      tmp[i,1] <-tmp[i-1,1]*(exp(inVec[5]*(1-(tmp[i-1,1]+inVec[1]*tmp[i-1,2])/micK)))
      
      tmp[i,2] <-tmp[i-1,2]*(exp(inVec[6]*(1-(tmp[i-1,2]+inVec[2]*tmp[i-1,1])/micK)))
    }
    return(tmp[10,])
  }
  
  
  
  for(k in 2:n.gen){
    
    oneGenPop <- apply(cbind(oneMat[,k-1], twoMat[,k-1]), 1, internalLV)
    
    nextPop <- apply(oneGenPop, 2, nextSeed, tProb1=inVec[3], tProb2=inVec[4])
    
    oneMat[,k] <- nextPop[1,]
    twoMat[,k] <- nextPop[2,]
    
  }
  
  return(cbind(mean(colMeans(oneMat)), mean(colMeans(twoMat))))
}

micMat <- apply(inMat, 1, paramTest, hostMat=hostMat)

tmat <- cbind(inMat, micMat[1,], micMat[2,])

plot(x=seq(min(tmat[,2]), max(tmat[,2]), by=0.1), y=seq(from=0, to=1000, by=(1000/(length(seq(min(tmat[,2]), max(tmat[,2]), by=0.1))-1))), type="n", xlab="competition coefficient", ylab="microbe 1 population")

points(x=tmat[,2],y=tmat[,7])

plot(x=seq(min(tmat[,3]), max(tmat[,3]), by=0.1), y=seq(from=0, to=1000, by=(1000/(length(seq(min(tmat[,3]), max(tmat[,3]), by=0.1))-1))), type="n", xlab="transition probability", ylab="microbe 1 population")

points(x=tmat[,3],y=tmat[,7])

plot(x=seq(min(tmat[,6]), max(tmat[,6])), y=seq(from=0, to=1000, by=(1000/(length(seq(min(tmat[,6]), max(tmat[,6])))-1))), type="n", xlab="intrinsic rate of increase", ylab="microbe 1 population")

points(x=tmat[,5],y=tmat[,7])


pairwiseDist <- function(testPoint, outVec)(sqrt((testPoint-outVec[1])^2+((micK-testPoint)-outVec[2])^2))

checkPt <- function(testPt, k){
  distances <- apply(tmat[,6:7], 1, pairwiseDist, testPoint=testPt)
  
  kNeighbors <- order(distances)[1:k]
  
  kParams <- inMat[kNeighbors,]
  
  weightedDist <- (1/distances[kNeighbors])^2
  
  par1est <- sum(kParams[,1]*(weightedDist/sum(weightedDist)))
  par2est <- sum(kParams[,2]*(weightedDist/sum(weightedDist)))
  par3est <- sum(kParams[,3]*(weightedDist/sum(weightedDist))) 
  par4est <- sum(kParams[,4]*(weightedDist/sum(weightedDist)))
  
  return(c(par1est, par2est, par3est, par4est))
}

checkPt(300, 4)

