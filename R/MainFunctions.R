isValidRIMNodesList <- function(l) {
  if(!is.list(l)) {
    return(F)
  }
  ranks = c()
  hasParent = numeric(length(l))
  for(i in 1:length(l)) {
    vec = l[[i]]
    if(!is.numeric(vec) || length(vec) != 5 || !(vec[3] %in% c(0,1)) ||
         (vec[3] == 0 && (vec[1] <= i || vec[2] <= i)) ||
         (vec[3] == 0 && (vec[1] > length(l) || vec[2] > length(l))) ||
         (vec[3] == 0 && (vec[1]%%1 != 0 || vec[2]%%1 != 0)) ||
         (vec[3] == 1 && vec[5]%%1 != 0)) {
      return(F)
    }
    if(vec[3] == 1) {
      ranks = c(ranks, vec[5])
    } else {
      hasParent[vec[1]] = hasParent[vec[1]] + 1
      hasParent[vec[2]] = hasParent[vec[2]] + 1
    }
  }
  numLeaves = (length(l)+1)/2
  if(length(ranks) != numLeaves || !all(sort(ranks) == 0:(numLeaves-1)) ||
       !all(hasParent == c(0,rep(1, length(l)-1)))) {
    return(F)
  }
  return(T)
}

isValidRIMSamples <- function(samples) {
  if(!is.matrix(samples) || nrow(samples) == 0 || ncol(samples) == 0 ||
       !all(sort(unique(as.numeric(samples))) == 0:(ncol(samples)-1))) {
    return(F)
  }
  if(!all(apply(samples, 1, function(x){length(unique(x)) == ncol(samples)}))) {
    return(F)
  }
  return(T)
}

isValidAveDiscMatrix = function(aveDiscMatrix) {
  if(!is.matrix(aveDiscMatrix) || ncol(aveDiscMatrix) != nrow(aveDiscMatrix) ||
       !all(aveDiscMatrix*lower.tri(aveDiscMatrix) == aveDiscMatrix) ||
       !all(aveDiscMatrix >= 0) || !all(aveDiscMatrix <= 1)) {
    return(F)
  }
  return(T)
}

isValidRanking = function(ranking) {
  if(!is.numeric(ranking) ||
       !all(sort(ranking) == 0:(length(ranking)-1))) {
    return(F)
  }
  return(T)
}

sampleFromRIM <- function(n, rimNodesList) {
  if(!is.numeric(n) || length(n) != 1 || n%%1 != 0 || n <= 0) {
    stop("ERROR: sampleFromRIM requires n be a positive integer.\n")
  }
  if(!isValidRIMNodesList(rimNodesList)) {
    stop("ERROR: incorrectly formatted rimNodesList in sampleFromRIM.\n")
  }
  return(RCPPSampleFromRIM(n, rimNodesList))
}

averageDiscMatrix <- function(samples) {
  if(!isValidRIMSamples(samples)) {
    stop("ERROR: input samples to averageDiscMatrix are incorrectly formatted.\n")
  }
  return(RCPPAverageDiscMatrix(samples))
}

structByDP <- function(aveDiscMatrix, refRanking, makeCanonical) {
  if(!isValidAveDiscMatrix(aveDiscMatrix)) {
    stop("ERROR: incorrectly formatted aveDiscMatrix inputted to structByDP.\n")
  }
  if(!isValidRanking(refRanking)) {
    stop("ERROR: incorrectly formatted refRanking inputted to structByDP.\n")
  }
  if(!is.logical(makeCanonical)) {
    stop("ERROR: input makeCanonical for structByDP is not a logical value.\n")
  }
  return(RCPPStructByDP(aveDiscMatrix, refRanking, makeCanonical))
}

SASearch <- function(aveDiscMatrix, refRanking, inverseTemp, maxIter, makeCanonical, verbose=F) {
  if(!isValidAveDiscMatrix(aveDiscMatrix)) {
    stop("ERROR: incorrectly formatted aveDiscMatrix inputted to SASearch\n")
  }
  if(!isValidRanking(refRanking)) {
    stop("ERROR: incorrectly formatted refRanking inputted to SASearch\n")
  }
  if(!is.numeric(inverseTemp) || inverseTemp < 0) {
    stop("ERROR: incorrectly formatted inverseTemp inputted to SASearch\n")
  }
  if(!is.numeric(maxIter) || maxIter%%1 != 0 || maxIter <= 0) {
    stop("ERROR: incorrectly formatted maxIter inputted to SASearch\n")
  }
  if(!is.logical(makeCanonical)) {
    stop("ERROR: input makeCanonical for SASearch is not a logical value.\n")
  }
  return(RCPPSASearch(aveDiscMatrix, refRanking, inverseTemp, maxIter, makeCanonical, verbose))
}

logProbRIM <- function(rimNodesList, aveDiscMatrix) {
  if(!isValidAveDiscMatrix(aveDiscMatrix)) {
    stop("ERROR: incorrectly formatted aveDiscMatrix inputted to logProbRIM\n")
  }
  if(!isValidRIMNodesList(rimNodesList)) {
    stop("ERROR: incorrectly formatted rimNodesList inputted to logProbRIM\n")
  }
  return(RCPPLogProbRIM(rimNodesList, aveDiscMatrix))
}

thetaMLERIM <- function(rimNodesList, aveDiscMatrix) {
  if(!isValidAveDiscMatrix(aveDiscMatrix)) {
    stop("ERROR: incorrectly formatted aveDiscMatrix inputted to thetaMLERIM\n")
  }
  if(!isValidRIMNodesList(rimNodesList)) {
    stop("ERROR: incorrectly formatted rimNodesList inputted to thetaMLERIM\n")
  }
  return(RCPPthetaMLERIM(rimNodesList, aveDiscMatrix))
}
