isValidRIMNodesMat <- function(l) {
  if (!is.matrix(l) || !is.numeric(l) || ncol(l) != 5 || nrow(l)  == 0) {
    return(F)
  }
  ranks = c()
  hasParent = numeric(nrow(l))
  for (i in 1:nrow(l)) {
    vec = l[i,]
    if (!(vec[3] %in% c(0,1)) ||
       (vec[3] == 0 && (vec[1] < i || vec[2] < i)) ||
       (vec[3] == 0 && (vec[1] >= length(l) || vec[2] >= length(l))) ||
       (vec[3] == 0 && (vec[1] %% 1 != 0 || vec[2] %% 1 != 0)) ||
       (vec[3] == 1 && vec[5] %% 1 != 0)) {
      return(F)
    }
    if (vec[3] == 1) {
      ranks = c(ranks, vec[5])
    } else {
      hasParent[vec[1] + 1] = hasParent[vec[1] + 1] + 1
      hasParent[vec[2] + 1] = hasParent[vec[2] + 1] + 1
    }
  }
  numLeaves = (nrow(l) + 1) / 2
  if (length(ranks) != numLeaves ||
      !all(sort(ranks) == 0:(numLeaves - 1)) ||
      !all(hasParent == c(0, rep(1, nrow(l) - 1)))) {
    return(F)
  }
  return(T)
}

isValidRIMSamples <- function(samples) {
  if (!is.matrix(samples) || nrow(samples) == 0 || ncol(samples) == 0 ||
       !all(sort(unique(as.numeric(samples))) == 0:(ncol(samples)-1))) {
    return(F)
  }
  if (!all(apply(samples, 1, function(x){length(unique(x)) == ncol(samples)}))) {
    return(F)
  }
  return(T)
}

isValidAveDiscMatrix = function(aveDiscMatrix) {
  if (!is.matrix(aveDiscMatrix) || ncol(aveDiscMatrix) != nrow(aveDiscMatrix) ||
       !all(aveDiscMatrix*lower.tri(aveDiscMatrix) == aveDiscMatrix) ||
       !all(aveDiscMatrix >= 0) || !all(aveDiscMatrix <= 1)) {
    return(F)
  }
  return(T)
}

isValidRanking = function(ranking) {
  if (!is.numeric(ranking) ||
       !all(sort(ranking) == 0:(length(ranking)-1))) {
    return(F)
  }
  return(T)
}

rRIM <- function(n, rimNodesMat) {
  if (!is.numeric(n) || length(n) != 1 || n %% 1 != 0 || n <= 0) {
    stop("ERROR: rRIM requires n be a positive integer.\n")
  }
  if (!isValidRIMNodesMat(rimNodesMat)) {
    stop("ERROR: incorrectly formatted rimNodesMat in rRIM.\n")
  }
  return(RCPPSampleFromRIM(n, rimNodesMat))
}

averageDiscMatrix <- function(samples) {
  if (!isValidRIMSamples(samples)) {
    stop("ERROR: input samples to averageDiscMatrix are incorrectly formatted.\n")
  }
  return(RCPPAverageDiscMatrix(samples))
}

structByDP <- function(aveDiscMatrix, refRanking, makeCanonical) {
  if (!isValidAveDiscMatrix(aveDiscMatrix)) {
    stop("ERROR: incorrectly formatted aveDiscMatrix inputted to structByDP.\n")
  }
  if (!isValidRanking(refRanking)) {
    stop("ERROR: incorrectly formatted refRanking inputted to structByDP.\n")
  }
  if (!is.logical(makeCanonical)) {
    stop("ERROR: input makeCanonical for structByDP is not a logical value.\n")
  }
  if (ncol(aveDiscMatrix) != length(refRanking)) {
    stop("ERROR: number of columns of discrepancy matrix not equal to length of reference ranking in structByDP.\n")
  }
  return(RCPPStructByDP(aveDiscMatrix, refRanking, makeCanonical))
}

SASearch <- function(aveDiscMatrix, refRanking, inverseTemp, maxIter, makeCanonical, verbose=F) {
  if (!isValidAveDiscMatrix(aveDiscMatrix)) {
    stop("ERROR: incorrectly formatted aveDiscMatrix inputted to SASearch.\n")
  }
  if (!isValidRanking(refRanking)) {
    stop("ERROR: incorrectly formatted refRanking inputted to SASearch.\n")
  }
  if (!is.numeric(inverseTemp) || inverseTemp < 0) {
    stop("ERROR: incorrectly formatted inverseTemp inputted to SASearch.\n")
  }
  if (!is.numeric(maxIter) || maxIter %% 1 != 0 || maxIter <= 0) {
    stop("ERROR: incorrectly formatted maxIter inputted to SASearch.\n")
  }
  if (!is.logical(makeCanonical)) {
    stop("ERROR: input makeCanonical for SASearch is not a logical value.\n")
  }
  if (!is.logical(verbose)) {
    stop("ERROR: input verbose for SASearch is not a logical value.\n")
  }
  if (ncol(aveDiscMatrix) != length(refRanking)) {
    stop("ERROR: number of columns of discrepancy matrix not equal to length of reference ranking in SASearch.\n")
  }
  return(RCPPSASearch(aveDiscMatrix, refRanking, inverseTemp, maxIter, makeCanonical, verbose))
}

pRIM <- function(rimNodesMat, aveDiscMatrix, log.p = FALSE) {
  if (!isValidAveDiscMatrix(aveDiscMatrix)) {
    stop("ERROR: incorrectly formatted aveDiscMatrix inputted to pRIM.\n")
  }
  if (!isValidRIMNodesMat(rimNodesMat)) {
    stop("ERROR: incorrectly formatted rimNodesMat inputted to pRIM.\n")
  }
  if ((nrow(rimNodesMat) + 1) / 2 != ncol(aveDiscMatrix)) {
    stop("ERROR: rimNodesMat and aveDiscMatrix are on different numbers of items in pRIM.\n")
  }
  if (log.p) {
    return(RCPPLogProbRIM(rimNodesMat, aveDiscMatrix))
  } else {
    return(exp(RCPPLogProbRIM(rimNodesMat, aveDiscMatrix)))
  }
}

thetaMLERIM <- function(rimNodesMat, aveDiscMatrix) {
  if (!isValidAveDiscMatrix(aveDiscMatrix)) {
    stop("ERROR: incorrectly formatted aveDiscMatrix inputted to thetaMLERIM\n")
  }
  if (!isValidRIMNodesMat(rimNodesMat)) {
    stop("ERROR: incorrectly formatted rimNodesMat inputted to thetaMLERIM\n")
  }
  if ((nrow(rimNodesMat) + 1) / 2 != ncol(aveDiscMatrix)) {
    stop("ERROR: rimNodesMat and aveDiscMatrix are on different numbers of items in thetaMLERIM\n")
  }
  return(RCPPthetaMLERIM(rimNodesMat, aveDiscMatrix))
}

treeMatrixToGraph = function(treeMatrix, colorLeftRight=T) {
  g = graph.empty(nrow(treeMatrix))
  V(g)$color = "#DDFFFF"
  V(g)$theta = 0
  V(g)$rank = 0

  numInternalNodes = (nrow(treeMatrix) + 1)/2 - 1
  maxAbsTheta = max(abs(treeMatrix[,4])) + .0000001
  colors = colorRampPalette(c("red", "white", "green"))(2*numInternalNodes + 1)
  for (i in 1:nrow(treeMatrix)) {
    if (treeMatrix[i,3] == 0) {
      g = add.edges(g, c(i, treeMatrix[i,1] + 1, i, treeMatrix[i,2] + 1))
      if(colorLeftRight) {
        V(g)$color[treeMatrix[i,1] + 1] = "green"
        V(g)$color[treeMatrix[i,2] + 1] = "yellow"
      } else {
        V(g)$color[i] = colors[floor(numInternalNodes*treeMatrix[i,4]/maxAbsTheta) + numInternalNodes + 1]
      }
      V(g)$theta[i] = treeMatrix[i,4]
    } else {
      V(g)$rank[i] = treeMatrix[i,5]
    }
  }
  V(g)$name = round(V(g)$theta + V(g)$rank, 3)
  return(g)
}


plotTreeMatrix = function(treeMatrix,
                          leafNames=0:((nrow(treeMatrix) + 1)/2 - 1),
                          roundDecimals=1,
                          leafSize=20,
                          internalNodeSize=15,
                          colorLeftRight=T,
                          arrow.size=1,
                          edge.color="black",
                          vertex.label.cex=1,
                          RIMLayout=T) {
  g = treeMatrixToGraph(treeMatrix, colorLeftRight)
  toRename = which(treeMatrix[,3] == 1)
  V(g)$name[toRename] = leafNames[V(g)$name[toRename] + 1]
  V(g)$name[-toRename] = round(as.numeric(V(g)$name[-toRename]), roundDecimals)

  vertexShapes = rep("circle", length(V(g)))
  vertexShapes[toRename] = "rectangle"
  vertexSizes = rep(internalNodeSize, length(V(g)))
  vertexSizes[toRename] = leafSize

  if(RIMLayout) {
    layout = layout.RIM(g, root=1)
  } else {
    layout = layout.reingold.tilford(g, root=1)
  }

  oldMar = par("mar")
  par(mar = rep(0,4))
  plot.igraph(g,
              layout = layout,
              edge.color="black",
              vertex.size=vertexSizes,
              vertex.shape=vertexShapes,
              vertex.label.cex=vertex.label.cex,
              edge.arrow.mode = 0)
  par(mar = oldMar)
}

leafOrdering = function(g, node=1) {
  outList = get.adjlist(g, mode = "out")
  if(length(outList[[node]]) == 0) {
    return(node)
  }
  ordering = c()
  for(i in outList[[node]]) {
    ordering = c(ordering, leafOrdering(g, i))
  }
  return(ordering)
}

layout.RIM = function(g, root=1) {
  outList = get.adjlist(g,mode="out")
  layoutMat = layout.reingold.tilford(g, root=root)
  layoutMat[,1] = layoutMat[,1]/2

  leafNodesInOrder = leafOrdering(g)
  layoutMat[leafNodesInOrder,2] = min(layoutMat[leafNodesInOrder,2]) + rep(c(-1,0), ceiling(length(leafNodesInOrder)/2))[1:length(leafNodesInOrder)]
  a = min(layoutMat[leafNodesInOrder,1])
  b = max(layoutMat[leafNodesInOrder,1])
  c = (b-a)/length(leafNodesInOrder)
  layoutMat[leafNodesInOrder,1] = seq(a-10*c, b+10*c, length=length(leafNodesInOrder))

  return(layoutMat)
}

rRimStructureHelper = function(range, thetaRange = c(.1, 1),
                             distinctChildren = T) {
  if (thetaRange[2] - thetaRange[1] < .25 && distinctChildren) {
    stop("distinctChildren is true but theta range is smaller than .25")
  }
  n = length(range)
  if (n == 1) {
    return(t(c(0,0,1,0,range[1])))
  }

  k = sample(1:(length(range) - 1),1)
  leftMat = rRimStructureHelper(range[1:k])
  rightMat = rRimStructureHelper(range[(k + 1):n])

  leftMat[,1:2] = leftMat[,1:2] + 1
  rightMat[,1:2] = rightMat[,1:2] + 1 + nrow(leftMat)

  rtheta = runif(1, thetaRange[1], thetaRange[2])
  while (distinctChildren && (abs(rtheta - leftMat[,4]) < .1 ||
                                   abs(rtheta - rightMat[,4]) < .1)) {
    rtheta = runif(1, thetaRange[1], thetaRange[2])
  }

  mat = rbind(t(c(1, nrow(leftMat) + 1, 0, rtheta, 0)), leftMat, rightMat)
  return(mat)
}

rRIMStructure = function(numLeaves) {
  range = 0:(numLeaves - 1)
  mat = rRimStructureHelper(range)
  for(i in 1:nrow(mat)) {
    if(mat[i, 3] == 1) {
      mat[i, 1:2] = 0
    }
  }
  return(mat)
}
