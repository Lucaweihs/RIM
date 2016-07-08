treeMatToList = function(treeMat) {
  treeList = vector("list", length=nrow(treeMat))
  for(i in 1:nrow(treeMat)) {
    treeList[[i]] = treeMat[i,]
  }
  return(treeList)
}

treeListToMat = function(treeList) {
  treeMat = matrix(numeric(5 * length(treeList)), ncol = 5)
  for (i in 1:nrow(treeMat)) {
    treeMat[i,] = treeList[[i]]
  }
  return(treeMat)
}

treeMatrixToStringHelper = function(treeMatrix, ind) {
  if(treeMatrix[ind,3] == 1) {
    return(as.character(treeMatrix[ind,5]))
  }
  left = treeMatrixToStringHelper(treeMatrix, treeMatrix[ind,1] + 1)
  right = treeMatrixToStringHelper(treeMatrix, treeMatrix[ind,2] + 1)
  return(paste("(", left, ")(", right, ")", sep=""))
}

treeMatrixToString = function(treeMatrix) {
  return(treeMatrixToStringHelper(treeMatrix, 1))
}
