#include <Rcpp.h>
#include"RIMTree.cpp"
#include"List.cpp"
using namespace Rcpp;

void printIntMatrix(int* mat, int nrow, int ncol) {
  for(int i=0; i<nrow; i++) {
    for(int j=0; j<ncol; j++) {
      printf("%d ", mat[i*ncol + j]);
    }
    printf("\n");
  }
}
void printDoubleMatrix(double* mat, int nrow, int ncol) {
  for(int i=0; i<nrow; i++) {
    for(int j=0; j<ncol; j++) {
      printf("%f ", mat[i*ncol + j]);
    }
    printf("\n");
  }
}

RIM::RIMTree* RIMTreeFromList(List rimNodesList) {
  int numNodesTotal = rimNodesList.length();
  RIM::RIMNode* rimNodes[numNodesTotal];
  bool isLeaf;
  int leftChildIndex, rightChildIndex, rank;
  double theta;
  int numLeafNodes = (numNodesTotal+1)/2;
  NumericVector y;
  RIM::RIMNode* newNode;
  
  for(int i=numNodesTotal-1; i >= 0; i--) {
    y = rimNodesList[i];
    leftChildIndex = y[0] - 1;
    rightChildIndex = y[1] - 1;
    isLeaf = (y[2] == 1);
    theta = y[3];
    rank = y[4];

    newNode = new RIM::RIMNode(theta, rank);
    rimNodes[i] = newNode;
    if(!isLeaf) {
      newNode->attachLeft(rimNodes[leftChildIndex]);
      newNode->attachRight(rimNodes[rightChildIndex]);
    }
  }
  
  return(new RIM::RIMTree(rimNodes[0]));
}

//RIM::RIMTree* RIMTreeFromMatrix(NumericMatrix rimNodesMatrix) {
//  int numNodesTotal = rimNodesMatrix.nrow();
//  RIM::RIMNode* rimNodes[numNodesTotal];
//  bool isLeaf;
//  int leftChildIndex, rightChildIndex, rank;
//  double theta;
//  int numLeafNodes = (numNodesTotal+1)/2;
//  RIM::RIMNode* newNode;
//  
//  for(int i=numNodesTotal-1; i >= 0; i--) {
//    leftChildIndex = rimNodesMatrix(i,0) - 1;
//    rightChildIndex = rimNodesMatrix(i,1) - 1;
//    isLeaf = (rimNodesMatrix(i,2) == 1);
//    theta = rimNodesMatrix(i,3);
//    rank = rimNodesMatrix(i,0);
//
//    newNode = new RIM::RIMNode(theta, rank);
//    rimNodes[i] = newNode;
//    if(!isLeaf) {
//      newNode->attachLeft(rimNodes[leftChildIndex]);
//      newNode->attachRight(rimNodes[rightChildIndex]);
//    }
//  }
//  
//  return(new RIM::RIMTree(rimNodes[0]));
//}

// [[Rcpp::export]]
NumericMatrix sampleFromRIM(NumericVector numSamplesVec, List rimNodesList) {
  int numNodesTotal = rimNodesList.length();
  int numLeafNodes = (numNodesTotal+1)/2;
  
  RIM::RIMTree* tree = RIMTreeFromList(rimNodesList);
  
  int numSamples = numSamplesVec[0];

  NumericVector randsVector;
  RIM::List<double>* rands;
  NumericMatrix rankings(numSamples, numLeafNodes);
  RIM::List<int>* ranking;
  
  // New way of doing it
  for(int i=0; i < numSamples; i++) {
    randsVector = runif(numLeafNodes*(numLeafNodes - 1));
    rands = new RIM::List<double>();
    for(int j=0; j < numLeafNodes*(numLeafNodes - 1); j++) {
      rands->appendValue(randsVector[j]);
    }
    ranking = tree->randomRanking(rands);
    
    ranking->restart();
    for(int j=0; j<numLeafNodes; j++) {
      rankings(i,j) = ranking->currentValue();
      ranking->next();
    }
    delete ranking;
    delete rands;
  }
  
  /* // Old way of doing it
  for(int i=0; i < numSamples; i++) {
    randsVector = runif(2*(numLeafNodes - 1));
    rands = new RIM::List<double>();
    for(int j=0; j < 2*(numLeafNodes - 1); j++) {
      rands->appendValue(randsVector[j]);
    }
    ranking = tree->randomRanking(rands);
    
    ranking->restart();
    for(int j=0; j<numLeafNodes; j++) {
      rankings(i,j) = ranking->currentValue();
      ranking->next();
    }
  }*/
  delete tree;
  return rankings;
}

// [[Rcpp::export]]
NumericVector RIMThetaMLEs(NumericMatrix samples, List rimNodesList) {
  RIM::RIMTree* tree = RIMTreeFromList(rimNodesList);
  int* samplesMat = (int*) malloc(sizeof(int)*samples.nrow()*samples.ncol());
  for(int i=0; i<samples.nrow(); i++) {
    for(int j=0; j<samples.ncol(); j++) {
      samplesMat[i*samples.ncol() + j] = samples(i,j);
    }
  }
  tree->mlThetaTree(samplesMat, samples.nrow()) ;
  RIM::List<double>* preOrderThetasList = tree->preOrderThetasList();
  NumericVector preOrderThetasVector(preOrderThetasList->length());
  preOrderThetasList->restart();
  for(int i=0; i<preOrderThetasList->length(); i++) {
    preOrderThetasVector[i] = preOrderThetasList->currentValue();
    preOrderThetasList->next();
  }
  
  free(samplesMat);
  delete preOrderThetasList;
  delete tree;
  return(preOrderThetasVector);
}

// [[Rcpp::export]]
NumericMatrix averageDiscMatrix(NumericMatrix samples) {
  int numLeaves = samples.ncol();
  int* samplesMat = (int*) malloc(sizeof(int)*samples.nrow()*numLeaves);
  for(int i=0; i<samples.nrow(); i++) {
    for(int j=0; j<numLeaves; j++) {
      samplesMat[i*numLeaves + j] = samples(i,j);
    }
  }
  
  int* discMat = RIM::RIMTree::discrepancyMatix(samplesMat, samples.nrow(), numLeaves);
  NumericMatrix aveDiscMat(numLeaves,numLeaves);
  for(int i=0; i < numLeaves; i++) {
    for(int j=0; j < numLeaves; j++) {
      aveDiscMat(i,j) = (1.0*discMat[i*numLeaves + j])/samples.nrow();
    }
  }
  free(discMat);
  free(samplesMat);
  return(aveDiscMat);
}

void costAndBackPointersToRIMTreeHelper(int n, RIM::RIMNode* curNode, int i, int j, int* backPointers, double* thetas) {
  if(i == j) {
    curNode->rank = i;
    return;
  }
  curNode->theta = thetas[i*n + j];
  int k = backPointers[i*n + j];
  RIM::RIMNode* leftNode = new RIM::RIMNode(0,0);
  RIM::RIMNode* rightNode = new RIM::RIMNode(0,0);
  costAndBackPointersToRIMTreeHelper(n, leftNode, i, k, backPointers, thetas);
  costAndBackPointersToRIMTreeHelper(n, rightNode, k + 1, j, backPointers, thetas);
  curNode->attachLeft(leftNode);
  curNode->attachRight(rightNode);
}

RIM::RIMTree* backPointersAndThetasToRIMTree(int n, int* backPointers, double* thetas) {
  RIM::RIMNode* root = new RIM::RIMNode(0,0);
  costAndBackPointersToRIMTreeHelper(n, root, 0, n - 1, backPointers, thetas);
  return(new RIM::RIMTree(root));
}

NumericMatrix RIMTreeToMatrix(RIM::RIMTree* tree) {
  double* treeMatrix = tree->treeToMatrix();
  int numLeaves = tree->getNumLeaves();
  NumericMatrix treeMatrixNumeric(2*numLeaves - 1, 5);
  for(int i=0; i < 2*numLeaves - 1; i++) {
    for(int j=0; j < 5; j++) {
      treeMatrixNumeric(i,j) = treeMatrix[i*5 + j];
    }
  }
  return(treeMatrixNumeric);
}

// [[Rcpp::export]]
NumericMatrix StructByDP(NumericMatrix aveDiscMatrix, bool makeCanonical) {
  int n = aveDiscMatrix.ncol();
  double* costMatrix = (double*) calloc(n*n, sizeof(double));
  int* backPointers = (int*) calloc(n*n, sizeof(int));
  double* thetas = (double*) calloc(n*n, sizeof(double));
  
  int m, L, R;
  double aveV = 0;
  double s;
  double theta;
  for(int l=2; l <= n; l++) {
    for(int j=1; j <= n - l + 1; j++) {
      m = j + l - 1;
      costMatrix[(j-1)*n + (m-1)] = std::numeric_limits<double>::max();
      for(int k=j; k <= m-1; k++) {
        // Calculating average disc
        aveV = 0;
        for(int j1=j; j1 <= k; j1++) {
          for(int m1=k+1; m1 <= m; m1++) {
            aveV += aveDiscMatrix(m1-1, j1-1);
          }
        }
        L = k - j + 1;
        R = m - k;
        theta = RIM::RIMTree::mlTheta(L, R, aveV, 0.1, 1000, .0000001);
        s = costMatrix[(j-1)*n + (k-1)] + costMatrix[k*n + (m-1)] + RIM::RIMTree::score(L, R, aveV, theta);
        //printf("k=%d, L=%d, R=%d, aveV=%f, theta=%f, score=%f, gaussianPoly=%f, s=%f \n", k, L, R, aveV, theta, RIM::RIMTree::score(L, R, aveV, theta), RIM::RIMTree::gaussianPoly(L, R, theta), s);
        if(s < costMatrix[(j-1)*n+(m-1)]) {
          costMatrix[(j-1)*n + (m-1)] = s;
          backPointers[(j-1)*n + (m-1)] = k - 1; // Added -1 here
          thetas[(j-1)*n + (m-1)] = theta;
        }
      }
    }
  }
    
  RIM::RIMTree* tree = backPointersAndThetasToRIMTree(n, backPointers, thetas);
  
  if(makeCanonical) {
    tree->transformToCanonical();
  }
  NumericMatrix treeMatrix = RIMTreeToMatrix(tree);
  delete tree;
  return(treeMatrix);
}

//NumericMatrix SASearch(NumericMatrix aveDiscMatrix, double inverseTemp, bool makeCanonical) {
  
//}





