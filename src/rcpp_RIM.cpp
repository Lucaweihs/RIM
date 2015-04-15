#include <Rcpp.h>
#include"RIMTree.cpp"
#include"List.cpp"
using namespace Rcpp;

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
    //printf("bar?\n");
    ranking = tree->randomRanking(rands);
    //printf("huh?\n");
    
    ranking->restart();
    for(int j=0; j<numLeafNodes; j++) {
      rankings(i,j) = ranking->currentValue();
      ranking->next();
    }
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
  return(preOrderThetasVector);
}






