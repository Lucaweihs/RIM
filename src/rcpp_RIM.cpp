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

/***
 * Generates a RIMTree from an input list.
 *
 * @param rimNodeList a Rcpp list that encodes the tree
 *        structure. In particular the list contains numeric
 *        vectors all of length 5. Each vector corresponds to
 *        a tree node and it's position in the list corresponds
 *        to the node's numbering. Each vector should have the form
 *        c(left child # + 1, right child # + 1, is leaf boolean, theta, rank)
 *        where
 *        - left child # = the number correspond to the left child of the node
 *        - right child # = the number correspond to the right child of the node
 *        - is leaf boolean = 1 if the node is a leaf, 0 otherwise
 *        - theta = the theta value corresponding to the node (for internal
 *          nodes)
 *        - rank = the rank of the node (for leaves)
 *        left/right child #s and theta are ignored if the node is a leaf and
 *        rank is ignored if the node is not a leaf. Note that the list should
 *        be  ordered topologically so that PARENTS ALWAYS COME BEFORE CHILDREN.
 *        As an example the list:
 *        list(
 *          c(2, 3, 0, -.1, 0),
 *          c(4, 5, 0, .8, 0),
 *          c(6, 7, 0, 1.6, 0),
 *          c(0, 0, 1, 0, 0),
 *          c(0, 0, 1, 0, 1),
 *          c(0, 0, 1, 0, 2),
 *          c(0, 0, 1, 0, 3))
 *        corresponds to a tree with a root node 0, the root having left and right
 *        children as 2,3 respectively. Moreover node 2 has children 3,4 which are
 *        both leaves and node 4 has children 5,6 which are also both leaves. If
 *        traversing the tree in preorder one visits the leaves in the order 3,4,5,6
 *        which corresponds to the ranking 0,1,2,3.
 * @return a pointer to a RIM::RIMTree corresponding to rimNodeList.
 ***/
RIM::RIMTree* RIMTreeFromList(List rimNodesList) {
  int numNodesTotal = rimNodesList.length();
  RIM::RIMNode* rimNodes[numNodesTotal];
  bool isLeaf;
  int leftChildIndex, rightChildIndex, rank;
  double theta;
  int numLeafNodes = (numNodesTotal+1)/2; // Since #nodes in binary tree is 2*numLeaves-1
  NumericVector y;
  RIM::RIMNode* newNode;

  // A loop that goes through the list and constructs the
  // RIM iteratively. Starts from the end of the list and moves
  // to the front so that children are visited before parents.
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
      if(leftChildIndex <= i || rightChildIndex <= i) {
        printf("ERROR: RIMTreeFromList expects children to come after parents in input list.\n");
        std::exit(1);
      }
      newNode->attachLeft(rimNodes[leftChildIndex]);
      newNode->attachRight(rimNodes[rightChildIndex]);
    }
  }

  return(new RIM::RIMTree(rimNodes[0]));
}

/***
 * Samples from the RIM represented by rimNodesList.
 *
 * @param numSamplesVec the number of samples to draw.
 * @param rimNodesList see documentation for
 *        RIMTreeFromList for the expected form of rimNodesList
 * @return a NumericMatrix with numSamplesVec rows each row representing
 *         a single draw from the RIM.
 ***/
// [[Rcpp::export]]
NumericMatrix RCPPSampleFromRIM(NumericVector numSamplesVec, List rimNodesList) {
  int numNodesTotal = rimNodesList.length();
  int numLeafNodes = (numNodesTotal+1)/2; // Since #nodes in binary tree is 2*numLeaves-1

  RIM::RIMTree* tree = RIMTreeFromList(rimNodesList); // The tree to sample from

  int numSamples = numSamplesVec[0];

  NumericVector randsVector;
  RIM::List<double>* rands;
  NumericMatrix rankings(numSamples, numLeafNodes); // This will be the matrix of rankings returned
  RIM::List<int>* ranking;

  // Run a loop numSamples times and generate a single
  // sample on each iteration.
  for(int i=0; i < numSamples; i++) {
    // Generate the random numbers that will be used to
    // create the randon sample and save these numbers in
    // a RIM::List.
    randsVector = runif(numLeafNodes*(numLeafNodes - 1));
    rands = new RIM::List<double>();
    for(int j=0; j < numLeafNodes*(numLeafNodes - 1); j++) {
      rands->appendValue(randsVector[j]);
    }
    // Using these random numbers generate a random ranking
    // from the RIM.
    ranking = tree->randomRanking(rands);

    // Add the random ranking to the matrix to be returned.
    ranking->restart();
    for(int j=0; j<numLeafNodes; j++) {
      rankings(i,j) = ranking->currentValue();
      ranking->next();
    }
    delete ranking;
    delete rands;
  }

  delete tree;
  return rankings;
}

/***
 * TO BE DETERMINED IF NECESSARY
 ***/
//// [[Rcpp::RCPPStructByDP]]
//NumericVector RIMThetaMLEs(NumericMatrix samples, List rimNodesList) {
//  RIM::RIMTree* tree = RIMTreeFromList(rimNodesList);
//  int* samplesMat = (int*) malloc(sizeof(int)*samples.nrow()*samples.ncol());
//  for(int i=0; i<samples.nrow(); i++) {
//    for(int j=0; j<samples.ncol(); j++) {
//      samplesMat[i*samples.ncol() + j] = samples(i,j);
//    }
//  }
//  tree->mlThetaTree(samplesMat, samples.nrow()) ;
//  RIM::List<double>* preOrderThetasList = tree->preOrderThetasList();
//  NumericVector preOrderThetasVector(preOrderThetasList->length());
//  preOrderThetasList->restart();
//  for(int i=0; i<preOrderThetasList->length(); i++) {
//    preOrderThetasVector[i] = preOrderThetasList->currentValue();
//    preOrderThetasList->next();
//  }
//
//  free(samplesMat);
//  delete preOrderThetasList;
//  delete tree;
//  return(preOrderThetasVector);
//}

/***
 * Creates the averaged discrepancy matrix correpsonding to the input samples
 * with respect to the reference permutation being the identity ranking.
 *
 * @param samples a NumericMatrix of rankings, each row
 *        should be a complete ranking.
 * @return the discrepancy matrix as a NumericMatrix.
 ***/
// [[Rcpp::export]]
NumericMatrix RCPPAverageDiscMatrix(NumericMatrix samples) {
  int numLeaves = samples.ncol();

  // A matrix of integers that we will copy samples into, this is necessary
  // as RIM::RIMTree::discrepancyMatix expects this format.
  int* samplesMat = (int*) malloc(sizeof(int)*samples.nrow()*numLeaves);
  for(int i=0; i<samples.nrow(); i++) {
    for(int j=0; j<numLeaves; j++) {
      samplesMat[i*numLeaves + j] = samples(i,j);
    }
  }

  // Create the matrix and then put it back into a NumericMatrix form, also
  // need to divide each entry by the number of samples to average them.
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

/***
 * This is a helper function for refRankingBackPointersAndThetasToRIMTree
 * which is itself a helper function for StructByDPRIMTree. Certain quantities
 * are created by the dynamic programming algorithm StructByDPRIMTree and this
 * function helps translate this quantities into a RIM::RIMTree. The tree is
 * created recursively.
 *
 * @param n # leaves (i.e. the number of items being ranked)
 * @param refRanking an array of length n containing the reference ranking
 * @param curNode the current node whose children need to be created.
 * @param i along with j determine the current lower/upper bounds (respectively)
 *        of refRanking that are being considered.
 * @param j see i.
 * @param backPointers an nxn matrix of of backPointers created by
 *        StructByDPRIMTree.
 * @param thetas an nxn matrix of ML theta values created by StructByDPRIMTree.
 ***/
void refRankingBackPointersAndThetasToRIMTreeHelper(int n, int* refRanking, RIM::RIMNode* curNode, int i, int j, int* backPointers, double* thetas) {
  if(i > j) {
    printf("ERROR: i must be <= j in refRankingBackPointersAndThetasToRIMTreeHelper.");
    std::exit(1);
  } else if(i == j) {
    // Current node is a leaf, no more recursion to be done. Set the rank
    // of the current node to the appropriate value.
    curNode->rank = refRanking[i];
    return;
  }
  // Current node is an internal node. Set the theta value of the node,
  // use backPointers to determine the correct splitting of [i,j], and then
  // recurse.
  curNode->theta = thetas[i*n + j];
  int k = backPointers[i*n + j];
  RIM::RIMNode* leftNode = new RIM::RIMNode(0,0);
  RIM::RIMNode* rightNode = new RIM::RIMNode(0,0);
  refRankingBackPointersAndThetasToRIMTreeHelper(n, refRanking, leftNode, i, k, backPointers, thetas);
  refRankingBackPointersAndThetasToRIMTreeHelper(n, refRanking, rightNode, k + 1, j, backPointers, thetas);
  curNode->attachLeft(leftNode);
  curNode->attachRight(rightNode);
}

/***
 * This is a helper function for StructByDPRIMTree. Certain quantities
 * are created by the dynamic programming algorithm StructByDPRIMTree and this
 * function helps translate this quantities into a RIM::RIMTree. The tree is
 * created recursively.
 *
 * @param n # leaves (i.e. the number of items being ranked)
 * @param refRanking an array of length n containing the reference ranking
 * @param backPointers an nxn matrix of of backPointers created by
 *        StructByDPRIMTree.
 * @param thetas an nxn matrix of ML theta values created by StructByDPRIMTree.
 * @return a RIM::RIMTree corresponding to the ML tree found by
 *         StructByDPRIMTree.
 ***/
RIM::RIMTree* refRankingBackPointersAndThetasToRIMTree(int n, int* refRanking, int* backPointers, double* thetas) {
  // Set up the root node and then use the helper function to construct the tree.
  RIM::RIMNode* root = new RIM::RIMNode(0,0);
  refRankingBackPointersAndThetasToRIMTreeHelper(n, refRanking, root, 0, n - 1, backPointers, thetas);
  return(new RIM::RIMTree(root));
}

/***
 * Takes an input RIM tree and outputs a NumericMatrix representing the tree.
 * The formatting is exactly the same as the input list to RIMTreeFromList
 * except that the rows of the outputted matrix correspond to the entries of the
 * list from RIMTreeFromList.
 *
 * @param tree the input RIM tree to convert into a NumericMatrix.
 * @return a NumericMatrix representing the input tree.
 ***/
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

/***
 * An implementation of the StructByDP algorithm presented by Meek and Meila
 * (2014). Creates an ML RIM tree structure given the reference ranking (i.e.
 * the ordering of the leaves of the tree).
 *
 * @param aveDiscMatrix the average discrepancy matrix corresponding to data.
 *        This should be created by the function RCPPAverageDiscMatrix.
 * @param refRanking an integer array of length equaling the number of items
 *        being ranked (i.e. the number of rows/columns of aveDiscMatrix).
 * @param makeCanonical if true then the algorithmenforces that the output
 *        tree is in canonical form (no internal nodes are negative).
 * @return the ML RIM tree corresponding to the input discrepancy matrix given
 *         the known reference ranking.
 ***/
RIM::RIMTree* StructByDPRIMTree(NumericMatrix aveDiscMatrix, int* refRanking, bool makeCanonical) {
  int n = aveDiscMatrix.ncol();
  double* costMatrix = (double*) calloc(n*n, sizeof(double));
  int* backPointers = (int*) calloc(n*n, sizeof(int));
  double* thetas = (double*) calloc(n*n, sizeof(double));

  int m, L, R;
  double aveV = 0;
  double s;
  double theta;
  int indMin, indMax;
  // A direct translation of the StructByDP algorithm from the paper of
  // Meek and Meila (2014). See the paper for justification and explaination.
  for(int l=2; l <= n; l++) {
    for(int j=1; j <= n - l + 1; j++) {
      m = j + l - 1;
      costMatrix[(j-1)*n + (m-1)] = std::numeric_limits<double>::max();
      for(int k=j; k <= m-1; k++) {
        // Calculating average disc
        aveV = 0;
        for(int j1=j; j1 <= k; j1++) {
          for(int m1=k+1; m1 <= m; m1++) {
            if(refRanking[j1-1] < refRanking[m1-1]) {
              aveV += aveDiscMatrix(refRanking[m1-1], refRanking[j1-1]);
            } else {
              aveV += 1 - aveDiscMatrix(refRanking[j1-1], refRanking[m1-1]);
            }
          }
        }
        L = k - j + 1;
        R = m - k;
        // These values for theta should probably not be hard coded.
        theta = RIM::RIMTree::mlTheta(L, R, aveV, 0.1, 2000, .0001);
        s = costMatrix[(j-1)*n + (k-1)] + costMatrix[k*n + (m-1)] + RIM::RIMTree::score(L, R, aveV, theta);
        if(s < costMatrix[(j-1)*n+(m-1)]) {
          costMatrix[(j-1)*n + (m-1)] = s;
          backPointers[(j-1)*n + (m-1)] = k - 1; // Added -1 here
          thetas[(j-1)*n + (m-1)] = theta;
        }
      }
    }
  }

  // Convert the output of the above loop into a RIM tree.
  RIM::RIMTree* tree = refRankingBackPointersAndThetasToRIMTree(n, refRanking, backPointers, thetas);
  if(makeCanonical) {
    tree->transformToCanonical();
  }
  return(tree);
}

/***
 * An implementation of the StructByDP algorithm presented by Meek and Meila
 * (2014). Creates a NumericMatrix representing the ML RIM tree structure
 * given the reference ranking (i.e. the ordering of the leaves of the tree).
 * For how this matrix is formatted see RIMTreeToMatrix.
 *
 * @param aveDiscMatrix the average discrepancy matrix corresponding to data.
 *        This should be created by the function RCPPAverageDiscMatrix.
 * @param refRanking an integer array of length equaling the number of items
 *        being ranked (i.e. the number of rows/columns of aveDiscMatrix).
 * @param makeCanonical if true then the algorithmenforces that the output
 *        tree is in canonical form (no internal nodes are negative).
 * @return a NumericMatrix representing ML RIM tree corresponding to the input
 *         discrepancy matrix given the known reference ranking. See
 *         the function RIMTreeToMatrix for how this matrix is formatted.
 ***/
// [[Rcpp::export]]
NumericMatrix RCPPStructByDP(NumericMatrix aveDiscMatrix, NumericVector refRanking, bool makeCanonical) {
  // Convert in input NumericVector refRanking to an int* as this is what
  // is expected by StructByDPRIMTree.
  int* refRankingArray = (int*) malloc(sizeof(int)*aveDiscMatrix.ncol());
  for(int i=0; i < aveDiscMatrix.ncol(); i++) {
    refRankingArray[i] = refRanking[i];
  }
  RIM::RIMTree* tree = StructByDPRIMTree(aveDiscMatrix, refRankingArray, makeCanonical);
  NumericMatrix treeMatrix = RIMTreeToMatrix(tree);
  delete tree;
  return(treeMatrix);
}

/***
 * Samples one ranking from the input RIM tree and inserts this ranking
 * into the input array rankingArray.
 *
 * @param tree the RIM tree to sample from.
 * @param rankingArray an int array of size the number of items being ranked
 *        (i.e. the number of leaves of the tree) which will be modified so that
 *        it corresponds to the ranking samples from the tree.
 ***/
void sampleOneFromRIMTree(RIM::RIMTree* tree, int* rankingArray) {
  int numLeafNodes = tree->getNumLeaves();
  RIM::List<int>* ranking;
  // In order to get a random ranking we must pass in an array of random values,
  // we create this array here.
  RIM::List<double>* rands = new RIM::List<double>();
  NumericVector randsVector = runif(numLeafNodes*(numLeafNodes - 1));
  for(int j=0; j < numLeafNodes*(numLeafNodes - 1); j++) {
    rands->appendValue(randsVector[j]);
  }
  // Create the random ranking.
  ranking = tree->randomRanking(rands);
  ranking->restart();

  // Insert the ranking into the array.
  for(int i=0; i < ranking->length(); i++) {
    rankingArray[i] = ranking->currentValue();
    ranking->next();
  }
  delete rands;
  delete ranking;
}

/***
 * Runs the SASearch algorithm of Meila and Meek (2004). Returns a NumericMatrix
 * representing the best tree found during the random search, (see the function
 * RIMTreeToMatrix for thow this matrix is formatted).
 *
 * @param aveDiscMatrix the average discrepancy matrix corresponding to data.
 *        This should be created by the function RCPPAverageDiscMatrix.
 * @param refRanking an integer array of length equaling the number of items
 *        being ranked (i.e. the number of rows/columns of aveDiscMatrix). This
 *        acts as the initial starting point of the search.
 * @param inverseTemp a parameter controlling the acceptance probability of new
 *        candidate rankings. Should be >=0 and a bigger value means less
 *        change of accepting candidate rankings.
 * @param maxIter the number of iterations to run the search.
 * @param makeCanonical if true then the algorithm enforces that all trees in
 *        the search are in canonical form (no internal nodes are negative).
 * @return a NumericMatrix representing the best RIM tree found during the
 *         found during the search and corresponding to the input discrepancy
 *         matrix See the function RIMTreeToMatrix for how this matrix is formatted.
 ***/
// [[Rcpp::export]]
NumericMatrix RCPPSASearch(NumericMatrix aveDiscMatrix, NumericVector refRanking, double inverseTemp, int maxIter, bool makeCanonical) {
  int numLeafNodes = aveDiscMatrix.nrow();
  // We will require the reference ranking to be a int vector so we put it in
  // that form now.
  int* ranking = (int*) malloc(sizeof(int)*numLeafNodes);
  for(int i=0; i < numLeafNodes; i++) {
     ranking[i] = refRanking[i];
  }

  // We will require that the aveDiscMatrix is a matrix of doubles, we transform
  // it into this form.
  double* aveDiscMatrixAsDouble = (double*) malloc(sizeof(double)*numLeafNodes*numLeafNodes);
  for(int i=0; i < numLeafNodes; i++) {
    for(int j=0; j < numLeafNodes; j++) {
      aveDiscMatrixAsDouble[i*numLeafNodes + j] = aveDiscMatrix(i,j);
    }
  }

// Create the first tree from the input ranking.
  RIM::RIMTree* curTree = StructByDPRIMTree(aveDiscMatrix, ranking, makeCanonical);
  RIM::RIMTree* tmpTree;
  RIM::RIMTree* bestTree = curTree;

  double lastLogProb = curTree->logProbability(aveDiscMatrixAsDouble);
  double bestLogProb = lastLogProb;
  double tmpLogProb;
  double rand;
  // The loop here is a direct translation of the algorithm SASearch of Meek and
  // Meila (2014), see the paper for more detail.
  for(int t=1; t <= maxIter; t++) {
    while(true) {
      sampleOneFromRIMTree(curTree, ranking);
      tmpTree = StructByDPRIMTree(aveDiscMatrix, ranking, makeCanonical);
      free(ranking);
      ranking = tmpTree->refRankingArray();
      delete tmpTree;
      tmpTree = StructByDPRIMTree(aveDiscMatrix, ranking, makeCanonical);
      tmpLogProb = tmpTree->logProbability(aveDiscMatrixAsDouble);
      rand = (runif(1))[0];
      if(exp(-inverseTemp*(lastLogProb - tmpLogProb)) > rand) {
        if(bestTree != curTree) {
          // If curTree was the best tree then don't delete it now, it will
          // be deleted below.
          delete curTree;
        }
        curTree = tmpTree;
        lastLogProb = tmpLogProb;
        break;
      }
    }
    if(bestLogProb < lastLogProb) {
      delete bestTree;
      bestTree = curTree;
      bestLogProb = lastLogProb;
    }
  }
  free(ranking);
  free(aveDiscMatrixAsDouble);

  NumericMatrix treeMatrix = RIMTreeToMatrix(bestTree);
  if(bestTree == curTree) {
    delete bestTree;
  } else {
    delete curTree;
    delete bestTree;
  }
  return(treeMatrix);
}

/***
 * Gets the log probability of a tree given some average discrepancy matrix
 * corresponding to the data.
 *
 * @param aveDiscMatrix the average discrepancy matrix corresponding to data.
 *        This should be created by the function RCPPAverageDiscMatrix.
 * @param rimNodeList see the function RIMTreeFromList for the appropriate
 *        formatting of this list.
 * @return the log probability of the tree given the data.
 ***/
// [[Rcpp::export]]
double RCPPLogProbOfTree(List rimNodesList, NumericMatrix aveDiscMatrix) {
  // We require that the aveDiscMatrix is a matrix of doubles, we transform
  // it into this form.
  int numLeafNodes = aveDiscMatrix.nrow();
  double* aveDiscMatrixAsDouble = (double*) malloc(sizeof(double)*numLeafNodes*numLeafNodes);
  for(int i=0; i < numLeafNodes; i++) {
    for(int j=0; j < numLeafNodes; j++) {
      aveDiscMatrixAsDouble[i*numLeafNodes + j] = aveDiscMatrix(i,j);
    }
  }
  RIM::RIMTree* tree = RIMTreeFromList(rimNodesList);
  double logProb = tree->logProbability(aveDiscMatrixAsDouble);

  delete tree;
  free(aveDiscMatrixAsDouble);
  return(logProb);
}




