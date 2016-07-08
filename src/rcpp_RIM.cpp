/***
 * Copyright (C) 2016 Luca Weihs
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#include "List.h"
#include "RIMTree.h"
#include <RcppArmadillo.h>

using namespace Rcpp;

/***
 * Generates a RIMTree from an input list.
 *
 * @param rimNodeMat a matrix that encodes the tree
 *        structure. In particular the matrix contains must be of dimension
 *        (number of nodes in RIM) x 5. Each row of the matrix corresponds to
 *        a tree node and it's position in the list corresponds
 *        to the node's numbering. Each vector should have the form
 *        c(left child #, right child #, is leaf boolean, theta, rank)
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
 *        matrix(c(1, 2, 0, -.1, 0,
 *                 3, 4, 0, 0.8, 0,
 *                 5, 6, 0, 1.6, 0,
 *                 0, 0, 1, 0.0, 0,
 *                 0, 0, 1, 0.0, 1,
 *                 0, 0, 1, 0.0, 2,
 *                 0, 0, 1, 0.0, 3), ncol = 5)
 *        corresponds to a tree with a root node 0, the root having left and right
 *        children as 1,2 respectively. Moreover node 1 has children 3,4 which are
 *        both leaves and node 3 has children 5,6 which are also both leaves. If
 *        traversing the tree in preorder one visits the leaves in the order
 *        3,4,5,6 which corresponds to the ranking 0,1,2,3.
 * @return a pointer to a RIM::RIMTree corresponding to rimNodeMat.
 ***/
RIM::RIMTree* RIMTreeFromMat(arma::mat rimNodesMat) {
  int numNodesTotal = rimNodesMat.n_rows;
  RIM::RIMNode* rimNodes[numNodesTotal];
  bool isLeaf;
  int leftChildIndex, rightChildIndex, rank;
  double theta;
  int numLeafNodes = (numNodesTotal + 1) / 2; // Since #nodes in binary tree is 2*numLeaves-1
  arma::vec y;
  RIM::RIMNode* newNode;

  // A loop that goes through the list and constructs the
  // RIM iteratively. Starts from the end of the list and moves
  // to the front so that children are visited before parents.
  for(int i = numNodesTotal - 1; i >= 0; i--) {
    y = rimNodesMat.row(i).t();
    leftChildIndex = y[0];
    rightChildIndex = y[1];
    isLeaf = (y[2] == 1);
    theta = y[3];
    rank = y[4];

    newNode = new RIM::RIMNode(theta, rank);
    rimNodes[i] = newNode;
    if(!isLeaf) {
      if(leftChildIndex <= i || rightChildIndex <= i) {
        Rcpp::stop("ERROR: RIMTreeFromMat expects children to come after "
                     "parents in input list.\n");
      }
      newNode->attachLeft(rimNodes[leftChildIndex]);
      newNode->attachRight(rimNodes[rightChildIndex]);
    }
  }

  return(new RIM::RIMTree(rimNodes[0]));
}

/***
 * Samples from the RIM represented by rimNodesMat.
 *
 * @param numSamplesVec the number of samples to draw.
 * @param rimNodesList see documentation for
 *        RIMTreeFromList for the expected form of rimNodesList
 * @return a NumericMatrix with numSamplesVec rows each row representing
 *         a single draw from the RIM.
 ***/
// [[Rcpp::export]]
arma::mat RCPPSampleFromRIM(NumericVector numSamplesVec, arma::mat rimNodesMat) {
  int numNodesTotal = rimNodesMat.n_rows;
  int numLeafNodes = (numNodesTotal + 1) / 2; // Since #nodes in binary tree is 2*numLeaves-1
  RIM::RIMTree* tree = RIMTreeFromMat(rimNodesMat); // The tree to sample from

  int numSamples = numSamplesVec[0];
  arma::mat rankings(numSamples, numLeafNodes); // This will be the matrix of rankings returned
  arma::ivec ranking;

  // Run a loop numSamples times and generate a single
  // sample on each iteration.
  for(int i = 0; i < numSamples; i++) {
    // Generate a random ranking from the RIM.
    ranking = tree->randomRanking();
    // Add the random ranking to the matrix to be returned.
    for(int j = 0; j < numLeafNodes; j++) {
      rankings(i,j) = ranking(j);
    }
  }

  delete tree;
  return rankings;
}

/***
 * Creates the averaged discrepancy matrix correpsonding to the input samples
 * with respect to the reference permutation being the identity ranking.
 *
 * @param samples a NumericMatrix of rankings, each row should be a complete
 *        ranking.
 * @return the discrepancy matrix as a arma::mat.
 ***/
// [[Rcpp::export]]
arma::mat RCPPAverageDiscMatrix(arma::imat samples) {
  int numLeaves = samples.n_cols;

  arma::imat discMat = RIM::RIMTree::discrepancyMatix(samples);
  arma::mat aveDiscMat(numLeaves, numLeaves);
  for(int i = 0; i < numLeaves; i++) {
    for(int j = 0; j < numLeaves; j++) {
      aveDiscMat(i,j) = (1.0 * discMat(i, j)) / samples.n_rows;
    }
  }
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
void refRankingBackPointersAndThetasToRIMTreeHelper(const arma::ivec& refRanking,
                                                    RIM::RIMNode* curNode,
                                                    int i,
                                                    int j,
                                                    const arma::imat& backPointers,
                                                    const arma::mat& thetas) {
  if(i > j) {
    Rcpp::stop("ERROR: i must be <= j in refRankingBackPointersAndThetasToRIMTreeHelper.");
  } else if(i == j) {
    // Current node is a leaf, no more recursion to be done. Set the rank
    // of the current node to the appropriate value.
    curNode->rank = refRanking[i];
    return;
  }
  // Current node is an internal node. Set the theta value of the node,
  // use backPointers to determine the correct splitting of [i,j], and then
  // recurse.
  curNode->theta = thetas(i, j);
  int k = backPointers(i, j);
  RIM::RIMNode* leftNode = new RIM::RIMNode(0,0);
  RIM::RIMNode* rightNode = new RIM::RIMNode(0,0);
  refRankingBackPointersAndThetasToRIMTreeHelper(refRanking, leftNode, i, k, backPointers, thetas);
  refRankingBackPointersAndThetasToRIMTreeHelper(refRanking, rightNode, k + 1, j, backPointers, thetas);
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
RIM::RIMTree* refRankingBackPointersAndThetasToRIMTree(
    int n, const arma::ivec& refRanking, const arma::imat& backPointers,
    const arma::mat& thetas) {
  // Set up the root node and then use the helper function to construct the tree.
  RIM::RIMNode* root = new RIM::RIMNode(0,0);
  refRankingBackPointersAndThetasToRIMTreeHelper(refRanking, root, 0, n - 1, backPointers, thetas);
  return(new RIM::RIMTree(root));
}

/***
 * Takes an input RIM tree and outputs a matrix representing the tree.
 * The formatting is exactly the same as the input list to RIMTreeFromMat.
 *
 * @param tree the input RIM tree to convert into a matrix.
 * @return a matrix representing the input tree.
 ***/
arma::mat RIMTreeToMatrix(RIM::RIMTree* tree) {
  return(tree->treeToMatrix());
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
 * @param makeCanonical if true then the algorithm forces that the output
 *        tree is in canonical form (no internal nodes are negative).
 * @return the ML RIM tree corresponding to the input discrepancy matrix given
 *         the known reference ranking.
 ***/
RIM::RIMTree* StructByDPRIMTree(
    arma::mat aveDiscMatrix, const arma::ivec& refRanking, bool makeCanonical) {

  int n = aveDiscMatrix.n_cols;

  arma::mat costMatrix = arma::zeros<arma::mat>(n, n);
  arma::imat backPointers = arma::zeros<arma::imat>(n, n);
  arma::mat thetas = arma::zeros<arma::mat>(n, n);

  int m, L, R;
  double aveV = 0;
  double s;
  double theta;
  int indMin, indMax;
  // A direct translation of the StructByDP algorithm from the paper of
  // Meek and Meila (2014). See the paper for justification and explaination.
  for (int l = 2; l <= n; l++) {
    for (int j = 1; j <= n - l + 1; j++) {
      m = j + l - 1;
      costMatrix(j - 1, m - 1) = std::numeric_limits<double>::max();
      for (int k = j; k <= m - 1; k++) {
        // Calculating average disc
        aveV = 0;
        for(int j1 = j; j1 <= k; j1++) {
          for(int m1 = k + 1; m1 <= m; m1++) {
            if(refRanking[j1 - 1] < refRanking[m1 - 1]) {
              aveV += aveDiscMatrix(refRanking[m1 - 1], refRanking[j1 - 1]);
            } else {
              aveV += 1 - aveDiscMatrix(refRanking[j1 - 1], refRanking[m1 - 1]);
            }
          }
        }
        L = k - j + 1;
        R = m - k;
        // These values for theta should probably not be hard coded.
        theta = RIM::RIMTree::mlTheta(L, R, aveV, 0, 200, .00001);
        s = costMatrix(j - 1, k - 1) + costMatrix(k, m - 1) +
          RIM::RIMTree::score(L, R, aveV, theta);
        if(s < costMatrix(j - 1, m - 1)) {
          costMatrix(j - 1, m - 1) = s;
          backPointers(j - 1, m - 1) = k - 1; // Added -1 here
          thetas(j - 1, m - 1) = theta;
        }
      }
    }
  }

  // Convert the output of the above loop into a RIM tree.
  RIM::RIMTree* tree = refRankingBackPointersAndThetasToRIMTree(n, refRanking,
                                                                backPointers,
                                                                thetas);
  if(makeCanonical) {
    tree->transformToCanonical();
  }
  return(tree);
}

/***
 * An implementation of the StructByDP algorithm presented by Meek and Meila
 * (2014). Creates a matrix representing the ML RIM tree structure
 * given the reference ranking (i.e. the ordering of the leaves of the tree).
 * For how this matrix is formatted see RIMTreeToMatrix.
 *
 * @param aveDiscMatrix the average discrepancy matrix corresponding to data.
 *        This should be created by the function RCPPAverageDiscMatrix.
 * @param refRanking an integer array of length equaling the number of items
 *        being ranked (i.e. the number of rows/columns of aveDiscMatrix).
 * @param makeCanonical if true then the algorithmenforces that the output
 *        tree is in canonical form (no internal nodes are negative).
 * @return a matrix representing ML RIM tree corresponding to the input
 *         discrepancy matrix given the known reference ranking. See
 *         the function RIMTreeToMatrix for how this matrix is formatted.
 ***/
// [[Rcpp::export]]
arma::mat RCPPStructByDP(arma::mat aveDiscMatrix, arma::ivec refRanking,
                             bool makeCanonical) {
  RIM::RIMTree* tree = StructByDPRIMTree(aveDiscMatrix, refRanking, makeCanonical);
  return(RIMTreeToMatrix(tree));
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
arma::mat RCPPSASearch(arma::mat aveDiscMatrix, arma::ivec refRanking,
                           double inverseTemp, int maxIter, bool makeCanonical,
                           bool verbose) {
  int numLeafNodes = aveDiscMatrix.n_rows;

  // Create the first tree from the input ranking.
  RIM::RIMTree* curTree = StructByDPRIMTree(aveDiscMatrix, refRanking, makeCanonical);
  RIM::RIMTree* tmpTree;
  RIM::RIMTree* bestTree = curTree;

  double lastLogProb = curTree->logProbability(aveDiscMatrix);
  double bestLogProb = lastLogProb;
  double tmpLogProb;
  double rand;
  int runsBeforeAcceptance = 1;
  // The loop here is a direct translation of the algorithm SASearch of Meek and
  // Meila (2014), see the paper for more detail.
  for (int t = 1; t <= maxIter; t++) {
    if (verbose && (t % 100) == 0) {
      Rprintf("t=%d\n", t);
      Rprintf("Average runs before acceptance in last 100=%f\n", runsBeforeAcceptance/100.0);
      runsBeforeAcceptance = 0;
    }
    while (true) {
      runsBeforeAcceptance++;
      arma::ivec ranking = curTree->randomRanking();
      tmpTree = StructByDPRIMTree(aveDiscMatrix, ranking, makeCanonical);
      ranking = tmpTree->refRankingArray();
      tmpTree = StructByDPRIMTree(aveDiscMatrix, ranking, makeCanonical);
      tmpLogProb = tmpTree->logProbability(aveDiscMatrix);
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
      } else {
        delete tmpTree;
      }
    }
    if(bestLogProb < lastLogProb) {
      delete bestTree;
      bestTree = curTree;
      bestLogProb = lastLogProb;
    }
  }

  arma::ivec ranking = bestTree->refRankingArray();
  tmpTree = StructByDPRIMTree(aveDiscMatrix, ranking, makeCanonical);

  arma::mat treeMatrix = RIMTreeToMatrix(tmpTree);
  delete tmpTree;
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
 * @param rimNodeMat see the function RIMTreeFromMat for the appropriate
 *        formatting of this matrix.
 * @return the log probability of the tree given the data.
 ***/
// [[Rcpp::export]]
double RCPPLogProbRIM(arma::mat rimNodesMat, arma::mat aveDiscMatrix) {
  // We require that the aveDiscMatrix is a matrix of doubles, we transform
  // it into this form.
  int numLeafNodes = aveDiscMatrix.n_rows;
  RIM::RIMTree* tree = RIMTreeFromMat(rimNodesMat);
  double logProb = tree->logProbability(aveDiscMatrix);

  delete tree;
  return(logProb);
}

/***
 * Takes an input averaged discrepancy matrix (from the data), and a RIM as a matrix
 * (see RIMTreeFromMat for how this matrix should be formatted) and outputs a
 * matrix of the tree with theta values equalling the MLE theta values
 * for the data.
 *
 * @param rimNodesMat a matrix representing the RIM tree.
 * @param aveDiscMatrix an averaged discrepancy matrix from the data.
 * @return a NumericMatrix representing the RIM tree with optimized theta values.
 ***/
// [[Rcpp::export]]
arma::mat RCPPthetaMLERIM(arma::mat rimNodesMat, arma::mat aveDiscMatrix) {

  RIM::RIMTree* tree = RIMTreeFromMat(rimNodesMat);
  tree->mlThetaTree(aveDiscMatrix);

  arma::mat treeMat = RIMTreeToMatrix(tree);
  delete tree;
  return(treeMat);
}
