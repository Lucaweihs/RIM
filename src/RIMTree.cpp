/*****
 * The main class that encodes all the necessary structure of a RIM tree model.
 * Includes all of the necessary functions to initialize, transform to
 * canonical form, find MLEs of theta values, generate samples, and query for
 * the log probability of the model given data.
 *****/

#ifndef RIM_RIMTree
#define RIM_RIMTree

#include<stdio.h>
#include<stdlib.h>
#include<cmath>
#include"RIMNode.cpp"

namespace RIM {
  class RIMTree {
    private:
      static const int MAX_ITER_DEFAULT = 5000; // Iterations when finding MLE for theta values
      static const double TOL_DEFAULT = .0001; // Tolerance in grad descent for finding MLE
                                                // theta values.
      RIMNode* root; // Root of the tree
      int* discMat; // Discrepancy matrix of the tree's reference permutation
                    // with respect to the identity permutation.
      int numLeaves; // Number of leaves in the tree.


      /***
       * A helper function for randomRanking. Randomly creates a partial ranking
       * at a current node given the partitial rankings recursively created at
       * the current nodes left and right children.
       *
       * @param currentNode the current node where the partial ranking should be
       *        constructed.
       * @param rands a list of random uniform(0,1) numbers which should have
       *        it's current pointer set to the next random number to be used.
       * @return a RIM::List which contains the randomly generated partial
       *         ranking at the current node.
       ***/
      RIM::List<int>* randomRankingHelper(RIM::RIMNode* currentNode, RIM::List<double>* rands) {
        if(currentNode->left == NULL && currentNode->right == NULL) {
          // If the current node is a leaf then the partial ranking is just
          // the current nodes rank. Return it as a single element in a list.
          RIM::List<int>* l = new List<int>();
          l->appendValue(currentNode->rank);
          return(l);
        }

        // Generate random partial rankings associated with the current nodes'
        // left and right children.
        RIM::List<int>* leftRanking = randomRankingHelper(currentNode->left, rands);
        RIM::List<int>* rightRanking = randomRankingHelper(currentNode->right, rands);
        int L = leftRanking->length();
        int R = rightRanking->length();

        // Generates a partial ranking by merging leftRanking and rightRanking
        // together in a generative fashion. See TODO for reference on why
        // this is appropriate.
        double* zMatrix = getZMatrix(L, R, currentNode->theta);
        double rand;
        int last = L;
        RIM::List<int>* partition = new RIM::List<int>();
        for(int i = R; i >= 1; i--) {
          rand = rands->currentValue();
          rands->next();
          for(int j = 0; j <= std::min(L,last); j++) {
            rand -= exp(-1.0*currentNode->theta*j)*zMatrix[j*(R+1) + (i-1)]/zMatrix[last*(R+1) + i];
            if(rand <= 0) {
              last = j;
              if(last != 0) {
                partition->appendValue(j);
              }
              break;
            }
          }
          if(last == 0) {
            break;
          }
        }
        free(zMatrix);

        leftRanking->joinWithPartition(rightRanking, partition);
        delete rightRanking;
        delete partition;
        return(leftRanking);
      }

      /***
       * A helper function for transformToCanonical. Recusively visits
       * a node and it's descendants negating theta and flipping subtrees
       * when theta is negative.
       *
       * @param node the current node that will have it's theta negated
       *        and subtrees flipped if it's theta is negative. The children
       *        of the node will then will have the procedure recursively applied
       *        to them.
       ***/
      void transformToCanonicalHelper(RIMNode* node) {
        if(node->left == NULL && node->right == NULL) {
          // If node is a leaf then there is nothing to do.
          return;
        }
        if(node->theta < 0) {
          node->flipSubTreesAndNegate();
        }
        transformToCanonicalHelper(node->left);
        transformToCanonicalHelper(node->right);
      }

      /***
       * A helper function for mlThetaTree. Recursively visits all internal
       * nodes of the tree and sets their theta value to be the MLE theta
       * corresponding to the discrepancy matrix aveDataDiscMat.
       *
       * @param curNode the current node that the ML theta value is being found
       *        for.
       * @param aveDataDiscMat the discrepancy matrix corresponding to the data.
       * @return a list of the reference ranking corresponding to the left and
       *         right subtrees of the curNode. This is used to recursively by
       *         the calls to mlThetaTreeHelper.
       ***/
      List<int>* mlThetaTreeHelper(RIMNode* curNode, double* aveDataDiscMat) {
        if(curNode->left == NULL && curNode->right == NULL) {
          // If curNode is a leaf then there is no theta to optimize instead,
          // just return the partial ranking corresponding to the leaf (i.e.
          // a list containing the rank of the leaf).
          RIM::List<int>* l = new List<int>();
          l->appendValue(curNode->rank);
          return(l);
        }

        // Recursively generate and set the ML theta values for the curNode's
        // left and right subtrees. Also get the partial rankings corresponding
        // to the left and right subtrees.
        RIM::List<int>* leftRanking = mlThetaTreeHelper(curNode->left, aveDataDiscMat);
        RIM::List<int>* rightRanking = mlThetaTreeHelper(curNode->right, aveDataDiscMat);

        // Loop over leftRanking and rightRanking to compute the sufficient
        // statistic aveDisc needed for the optimization. This is descibed in
        // Meek and Meila (2014).
        double aveDisc = 0;
        int rankLeft, rankRight;
        int minRankIndex, maxRankIndex, index;
        leftRanking->restart();
        for(int i=0; i<leftRanking->length(); i++) {
          rankLeft = leftRanking->currentValue();
          leftRanking->next();

          rightRanking->restart();
          for(int j=0; j<rightRanking->length(); j++) {
            rankRight = rightRanking->currentValue();
            rightRanking->next();

            minRankIndex = std::min(rankLeft,rankRight);
            maxRankIndex = std::max(rankLeft,rankRight);
            index = maxRankIndex*this->numLeaves + minRankIndex;

            if(this->discMat[index] != 0) {
              aveDisc += 1 - aveDataDiscMat[index];
            } else {
              aveDisc += aveDataDiscMat[index];
            }
          }
        }

        // Find the ML theta value and set curNode's theta to be the value.
        curNode->theta = mlTheta(leftRanking->length(), rightRanking->length(), aveDisc, 0, MAX_ITER_DEFAULT, TOL_DEFAULT);

        // Append the rightRanking to the end of leftRanking and then
        // return the left ranking.
        RIM::List<int> tmp = RIM::List<int>();
        leftRanking->joinWithPartition(rightRanking, &tmp);
        delete rightRanking;
        return(leftRanking);
      }

    public:
      /***
       * An initializer of RIMTree. Takes as input a root node which should
       * have the entire tree structure already associated with it.
       *
       * @param newRoot a RIM tree node that should have a tree structure
       *        already attached.
       * @return a RIM tree with newRoot as the root and with the tree structure
       *         attached to to newRoot.
       ***/
      RIMTree(RIMNode* newRoot) {
        root = newRoot;
        // Just a little work to set up the discrepancy matrix associated
        // to the tree.
        RIM::List<int> rankingList = RIM::List<int>();
        root->refRanking(&rankingList);
        int* rankingArray = refRankingListToArray(&rankingList);
        numLeaves = rankingList.length();
        discMat = discrepancyMatix(rankingArray, 1, numLeaves);
        free(rankingArray);
      }

      /***
       * Deletes the tree.
       ***/
      ~RIMTree() {
        root->deleteAncestors();
        delete root;
        free(discMat);
      }

      /***
       * Returns the number of leaves in the tree.
       *
       * @return number of leaves in the tree.
       ***/
      int getNumLeaves() {
        return(numLeaves);
      }

      /***
       * Returns the theta values of the internal nodes in preorder.
       *
       * @return a RIM::List of the theta values.
       ***/
      RIM::List<double>* preOrderThetasList() {
        RIM::List<double>* l = new RIM::List<double>();
        preOrderThetasListHelper(l, root);
        return(l);
      }

      /***
       * A helper function for preOrderThetasList.
       *
       * @param l a RIM::List of the theta values recorded so far.
       * @param curNode the current node of the tree.
       ***/
      void preOrderThetasListHelper(RIM::List<double>* l, RIM::RIMNode* curNode) {
        if(curNode->left == NULL && curNode->right == NULL) {
          return;
        }
        preOrderThetasListHelper(l, curNode->left);
        l->appendValue(curNode->theta);
        preOrderThetasListHelper(l, curNode->right);
      }

      /***
       * A function that, given the averaged discrepancy matrix of the data,
       * updates all of the theta values of the tree to be their MLE.
       *
       * @param aveDiscMat an averaged discrepancy matrix over data.
       ***/
      void mlThetaTree(double* aveDiscMat) {
        mlThetaTreeHelper(root, aveDiscMat);
      }

      /***
       * Given n*(n-1) random numbers in a RIM::List (n=number of leaves in the
       * RIM tree) returns a random sample from the RIM as a RIM::List.
       *
       * @param rands a list of n*(n-1) random uniform(0,1) numbers (n=number of
       *         leaves in the RIM tree).
       * @returns a RIM:List containing the random ranking.
       ***/
      RIM::List<int>* randomRanking(RIM::List<double>* rands) {
        rands->restart();
        if(this->numLeaves*(this->numLeaves-1) != rands->length()) {
          printf("\nERROR: Not enough random numbers inputted to randomRanking.");
          std::exit(1);
        }
        RIM::List<int>* l = randomRankingHelper(root, rands);
        return(l);
      }

      /***
       * Returns the reference ranking of the RIM tree as an int array.
       *
       * @return the reference ranking as an length numLeaves int array.
       ***/
      int* refRankingArray() {
        RIM::List<int> l = RIM::List<int>();
        root->refRanking(&l);
        return(refRankingListToArray(&l));
      }

      /***
       * Creates an (# nodes in tree)x5 matrix representing the RIM tree. See
       * the documentation for RIMTreeFromList in rcpp_RIM.cpp for how this
       * matrix is formatted.
       *
       * @return a matrix of doubles representing the tree.
       ***/
      double* treeToMatrix() {
        RIM::List<RIM::RIMNode*> nodeList = List<RIM::RIMNode*>();

        // The matrix to be returned.
        double* treeMat = (double*) malloc(sizeof(double)*(2*numLeaves-1)*5);
        for(int i=0; i<(2*numLeaves-1)*5; i++) {
          treeMat[i] = 0;
        }

        int k = 0;
        int length = 1;
        RIM::RIMNode* curNode;

        // Performs a breadth-first search of the tree and constructs the matrix
        // along the way. This ensures that parents will always come before all
        // of their children in the matrix.
        nodeList.appendValue(root);
        nodeList.restart();
        for(int k = 0; k < 2*numLeaves-1; k++) {
          curNode = nodeList.currentValue();
          if(curNode->left == NULL && curNode->right == NULL) {
            // If node is a leaf then things are simple.
            treeMat[k*5 + 0] = 0;
            treeMat[k*5 + 1] = 0;
            treeMat[k*5 + 2] = 1;
            treeMat[k*5 + 3] = 0;
            treeMat[k*5 + 4] = curNode->rank;
          } else {
            // If node is internal then make sure to set which nodes are it's
            // children carefully. Everything else is simple.
            treeMat[k*5 + 0] = length + 1;
            treeMat[k*5 + 1] = length + 2;
            treeMat[k*5 + 2] = 0;
            treeMat[k*5 + 3] = curNode->theta;
            treeMat[k*5 + 4] = 0;
            nodeList.appendValue(curNode->left);
            nodeList.appendValue(curNode->right);
            length += 2;
          }
          nodeList.next();
        }
        return(treeMat);
      }

      /***
       * Takes a reference ranking as a RIM list and outputs the reference
       * ranking as an int array.
       *
       * @param l the reference ranking as an array.
       * @return an int array version of the reference ranking.
       ***/
      static int* refRankingListToArray(RIM::List<int>* l) {
        int* refArray = (int*) malloc(sizeof(int)*l->length());
        l->restart();
        for(int i=0; i<l->length(); i++) {
          refArray[i] = l->currentValue();
          l->next();
        }
        return(refArray);
      }

      /***
       * Returns the reference ranking of the tree as a RIM list.
       *
       * @return a RIM list that has the reference ranking in order.
       ***/
      RIM::List<int>* refRankingList() {
        RIM::List<int>* l = new RIM::List<int>();
        root->refRanking(l);
        return(l);
      }

      /***
       * Transforms the RIM tree into it's canonical form (no negative thetas).
       ***/
      void transformToCanonical() {
        // Transform the tree structure.
        transformToCanonicalHelper(root);
        // Ensure that the discrepancy matrix is changed accordingly.
        RIM::List<int> rankingList = RIM::List<int>();
        root->refRanking(&rankingList);
        int* rankingArray = refRankingListToArray(&rankingList);
        numLeaves = rankingList.length();
        free(discMat);
        discMat = discrepancyMatix(rankingArray, 1, numLeaves);
        free(rankingArray);
      }

      /***
       * A helper function for logProbability. Recurses down the tree and
       * updates the log probabability of the model with respect to the input
       * average discrepancy matrix as it goes.
       *
       * @param curNode the current node it's computing the log probability for.
       * @param logProb a pointer to a double which contains the currently
       *        computed log probability. This value is updated on each
       *        recursion.
       * @param aveDiscMatrix the averaged discrepancy matrix of the data.
       * @return a RIM list containing the part of the reference ranking
       *         corresponding to the leaves of curNode.
       ***/
      RIM::List<int>* logProbHelper(RIM::RIMNode* curNode, double* logProb, double* aveDiscMatrix) {
        if(curNode->left == NULL && curNode->right == NULL) {
          // If curNode is a leaf then there is nothing to do to logProb as the
          // contribution is 0. Just insert the curNode's rank into a list and
          // return it.
          RIM::List<int>* l = new RIM::List<int>();
          l->appendValue(curNode->rank);
          return(l);
        }
        // Find the contributions of the left and right subtrees of curNode to
        // the log probability, also find the portions of the reference ranking
        // at the left and right subtrees of the curNode.
        RIM::List<int>* leftRanking = logProbHelper(curNode->left, logProb, aveDiscMatrix);
        RIM::List<int>* rightRanking = logProbHelper(curNode->right, logProb, aveDiscMatrix);

        double aveDisc = 0;
        int L = leftRanking->length();
        int R = rightRanking->length();
        leftRanking->restart();
        rightRanking->restart();
        int ind1, ind2;
        // Loop through the right and left parts of the reference ranking
        // and compute the sufficient statistic used to compute the log
        // probability at this node. This is discussed in detail in Meek and
        // Meila (2014).
        for(int i=0; i < L; i++) {
          ind1 = leftRanking->currentValue();
          leftRanking->next();

          rightRanking->restart();
          for(int j=0; j < R; j++) {
            ind2 = rightRanking->currentValue();

            rightRanking->next();
            if(ind1 < ind2) {
              aveDisc += aveDiscMatrix[ind2*this->numLeaves + ind1];
            } else {
              aveDisc += 1 - aveDiscMatrix[ind1*this->numLeaves + ind2];
            }
          }
        }
        // Here there is a negative since the score is for the negated log
        // likelihood and we wish to compute the log likelihood.
        (*logProb) -= score(L, R, aveDisc, curNode->theta);
        // Append rightRanking to leftRanking and return leftRanking.
        RIM::List<int> l = RIM::List<int>();
        leftRanking->joinWithPartition(rightRanking, &l);
        delete rightRanking;
        return(leftRanking);
      }

       /***
       * Computes the log probability of the model given the averaged
       * discrepancy matrix of the data.
       *
       * @param aveDiscMatrix the averaged discrepancy matrix of the data.
       * @return the log probability of the model.
       ***/
      double logProbability(double* aveDiscMatrix) {
        double logProb = 0;
        RIM::List<int>* l = logProbHelper(root, &logProb, aveDiscMatrix);
        delete l;
        return(logProb);
      }

      /***
       * Finds the ML theta value for the given parameters via gradient descent.
       * See Meek and Meila (2014) for details.
       *
       * @param L the number of left rankings.
       * @param R the number of right rankings.
       * @param aveDisc the sufficient statistic for the node.
       * @param theta the initializing theta value.
       * @param maxIter the maximum number of iterations to run the gradient
       *        descent.
       * @param tol the tolerance after which to stop the graident descent.
       *        should usually be around 10^(-5).
       * @return the ML theta value
       ***/
      static double mlTheta(int L, int R, double aveDisc, double theta, int maxIter, double tol) {
        int i = 0;
        int k = 0;
        double lastScore, grad;
        double curScore = score(L, R, aveDisc, theta);
        // A standard gradient descent algorithm with an Armijo line search.
        for(i = 0; i < maxIter; i++) {
          lastScore = curScore;
          grad = gradScore(L, R, aveDisc, theta);
          curScore = score(L, R, aveDisc, theta - grad);
          k = 0;
          while(curScore > lastScore && ++k < 100) {
            grad = grad/2;
            curScore = score(L, R, aveDisc, theta - grad);
          }
          theta = theta - grad;
          if(1 - curScore/lastScore < tol) {
            return(theta);
          }
        }
        //printf("WARNING: Maximum Iterations reached!\n");
        return(theta);
      }

      /***
       * Returns the gradient of the score of a node with the given parameters.
       *
       * @param L the number of left rankings.
       * @param R the number of right rankings.
       * @param aveDisc the sufficient statistic associated with the node.
       * @param theta the theta value at which to take the gradient.
       * @return the gradient of the score function at theta.
       ***/
      static double gradScore(int L, int R, double aveDisc, double theta) {
        double gScore = aveDisc;
        int n = std::max(L, R);
        if(theta != 0) {
          for(int i=n+1; i <= L+R; i++) {
            gScore += i*exp(-theta*i)/(1-exp(-theta*i)) - (i-n)*exp(-theta*(i-n))/(1-exp(-theta*(i-n)));
          }
          return(gScore);
        } else {
          return(gScore - R*L/2.0);
        }
      }

      /***
       * Returns the score of a node with the given parameters.
       *
       * @param L the number of left rankings.
       * @param R the number of right rankings.
       * @param aveDisc the sufficient statistic associated with the node.
       * @param theta the theta value at which to take the gradient.
       * @return the score function evaluated at theta.
       ***/
      static double score(int L, int R, double aveDisc, double theta) {
        return(theta*aveDisc + log(gaussianPoly(L, R, theta)));
      }

      /***
       * A matrix of normalization constants for values of L and R ranging from
       * 0,...,L and 0,...,R (sloppy notation since L,R are used in two
       * different contexts) and with parameter theta.
       *
       * @param L the number of left rankings.
       * @param R the number of right rankings.
       * @param theta the theta value.
       * @return a matrix of normalization constants.
       ***/
      static double* getZMatrix(int L, int R, double theta) {
        double* zMatrix = (double*) malloc(sizeof(double)*(L+1)*(R+1));
        double q = exp(-theta);
        // Compute the Z values for all different pairs i=1,...,L and j=1,...,R,
        // and save these into a matrix.
        if(theta != 0) {
          for(int i=0; i <= L; i++) {
            for(int j=0; j <= R; j++) {
              if(i == 0 && j == 0) {
                zMatrix[i*(R+1) + j] = 1;
              } else if(j == 0) {
                zMatrix[i*(R+1) + j] = zMatrix[(i-1)*(R+1) + j]*(1-pow(q, i+j))/(1-pow(q, i));
              } else {
                zMatrix[i*(R+1) + j] = zMatrix[i*(R+1) + (j-1)]*(1-pow(q, i+j))/(1-pow(q, j));
              }
            }
          }
        } else {
          // The normalization constant has a distinctly different form when
          // theta=0 as it is not techinically defined at theta=0 so we instead
          // use it's limit.
          for(int i=0; i <= L; i++) {
            for(int j=0; j <= R; j++) {
              if(i == 0 && j == 0) {
                zMatrix[i*(R+1) + j] = 1;
              } else if(j == 0) {
                zMatrix[i*(R+1) + j] = zMatrix[(i-1)*(R+1) + j]*(i+j)/(1.0*i);
              } else {
                zMatrix[i*(R+1) + j] = zMatrix[i*(R+1) + (j-1)]*(i+j)/(1.0*j);
              }
            }
          }
        }
        return(zMatrix);
      }

      /***
       * Returns the gaussian polynomial with parameters n=L,m=R,q=exp(-theta).
       * See Meek and Meila (2014) for more details.
       *
       * @param L the number of left rankings.
       * @param R the number of right rankings.
       * @param theta the theta value.
       ***/
      static double gaussianPoly(int L, int R, double theta) {
        if(L < 0 || R < 0) {
          printf("ERROR: One of L,R<0 in gaussianPoly.\n");
          std::exit(1);
        } else if(L == 0 || R == 0) {
          return(1.0);
        }
        double q, logG;
        int n = std::max(L,R);
        logG = 0;

        if(theta == 0) {
          for(int i=n+1; i<=R+L; i++) {
            logG += log(i);
            logG -= log(i-n);
          }
        } else {
          q = exp(-theta);
          for(int i=n+1; i<=R+L; i++) {
            logG += log(fabs(1.0-pow(q, i)));
            logG -= log(fabs(1-pow(q, i-n)));
          }
        }
        return(exp(logG));
      }

      /***
       * Computes the (nCol)x(nCol) discrepancy matrix corresponding to the
       * (nRow)x(nCol) matrix of rankings (one ranking per row) with respect
       * to the identity ranking. Rankings must permutations of the numbers
       * 0,1,....,nCol-1.
       *
       * @param rankings (nRow)x(nCol) matrix of rankings.
       * @param nRow the number of rows of rankings.
       * @param nCol the number of columns of rankings.
       * @return an (nCol)x(nCol) discrepancy matrix.
       ***/
      static int* discrepancyMatix(int* rankings, int nRow, int nCol) {
        int* discMat = (int*) malloc(sizeof(int)*nCol*nCol);
        int i,j,k;
        int r1, r2;

        // Fill discMat with 0s
        for(j=0; j<nCol; j++) {
          for(k=0; k<nCol; k++) {
            discMat[j*nCol + k] = 0;
          }
        }

        // Compute the discMat sums
        for(i=0; i<nRow; i++) {
          for(j=0; j<nCol; j++) {
            r1 = rankings[i*nCol + j];
            for(k=j+1; k<nCol; k++) {
              r2 = rankings[i*nCol + k];
              if(r1 > r2) {
                discMat[r1*nCol + r2] += 1;
              }
            }
          }
        }
        return(discMat);
      }

      /***
       * Exactly the same as discrepancyMatix except each entry of the matrix
       * is divided by nRow so that it an average of discrepancy matrices.
       *
       * @param rankings (nRow)x(nCol) matrix of rankings.
       * @param nRow the number of rows of rankings.
       * @param nCol the number of columns of rankings.
       * @return an (nCol)x(nCol) averaged discrepancy matrix.
       ***/
      static double* averagedDiscrepancyMatrix(int* rankings, int nRow, int nCol) {
        int* discMatrix = discrepancyMatix(rankings, nRow, nCol);
        double* aveDataDiscMat = (double*) malloc(sizeof(double)*nCol*nCol);
        for(int i=0; i<nCol; i++) {
          for(int j=0; j<nCol; j++) {
            aveDataDiscMat[i*nCol + j] = discMatrix[i*nCol + j] / (1.0*nRow);
          }
        }
        free(discMatrix);
        return(aveDataDiscMat);
      }
  };
}

#endif
