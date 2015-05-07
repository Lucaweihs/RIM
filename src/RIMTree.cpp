#ifndef RIM_RIMTree
#define RIM_RIMTree

//#include<Rcpp.h>
#include<stdio.h>
#include<stdlib.h>
#include<cmath>
#include"RIMNode.cpp"
#include"PartitionCache.cpp"

namespace RIM {
  class RIMTree {
    private:
      static const int MAX_ITER_DEFAULT = 100;
      static const double TOL_DEFAULT = .000001;
      RIMNode* root;
      int* discMat;
      PartitionCache* pc;
      int numLeaves;
      
      static int getMaxNLimitParts(RIMNode* node, int* n, int* limit, int* parts) {
        if(node->left == NULL && node->right==NULL) {
          *n = std::max(*n, 1);
          *limit = std::max(*limit, 1);
          *parts = std::max(*parts, 1);
          return(1);
        }
        
        int leftSize = getMaxNLimitParts(node->left, n, limit, parts);
        int rightSize = getMaxNLimitParts(node->right, n, limit, parts);
        
        *n = std::max(*n, leftSize*rightSize);
        *limit = std::max(*limit, leftSize);
        *parts = std::max(*parts, rightSize);
        return(rightSize+leftSize+1);
      }
      
      RIM::List<int>* randomRankingHelper(RIM::RIMNode* currentNode, RIM::List<double>* rands) {
        //printf("Confused2\n");
        //if(currentNode == NULL) {
        //  printf("What!!!!?\n");
        //}
        //printf("What?\n");
        if(currentNode->left == NULL && currentNode->right == NULL) {
          //printf("Confused1\n");
          RIM::List<int>* l = new List<int>();
          l->appendValue(currentNode->rank);
          return(l);
        }
        //printf("Confused0\n");
        //if((currentNode->left == NULL) || (currentNode->right == NULL)) {
        //  printf("FUUUUUUCK\n");
        //}
        //printf("rawr\n");
        RIM::List<int>* leftRanking = randomRankingHelper(currentNode->left, rands);
        RIM::List<int>* rightRanking = randomRankingHelper(currentNode->right, rands);
        //printf("rawr2\n");
        int L = leftRanking->length();
        int R = rightRanking->length();
        
        // New way of doing it
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
        
        /*
        // Old way of doing it
        int maxInversions = leftRanking->length()*rightRanking->length();
        double Z = gaussianPoly(leftRanking->length(), rightRanking->length(), currentNode->theta);
        int numInversions = -1;
        double rand = rands->currentValue();
        rands->next();
        
        while(rand > 0) {
          numInversions += 1;
          rand -= exp(-1.0*(currentNode->theta)*numInversions)*pc->countPartitions(numInversions, leftRanking->length(), rightRanking->length())/Z;
        }

        rand = rands->currentValue();
        rands->next();
        RIM::List<int>* partition = pc->randomPartition(rand, numInversions, leftRanking->length(), rightRanking->length());
        */
                
        leftRanking->joinWithPartition(rightRanking, partition);
        delete rightRanking;
        delete partition;
        return(leftRanking);
      }
      
      void transformToCanonicalHelper(RIMNode* node) {
        if(node->left == NULL && node->right == NULL) {
          return;
        }
        if(node->theta < 0) {
          node->flipSubTreesAndNegate();
        }
        transformToCanonicalHelper(node->left);
        transformToCanonicalHelper(node->right);
      }
      
      List<int>* mlThetaTreeHelper(RIMNode* curNode, double* aveDataDiscMat) {
        if(curNode->left == NULL && curNode->right == NULL) {
          RIM::List<int>* l = new List<int>();
          l->appendValue(curNode->rank);
          return(l);
        }
        
        RIM::List<int>* leftRanking = mlThetaTreeHelper(curNode->left, aveDataDiscMat);
        RIM::List<int>* rightRanking = mlThetaTreeHelper(curNode->right, aveDataDiscMat);
        
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
        
        curNode->theta = mlTheta(leftRanking->length(), rightRanking->length(), aveDisc, curNode->theta, MAX_ITER_DEFAULT, TOL_DEFAULT);
        
        RIM::List<int> tmp = RIM::List<int>();
        leftRanking->joinWithPartition(rightRanking, &tmp);
        delete rightRanking;
        return(leftRanking);
      }
    
    public:
      RIMTree(RIMNode* newRoot) {
        root = newRoot;
        RIM::List<int> rankingList = RIM::List<int>();
        root->refRanking(&rankingList);
        int* rankingArray = refRankingListToArray(&rankingList);
        
        numLeaves = rankingList.length();
        discMat = discrepancyMatix(rankingArray, 1, numLeaves);
        free(rankingArray);
        
        int n, limit, parts;
        n = 0; limit = 0; parts = 0;
        getMaxNLimitParts(root, &n, &limit, &parts);
        pc = new PartitionCache(n, limit, parts);
      }
      
      ~RIMTree() {
        root->deleteAncestors();
        delete root;
        free(discMat);
        delete pc;
      }
      
      RIMNode* getRoot() {
        return(root);
      }
      
      int getNumLeaves() {
        return(numLeaves);
      }
      
      RIM::List<double>* preOrderThetasList() {
        RIM::List<double>* l = new RIM::List<double>();
        preOrderThetasListHelper(l, root);
        return(l);
      }
      
      void preOrderThetasListHelper(RIM::List<double>* l, RIM::RIMNode* curNode) {
        if(curNode->left == NULL && curNode->right == NULL) {
          return;
        }
        preOrderThetasListHelper(l, curNode->left);
        l->appendValue(curNode->theta);
        preOrderThetasListHelper(l, curNode->right);
      }
      
      void mlThetaTree(int* rankings, int nRow) {
        int* dataDiscMat = discrepancyMatix(rankings, nRow, this->numLeaves);
        double* aveDataDiscMat = (double*) malloc(sizeof(double)*this->numLeaves*this->numLeaves);
        for(int i=0; i<this->numLeaves; i++) {
          for(int j=0; j<this->numLeaves; j++) {
            aveDataDiscMat[i*this->numLeaves + j] = dataDiscMat[i*this->numLeaves + j] / (1.0*nRow); 
          }
        }
        mlThetaTreeHelper(root, aveDataDiscMat);
        free(dataDiscMat);
        free(aveDataDiscMat);
      }
      
      RIM::List<int>* randomRanking(RIM::List<double>* rands) {
        rands->restart();
        RIM::List<int>* l = randomRankingHelper(root, rands);
        if(l->length()*(l->length()-1) != rands->length()) {
          printf("\nERROR: Not enough random numbers inputted to randomRanking.");
          std::exit(1);
        }
        // Old Way
        //if(2*(l->length()-1) != rands->length()) {
        //  std::exit(1);
        //}
        return(l);
      }
      
      int* refRankingArray() {
        RIM::List<int> l = RIM::List<int>();
        root->refRanking(&l);
        return(refRankingListToArray(&l));
      }
      
      double* treeToMatrix() {
        RIM::List<RIM::RIMNode*> nodeList = List<RIM::RIMNode*>();
        
        double* treeMat = (double*) malloc(sizeof(double)*(2*numLeaves-1)*5);
        for(int i=0; i<(2*numLeaves-1)*5; i++) {
          treeMat[i] = 0;
        }
        
        int k = 0;
        int length = 1;
        RIM::RIMNode* curNode;
        
        nodeList.appendValue(root);
        nodeList.restart();
        for(int k = 0; k < 2*numLeaves-1; k++) {
          curNode = nodeList.currentValue();
          if(curNode->left == NULL && curNode->right == NULL) {
            treeMat[k*5 + 0] = 0;
            treeMat[k*5 + 1] = 0;
            treeMat[k*5 + 2] = 1;
            treeMat[k*5 + 3] = 0;
            treeMat[k*5 + 4] = curNode->rank;
          } else {
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
      
      static int* refRankingListToArray(RIM::List<int>* l) {
        int* refArray = (int*) malloc(sizeof(int)*l->length());
        l->restart();
        for(int i=0; i<l->length(); i++) {
          refArray[i] = l->currentValue();
          l->next();
        }
        return(refArray);
      }
      
      RIM::List<int>* refRankingList() {
        RIM::List<int>* l = new RIM::List<int>();
        root->refRanking(l);
        return(l);
      }
      
      void transformToCanonical() {
        transformToCanonicalHelper(root);
      }
      
      RIM::List<int>* logProbHelper(RIM::RIMNode* curNode, double* logProb, double* aveDiscMatrix) {
        if(curNode->left == NULL && curNode->right == NULL) {
          RIM::List<int>* l = new RIM::List<int>();
          l->appendValue(curNode->rank);
          return(l);
        }
        //printf("yar\n");
        RIM::List<int>* leftRanking = logProbHelper(curNode->left, logProb, aveDiscMatrix);
        RIM::List<int>* rightRanking = logProbHelper(curNode->right, logProb, aveDiscMatrix);
        
        double aveDisc = 0;
        int L = leftRanking->length();
        int R = rightRanking->length();
        
        //printf("L=%d, R=%d, numLeaves=%d\n", L, R, this->numLeaves);
        
        leftRanking->restart();
        rightRanking->restart();
        int ind1, ind2;
        for(int i=0; i < L; i++) {
          ind1 = leftRanking->currentValue();
          leftRanking->next();
          for(int j=0; j < R; j++) {
            ind2 = rightRanking->currentValue();
            rightRanking->next();
            //printf("%d %d \n", ind1, ind2);
            if(ind1 < ind2) {
              aveDisc += aveDiscMatrix[ind2*this->numLeaves + ind1];
            } else {
              aveDisc += 1 - aveDiscMatrix[ind1*this->numLeaves + ind2];
            }
          }
        }
        //printf("bar\n");
        (*logProb) -= score(L, R, aveDisc, curNode->theta);
        //printf("tar");
        RIM::List<int> l = RIM::List<int>();
        leftRanking->joinWithPartition(rightRanking, &l);
        delete rightRanking;
        return(leftRanking);
      }
      
      double logProbability(double* aveDiscMatrix) {
        double logProb = 0;
        RIM::List<int>* l = logProbHelper(root, &logProb, aveDiscMatrix);
        delete l;
        return(logProb);
      }
      
      static double mlTheta(int L, int R, double aveDisc, double theta, int maxIter, double tol) {
        int i = 0;
        double lastScore = score(L, R, aveDisc, theta);
        double curScore, grad;
        while(true) {
          grad = gradScore(L, R, aveDisc, theta);
          curScore = score(L, R, aveDisc, theta - grad);
          while(curScore > lastScore) {
            grad = grad/2;
            curScore = score(L, R, aveDisc, theta - grad);
          }
          theta = theta - grad;
          i++;
          if(i != 1 && 1 - curScore/lastScore < tol) {
            return(theta);
          } else {
            lastScore = curScore;
          }
          if(i >= maxIter) {
            printf("WARNING: Maximum Iterations reached!\n");
            return(theta);
          }
        }
      }
      
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
      
      static double score(int L, int R, double aveDisc, double theta) {
        return(theta*aveDisc + log(gaussianPoly(L, R, theta)));
      }
      
      static double* getZMatrix(int L, int R, double theta) {
        double* zMatrix = (double*) malloc(sizeof(double)*(L+1)*(R+1));
        double q = exp(-theta);
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
      
      static double gaussianPoly(int L, int R, double theta) {
        if(L < 0 || R < 0) {
          printf("\nERROR: One of L,R<0 in gaussianPoly.");
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
  };
}

#endif