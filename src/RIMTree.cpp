#ifndef RIM_RIMTree
#define RIM_RIMTree

#include<stdio.h>
#include<stdlib.h>
#include<cmath>
#include"RIMNode.cpp"
#include"PartitionCache.cpp"

namespace RIM {
  class RIMTree {
    private:
      RIMNode* root;
      int* discMat;
      PartitionCache* pc;
      
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
      
      RIM::List<int>* randomRankingHelper(RIM::RIMNode* currentNode, RIM::List<int>* rands) {
        RIM::List<int>* l = new List<int>();
        if(currentNode->left == NULL && currentNode->right == NULL) {
          l->appendValue(currentNode->rank);
          return(l);
        }
        RIM::List<int>* leftRanking = randomRankingHelper(currentNode->left, rands);
        RIM::List<int>* rightRanking = randomRankingHelper(currentNode->right, rands);
        
        int maxInversions = leftRanking->length()*rightRanking->length();
        double Z = gaussianPoly(leftRanking->length(), rightRanking->length(), currentNode->theta);
        int numInversions = -1;
        double rand = rands->currentValue();
        rands->next();
        
        while(rand > 0) {
          numInversions += 1;
          rand -= exp(-1.0*(currentNode->theta)*numInversions)/Z;
        }
        
        rand = rands->currentValue();
        rands->next();
        RIM::List<int>* partition = pc->randomPartition(rand, numInversions, leftRanking->length(), rightRanking->length());
        leftRanking->joinWithPartition(rightRanking, partition);
      }
    
    public:
      RIMTree(RIMNode* newRoot) {
        root = newRoot;
        discMat = NULL;
        int n, limit, parts;
        getMaxNLimitParts(root, &n, &limit, &parts);
        pc = new PartitionCache(n, limit, parts);
      }
      
      RIM::List<int>* randomRanking(RIM::List<int>* rands) {
        rands->restart();
        return(randomRankingHelper(root, rands));
      }
      
      static double gaussianPoly(int L, int R, double theta) {
        if(L <= 0 || R <= 0) {
          std::exit(1);
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