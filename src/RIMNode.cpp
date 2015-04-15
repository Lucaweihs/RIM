#ifndef RIM_RIMNode
#define RIM_RIMNode

#include<stdio.h>
#include<stdlib.h>
#include<iostream>
#include<string>
#include"List.cpp"

namespace RIM {
  class RIMNode {
    private:
      int attachRankingToLeavesHelper(List<int>* ranking) {
        int treeSize = 1;
        if(left == NULL && right == NULL) {
          rank = ranking->currentValue();
          ranking->next();
          return(treeSize);
        }
        if(left != NULL) {
          treeSize += left->attachRankingToLeavesHelper(ranking);
        }
        if(right != NULL) {
          treeSize += right->attachRankingToLeavesHelper(ranking);
        }
        return(treeSize);
      }
    
    public:
      double theta;
      int rank;
      RIMNode* left;
      RIMNode* right;
      RIMNode* parent;
      
      RIMNode(double initTheta, int initRank) {
        left = NULL;
        right = NULL;
        parent = NULL;
        theta = initTheta;
        rank = initRank;
      }
      
      void attachLeft(RIMNode* nodeForLeft) {
        left = nodeForLeft;
        nodeForLeft->parent = this;
      }
      
      void attachRight(RIMNode* nodeForRight) {
        right = nodeForRight;
        nodeForRight->parent = this;
      }
      
      void attachRankingToLeaves(RIM::List<int>* ranking) {
        if(ranking->length() == 0) {
          return;
        }
        ranking->restart();
        int treeSize = attachRankingToLeavesHelper(ranking);
        if(treeSize != 2*ranking->length() - 1) {
          printf("\nERROR: Ranking too attach was too short or too long for RIM tree.\n");
          std::exit(1);
        }
      }
      
      void flipSubTreesAndNegate() {
        theta = -1.0*theta;
        RIMNode *tmp;
        tmp = left;
        left = right;
        right = tmp;
      }
      
      int numLeaves() {
        RIM::List<int>* l = new RIM::List<int>();
        refRanking(l);
        return(l->length());
      }
      
      void refRanking(RIM::List<int>* l) {
        if(left == NULL && right == NULL) {
          l->appendValue(rank);
        }
        if(left != NULL) {
          left->refRanking(l);
        }
        if(right != NULL) {
          right->refRanking(l);
        }
      }
      
      void preOrderPrint() {
        if(left == NULL && right == NULL) {
          printf("%d", rank);
        }
        if(left != NULL) {
          left->preOrderPrint();
        }
        if(right != NULL) {
          right->preOrderPrint();
        }
      }
  };
}

#endif