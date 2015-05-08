/*****
 * The nodes used by the RIMTree class.
 *****/

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

      /***
       * A helper function for attachRankingToLeaves. Recursively attaches the
       * input rankings to the leaves of the node.
       *
       * @param ranking a RIM::List of rankings with internal current pointer
       *        set appropriately.
       * @return the size of the tree rooted at the current node.
       ***/
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
      double theta; // The theta value, only for internal nodes.
      int rank; // The rank, only for leaves.
      RIMNode* left; // Left child of node.
      RIMNode* right; // Right child of node.
      RIMNode* parent; // Parent of the node.

      /***
       * Initializes a RIMNode with given theta and rank values. One of these
       * will be ignored if the node is a leaf or internal node, the ignored
       * value may be set to anything.
       *
       * @param initTheta the initalizing theta value.
       * @param initRank the initalizing rank value.
       ***/
      RIMNode(double initTheta, int initRank) {
        left = NULL;
        right = NULL;
        parent = NULL;
        theta = initTheta;
        rank = initRank;
      }

      /***
       * Deletes all the ancestors of the node. Used when desiring to destroy
       * the entire subtree rooted at a node (except for the node itself).
       ***/
      void deleteAncestors() {
        if(left == NULL && right == NULL) {
          // If the node is a leaf there is nothing to do.
          return;
        }
        left->deleteAncestors();
        right->deleteAncestors();
        delete left;
        delete right;
        left = NULL;
        right = NULL;
      }

      /***
       * Attaches the input node as a left child of the current node.
       *
       * @param nodeForLeft the node to attach as the left child.
       ***/
      void attachLeft(RIMNode* nodeForLeft) {
        left = nodeForLeft;
        nodeForLeft->parent = this;
      }

      /***
       * Attaches the input node as a right child of the current node.
       *
       * @param nodeForRight the node to attach as the right child.
       ***/
      void attachRight(RIMNode* nodeForRight) {
        right = nodeForRight;
        nodeForRight->parent = this;
      }

      /***
       * Attaches the input ranking to the leaves of the node (i.e. replaces the
       * ranks current at the leaf nodes with those from the input ranking).
       *
       * @param ranking a RIM::List of rankings to attach.
       ***/
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

      /***
       * Flips the left/right children of the current node and negates the theta
       * value. Useful for putting a tree into canonical form.
       ***/
      void flipSubTreesAndNegate() {
        theta = -1.0*theta;
        RIMNode *tmp;
        tmp = left;
        left = right;
        right = tmp;
      }

      /***
       * Returns the number of leaves of the subtree rooted at the node.
       *
       * @return the number of leaves.
       ***/
      int numLeaves() {
        RIM::List<int> l = RIM::List<int>();
        refRanking(&l);
        return(l.length());
      }

      /***
       * Modifies the input list to contain the reference ranking of the
       * leaves of the current node.
       ***/
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

      /***
       * Prints the ranks of the leaves of the current node in preorder.
       ***/
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
