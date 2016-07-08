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

#include "RIMNode.h"

using namespace RIM;

/***
 * A helper function for attachRankingToLeaves. Recursively attaches the
 * input rankings to the leaves of the node.
 *
 * @param ranking a RIM::List of rankings with internal current pointer
 *        set appropriately.
 * @return the size of the tree rooted at the current node.
 ***/
int RIMNode::attachRankingToLeavesHelper(List<int>* ranking) {
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

/***
 * Initializes a RIMNode with given theta and rank values. One of these
 * will be ignored if the node is a leaf or internal node, the ignored
 * value may be set to anything.
 *
 * @param initTheta the initalizing theta value.
 * @param initRank the initalizing rank value.
 ***/
RIMNode::RIMNode(double initTheta, int initRank) {
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
void RIMNode::deleteAncestors() {
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
void RIMNode::attachLeft(RIMNode* nodeForLeft) {
  left = nodeForLeft;
  nodeForLeft->parent = this;
}

/***
 * Attaches the input node as a right child of the current node.
 *
 * @param nodeForRight the node to attach as the right child.
 ***/
void RIMNode::attachRight(RIMNode* nodeForRight) {
  right = nodeForRight;
  nodeForRight->parent = this;
}

/***
 * Attaches the input ranking to the leaves of the node (i.e. replaces the
 * ranks current at the leaf nodes with those from the input ranking).
 *
 * @param ranking a RIM::List of rankings to attach.
 ***/
void RIMNode::attachRankingToLeaves(RIM::List<int>* ranking) {
  if(ranking->length() == 0) {
    return;
  }
  ranking->restart();
  int treeSize = attachRankingToLeavesHelper(ranking);
  if(treeSize != 2 * ranking->length() - 1) {
    Rcpp::stop("\nERROR: Ranking too attach was too short or too long for RIM tree.\n");
  }
}

/***
 * Flips the left/right children of the current node and negates the theta
 * value. Useful for putting a tree into canonical form.
 ***/
void RIMNode::flipSubTreesAndNegate() {
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
int RIMNode::numLeaves() {
  RIM::List<int> l;
  refRanking(&l);
  return(l.length());
}

/***
 * Modifies the input list to contain the reference ranking of the
 * leaves of the current node.
 ***/
void RIMNode::refRanking(RIM::List<int>* l) {
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
void RIMNode::preOrderPrint() {
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

/***
 * Makes a copy of the current RIM node and its entire subtree
 */
RIMNode* RIMNode::deepCopy() {
  RIMNode* newRoot = new RIMNode(theta, rank);
  if (left != NULL) {
    newRoot->attachLeft(left->deepCopy());
  }
  if (right != NULL) {
    newRoot->attachRight(right->deepCopy());
  }
  return newRoot;
}

/***
 * Returns true if node is a leaf and false otherwise.
 */
bool RIMNode::isLeaf() {
  if (left == NULL && right == NULL) {
    return true;
  } else if (left != NULL && right != NULL) {
    return false;
  } else {
    Rcpp::stop("RIMNode in invalid state.");
    return false;
  }
}
