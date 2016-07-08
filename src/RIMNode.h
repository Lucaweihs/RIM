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

/*****
 * The nodes used by the RIMTree class.
 *****/

#ifndef RIM_RIMNode
#define RIM_RIMNode

#include<stdio.h>
#include<stdlib.h>
#include<iostream>
#include<string>
#include"List.h"

namespace RIM {
class RIMNode {
private:
  int attachRankingToLeavesHelper(List<int>* ranking);

  // We disable copying by making this function private and providing
  // no implementation
  RIMNode(const RIMNode& that);

public:
  double theta; // The theta value, only for internal nodes.
  int rank; // The rank, only for leaves.
  RIMNode* left; // Left child of node.
  RIMNode* right; // Right child of node.
  RIMNode* parent; // Parent of the node.

  RIMNode(double initTheta, int initRank);

  bool isLeaf();
  RIMNode* deepCopy();
  void deleteAncestors();
  void attachLeft(RIMNode* nodeForLeft);
  void attachRight(RIMNode* nodeForRight);
  void attachRankingToLeaves(RIM::List<int>* ranking);
  void flipSubTreesAndNegate();
  int numLeaves();
  void refRanking(RIM::List<int>* l);
  void preOrderPrint();
};
}

#endif
