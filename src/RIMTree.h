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
#include"RIMNode.h"
#include"RcppArmadillo.h"

namespace RIM {
class RIMTree {
private:
  static const int MAX_ITER_DEFAULT = 200; // Iterations when finding MLE for theta values
  static const double TOL_DEFAULT = .00001; // Tolerance in grad descent for finding MLE
  // theta values.
  RIMNode* root; // Root of the tree
  arma::imat discMat; // Discrepancy matrix of the tree's reference permutation
  // with respect to the identity permutation.
  int numLeaves; // Number of leaves in the tree.

  RIM::List<int>* randomRankingHelper(RIM::RIMNode* currentNode) const;
  void transformToCanonicalHelper(RIMNode* node);
  List<int>* mlThetaTreeHelper(RIMNode* curNode, const arma::mat& aveDataDiscMat);

  // We disable copying by making this function private and providing
  // no implementation
  RIMTree(const RIMTree& that);

public:
  RIMTree(RIMNode* newRoot);
  ~RIMTree();

  int getNumLeaves() const;
  RIM::List<double>* preOrderThetasList() const;
  void preOrderThetasListHelper(RIM::List<double>* l, RIM::RIMNode* curNode) const;
  void mlThetaTree(const arma::mat& aveDiscMat);
  arma::ivec randomRanking() const;
  arma::ivec refRankingArray() const;
  arma::mat treeToMatrix() const;
  RIM::List<int>* refRankingList() const;
  void transformToCanonical();
  RIM::List<int>* logProbHelper(RIM::RIMNode* curNode, double& logProb,
                                const arma::mat& aveDiscMatrix) const;
  double logProbability(const arma::mat& aveDiscMatrix) const;

  static arma::ivec refRankingListToArray(RIM::List<int>* l);
  static double mlThetaOld(int L, int R, double aveDisc, double theta, int maxIter, double tol);
  static double mlTheta(int L, int R, double aveDisc, double theta, int maxIter, double tol);
  static double gradScore(int L, int R, double aveDisc, double theta);
  static double gradScoreExpNegTheta(int L, int R, double aveDisc, double expNegTheta);
  static double score(int L, int R, double aveDisc, double theta);
  static arma::mat getZMatrix(int L, int R, double theta);
  static double logGaussianPoly(int L, int R, double theta);
  static double gaussianPoly(int L, int R, double theta);
  static arma::imat discrepancyMatix(const arma::imat& rankings);
  static arma::mat averagedDiscrepancyMatrix(const arma::imat& rankings);
};
}

#endif
