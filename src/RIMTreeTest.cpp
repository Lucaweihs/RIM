#include"RIMTree.cpp"

#include<stdio.h>
#include<ctype.h>
#include<assert.h>

int main() {
  
  // Testing discrepency matrices
  int ranking[16] = {3,2,1,0, 0,1,2,3, 0,1,3,2, 1,0,3,2};
  int whatMatShouldBe[16] = {0,0,0,0, 2,0,0,0, 1,1,0,0, 1,1,3,0};
  int* mat = RIM::RIMTree::discrepancyMatix(ranking, 4, 4);
  for(int i=0; i<4; i++) {
    for(int j=0; j<4; j++) {
      //printf("%d ", (mat[i*4+j]));
      assert(mat[i*4+j] == whatMatShouldBe[i*4+j]);
    }
    //printf("\n");
  }

  // Testing Tree Structure
  RIM::RIMNode* root = new RIM::RIMNode(.12372, 1237);
  root->attachLeft(new RIM::RIMNode(12312.1, 233));
  root->left->attachLeft(new RIM::RIMNode(.5984, 1));
  root->left->attachRight(new RIM::RIMNode(.8283, 2));
  root->attachRight(new RIM::RIMNode(-.1232, 5));
  root->right->attachLeft(new RIM::RIMNode(.488, 3));
  root->right->attachRight(new RIM::RIMNode(.82819, 4));
  
  RIM::List<int> l = RIM::List<int>();
  root->refRanking(&l);
  int i = 1;
  l.restart();
  assert(l.currentValue() == i);
  while(!l.atEnd()) {
    i++;
    l.next();
    assert(l.currentValue() == i);
  }
  
  // Testing Gaussian Poly Computation
  assert(RIM::RIMTree::gaussianPoly(1, 1, 0) == 2.0);
  assert(fabs(RIM::RIMTree::gaussianPoly(5, 3, 0) - 56.0) < .00001);
  assert(fabs(RIM::RIMTree::gaussianPoly(5, 8, -.1) - 11976.49) < .1);
    
  return(0);
}