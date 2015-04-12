#include"RIMTree.cpp"

#include<stdio.h>
#include<ctype.h>
#include<assert.h>

void printList(RIM::List<int>* l) {
    l->restart();
    printf("%d,", l->currentValue());
    while(!l->atEnd()) {
      l->next();
      printf("%d,", l->currentValue());
    }
    printf("\n");
}

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
  RIM::RIMNode* root = new RIM::RIMNode(-.1, -1);
  root->attachLeft(new RIM::RIMNode(.8, -1));
  root->left->attachLeft(new RIM::RIMNode(0, 1));
  root->left->attachRight(new RIM::RIMNode(0, 2));
  root->attachRight(new RIM::RIMNode(1.6, -1));
  root->right->attachLeft(new RIM::RIMNode(0, 3));
  root->right->attachRight(new RIM::RIMNode(0, 4));
  
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
  
  // Testing random sample
  RIM::RIMTree rim = RIM::RIMTree(root);
  RIM::List<double>* rands = new RIM::List<double>();
  rands->appendValue(0.7); rands->appendValue(0.8940123);
  rands->appendValue(0.1154309); rands->appendValue(0.6486586);
  rands->appendValue(0.9428670); rands->appendValue(0.2261249);
  
  printList(rim.randomRanking(rands));
  
  printf("Complete\n");
  return(0);
}