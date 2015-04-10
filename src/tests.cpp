#include"rim.h"
#include<stdio.h>
#include<ctype.h>
#include<assert.h>
  
void printMat(int* mat, int nRow, int nCol) {
  for(int i=0; i<nRow; i++) {
      for(int j=0; j<nCol; j++) {
          printf("%d ", mat[i*nCol + j]);
      }
    printf("\n");
  }
}

void preOrderPrintHelper(RimNode* current) {
  if(current->left == NULL && current->right == NULL) {
    printf("%d ", *((int*) current->element->value));
  }
  if(current->left != NULL) {
    preOrderPrintHelper(current->left);
  }
  if(current->right != NULL) {
    preOrderPrintHelper(current->right);
  }
}
void preOrderPrint(RimNode* current) {
  preOrderPrintHelper(current);
  printf("\n");
}

int main() {
  List* l = createEmptyList();

  // Testing empty list properties
  assert(length(l) == 0);
  assert(atEnd(l) == true);
  assert(next(l) == false);
  assert(firstValue(l) == NULL);
  assert(lastValue(l) == NULL);
  assert(currentValue(l) == NULL);
  assert(firstNode(l) == NULL);
  assert(lastNode(l) == NULL);
  assert(currentNode(l) == NULL);

  // Add some things to the list and see if it still works
  int a[8] = {0,1,2,3,4,5,6,7};
  appendValue(l, &(a[0]));
  
  assert(length(l) == 1);
  assert(atEnd(l) == false);
  assert(next(l) == true);
  assert(firstValue(l) == &(a[0]));
  assert(lastValue(l) == &(a[0]));
  restart(l);
  assert(currentValue(l) == &(a[0]));
  assert(firstNode(l) != NULL);
  assert(lastNode(l) != NULL);
  assert(currentNode(l) != NULL);
  assert(next(l) == false);
  assert(atEnd(l) == true);
  
  for(int i=1; i<8; i++) {
    appendValue(l, &(a[i]));
  }
  assert(length(l) == 8);
  assert(atEnd(l) == false);
  assert(next(l) == true);
  assert(firstValue(l) == &(a[0]));
  assert(lastValue(l) == &(a[7]));
  assert(currentValue(l) == &(a[1]));
  next(l);
  next(l);
  next(l);
  assert(firstValue(l) == &(a[0]));
  assert(lastValue(l) == &(a[7]));
  assert(currentValue(l) == &(a[4]));
  assert(firstNode(l)->value == &(a[0]));
  assert(lastNode(l)->value == &(a[7]));
  assert(currentNode(l)->value == &(a[4]));
  
  // Try adding a node
  int* c = (int*) malloc(sizeof(int));
  *c = 8;
  ListNode* listNode = (ListNode*) malloc(sizeof(ListNode));
  listNode->value = c;
  appendNode(l, listNode);
  
  assert(length(l) == 9);
  assert(atEnd(l) == false);
  assert(next(l) == true);
  assert(firstValue(l) == &(a[0]));
  assert(lastValue(l) == c);
  assert(currentValue(l) == &(a[5]));
  while(next(l) == true);
  assert(currentValue(l) == c);
  assert(lastNode(l) == listNode);
  assert(currentNode(l) == listNode);
  
  // Testing discrepency matrices
  int ranking[16] = {3,2,1,0, 0,1,2,3, 0,1,3,2, 1,0,3,2};
  int whatMatShouldBe[16] = {0,0,0,0, 2,0,0,0, 1,1,0,0, 1,1,3,0};
  int* mat = discrepancyMat(ranking, 4, 4);
  for(int i=0; i<4; i++) {
    for(int j=0; j<4; j++) {
      assert(mat[i*4+j] == whatMatShouldBe[i*4+j]);
    }
  }
  
  // Testing rankingFromList
  int* ranking1 = rankingFromList(l);
  for(int i=0; i<9; i++) {
    assert(i == ranking1[i]);
  }
  
  // Testing Log-Likelihood computaton
  // todo: destroy old list
  int ranking2[4] = {1,0,3,2};
  l = createEmptyList();
  for(int i=0; i<4; i++) {
    appendValue(l, &(ranking2[i]));
  }
  RimNode* root = createEmptyRimNode();
  root->theta = -.1;
  attachLeft(root, createEmptyRimNode());
  attachLeft(root->left, createEmptyRimNode());
  attachRight(root->left, createEmptyRimNode());
  attachRight(root, createEmptyRimNode());
  attachLeft(root->right, createEmptyRimNode());
  attachRight(root->right, createEmptyRimNode());
  
  root->left->theta = .8;
  root->right->theta = 1.6;
  attachRankingToLeaves(root, l);
  preOrderPrint(root);
  
  return(0);
}

