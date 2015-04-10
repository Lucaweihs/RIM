#include"rim.h"
#include<stdio.h>
#include<stdlib.h>

List* createEmptyList() {
  List* l = (List*) malloc(sizeof(List));
  l->sentinal = (ListNode*) malloc(sizeof(ListNode));
  l->sentinal->next = l->sentinal;
  l->sentinal->prev = l->sentinal;
  l->sentinal->value = NULL;
  l->current = l->sentinal;
  l->length = 0;
  return(l);
}
int length(List* l) {
  return(l->length);
}
bool atEnd(List* l) {
  if(l->sentinal == l->current->next) {
    return(true);
  } else {
    return(false);
  }
}
bool next(List* l) {
  if(true == atEnd(l)) {
    return(false);
  } else {
    l->current = l->current->next;
    return(true);
  }
}
void restart(List* l) {
    l->current = l->sentinal->next;
}
ListNode* firstNode(List* l) {
  if(l->length == 0) {
    return(NULL);
  }
  return(l->sentinal->next);
}
void* firstValue(List* l) {
  if(l->length == 0) {
    return(NULL);
  }
  return(l->sentinal->next->value);
}
ListNode* lastNode(List* l) {
  if(l->length == 0) {
    return(NULL);
  }
  return(l->sentinal->prev);
}
void* lastValue(List* l) {
  if(l->length == 0) {
    return(NULL);
  }
  return(l->sentinal->prev->value);
}
ListNode* currentNode(List* l) {
  if(l->length == 0) {
    return(NULL);
  }
  return(l->current);
}
void* currentValue(List* l) {
  if(l->length == 0) {
    return(NULL);
  }
  return(l->current->value);
}
void appendNode(List* l, ListNode* listNode) {
  listNode->next = l->sentinal;
  listNode->prev = l->sentinal->prev;
  l->sentinal->prev->next = listNode;
  l->sentinal->prev = listNode;
  l->length += 1;
}
void appendValue(List* l, void* value) {
  ListNode* listNode = (ListNode*) malloc(sizeof(ListNode));
  listNode->value = value;
  appendNode(l, listNode);
}


int* discrepancyMat(int* rankings, int nRow, int nCol) {
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

int* rankingFromList(List* rList) {
  int n = length(rList);
  if(n == 0) {
    return NULL;
  }
  int* ranking = (int*) malloc(sizeof(int)*n);
  int i;
  
  restart(rList);
  for(i=0; i<n; i++) {
    ranking[i] = *((int*) currentValue(rList));
    next(rList);
  }
  return(ranking);
}

RimNode* createEmptyRimNode() {
  RimNode* rimNode = (RimNode*) malloc(sizeof(RimNode));
  rimNode->left = NULL;
  rimNode->right = NULL;
  rimNode->parent = NULL;
  rimNode->theta = 0;
  rimNode->element = NULL;
  return(rimNode);
}

void attachLeft(RimNode* node, RimNode* nodeForLeft) {
  node->left = nodeForLeft;
  nodeForLeft->parent = node;
}

void attachRight(RimNode* node, RimNode* nodeForRight) {
  node->right = nodeForRight;
  nodeForRight->parent = node;
}

void attachRankingToLeaves(RimNode* root, List* ranking) {
  if(length(ranking) == 0) {
    return;
  }
  restart(ranking);
  int treeSize = attachRankingToLeavesHelper(root, ranking);
  if(treeSize != 2*length(ranking) - 1) {
    exit(1);
  }
}

int attachRankingToLeavesHelper(RimNode* current, List* ranking) {
  int treeSize = 0;
  if(current->left == NULL && current->right == NULL) {
    current->element = currentNode(ranking);
    next(ranking);
    return(1);
  }
  printf("%f\n", current->theta);
  if(current->left != NULL) {
    treeSize += attachRankingToLeavesHelper(current->left, ranking);
  }
  
  if(current->right != NULL) {
    treeSize += attachRankingToLeavesHelper(current->right, ranking);
  }
  
  return(treeSize + 1);
}

/*
RimNode* minNode(RimNode* current) {
  if(current->left == NULL && current->right == NULL) {
    return(current);
  }
  return(minNode(current->left));
}

RimNode* maxNode(RimNode* current) {
  if(current->left == NULL && current->right == NULL) {
    return(current);
  }
  return(maxNode(current->right));
}
*/

List* randomPartition(int n, int limit, int maxParts) {
  int* countMat = (int *) malloc(sizeof(int)*n*limit*maxParts);
  int total = countPartitions(countMat, n, limit, maxParts, n, limit, maxParts);
  int which = rand() % total;
  List* partition = createEmptyList();
  int i, currentN, currentLimit, currentMaxParts, ind, count;
  currentN = n;
  currentLimit = limit;
  currentMaxParts = maxParts;
  int* tmp;
  
  while(currentN != 0) {
    for(i=1; i<min(currentLimit, currentN); i++) {
      ind = (currentN-i)*(limit*maxParts) + i*maxParts + (currentMaxParts-1);
      count = countMat[ind];
      if(count >= which) {
        break;
      }
      which -= count;
    }
    tmp = (int*) malloc(sizeof(int));
    tmp = i;
    appendValue(partition, &tmp);
    currentN -= i;
    currentLimit = i;
  }
  return(partition);
}

int countPartitions(int* countMat, int n, int limit, int maxParts, int currentN, int currentLimit, int currentMaxParts) {
  int ind = currentN*(limit*maxParts) + currentLimit*maxParts + currentMaxParts;
  if(countMat[ind] != 0) {
    return(countMat[ind]);
  }
  else if(currentN == 0) {
    countMat[ind] = 1;
    return(1);
  }
  else if(currentN > currentLimit*currentMaxParts) {
    return(0);
  } 
  else {
    for(int newLimit=1; newLimit <= min(currentLimit, currentN); newLimit++) {
      countMat[ind] += countPartitions(countMat, n, limit, maxParts, currentN-newLimit, newLimit, currentMaxParts-1);
    }
    return(countMat[ind]);
  }
}

void flipSubTrees(RimNode* root) {
  RimNode *tmp;
  tmp = root->left;
  root->left = root->right;
  root->right = tmp;
}

RimTree* createRimTree(RimNode* root, List* refList) {
  RimTree* rimTree = (RimTree*) malloc(sizeof(RimTree));
  rimTree->root = root;
  rimTree->refList = refList;
  rimTree->refArray = rankingFromList(refList);
  rimTree->discMat = discrepancyMat(rimTree->refArray, 1, length(refList));
  return(rimTree); 
}