#include"List.cpp"
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
  RIM::List<int> l; //= new List<int>();
  
  // Testing empty list properties
  assert(l.length() == 0);
  
  assert(l.atEnd() == true);
  assert(l.next() == false);

  // Add some things to the list and see if it still works
  l.appendValue(0);
  assert(l.length() == 1);
  assert(l.atEnd() == false);
  assert(l.next() == true);
  assert(l.firstValue() == 0);
  assert(l.lastValue() == 0);
  l.restart();
  assert(l.currentValue() == 0);
  assert(l.firstNode() != NULL);
  assert(l.lastNode() != NULL);
  assert(l.currentNode() != NULL);
  assert(l.next() == false);
  assert(l.atEnd() == true);
  
  for(int i=1; i<8; i++) {
    l.appendValue(i);
  }
  assert(l.length() == 8);
  assert(l.atEnd() == false);
  assert(l.next() == true);
  assert(l.firstValue() == 0);
  assert(l.lastValue() == 7);
  assert(l.currentValue() == 1);
  l.next();
  l.next();
  l.next();
  assert(l.firstValue() == 0);
  assert(l.lastValue() == 7);
  assert(l.currentValue() == 4);
  assert(l.firstNode()->value == 0);
  assert(l.lastNode()->value == 7);
  assert(l.currentNode()->value == 4);
  
  
  // Try adding a node
  RIM::ListNode<int>* listNode = new RIM::ListNode<int>(8);
  l.appendNode(listNode);
  
  assert(l.length() == 9);
  assert(l.atEnd() == false);
  assert(l.next() == true);
  assert(l.firstValue() == 0);
  assert(l.lastValue() == 8);
  assert(l.currentValue() == 5);
  while(l.next() == true);
  assert(l.currentValue() == 8);
  assert(l.lastNode() == listNode);
  assert(l.currentNode() == listNode);
  
  // Testing joinWithPartition
  RIM::List<int> l1;
  RIM::List<int> l2;
  RIM::List<int> partition;

  // joinWithPartition 1
  l1 = RIM::List<int>();
  l2 = RIM::List<int>();
  l1.appendValue(1); l1.appendValue(2);
  l1.appendValue(3); l1.appendValue(4);
  l2.appendValue(5); l2.appendValue(6);
  l2.appendValue(7); l2.appendValue(8);
  
  partition = RIM::List<int>();
  l1.joinWithPartition(&l2, &partition);
  l1.restart();
  for(int i=0; i<8; i++) {
    assert(l1.currentValue() == i+1);
    l1.next();
  }
  assert(l1.atEnd());
  
  // joinWithPartition 2
  int a2[8] = {5,6,7,8,1,2,3,4};

  l1 = RIM::List<int>();
  l2 = RIM::List<int>();
  l1.appendValue(1); l1.appendValue(2);
  l1.appendValue(3); l1.appendValue(4);
  l2.appendValue(5); l2.appendValue(6);
  l2.appendValue(7); l2.appendValue(8);
  
  partition = RIM::List<int>();
  partition.appendValue(4); partition.appendValue(4);
  partition.appendValue(4); partition.appendValue(4);
  l1.joinWithPartition(&l2, &partition);
  
  l1.restart();
  for(int i=0; i<8; i++) {
    assert(l1.currentValue() == a2[i]);
    l1.next();
  }
  assert(l1.atEnd());

  // joinWithPartition 3
  int a3[8] = {5,1,6,2,7,3,8,4};

  l1 = RIM::List<int>();
  l2 = RIM::List<int>();
  l1.appendValue(1); l1.appendValue(2);
  l1.appendValue(3); l1.appendValue(4);
  l2.appendValue(5); l2.appendValue(6);
  l2.appendValue(7); l2.appendValue(8);
  
  partition = RIM::List<int>();
  partition.appendValue(4); partition.appendValue(3);
  partition.appendValue(2); partition.appendValue(1);
  l1.joinWithPartition(&l2, &partition);
  l1.restart();
  for(int i=0; i<8; i++) {
    assert(l1.currentValue() == a3[i]);
    l1.next();
  }
  assert(l1.atEnd());
  
  // joinWithPartition 4
  int a4[8] = {1,5,2,6,3,7,4,8};

  l1 = RIM::List<int>();
  l2 = RIM::List<int>();
  l1.appendValue(1); l1.appendValue(2);
  l1.appendValue(3); l1.appendValue(4);
  l2.appendValue(5); l2.appendValue(6);
  l2.appendValue(7); l2.appendValue(8);
  
  partition = RIM::List<int>();
  partition.appendValue(3); partition.appendValue(2);
  partition.appendValue(1);
  l1.joinWithPartition(&l2, &partition);
  l1.restart();
  for(int i=0; i<8; i++) {
    assert(l1.currentValue() == a4[i]);
    l1.next();
  }
  assert(l1.atEnd());
  
  // joinWithPartition 5
  int a5[8] = {5,1,2,6,7,3,4,8};

  l1 = RIM::List<int>();
  l2 = RIM::List<int>();
  l1.appendValue(1); l1.appendValue(2);
  l1.appendValue(3); l1.appendValue(4);
  l2.appendValue(5); l2.appendValue(6);
  l2.appendValue(7); l2.appendValue(8);
  
  partition = RIM::List<int>();
  partition.appendValue(4); partition.appendValue(2);
  partition.appendValue(2);
  l1.joinWithPartition(&l2, &partition);
  l1.restart();
  for(int i=0; i<8; i++) {
    assert(l1.currentValue() == a5[i]);
    l1.next();
  }
  assert(l1.atEnd());
  
  // joinWithPartition 6
  int a6[8] = {1,5,2,6,7,3,8,4};

  l1 = RIM::List<int>();
  l2 = RIM::List<int>();
  l1.appendValue(1); l1.appendValue(2);
  l1.appendValue(3); l1.appendValue(4);
  l2.appendValue(5); l2.appendValue(6);
  l2.appendValue(7); l2.appendValue(8);
  
  partition = RIM::List<int>();
  partition.appendValue(3); partition.appendValue(2);
  partition.appendValue(2); partition.appendValue(1);
  l1.joinWithPartition(&l2, &partition);
  l1.restart();
  for(int i=0; i<8; i++) {
    assert(l1.currentValue() == a6[i]);
    l1.next();
  }
  assert(l1.atEnd());

  return(0);
}