#include"List.cpp"
#include<stdio.h>
#include<ctype.h>
#include<assert.h>


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
  
  return(0);
}