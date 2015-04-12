#include"PartitionCache.cpp"

#include<stdio.h>
#include<ctype.h>
#include<assert.h>

void printList(RIM::List<int>* l) {
  if(l->length() == 0) {
    printf("\n");
    return;
  }
  l->restart();
  printf("%d,", l->currentValue());
  while(!l->atEnd()) {
    l->next();
    printf("%d,", l->currentValue());
  }
  printf("\n");
}

int main() {
  // Testing Parition Cache
  RIM::PartitionCache pc = RIM::PartitionCache(9, 3, 3);
  
  assert(pc.atIndex(1,1,1) == -1);
  
  assert(pc.countPartitions(0,0,0) == 1);
  assert(pc.countPartitions(1,0,0) == 0);
  assert(pc.countPartitions(0,1,0) == 1);
  assert(pc.countPartitions(0,0,1) == 1);
  assert(pc.countPartitions(1,1,0) == 0);
  assert(pc.countPartitions(1,0,1) == 0);
  assert(pc.countPartitions(0,1,1) == 1);
  
  assert(pc.countPartitions(4,3,3) == 3);
  assert(pc.countPartitions(4,2,2) == 1);
  assert(pc.countPartitions(4,3,2) == 2);
  
  assert(pc.countPartitions(3,3,3) == 3);
  assert(pc.countPartitions(4,3,3) == 3);
  assert(pc.countPartitions(5,3,3) == 3);
  assert(pc.countPartitions(6,3,3) == 3);
  assert(pc.countPartitions(7,3,3) == 2);
  assert(pc.countPartitions(8,3,3) == 1);
  assert(pc.countPartitions(9,3,3) == 1);
  
  //printList(pc.randomPartition(.999, 5, 3, 3));
  pc = RIM::PartitionCache(8, 6, 6);
  assert(pc.countPartitions(8,6,6) == 18);
  for(int i=1; i<=18; i++) {
    printList(pc.randomPartition(i/18.0, 8, 6, 6));
  }
  
  return(0);
}