#ifndef RIM_List
#define RIM_List

#include"ListNode.cpp"
#include<cstdlib>

namespace RIM {
  template <class T>
  class List {
    private:
      ListNode<T>* sentinal;
      ListNode<T>* current;
      int listLength;
    
    public:
      List() {
        sentinal = new ListNode<T>();
        current = sentinal;
        listLength = 0;
      }
      
      int length() {
        return(listLength);
      }
      
      bool atEnd() {
        return(current->next == sentinal);
      }
      
      bool next() {
        if(atEnd()) {
          return(false);
        } else{
          current = current->next;
          return(true);
        }
      }
      
      void restart() {
        current = sentinal->next;
      }
      
      ListNode<T>* firstNode() {
        if(length() == 0) {
          return(NULL);
        }
        return(sentinal->next);
      }
      
      T firstValue() {
        if(length() == 0) {
          std::exit(1);
        }
        return(sentinal->next->value);
      }
      
      ListNode<T>* lastNode() {
        if(length() == 0) {
          return(NULL);
        }
        return(sentinal->prev);
      }
      
      T lastValue() {
        if(length() == 0) {
          std::exit(1);
        }
        return(sentinal->prev->value);
      }
      
      ListNode<T>* currentNode() {
        if(length() == 0) {
          std::exit(1);
        }
        return(current);
      }
      
      T currentValue() {
        if(length() == 0) {
          std::exit(1);
        }
        return(current->value);
      }
      
      void appendNode(ListNode<T>* listNode) {
        listNode->next = sentinal;
        listNode->prev = sentinal->prev;
        sentinal->prev->next = listNode;
        sentinal->prev = listNode;
        listLength += 1;
      }
      
      void appendValue(T value) {
        ListNode<T>* listNode = new ListNode<T>(value);
        appendNode(listNode);
      }
      
      void insertNodeAfterCurrent(ListNode<T>* listNode) {
        listNode->next = current->next;
        listNode->prev = current;
        current->next->prev = listNode;
        current->next = listNode;
        listLength += 1;
      }
      
      void insertNodeBeforeCurrent(ListNode<T>* listNode) {
        listNode->next = current;
        listNode->prev = current->prev;
        current->prev->next = listNode;
        current->prev = listNode;
        listLength += 1;
      }
      
      void joinWithPartition(List<T>* listToMerge, List<int>* partition) {
        if(partition->length() > listToMerge->length()) {
          std::exit(1);
        }
        partition->restart();
        this->restart();
        listToMerge->restart();
        
        int initialLength = this->length();
        int lastPartitionNum = initialLength;
        int curPartitionNum = initialLength;
        //int lastPartitionNum = partition->currentValue();
        //int curPartitionNum = lastPartitionNum;
        /*partition->next();
        for(int i=0; i<initialLength-lastPartitionNum-1; i++) {
          this->next();
        }*/
        
        ListNode<int>* curNodeToMerge = listToMerge->sentinal->next;
        ListNode<int>* nextNodeToMerge;
        int k = 1;
        
        while(true) {
          nextNodeToMerge = curNodeToMerge->next;
          
          if(k > partition->length()) {
            curPartitionNum = -1;
            while(!this->atEnd()) {
              this->next();
            }
          } else {
            curPartitionNum = partition->currentValue();
            partition->next();
            if(lastPartitionNum == initialLength) {
              for(int i=0; i<lastPartitionNum - curPartitionNum - 1; i++) {
                this->next();
              }
            } else {
              for(int i=0; i<lastPartitionNum - curPartitionNum; i++) {
                this->next();
              }
            }
            lastPartitionNum = curPartitionNum;
            k++;
          }
          
          if(curPartitionNum == initialLength) {
            this->insertNodeBeforeCurrent(curNodeToMerge);
          } else {
            this->insertNodeAfterCurrent(curNodeToMerge);
            this->next();
          }
          
          if(nextNodeToMerge == listToMerge->sentinal) {
            break;
          }
          curNodeToMerge = nextNodeToMerge;
        }
        /*
        if(curPartitionNum == initialLength) {
          this->insertNodeBeforeCurrent(curNodeToMerge);
        } else {
          this->insertNodeAfterCurrent(curNodeToMerge);
          this->next();
        }
        curNodeToMerge = nextNodeToMerge;
        nextNodeToMerge = curNodeToMerge->next;
        
        while(curNodeToMerge != listToMerge->sentinal) {
          curPartitionNum = partition->currentValue();
          partition->next();
          for(int i=0; i<lastPartitionNum - curPartitionNum; i++) {
            this->next();
          }
          
          if(curPartitionNum == initialLength) {
            this->insertNodeBeforeCurrent(curNodeToMerge);
          } else {
            this->insertNodeAfterCurrent(curNodeToMerge);
            this->next();
          }
          curNodeToMerge = nextNodeToMerge;
          nextNodeToMerge = curNodeToMerge->next;
          
          lastPartitionNum = curPartitionNum;
        }
        */
        listToMerge->sentinal->next = listToMerge->sentinal;
        listToMerge->sentinal->prev = listToMerge->sentinal;
        listToMerge->listLength = 0;
      }
    
  };
}

#endif