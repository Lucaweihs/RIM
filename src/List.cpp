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
  };
}