#ifndef RIM_ListNode
#define RIM_ListNode

#include<stdio.h>

namespace RIM {
  template <class T>
  class ListNode {
    public:
      ListNode<T>* prev;
      ListNode<T>* next;
      T value;

      ListNode() {
        prev = this;
        next = this;
      }

      ListNode(T val) {
        prev = this;
        next = this;
        value = val;
      }
  };
}

#endif
