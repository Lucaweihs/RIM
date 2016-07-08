/***
 * Copyright (C) 2016 Luca Weihs
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef RIM_List
#define RIM_List

#include "RcppArmadillo.h"
#include"ListNode.h"
#include<cstdlib>

namespace RIM {
  template <class T>
  class List {
    private:
      ListNode<T>* sentinal;
      ListNode<T>* current;
      int listLength;

      // We disable copying by making this function private and providing
      // no implementation
      List(const List& that);

    public:
      List() {
        sentinal = new ListNode<T>();
        current = sentinal;
        listLength = 0;
      }

      ~List() {
        RIM::ListNode<T>* nextNode = sentinal->next->next;
        RIM::ListNode<T>* curNode = sentinal->next;
        while(curNode != sentinal) {
          delete curNode;
          curNode = nextNode;
          nextNode = nextNode->next;
        }
        delete sentinal;
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
          Rcpp::stop("Attempting to get first value of empty list.\n");
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
          Rcpp::stop("Attempting to get last value of empty list.\n");
        }
        return(sentinal->prev->value);
      }

      ListNode<T>* currentNode() {
        if (length() == 0) {
          Rcpp::stop("Attempting to get current node of empty list.\n");
        }
        return(current);
      }

      T currentValue() {
        if (length() == 0) {
          Rcpp::stop("Attempting to get current value of empty list.\n");
        }
        return(current->value);
      }

      void removeCurrent() {
        if (length() == 0) {
          Rcpp::stop("Attempting to remove current node of empty list.\n");
        }
        if (current == sentinal) {
          Rcpp::stop("Attempting to remove sentinal, this shouldn't be possible.\n");
        }
        ListNode<T>* prev = current->prev;
        ListNode<T>* next = current->next;
        delete current;
        prev->next = next;
        next->prev = prev;
        current = prev;
        listLength--;
      }

      void appendNode(ListNode<T>* listNode) {
        listNode->next = sentinal;
        listNode->prev = sentinal->prev;
        sentinal->prev->next = listNode;
        sentinal->prev = listNode;
        listLength += 1;
      }

      void appendValue(T value) {
        appendNode(new ListNode<T>(value));
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

      void extend(List<T>* listToAppend) {
        ListNode<T>* otherSentinal = listToAppend->sentinal;
        ListNode<T>* otherCurrent = otherSentinal->next;
        ListNode<T>* otherNext = otherSentinal->next;
        while (otherCurrent != otherSentinal) {
          otherNext = otherCurrent->next;
          appendNode(otherCurrent);
          otherCurrent = otherNext;
        }

        otherSentinal->next = otherSentinal;
        otherSentinal->prev = otherSentinal;
        listToAppend->listLength = 0;
      }

      void joinWithPartition(List<T>* listToMerge, List<int>* partition) {
        if(partition->length() > listToMerge->length()) {
          Rcpp::stop("Attempting merge two lists with a partition that is too long.\n");
        }
        partition->restart();
        this->restart();
        listToMerge->restart();

        int initialLength = this->length();
        int lastPartitionNum = initialLength;
        int curPartitionNum = initialLength;

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

        listToMerge->sentinal->next = listToMerge->sentinal;
        listToMerge->sentinal->prev = listToMerge->sentinal;
        listToMerge->listLength = 0;
      }

  };
}

#endif
