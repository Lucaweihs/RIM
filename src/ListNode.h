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

#ifndef RIM_ListNode
#define RIM_ListNode

#include<stdio.h>

namespace RIM {
  template <class T>
  class ListNode {
    private:
      // We disable copying by making this function private and providing
      // no implementation
      ListNode(const ListNode& that);

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
