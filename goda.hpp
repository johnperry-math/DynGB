#ifndef __GODA_HPP_
#define __GODA_HPP_

/*****************************************************************************\
* This file is part of DynGB.                                                 *
*                                                                             *
* DynGB is free software: you can redistribute it and/or modify               *
* it under the terms of the GNU General Public License as published by        *
* the Free Software Foundation, either version 2 of the License, or           *
* (at your option) any later version.                                         *
*                                                                             *
* Foobar is distributed in the hope that it will be useful,                   *
* but WITHOUT ANY WARRANTY; without even the implied warranty of              *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the               *
* GNU General Public License for more details.                                *
*                                                                             *
* You should have received a copy of the GNU General Public License           *
* along with DynGB. If not, see <http://www.gnu.org/licenses/>.               *
\*****************************************************************************/

#include <cstdlib>
#include <iostream>

using std::cout; using std::endl;

#include "system_constants.hpp"

/**
  @defgroup memorygroup Memory Management
  @brief Classes and other structures related to memory management
  @details These classes create memory pools that allocate large blocks of
    memory, each of which contains an enormous number of subblocks that can
    be obtained and released in \f$O(1)\f$ time
    (or at least much more efficiently than \c new() and \c delete()).
*/

/**
  @ingroup memorygroup
  @union goda_block
  @author John Perry
  @date 2016
  @brief heart of the memory pool allocator (Grading_Order_Data_Allocator)
*/
template <typename TYPE>
union goda_block {
  /** @brief pointer to the next free block */
  goda_block * next;
  /** @brief the data contained in this block */
  TYPE * data;
};

/**
  @ingroup memorygroup
  @class Grading_Order_Data_Allocator
  @author John Perry
  @date 2016
  @brief special memory pool allocator for Grevlex_Order_Data and
    WGrevlex_Order_Data
  @details This is a quick-n-dirty memory pool allocator for Grevlex_Order_Data
    and WGrevlex_Order_Data. It isn't quite \f$O(1)\f$ but the plan is to fix
    that eventually (I have to locate my notes from class).
  @warning This is initialized to a certain number of variables, which cannot
    change. If you want to do this for a different number of variables,
    create another allocator.
*/
template <typename TYPE>
class Grading_Order_Data_Allocator {
public:
  /**
    @brief allocates a new superblock of almost 10000 blocks
    @return the superblock
  */
  goda_block<TYPE> * allocate_new_block() {
    goda_block<TYPE> * result
        = (goda_block<TYPE> *)malloc(10000 * data_size*sizeof(TYPE)*sizeof(long)/sizeof(TYPE));
    goda_block<TYPE> * tmp = result;
    for (unsigned i=0; i < 9998; ++i) {
      tmp->next = &(tmp[data_size]);
      tmp += data_size;
    }
    tmp->next = nullptr;
    return result;
  }
  /**
    @brief sets allocator up for blocks of \f$n\f$ of type @c TYPE.
    @param n number of objects of type @c TYPE that should be allocated initially
    @details The Grading_Order_Data_Allocator can and will allocate new blocks
      of @p n objects when it runs out of room.
  */
  Grading_Order_Data_Allocator(NVAR_TYPE n)
  : data_size(n)
  {
    big_blocks = block = allocate_new_block();
    block->next = nullptr; block += data_size;
  }
  /** @brief releases all memory &mdash; you'd better have freed yours! */
  ~Grading_Order_Data_Allocator() {
    while (big_blocks != nullptr) {
      goda_block<TYPE> * next_block = big_blocks->next;
      free(big_blocks);
      big_blocks = next_block;
    }
  }
  /**
    @brief allocates and returns a block of memory
    @return block of memory of the size specified in constructor
  */
  TYPE * get_new_block() {
    if (block == nullptr) {
      block = allocate_new_block();
      block->next = big_blocks; big_blocks = block; block += data_size;
    }
    TYPE * result = (TYPE *)block;
    block = block->next;
    return result;
  }
  /**
    @brief returns a block of memory that is no longer needed to the pool
    @param freed_block block to return to the pool
  */
  void return_used_block(TYPE * freed_block) {
    goda_block<TYPE> * new_head = (goda_block<TYPE> *)freed_block;
    new_head->next = block;
    block = new_head;
  }
protected:
  /** @brief how many words to step from one block to the next */
  const unsigned data_size;
  /** @brief pointer to the next free block */
  goda_block<TYPE> * block;
  /** @brief pointer to the superblock of all blocks */
  goda_block<TYPE> * big_blocks;
};

#endif