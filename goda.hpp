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
* DynGB is distributed in the hope that it will be useful,                    *
* but WITHOUT ANY WARRANTY; without even the implied warranty of              *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the               *
* GNU General Public License for more details.                                *
*                                                                             *
* You should have received a copy of the GNU General Public License           *
* along with DynGB. If not, see <http://www.gnu.org/licenses/>.               *
\*****************************************************************************/

#include <cstdlib>
#include <iostream>
#include <string>
#include <atomic>

using std::cout; using std::endl; using std::size_t; using std::string;
using std::atomic_flag;

#include "system_constants.hpp"

const int BLOCK_COUNT = 2 << 16;

/**
  @defgroup memorygroup Memory Management
  @brief Classes and other structures related to memory management
  @details These classes create memory pools that allocate large blocks of
    memory, each of which contains an enormous number of subblocks that can
    be obtained and released in @f$O(1)@f$ time
    (or at least much more efficiently than @c new() and @c delete()).
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
    and WGrevlex_Order_Data. It isn't quite @f$O(1)@f$ but the plan is to fix
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
    if (elaborate) {
      cout << name << " allocating " << BLOCK_COUNT << " blocks\n";
      cout << num_blocks + BLOCK_COUNT << " blocks total\n";
    }
    goda_block<TYPE> * result
        = (goda_block<TYPE> *)malloc(BLOCK_COUNT * data_size*sizeof(size_t));
    goda_block<TYPE> * tmp = result;
    for (unsigned i=0; i < BLOCK_COUNT - 1; ++i) {
      tmp[i*data_size].next = &(tmp[(i+1)*data_size]);
    }
    tmp[(BLOCK_COUNT - 1)*data_size].next = nullptr;
    return result;
  }
  /**
    @brief sets allocator up for blocks of @f$n@f$ of type @c TYPE.
    @param n number of objects of type @c TYPE that should be allocated initially
    @param block_name an optional name, for debugging purposes
    @param verbose whether to udpate occasionally during use
    @details The Grading_Order_Data_Allocator can and will allocate new blocks
      of @p n objects when it runs out of room.
  */
  Grading_Order_Data_Allocator(
      NVAR_TYPE n, const string & block_name = "", bool verbose=false
  )
  : data_size(n), name(block_name), elaborate(verbose)
  {
    in_use.clear();
    big_blocks = block = allocate_new_block();
    block->next = nullptr; block += data_size;
  }
  /** @brief releases all memory &mdash; you'd better have freed yours! */
  ~Grading_Order_Data_Allocator() {
    cout << "at termination, " << name << " has " << num_blocks << " blocks\n";
    while (big_blocks != nullptr) {
      goda_block<TYPE> * next_block = big_blocks->next;
      cout << name << " freeing " << BLOCK_COUNT << " blocks\n";
      if (num_blocks < BLOCK_COUNT) num_blocks = 0;
      else num_blocks -= BLOCK_COUNT;
      free(big_blocks);
      big_blocks = next_block;
    }
    cout << name << " destructs with " << num_blocks << " blocks\n";
  }
  /**
    @brief allocates and returns a block of memory
    @return block of memory of the size specified in constructor
  */
  TYPE * get_new_block() {
    while (in_use.test_and_set()) { /* */ }
    if (block == nullptr) {
      block = allocate_new_block();
      block->next = big_blocks; big_blocks = block; block += data_size;
    }
    ++num_blocks;
    if (elaborate)
      cout << name << " increase to " << num_blocks << " blocks in use\n";
    TYPE * result = reinterpret_cast<TYPE *>(block);
    block = block->next;
    in_use.clear();
    return result;
  }
  /**
    @brief returns a block of memory that is no longer needed to the pool
    @param freed_block block to return to the pool
  */
  void return_used_block(TYPE * freed_block) {
    while (in_use.test_and_set()) { /* */ }
    goda_block<TYPE> * new_head = (goda_block<TYPE> *)freed_block;
    --num_blocks;
    if (elaborate)
      cout << name << " decrease to " << num_blocks << " blocks in use\n";
    new_head->next = block;
    block = new_head;
    in_use.clear();
  }
  /**
    @brief indicates how many blocks this has allocated
  */
  void report() {
    cout << name << "has allocated " << num_blocks << " blocks\n";
  }
protected:
  /** @brief how many words to step from one block to the next */
  const unsigned data_size;
  /** @brief pointer to the next free block */
  goda_block<TYPE> * block;
  /** @brief pointer to the superblock of all blocks */
  goda_block<TYPE> * big_blocks;
  /** @brief number of blocks allocated &mdash; mainly for debugging */
  unsigned int num_blocks = 0;
  /** @brief name of allocator &mdash; mainly for debugging */
  const string name = "";
  /** @brief whether to elaborate during certain activities &mdash; mainly for debugging */
  bool elaborate = false;
  /** @brief used for a spinlock for thread safety */
  atomic_flag in_use = ATOMIC_FLAG_INIT;
};

#endif