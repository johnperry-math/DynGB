#ifndef __ALGORITHM_BUCHBERGER_EXPLORER_HPP_
#define __ALGORITHM_BUCHBERGER_EXPLORER_HPP_

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

#include <list>
#include <vector>
#include <iostream>
#include <iterator>

using std::list;
using std::vector;
using std::cout; using std::endl;

#include "system_constants.hpp"

#include "fields.hpp"
#include "monomial.hpp"
#include "polynomial.hpp"
#include "critical_pair.hpp"
#include "normal_strategy.hpp"
#include "polynomial_array.hpp"
#include "polynomial_geobucket.hpp"
#include "polynomial_linked_list.hpp"
#include "polynomial_double_buffered.hpp"

#include "sugar_strategy.hpp"
#include "weighted_sugar_strategy.hpp"

#include "algorithm_buchberger_basic.hpp"

/**
  @brief Alternate implementation of Buchberger&rsquo;s algorithm,
    for parallelization.
  @ingroup GBComputation
  @details Computes a Gr&ouml;bner basis by selecting \c number_to_advance
    pairs, reducing them as completely as possible, comparing their Hilbert
    functions, accepting the one with the smallest Hilbert function, and
    returning the others to the list of pairs.
*/
list<Constant_Polynomial *> buchberger_explorer(
    const vector<Abstract_Polynomial *> &F,
    SPolyCreationFlags method = GEOBUCKETS,
    StrategyFlags strategy = NORMAL_STRATEGY,
    WT_TYPE * strategy_weights = nullptr,
    const int comm_id = 0,
    const int comm_size = 1
);

#endif