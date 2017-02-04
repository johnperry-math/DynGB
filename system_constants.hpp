#ifndef __CONSTANTS_HPP_
#define __CONSTANTS_HPP_

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

/** Constants for the system */

#include <cstdint>
#include <gmpxx.h>

// for coefficients

#define COEF_TYPE int64_t
#define UCOEF_TYPE uint64_t
#define MPQCOEF_TYPE mpq_class
#define MPZCOEF_TYPE mpz_class

// for monomials
#define EXP_TYPE uint32_t
#define NVAR_TYPE uint64_t

// size of a degree
#define DEG_TYPE uint64_t

// entries in a ray
#define RAYENT_TYPE uint64_t

// entries in a constraint
#define CONSTR_TYPE int32_t

// size of dot product
#define DOTPROD_TYPE int64_t

// entries in a weighted monomial ordering -- needs to be the same as RAY_TYPE
#define WT_TYPE RAYENT_TYPE

#endif