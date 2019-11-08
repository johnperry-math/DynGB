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
* DynGB is distributed in the hope that it will be useful,                    *
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

#define COEF_TYPE int32_t
#define UCOEF_TYPE uint32_t // see OVERFLOW_MASK below
#define MPQCOEF_TYPE mpq_class
#define MPZCOEF_TYPE mpz_class

// to catch overflow; this should be the same size as UCOEF_TYPE
const UCOEF_TYPE OVERFLOW_MASK = ((UCOEF_TYPE )1) << 31;

// for monomials
#define EXP_TYPE uint16_t
#define NVAR_TYPE uint16_t
#define MASK_SIZE 256

// size of a degree
#define DEG_TYPE uint32_t

// entries in a ray
#define RAYENT_TYPE uint32_t

// entries in a constraint
#define CONSTR_TYPE int32_t

// size of dot product
#define DOTPROD_TYPE int32_t

// entries in a weighted monomial ordering -- needs to be the same as RAY_TYPE
#define WT_TYPE RAYENT_TYPE

// ordering for dynamic computation
#define ORDERING_TYPE WGrevlex

#endif