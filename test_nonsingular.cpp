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

#include <iostream>
#include "monomial_ordering.hpp"

#define ROWDIM 4
#define COLDIM 4

int main() {
  unsigned long **A = new unsigned long * [ROWDIM];
  for (unsigned i = 0; i < ROWDIM; ++i)
    A[i] = new unsigned long [COLDIM];
  /* A[0][0] = A[0][1] = A[0][2] = A[0][3] = 1;
  A[1][0] = A[1][1] = A[1][2] = 1; A[1][3] = 0;
  A[2][0] = A[2][1] = 1; A[2][2] = A[2][3] = 0;
  A[3][0] = 1; A[3][1] = A[3][2] = A[3][3] = 0; */
  A[0][0] = A[0][1] = 0; A[0][2] = A[0][3] = 1;
  A[1][0] = A[1][1] = A[1][2] = 1; A[1][3] = 2;
  A[2][0] = 1; A[2][1] = 2; A[2][2] = A[2][3] = 3;
  A[3][0] = A[3][1] = A[3][2] = 3; A[3][3] = 6;
  std::cout << "Nonsingular? " << nonsingular(ROWDIM, COLDIM, const_cast<const unsigned long **>(A)) << std::endl;
}