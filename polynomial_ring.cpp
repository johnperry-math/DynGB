#ifndef __POLYNOMIAL_RING_CPP_
#define __POLYNOMIAL_RING_CPP_

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

#include <cstring>

#include "polynomial_ring.hpp"

Polynomial_Ring::Polynomial_Ring(NVAR_TYPE num_vars,
  Prime_Field & field,
  string * new_names
) : F(field) {
  n = num_vars;
  names = nullptr;
  if (new_names != nullptr)
    set_names(new_names, n);
}

Polynomial_Ring::~Polynomial_Ring() {
  if (names != nullptr)
    delete [] names;
}

bool Polynomial_Ring::set_names(string * new_names, NVAR_TYPE length) {
  bool result = false;
  if (new_names != nullptr) {
    if (names == nullptr)
      names = new string [n];
    for (NVAR_TYPE i = 0; i < n; ++i) {
      names[i] = new_names[i];
    }
    result = true;
  }
  return result;
}

NVAR_TYPE Polynomial_Ring::number_of_variables() const { return n; }

Indeterminate * Polynomial_Ring::indeterminates() {
  Indeterminate * result = static_cast<Indeterminate *>(
      malloc(sizeof(Indeterminate) * n)
  );
  for (NVAR_TYPE i = 0; i < n; ++i)
    result[i] = Indeterminate(*this, i);
  return result;
}

Prime_Field & Polynomial_Ring::ground_field() const { return F; }

const string Polynomial_Ring::name(NVAR_TYPE i) const {
  if (names == nullptr) {
    char my_name[2];
    my_name[0] = 'x';
    my_name[1] = i + char(i);
    return my_name;
  } else
    return names[i];
}

#endif