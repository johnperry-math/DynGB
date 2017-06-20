#include <iostream>
using std::cout; using std::endl;

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

#include "fields.hpp"

int main()
{
  Prime_Field F7(7);
  Prime_Field_Element one_mod_7(1, &F7);
  cout << one_mod_7 << ":\t" << one_mod_7.inverse() << endl;
  unsigned p = 11;
  Prime_Field Fp(p);
  Prime_Field_Element sixteen_mod_p(16, &Fp);
  Prime_Field_Element sixteen_mod_p_inverse(sixteen_mod_p.inverse(), &Fp);
  cout << sixteen_mod_p << ":\t" << sixteen_mod_p_inverse << endl;
  cout << "verification: " << sixteen_mod_p * sixteen_mod_p_inverse << endl;
  p = 21;
  Prime_Field Fq(p);
  Prime_Field_Element nineteen_mod_p(19, &Fq);
  Prime_Field_Element nineteen_mod_p_inverse(nineteen_mod_p.inverse(), &Fq);
  cout << nineteen_mod_p << ":\t" << nineteen_mod_p_inverse << endl;
  cout << "verification: " << nineteen_mod_p * nineteen_mod_p_inverse << endl;
}