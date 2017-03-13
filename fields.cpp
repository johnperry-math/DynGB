#ifndef __FIELDS_CPP_
#define __FIELDS_CPP_

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

#include "fields.hpp"

Prime_Field::Prime_Field(UCOEF_TYPE modulus, bool show_modulus) {
  m = modulus;
  Fi = new COEF_TYPE [m] {0};
  Fi[1] = 1;
  Fi[m-1] = m - 1;
  for (UCOEF_TYPE a = 2; a < m - 1; ++a) {
    bool negative = true;
    COEF_TYPE s = 1; COEF_TYPE t = 0;
    COEF_TYPE u = 0; COEF_TYPE v = 1;
    COEF_TYPE c = (COEF_TYPE )m; COEF_TYPE d = a;
    while (d != 0) {
      COEF_TYPE q = c / d;   COEF_TYPE r = c - q*d;
      COEF_TYPE w = u*q + s; COEF_TYPE x = v*q + t;
      s = u; t = v; u = w; v = x;
      c = d; d = r;
      negative = not negative;
    }
    if (negative) Fi[a] = m - t; else Fi[a] = t;
  }
  print_modulus = show_modulus;
}

Prime_Field::Prime_Field(const Prime_Field & F) {
  m = F.m;
  Fi = new COEF_TYPE [m];
  for (unsigned i = 0; i < m; ++i)
    Fi[i] = F.Fi[i];
  print_modulus = F.print_modulus;
}

Prime_Field::~Prime_Field() { delete [] Fi; }

COEF_TYPE Prime_Field::inverse(COEF_TYPE a) {
  return Fi[a];
}

void Prime_Field::set_print_modulus(bool b) { print_modulus = b; }

bool Prime_Field::get_print_modulus() { return print_modulus; }

Prime_Field_Element::Prime_Field_Element(Prime_Field *field)
    : F(field)
{ a = 0; m = F->modulus(); }

Prime_Field_Element::Prime_Field_Element(COEF_TYPE value, Prime_Field *field)
    : F(field)
{
  a = value;
  m = F->modulus();
  if (a > m) a %= m;
  while (a < 0) a += m;
}

COEF_TYPE Prime_Field_Element::value() const { return a; }

unsigned Prime_Field_Element::modulus() const { return m; }

Prime_Field * Prime_Field_Element::field() const { return F; }

bool Prime_Field_Element::like(const Prime_Field_Element & b) const { return m == b.m; }

COEF_TYPE Prime_Field_Element::inverse() const { return F->inverse(a); }

//bool Prime_Field_Element::is_zero() const { return a == 0; }

bool Prime_Field_Element::is_one() const { return a == 1; }

void Prime_Field_Element::negate() { a = m - a; }

void Prime_Field_Element::operator +=(const Prime_Field_Element &other) {
  a += other.a;
  if (a >= m) a -= m;
}

void Prime_Field_Element::operator -=(const Prime_Field_Element &other) {
  a -= other.a;
  if (a < 0) a += m;
}

void Prime_Field_Element::operator *=(const Prime_Field_Element &other) {
  a *= other.a;
  if (a >= m) a %= m;
}

void Prime_Field_Element::operator *=(const COEF_TYPE b) {
  a *= b;
  if (a >= m) a %= m;
}

void Prime_Field_Element::operator /=(const Prime_Field_Element &other) {
  a *= other.inverse();
  if (a >= m) a %= m;
}

/** Assume the terms have the same modulus. Behavior undefined if not! */
Prime_Field_Element operator+(const Prime_Field_Element &a,
                               const Prime_Field_Element &b)
{
  Prime_Field_Element c(a.a, a.F);
  c.a += b.a;
  c.a %= c.m;
  return c;
}

/** Assume the terms have the same modulus. Behavior undefined if not! */
Prime_Field_Element operator-(const Prime_Field_Element &a,
                               const Prime_Field_Element &b)
{
  Prime_Field_Element c(a.a, a.F);
  c.a -= b.a;
  while (c.a < 0) c.a += c.m;
  return c;
}

/** Assume the terms have the same modulus. Behavior undefined if not! */
Prime_Field_Element operator*(const Prime_Field_Element &a,
                               const Prime_Field_Element &b)
{
  Prime_Field_Element c(a.a, a.F);
  c.a *= b.a;
  c.a %= c.m;
  return c;
}

Prime_Field_Element operator+(const Prime_Field_Element &a,
                               const int b)
{
  Prime_Field_Element c(a.a, a.F);
  c.a += b;
  c.a %= c.m;
  return c;
}

Prime_Field_Element operator-(const Prime_Field_Element &a,
                               const int b)
{
  Prime_Field_Element c(a.a, a.F);
  c.a -= b;
  while (c.a < 0) c.a += c.m;
  return c;
}

Prime_Field_Element operator*(const Prime_Field_Element &a,
                               const int b)
{
  Prime_Field_Element c(a.a, a.F);
  c.a *= b;
  c.a %= c.m;
  return c;
}

Prime_Field_Element operator-(const Prime_Field_Element &a)
{
  Prime_Field_Element c(-a.a, a.F);
  return c;
}

ostream & operator <<(ostream & os, const Prime_Field_Element &a)
{
  os << a.a;
  if (a.F->get_print_modulus())
    os << '_' << a.m;
  return os;
}

Prime_Field_Element Prime_Field::unity() {
  return Prime_Field_Element(1,this);
}

Prime_Field_Element Prime_Field::zero() {
  return Prime_Field_Element(0,this);
}

#endif