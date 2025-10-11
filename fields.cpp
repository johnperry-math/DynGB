#ifndef FIELDS_CPP_
#define FIELDS_CPP_

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

#include <gsl/gsl>

#include "fields.hpp"
#include "system_constants.hpp"

Prime_Field::Prime_Field(UCOEF_TYPE modulus, bool show_modulus)
    : m{gsl::narrow<COEF_TYPE>(modulus)}, print_modulus{show_modulus} {
  Fi.resize(m);
  Fi[1] = 1;
  Fi[m - 1] = m - 1;
  for (COEF_TYPE a = 2; a < m - 1; ++a) {
    bool negative = true;
    COEF_TYPE s = 1;
    COEF_TYPE t = 0;
    COEF_TYPE u = 0;
    COEF_TYPE v = 1;
    COEF_TYPE c = m;
    COEF_TYPE d = a;
    while (d != 0) {
      COEF_TYPE q = c / d;
      COEF_TYPE r = c - q * d;
      COEF_TYPE w = u * q + s;
      COEF_TYPE x = v * q + t;
      s = u;
      t = v;
      u = w;
      v = x;
      c = d;
      d = r;
      negative = not negative;
    }
    if (negative)
      Fi[a] = m - t;
    else
      Fi[a] = t;
  }
}

Prime_Field::Prime_Field(const Prime_Field &F) : m{F.m} {
  Fi.resize(m);
  for (unsigned i = 0; i < m; ++i)
    Fi[i] = F.Fi[i];
  print_modulus = F.print_modulus;
}

void Prime_Field::set_print_modulus(bool b) { print_modulus = b; }

auto Prime_Field::get_print_modulus() const -> bool { return print_modulus; }

Prime_Field_Element::Prime_Field_Element(const Prime_Field *field)
    : F(field), a{0}, m{F->modulus()} {}

Prime_Field_Element::Prime_Field_Element(COEF_TYPE value,
                                         const Prime_Field *field)
    : F(field), a{value}, m{F->modulus()} {
  if (a > m)
    a %= m;
  while (a < 0)
    a += m;
}

auto Prime_Field_Element::value() const -> COEF_TYPE { return a; }

auto Prime_Field_Element::modulus() const -> COEF_TYPE { return m; }

auto Prime_Field_Element::field() const -> const Prime_Field * { return F; }

auto Prime_Field_Element::like(const Prime_Field_Element &b) const -> bool {
  return m == b.m;
}

auto Prime_Field_Element::inverse() const -> COEF_TYPE { return F->inverse(a); }

// bool Prime_Field_Element::is_zero() const { return a == 0; }

auto Prime_Field_Element::is_one() const -> bool { return a == 1; }

void Prime_Field_Element::negate() { a = m - a; }

void Prime_Field_Element::operator+=(const Prime_Field_Element &other) {
  a += other.a;
  if (a >= m)
    a -= m;
}

void Prime_Field_Element::operator-=(const Prime_Field_Element &other) {
  a -= other.a;
  if (a < 0)
    a += m;
}

void Prime_Field_Element::operator*=(const Prime_Field_Element &other) {
  a *= other.a;
  if (a >= m)
    a %= m;
}

void Prime_Field_Element::operator*=(const COEF_TYPE b) {
  a *= b;
  if (a >= m)
    a %= m;
}

void Prime_Field_Element::operator/=(const Prime_Field_Element &other) {
  a *= other.inverse();
  if (a >= m)
    a %= m;
}

/** Assume the terms have the same modulus. Behavior undefined if not! */
auto operator+(const Prime_Field_Element &a, const Prime_Field_Element &b)
    -> Prime_Field_Element {
  Prime_Field_Element c(a.a, a.F);
  c.a += b.a;
  c.a %= c.m;
  return c;
}

/** Assume the terms have the same modulus. Behavior undefined if not! */
auto operator-(const Prime_Field_Element &a, const Prime_Field_Element &b)
    -> Prime_Field_Element {
  Prime_Field_Element c(a.a, a.F);
  c.a -= b.a;
  while (c.a < 0)
    c.a += c.m;
  return c;
}

/** Assume the terms have the same modulus. Behavior undefined if not! */
auto operator*(const Prime_Field_Element &a, const Prime_Field_Element &b)
    -> Prime_Field_Element {
  Prime_Field_Element c(a.a, a.F);
  c.a *= b.a;
  c.a %= c.m;
  return c;
}

auto operator+(const Prime_Field_Element &a, const int b)
    -> Prime_Field_Element {
  Prime_Field_Element c(a.a, a.F);
  c.a += b;
  c.a %= c.m;
  return c;
}

auto operator-(const Prime_Field_Element &a, const int b)
    -> Prime_Field_Element {
  Prime_Field_Element c(a.a, a.F);
  c.a -= b;
  while (c.a < 0)
    c.a += c.m;
  return c;
}

auto operator*(const Prime_Field_Element &a, const int b)
    -> Prime_Field_Element {
  Prime_Field_Element c(a.a, a.F);
  c.a *= b;
  c.a %= c.m;
  return c;
}

auto operator-(const Prime_Field_Element &a) -> Prime_Field_Element {
  Prime_Field_Element c(-a.a, a.F);
  return c;
}

auto operator<<(ostream &os, const Prime_Field_Element &a) -> ostream & {
  os << a.a;
  if (a.F->get_print_modulus())
    os << '_' << a.m;
  return os;
}

auto Prime_Field::unity() const -> Prime_Field_Element { return {1, this}; }

auto Prime_Field::zero() const -> Prime_Field_Element { return {0, this}; }

#endif