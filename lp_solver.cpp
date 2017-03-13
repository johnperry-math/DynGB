#ifndef __LP_SOLVER_CPP_
#define __LP_SOLVER_CPP_

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

#include "lp_solver.hpp"

constraint::constraint(NVAR_TYPE num_variables, CONSTR_TYPE coeffs [])
{
  nvars = num_variables;
  coefficients = new CONSTR_TYPE[nvars];
  for (NVAR_TYPE i = 0; i < nvars; ++i)
    coefficients[i] = coeffs[i];
}

constraint::constraint(vector<CONSTR_TYPE> &coeffs)
{
  nvars = coeffs.size();
  coefficients = new CONSTR_TYPE[nvars];
  for (NVAR_TYPE i = 0; i < nvars; ++i)
    coefficients[i] = coeffs[i];
}

constraint::constraint(const constraint &old_constraint)
{
  nvars = old_constraint.nvars;
  coefficients = new CONSTR_TYPE[nvars];
  for (NVAR_TYPE i = 0; i < nvars; ++i)
    coefficients[i] = old_constraint.coefficients[i];
}

// ordering is lexicographic
bool operator < (const constraint &first, const constraint &second)
{
  bool result = !(first == second);
  bool checking = result;
  for (NVAR_TYPE i = 0; checking and i < first.nvars; ++i)
    if (first[i] != second[i])
    {
      checking = false;
      result = first[i] < second[i];
    }
  return result;
}

bool operator == (const constraint &first, const constraint &second)
{
  bool result = true;
  for (NVAR_TYPE i = 0; result and i < first.get_number_of_variables(); ++i) {
    if (first[i] != second[i])
      result = false;
  }
  return result;
}

bool operator != (constraint &first, constraint &second)
{
  bool result = true;
  for (NVAR_TYPE i = 0; result and i < first.get_number_of_variables(); ++i)
    if (first[i] != second[i])
      result = false;
  return !result;
}

// formatted as sum of products of coefficients and variables
ostream & operator<<(ostream & ostr, const constraint &c)
{
  bool first = true;
  ostr << "0 ≤ ";
  for (NVAR_TYPE i = 0; i < c.nvars; ++i)
  {
    if (c[i] < 0)
    {
      if (first) { ostr << '-'; first = false; } else ostr << "- ";
      if (c[i] != -1)
        ostr << -c[i] << 'x' << i << ' ';
      else
        ostr << 'x' << i << ' ';
    }
    else if (c[i] > 0)
    {
      if (first) { first = false; } else ostr << "+ ";
      if (c[i] != 1)
        ostr << c[i] << 'x' << i << ' ';
      else
        ostr << 'x' << i << ' ';
    }
    else
    {
      // do nothing when c[i] == 0
    }
  }
  return ostr;
}

constraint::~constraint()
{
  delete [] coefficients;
}

/**
  @brief memory manager for ray entries
  @ingroup memorygroup
  @details Automatically initialized, but clients need to call the destructor
    when finished.
*/
Grading_Order_Data_Allocator<DEG_TYPE> * doda = nullptr;

unsigned invocations = 0;

inline DEG_TYPE * ray_data_allocation(NVAR_TYPE n) {
  invocations++;
  if (doda == nullptr) doda = new Grading_Order_Data_Allocator<DEG_TYPE>(n);
  return doda->get_new_block();
}

ray::ray(NVAR_TYPE dimension, long direction)
{
  dim = dimension;
  //coords = new RAYENT_TYPE[dim];
  coords = ray_data_allocation(dim);
  for (NVAR_TYPE i = 0; i < dim; ++i) coords[i] = 0;
  if (direction >= 0) {
    coords[direction] = 1;
  }
}

ray::ray(NVAR_TYPE dimension, const RAYENT_TYPE entries [])
{
  dim = dimension;
  //coords = new RAYENT_TYPE[dim];
  coords = ray_data_allocation(dim);
  for (NVAR_TYPE i = 0; i < dim; ++i)
    coords[i] = entries[i];
}

ray::ray(NVAR_TYPE dimension, const EXP_TYPE entries [])
{
  dim = dimension;
  //coords = new RAYENT_TYPE[dim];
  coords = ray_data_allocation(dim);
  for (NVAR_TYPE i = 0; i < dim; ++i)
    coords[i] = entries[i];
}

ray::ray(const vector<RAYENT_TYPE> &entries)
{
  dim = entries.size();
  //coords = new RAYENT_TYPE[dim];
  coords = ray_data_allocation(dim);
  for (NVAR_TYPE i = 0; i < dim; ++i)
    coords[i] = entries[i];
}

ray::ray(const ray &old_ray)
        : dim(old_ray.dim)
{
  //coords = new RAYENT_TYPE[dim];
  static unsigned long invocations;
  coords = ray_data_allocation(dim);
  for (RAYENT_TYPE i = 0; i < dim; ++i)
    coords[i] = old_ray.coords[i];
}

ray::~ray()
{
  //delete [] coords;
  doda->return_used_block(coords);
}

void ray::simplify_ray()
{
  RAYENT_TYPE * w = coords;
  RAYENT_TYPE gcd = 0;
  for (NVAR_TYPE i = 0; i < get_dimension(); ++i)
  {
    if (gcd == 0)
      gcd = w[i];
    else if (w[i] != 0)
    {
      RAYENT_TYPE r = (gcd < w[i] ? gcd : w[i]);
      RAYENT_TYPE s = (gcd < w[i] ? w[i] : gcd);
      while (r != 0)
      {
        RAYENT_TYPE t = s % r;
        s = r;
        r = t;
      }
      gcd = s;
    }
  }
  if (gcd > 1)
    for (NVAR_TYPE i = 0; i < get_dimension(); ++i)
      w[i] /= gcd;
}

DOTPROD_TYPE ray::obtain_dot_product(const constraint &hyperplane) const
{
  DOTPROD_TYPE result = 0;
  const CONSTR_TYPE * coeffs = hyperplane.coeffs();
  for (NVAR_TYPE i = 0; i < dim; ++i)
    //result += coords[i] * hyperplane[i];
    result += coords[i] * coeffs[i];
  return result;
}

ray operator*(RAYENT_TYPE a, ray &r)
{
  NVAR_TYPE d = r.get_dimension();
  //RAYENT_TYPE *coords = new RAYENT_TYPE[d];
  RAYENT_TYPE *coords = ray_data_allocation(d);
  for (NVAR_TYPE i = 0; i < d; ++i)
    coords[i] = a*r[i];
  ray result(d, coords);
  //delete [] coords;
  doda->return_used_block(coords);
  return result;
}

RAYENT_TYPE operator*(const ray &r1, const ray &r2)
{
  RAYENT_TYPE result = 0;
  for (NVAR_TYPE i = 0; i < r1.get_dimension(); ++i)
    result += r1[i]*r2[i];
  return result;
}

RAYENT_TYPE operator*(ray &r1, vector<long> &r2)
{
  RAYENT_TYPE result = 0;
  for (NVAR_TYPE i = 0; i < r1.get_dimension(); ++i)
    result += r1[i]*r2[i];
  return result;
}

RAYENT_TYPE operator*( vector<long> &r1, ray &r2)
{
  RAYENT_TYPE result = 0;
  for (NVAR_TYPE i = 0; i < r1.size(); ++i)
    result += r1[i]*r2[i];
  return result;
}

ray operator+(ray &r1, ray &r2)
{
  NVAR_TYPE d = r1.get_dimension();
  RAYENT_TYPE * coords = ray_data_allocation(d);
  for (NVAR_TYPE i = 0; i < d; ++i)
    coords[i] = r1[i] + r2[i];
  ray result(d, coords);
  doda->return_used_block(coords);
  return result;
}

ray operator-(const ray &r1, const ray &r2)
{
  NVAR_TYPE d = r1.get_dimension();
  RAYENT_TYPE *coords = ray_data_allocation(d);
  for (NVAR_TYPE i = 0; i < d; ++i)
    coords[i] = r1[i] - r2[i];
  ray result(d, coords);
  doda->return_used_block(coords);
  return result;
}

ray ray_sum(const set<ray> &rs)
{
  NVAR_TYPE d = 0;
  RAYENT_TYPE *coords = nullptr;
  for (auto riter = rs.begin(); riter != rs.end(); ++riter)
  {
    ray r = *riter;
    if (coords == nullptr)
    {
      d = r.get_dimension();
      coords = ray_data_allocation(d);
      for (NVAR_TYPE i = 0; i < d; ++i) coords[i] = 0;
    }
    for (NVAR_TYPE i = 0; i < d; ++i)
    {
      coords[i] += r[i];
    }
  }
  ray result(d, coords);
  doda->return_used_block(coords);
  return result;
}

bool operator==(const ray &r1, const ray &r2)
{
  bool result = true;
  for (NVAR_TYPE i = 0; result and i < r1.get_dimension(); ++i)
    result = r1[i] == r2[i];
  return result;
}

bool operator!=(const ray &r1, const ray &r2)
{
  bool result = true;
  for (NVAR_TYPE i = 0; result and i < r1.get_dimension(); ++i)
    result = r1[i] == r2[i];
  return !result;
}

ostream & operator<<(ostream & ostr, const ray &r)
{
  ostr << "( ";
  NVAR_TYPE i;
  for (i = 0; i < r.dim - 1; ++i)
    ostr << r[i] << ", ";
  ostr << r[i] << " )";
  return ostr;
}

ray & ray::operator=(const ray &other)
{
  if (!(*this == other))
  {
    // resize coords if need be
    if (dim != other.dim)
    {
      dim = other.dim;
      //delete [] coords;
      //coords = new RAYENT_TYPE[dim];
      doda->return_used_block(coords);
      coords = ray_data_allocation(dim);
    }
    // copy coords
    for (NVAR_TYPE i = 0; i < dim; ++i)
      coords[i] = other.coords[i];
  }
  return *this;
}

void ray::swap(ray &other)
{
  RAYENT_TYPE tmpval;
  for (NVAR_TYPE i = 0; i < dim; ++i)
  {
    tmpval = coords[i];
    coords[i] = other.coords[i];
    other.coords[i] = tmpval;
  }
}

bool operator<(const ray &first_ray, const ray &second_ray)
{
  bool result = true;
  NVAR_TYPE i = 0;
  bool equal = true;
  while (equal and i < first_ray.dim)
  {
    if (first_ray[i] != second_ray[i])
    {
      equal = false;
      result = first_ray[i] < second_ray[i];
    }
    ++i;
  }
  return result and (not equal);
}

const set<ray> & LP_Solver::get_rays() { return rays; }

#endif