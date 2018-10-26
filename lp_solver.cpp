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
* DynGB is distributed in the hope that it will be useful,                    *
* but WITHOUT ANY WARRANTY; without even the implied warranty of              *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the               *
* GNU General Public License for more details.                                *
*                                                                             *
* You should have received a copy of the GNU General Public License           *
* along with DynGB. If not, see <http://www.gnu.org/licenses/>.               *
\*****************************************************************************/

#include "lp_solver.hpp"
#include "glpk_solver.hpp"
#include "ppl_solver.hpp"
#include "skeleton.hpp"

namespace LP_Solvers {

Constraint::Constraint(NVAR_TYPE num_variables, CONSTR_TYPE coeffs [])
{
  nvars = num_variables;
  coefficients = new CONSTR_TYPE[nvars];
  for (NVAR_TYPE i = 0; i < nvars; ++i)
    coefficients[i] = coeffs[i];
}

Constraint::Constraint(vector<CONSTR_TYPE> &coeffs)
{
  nvars = coeffs.size();
  coefficients = new CONSTR_TYPE[nvars];
  for (NVAR_TYPE i = 0; i < nvars; ++i)
    coefficients[i] = coeffs[i];
}

Constraint::Constraint(const Constraint &old_constraint)
{
  nvars = old_constraint.nvars;
  coefficients = new CONSTR_TYPE[nvars];
  for (NVAR_TYPE i = 0; i < nvars; ++i)
    coefficients[i] = old_constraint.coefficients[i];
}

Constraint & Constraint::operator = (const Constraint & other) {
  if (this != &other) {
    nvars = other.nvars;
    coefficients = new CONSTR_TYPE[nvars];
    for (NVAR_TYPE i = 0; i < nvars; ++i)
      coefficients[i] = other.coefficients[i];
  }
  return *this;
}

// ordering is lexicographic
bool operator < (const Constraint &first, const Constraint &second)
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

bool operator == (const Constraint &first, const Constraint &second)
{
  bool result = true;
  for (NVAR_TYPE i = 0; result and i < first.get_number_of_variables(); ++i) {
    if (first[i] != second[i])
      result = false;
  }
  return result;
}

bool operator != (const Constraint &first, const Constraint &second)
{
  bool result = true;
  for (NVAR_TYPE i = 0; result and i < first.get_number_of_variables(); ++i)
    if (first[i] != second[i])
      result = false;
  return !result;
}

// formatted as sum of products of coefficients and variables
ostream & operator<<(ostream & ostr, const Constraint &c)
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

Constraint::~Constraint()
{
  delete [] coefficients;
}

Grading_Order_Data_Allocator<DEG_TYPE> * doda = nullptr;

/** @brief used to count the number of invocations of ray_data_allocation() */
unsigned invocations = 0;

/**
  @brief allocates data for a ray
  @param n number of entries needed for the ray
  @return a new block of memory for a ray
*/
inline DEG_TYPE * ray_data_allocation(NVAR_TYPE n) {
  invocations++;
  if (doda == nullptr) doda = new Grading_Order_Data_Allocator<DEG_TYPE>(n, "doda");
  return doda->get_new_block();
}

Ray::Ray(NVAR_TYPE dimension, long direction)
{
  dim = dimension;
  //coords = new RAYENT_TYPE[dim];
  coords = ray_data_allocation(dim);
  for (NVAR_TYPE i = 0; i < dim; ++i) coords[i] = 0;
  if (direction >= 0) {
    coords[direction] = 1;
  }
}

Ray::Ray(NVAR_TYPE dimension, const RAYENT_TYPE entries [])
{
  dim = dimension;
  //coords = new RAYENT_TYPE[dim];
  coords = ray_data_allocation(dim);
  for (NVAR_TYPE i = 0; i < dim; ++i)
    coords[i] = entries[i];
}

Ray::Ray(NVAR_TYPE dimension, const EXP_TYPE entries [])
{
  dim = dimension;
  //coords = new RAYENT_TYPE[dim];
  coords = ray_data_allocation(dim);
  for (NVAR_TYPE i = 0; i < dim; ++i)
    coords[i] = entries[i];
}

Ray::Ray(const vector<RAYENT_TYPE> &entries)
{
  dim = entries.size();
  //coords = new RAYENT_TYPE[dim];
  coords = ray_data_allocation(dim);
  for (NVAR_TYPE i = 0; i < dim; ++i)
    coords[i] = entries[i];
}

Ray::Ray(const Ray &old_ray)
        : dim(old_ray.dim)
{
  //coords = new RAYENT_TYPE[dim];
  coords = ray_data_allocation(dim);
  for (RAYENT_TYPE i = 0; i < dim; ++i)
    coords[i] = old_ray.coords[i];
}

Ray::~Ray()
{
  //delete [] coords;
  doda->return_used_block(coords);
}

void Ray::simplify_ray()
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

DOTPROD_TYPE Ray::obtain_dot_product(const Constraint &hyperplane) const
{
  DOTPROD_TYPE result = 0;
  const CONSTR_TYPE * coeffs = hyperplane.coeffs();
  for (NVAR_TYPE i = 0; i < dim; ++i)
    //result += coords[i] * hyperplane[i];
    result += coords[i] * coeffs[i];
  return result;
}

Ray operator*(const RAYENT_TYPE a, const Ray &r)
{
  NVAR_TYPE d = r.get_dimension();
  //RAYENT_TYPE *coords = new RAYENT_TYPE[d];
  RAYENT_TYPE *coords = ray_data_allocation(d);
  for (NVAR_TYPE i = 0; i < d; ++i)
    coords[i] = a*r[i];
  Ray result(d, coords);
  //delete [] coords;
  doda->return_used_block(coords);
  return result;
}

RAYENT_TYPE operator*(const Ray &r1, const Ray &r2)
{
  RAYENT_TYPE result = 0;
  for (NVAR_TYPE i = 0; i < r1.get_dimension(); ++i)
    result += r1[i]*r2[i];
  return result;
}

RAYENT_TYPE operator*(const Ray &r1, const vector<long> &r2)
{
  RAYENT_TYPE result = 0;
  for (NVAR_TYPE i = 0; i < r1.get_dimension(); ++i)
    result += r1[i]*r2[i];
  return result;
}

RAYENT_TYPE operator*( const vector<long> & r1, const Ray &r2)
{
  RAYENT_TYPE result = 0;
  for (NVAR_TYPE i = 0; i < r1.size(); ++i)
    result += r1[i]*r2[i];
  return result;
}

Ray operator+(const Ray &r1, const Ray &r2)
{
  NVAR_TYPE d = r1.get_dimension();
  RAYENT_TYPE * coords = ray_data_allocation(d);
  for (NVAR_TYPE i = 0; i < d; ++i)
    coords[i] = r1[i] + r2[i];
  Ray result(d, coords);
  doda->return_used_block(coords);
  return result;
}

Ray operator-(const Ray &r1, const Ray &r2)
{
  NVAR_TYPE d = r1.get_dimension();
  RAYENT_TYPE *coords = ray_data_allocation(d);
  for (NVAR_TYPE i = 0; i < d; ++i)
    coords[i] = r1[i] - r2[i];
  Ray result(d, coords);
  doda->return_used_block(coords);
  return result;
}

Ray ray_sum(const set<Ray> &rs)
{
  NVAR_TYPE d = 0;
  RAYENT_TYPE *coords = nullptr;
  for (auto riter = rs.begin(); riter != rs.end(); ++riter)
  {
    Ray r = *riter;
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
  Ray result(d, coords);
  doda->return_used_block(coords);
  return result;
}

bool operator==(const Ray &r1, const Ray &r2)
{
  bool result = true;
  for (NVAR_TYPE i = 0; result and i < r1.get_dimension(); ++i)
    result = r1[i] == r2[i];
  return result;
}

bool operator!=(const Ray &r1, const Ray &r2)
{
  bool result = true;
  for (NVAR_TYPE i = 0; result and i < r1.get_dimension(); ++i)
    result = r1[i] == r2[i];
  return !result;
}

ostream & operator<<(ostream & ostr, const Ray &r)
{
  ostr << "( ";
  NVAR_TYPE i;
  for (i = 0; i < r.dim - 1; ++i)
    ostr << r[i] << ", ";
  ostr << r[i] << " )";
  return ostr;
}

Ray & Ray::operator=(const Ray &other)
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

void Ray::swap(Ray &other)
{
  for (NVAR_TYPE i = 0; i < dim; ++i)
  {
    RAYENT_TYPE tmpval = coords[i];
    coords[i] = other.coords[i];
    other.coords[i] = tmpval;
  }
}

bool operator<(const Ray &first_ray, const Ray &second_ray)
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

const set<Ray> & LP_Solver::get_rays() const { return rays; }

ostream & operator<<(ostream & os, const LP_Solver & s) {
  const Skeleton * Skel_src;
  const GLPK_Solver * GLPK_src;
  const PPL_Solver * PPL_src;
  if ((Skel_src = dynamic_cast<const Skeleton *>(&s)) != nullptr)
    os << *Skel_src;
  else if ((GLPK_src = dynamic_cast<const GLPK_Solver *>(&s)) != nullptr)
    os << *GLPK_src;
  else if ((PPL_src = dynamic_cast<const PPL_Solver *>(&s)) != nullptr)
    os << *PPL_src;
  return os;
}

}

#endif