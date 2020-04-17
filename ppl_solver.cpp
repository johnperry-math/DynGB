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

#include <gmp.h>

#include "ppl_solver.hpp"

#ifndef __PPL_SOLVER_CPP
#define __PPL_SOLVER_CPP

using PPL::Variable;
using PPL::Generator_System;
using PPL::Linear_Expression;
using PPL::Polyhedron;
using PPL::NNC_Polyhedron;
using PPL::Constraint_System;
using PPL::raw_value;

namespace LP_Solvers {

/** @brief shorthand for <c>PPL::Constraint</c> */
typedef PPL::Constraint PPL_Constraint;
/** @brief shorthand for <c>PPL::Generator</c> */
typedef PPL::Generator PPL_Generator;

unsigned PPL_Solver::instances = 0;

void PPL_Solver::setup_rays() {
  rays.clear();
  if (not lp->is_empty()) {
    const Generator_System & G = lp->generators();
    for (const PPL_Generator & g : G) {
      if (g.is_ray()) {
        for (NVAR_TYPE i = 0; i < n; ++i)
          ray_data[i] = g.coefficient(*(X[i])).get_ui();
        rays.emplace(n, ray_data);
      }
    }
  }
}

PPL_Solver::PPL_Solver(NVAR_TYPE num_vars) {
  n = num_vars;
  ray_data = new RAYENT_TYPE[n];
  if (instances == 0)
    Parma_Polyhedra_Library::initialize();
  ++instances;
  lp = new NNC_Polyhedron(n);
  X = (Variable **)malloc(sizeof(Variable *)*n);
  for (NVAR_TYPE i = 0; i < n; ++i)
    X[i] = new Variable(i);
  Constraint_System cs;
  for (NVAR_TYPE i = 0; i < n; ++i) {
    cs.insert(*(X[i]) >= 0);
  }
  lp->refine_with_constraints(cs);
  m = n;
  setup_rays();
}

PPL_Solver::PPL_Solver(const PPL_Solver & other) {
  m = 0;
  n = other.n;
  rays = other.rays;
  ray_data = new RAYENT_TYPE[n];
  lp = new NNC_Polyhedron(*(other.lp));
  X = (Variable **)malloc(sizeof(Variable *)*n);
  for (NVAR_TYPE i = 0; i < n; ++i)
    X[i] = new Variable(i);
  ++instances;
}

PPL_Solver & PPL_Solver::operator = (const PPL_Solver & other) {
  if (this != &other) copy(&other);
  return *this;
}

bool PPL_Solver::copy(const LP_Solver * old_solver) {
  const PPL_Solver * other = dynamic_cast<const PPL_Solver *>(old_solver);
  if (other != nullptr) {
    n = other->n;
    rays = other->rays;
    delete [] ray_data;
    delete lp;
    //for (NVAR_TYPE i = 0; i < n; ++i)
    //  delete X[i];
    //free(X);
    ray_data = new RAYENT_TYPE[n];
    lp = new NNC_Polyhedron(*(other->lp));
    //X = (Variable **)malloc(sizeof(Variable *)*n);
    //for (NVAR_TYPE i = 0; i < n; ++i)
    //  X[i] = new Variable(i);
    ++instances;
  }
  return (other != nullptr);
}

PPL_Solver::~PPL_Solver() {
  delete [] ray_data;
  delete lp;
  for (NVAR_TYPE i = 0; i < n; ++i)
    delete X[i];
  free(X);
  --instances;
  if (instances == 0)
    Parma_Polyhedra_Library::finalize();
}

bool PPL_Solver::solve(const Constraint & c) {
  Linear_Expression ineq;
  for (NVAR_TYPE i = 0; i < n; ++i) {
    if (c[i] != 0)
      ineq += Linear_Expression(c[i]*(*X[i]));
  }
  PPL_Constraint pc(ineq >= 1);
  lp->refine_with_constraint(pc);
  setup_rays();
  return (rays.size() > 0);
}

bool PPL_Solver::solve(const vector<Constraint> &C) {
  Constraint_System cs;
  for (const Constraint & c : C) {
    Linear_Expression ineq;
    for (NVAR_TYPE i = 0; i < n; ++i) {
      if (c[i] != 0)
        ineq += Linear_Expression(c[i]*(*X[i]));
    }
    PPL_Constraint pc(ineq >= 1);
    cs.insert(pc);
  }
  lp->refine_with_constraints(cs);
  setup_rays();
  return (rays.size() > 0);
}

ostream & operator << (ostream & ostr, const PPL_Solver &skel)
{
  // header, start constraints
  ostr << "Skeleton defined by constraints" << endl;
  /*for (auto & c : skel.lp->constraints()) {
    c.ascii_dump(ostr);
  }*/
  ostr << "has " << skel.rays.size() << " rays" << endl;
  for (auto & r : skel.rays) {
    ostr << '\t' << r << endl;
  }
  ostr << "connected in edges" << endl;
  for (auto & g : skel.lp->generators()) {
    if (g.is_line()) {
      g.ascii_dump(ostr);
    }
  }
  // footer
  ostr << "End of skeleton" << endl;
  return ostr;
}

}

#endif