#ifndef __GLPK_SOLVER_CPP_
#define __GLPK_SOLVER_CPP_

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

#include <cmath>

#include "glpk_solver.hpp"

namespace LP_Solvers {

#define MIN_X 0.01

GLPK_Solver::GLPK_Solver(NVAR_TYPE nx) {
  dirty = true;
  n = nx; m = 0;
  lp = glp_create_prob();
  glp_set_obj_dir(lp, GLP_MIN);
  glp_add_cols(lp, n);
  for (NVAR_TYPE i = 1; i <= n; ++i) {
    glp_set_col_bnds(lp, i, GLP_LO, MIN_X, 0.0);
    glp_set_obj_coef(lp, i, 1.0);
  }
  glp_init_smcp(&smcp);
  glp_term_out(GLP_OFF);
  smcp.msg_lev = GLP_MSG_OFF;
  smcp.presolve = GLP_ON;
  row_data = new double[1+n]; // GLPK wants to start from index 1
  row_indx = new int[1+n];
  ray_data = new RAYENT_TYPE[n];
}

GLPK_Solver::GLPK_Solver(const GLPK_Solver &other) {
  dirty = other.dirty;
  m = other.m; n = other.n;
  row_data = new double[1+n];
  row_indx = new int[1+n];
  ray_data = new RAYENT_TYPE[n];
  lp = glp_create_prob();
  glp_copy_prob(lp, other.lp, GLP_OFF);
  smcp = other.smcp;
  rays = other.rays;
}

bool GLPK_Solver::copy(const LP_Solver * old_solver) {
  const GLPK_Solver * other = dynamic_cast<const GLPK_Solver *>(old_solver);
  if (other != nullptr) {
    dirty = other->dirty;
    if (n != other->n) {
      delete [] row_data; delete [] row_indx; delete [] ray_data;
      row_data = new double[1+n];
      row_indx = new int[1+n];
      ray_data = new RAYENT_TYPE[n];
    }
    m = other->m; n = other->n;
    glp_delete_prob(lp);
    lp = glp_create_prob();
    glp_copy_prob(lp, other->lp, GLP_OFF);
    smcp = other->smcp;
    rays = other->rays;
  }
  return (other != nullptr);
}

GLPK_Solver::~GLPK_Solver() {
  delete [] row_indx;
  delete [] row_data;
  delete [] ray_data;
  glp_delete_prob(lp);
}

bool GLPK_Solver::solve(const vector<Constraint> & newvecs) {
  dirty = true;
  int glp_result = 0;
  int new_m = newvecs.size();
  if (new_m > 0) {
    glp_add_rows(lp, new_m);
    for (int i = 0; i < new_m; ++i) {
      int num_valid = 0;
      for (int k = 0; k < n; ++k) {
        if (newvecs[i][k] != 0) {
          ++num_valid;
          row_indx[num_valid] = 1 + k;
          row_data[num_valid] = newvecs[i][k];
        }
      }
      glp_set_mat_row(lp, 1 + m + i, num_valid, row_indx, row_data);
      glp_set_row_bnds(lp, 1 + m + i, GLP_LO, MIN_X, 0.0);
    }
    if (m != 0) {
      smcp.presolve = GLP_OFF;
      glp_warm_up(lp);
    }
    m += new_m;
    glp_result = glp_simplex(lp, &smcp);
    if (glp_result == 0)
      glp_result = glp_exact(lp, nullptr);
  }
  int status = glp_get_status(lp);
  return (glp_result == 0 and (status == GLP_OPT or status == GLP_FEAS));
}

bool GLPK_Solver::solve(const Constraint & newvec) {
  dirty = true;
  glp_add_rows(lp, 1);
  int num_valid = 0;
  for (int k = 0; k < n; ++k) {
    if (newvec[k] != 0) {
      ++num_valid;
      row_indx[num_valid] = 1 + k;
      row_data[num_valid] = newvec[k];
    }
  }
  glp_set_mat_row(lp, 1 + m, num_valid, row_indx, row_data);
  glp_set_row_bnds(lp, 1 + m, GLP_LO, MIN_X, 0.0);
  if (m != 0) {
    smcp.presolve = GLP_OFF;
    glp_warm_up(lp);
  }
  m += 1;
  int glp_result = glp_simplex(lp, &smcp);
  smcp.presolve = GLP_OFF;
  if (glp_result == 0)
    glp_result = glp_exact(lp, nullptr);
  int status = glp_get_status(lp);
  return (glp_result == 0 and (status == GLP_OPT or status == GLP_FEAS));
}

const set<Ray> & GLPK_Solver::get_rays() {
  if (dirty) {
    rays.clear();
    // the next few lines add a row
    // that pushes the solution beyond the actual minimum
    double curr_min = glp_get_obj_val(lp);
    glp_add_rows(lp, 1);
    for (NVAR_TYPE i = 0; i < n; ++i) {
      row_indx[1+i] = 1+i;
      row_data[1+i] = 1.0;
    }
    glp_set_mat_row(lp, 1 + m, n, row_indx, row_data);
    glp_set_row_bnds(lp, 1 + m, GLP_DB, curr_min + 1, curr_min + 2);
    // now we vary the objective function to minimize each variable
    for (int i = 0; i < n; ++i)
      glp_set_obj_coef(lp, 1 + i, 0.0);
    for (int i = 0; i < n; ++i) {
      glp_set_obj_coef(lp, 1 + i, 1.0);
      if (i > 0)
        glp_set_obj_coef(lp, i, 0.0);
      // minimize first
      glp_simplex(lp, &smcp);
      glp_get_status(lp);
      glp_exact(lp, nullptr);
      for (NVAR_TYPE j = 0; j < n; ++j)
        ray_data[j] = static_cast<RAYENT_TYPE>(
            round(1/MIN_X * glp_get_col_prim(lp, 1 + j))
        );
      rays.emplace(n, ray_data);
      //now maximize
      glp_set_obj_dir(lp, GLP_MAX);
      glp_simplex(lp, &smcp);
      glp_get_status(lp);
      glp_exact(lp, nullptr);
      for (NVAR_TYPE j = 0; j < n; ++j)
        ray_data[j] = static_cast<RAYENT_TYPE>(
            round(1/MIN_X * glp_get_col_prim(lp, 1 + j))
        );
      rays.emplace(n, ray_data);
      glp_set_obj_dir(lp, GLP_MIN);
    }
    // fix the objective function
    for (unsigned i = 0; i < n; ++i)
      glp_set_obj_coef(lp, 1 + i, 1.0);
    // delete the row we added
    row_indx[1] = 1 + m;
    glp_set_row_stat(lp, 1 + m, GLP_BS);
    glp_del_rows(lp, 1, row_indx);
    glp_std_basis(lp);
    dirty = false;
  }
  return rays;
}

}

#endif