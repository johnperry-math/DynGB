#ifndef __PPL_SOLVER_HPP_
#define __PPL_SOLVER_HPP_

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

#include <ppl.hh>
namespace PPL = Parma_Polyhedra_Library;

#include "system_constants.hpp"

#include "lp_solver.hpp"

namespace LP_Solvers {

/**
  @brief approximate skeleton of a polyhedral cone, using PPL linear solver
  @author John Perry
  \version 1.0
  @date January 2017
  @copyright The University of Southern Mississippi
  @ingroup CLSSolvers
  @details This class serves as an interface to PPL @cite BagnaraHZ08SCP,
      which we can use to find the skeleton to a polyhedral cone.
*/
class PPL_Solver : public LP_Solver {
public:
  /** @name Construction */
  ///@{
  /** @brief initializes solver for @f$ n @f$ variables */
  PPL_Solver(NVAR_TYPE n);
  /** @brief copy constructor (deep copy) */
  PPL_Solver(const PPL_Solver &);
  virtual bool copy(const LP_Solver *) override;
  ///@}
  /** @name Destruction */
  ///@{
  virtual ~PPL_Solver();
  ///@}
  /** @name Basic properties */
  ///@{
  virtual NVAR_TYPE get_dimension() const override { return n; }
  virtual unsigned long get_number_of_constraints() const override { return m; }
  ///@}
  /** @name Modification */
  ///@{
  virtual bool solve(const Constraint &) override;
  virtual bool solve(const vector<Constraint> &) override;
  /** @brief clear the current set of rays and extracts the ones contained in lp */
  virtual void setup_rays();
  ///@}
protected:
  PPL::NNC_Polyhedron * lp; /**< @brief PPL problem interface */
  NVAR_TYPE n; /**< @brief number of variables */
  unsigned m; /**< @brief number of constraints */
  static unsigned instances; /**< @brief number of PPL instances */
  PPL::Variable ** X; /**< @brief array of variables */
  RAYENT_TYPE * ray_data; /**< @brief used to retrieve rays */
};

}

#endif