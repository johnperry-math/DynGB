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

#include <set>
#include <cstring>
#include <iostream>

using std::set;
using std::cout; using std::endl;

#include "system_constants.hpp"

#include "goda.hpp"
#include "cyclic_n.hpp"
#include "polynomial.hpp"
#include "strategies.hpp"
#include "dynamic_engine.hpp"
using LP_Solvers::doda;
#include "monomial_ordering.hpp"
#include "particular_orderings.hpp"
#include "polynomial_linked_list.hpp"
#include "algorithm_buchberger_dynamic.hpp"

extern Monomial_Ordering * generic_grevlex_ptr;
extern Grading_Order_Data_Allocator<WT_TYPE> * goda;
extern Grading_Order_Data_Allocator<DEG_TYPE> * doda;
extern Grading_Order_Data_Allocator<EXP_TYPE> * moda;
extern Grading_Order_Data_Allocator<Monomial> * monoda;
extern Grading_Order_Data_Allocator<Monomial_Node> * monododa;
extern Grading_Order_Data_Allocator<CachedWGrevlex_Ordering> * woda;

// Forward declarations
bool meaningful_arguments(
  int, char **, bool &, int &, int &, SPolyCreationFlags &,
  StrategyFlags &, Dynamic_Heuristic &, DynamicSolver &
);

void give_help();

int main(int argc, char *argv[]) {
  bool homog;
  int modulus, numvars;
  SPolyCreationFlags method;
  StrategyFlags strategy = StrategyFlags::SUGAR_STRATEGY;
  Dynamic_Heuristic heuristic = Dynamic_Heuristic::ORD_HILBERT_THEN_DEG;
  DynamicSolver solver = SKELETON_SOLVER;
  WT_TYPE * grading = nullptr;
  CachedWGrevlex_Ordering * mord = nullptr;
  const CachedWGrevlex_Ordering * ford = nullptr; // final ordering
  if (not meaningful_arguments(
        argc, argv, homog, modulus, numvars, method, strategy, heuristic, solver
      )) {
    give_help();
  } else {
    int true_numvars = (homog) ? numvars + 1 : numvars;
    grading = new WT_TYPE [true_numvars];
    for (NVAR_TYPE i = 0; i < true_numvars; ++i) grading[i] = 1;
    mord = new CachedWGrevlex_Ordering(true_numvars, grading);
    Prime_Field FF(modulus);
    // set up the basis
    list<Abstract_Polynomial *> F = cyclic_n(numvars, FF, homog, mord);
    // message
    cout << "Computing a Groebner basis for:\n";
    for (Abstract_Polynomial * f : F)
      cout << '\t' << *f << endl;
    cout << "Using ";
    switch(solver) {
      case GLPK_SOLVER: cout << "GLPK\n"; break;
      case SKELETON_SOLVER: cout << "Skeleton\n"; break;
      case PPL_SOLVER: cout << "PPL\n"; break;
      default: cout << "Unknown solver\n"; break;
    }
    // compute basis
    list<Constant_Polynomial *> G = buchberger_dynamic(
        F, method, strategy, grading, (Dynamic_Heuristic )heuristic, solver
    );
    // display basis
    cout << G.size() << " polynomials in basis:\n";
    /*for (list<Constant_Polynomial *>::const_iterator g = G.begin(); g != G.end(); ++g)
      cout << '\t' << *(*g) << endl;*/
    Polynomial_Ring * R = & (G.front()->base_ring());
    ford = static_cast<const CachedWGrevlex_Ordering *>(
        G.front()->monomial_ordering()
    );
    cout << G.size() << " leading monomials:\n";
    for (Constant_Polynomial * g : G) {
      cout << g->leading_monomial() << ", ";
      //cout << *g << endl;
      delete g;
    }
    cout << endl;
    for (Abstract_Polynomial * f : F) delete f;
    delete R;
  }
  if (ford != nullptr) {
    delete [] ford->order_weights();
    delete ford;
  }
  if (mord != generic_grevlex_ptr) delete mord;
  if (goda != nullptr) delete goda;
  if (doda != nullptr) delete doda;
  if (moda != nullptr) delete moda;
  if (monoda != nullptr) delete monoda;
  if (monododa != nullptr) delete monododa;
  if (woda != nullptr) delete woda;
  if (grading != nullptr) delete [] grading;
  cout << "bye\n";
}

enum order_flags { GENERIC_GREVLEX = 0, GREVLEX, LEX, WGREVLEX };

bool meaningful_arguments(
    int argc, char *argv[],
    bool & homogeneous, int & modulus, int & numvars,
    SPolyCreationFlags & method, StrategyFlags & strategy,
    Dynamic_Heuristic & heuristic, DynamicSolver & solver
) {
  modulus = 43;
  method = SPolyCreationFlags::GEOBUCKETS;
  homogeneous = false;
  WT_TYPE * weights = nullptr;
  unsigned int order_flag = 0;
  bool good_args = (argc > 1);
  if (good_args) {
    for (int i = 1; good_args and i < argc; ++i) {
      if (!strcmp(argv[i],"hom") or !strcmp(argv[i],"homog")
          or !strcmp(argv[i],"homogeneous"))
        homogeneous = true;
      else {
        int j = 0;
        for (/* */; argv[i][j] != '=' and argv[i][j] != '\0'; ++j) { /* */ }
        if (argv[i][j] != '=') {
          good_args = false;
          cout << "Options must have form <option>=<value>.\n";
        }
        else {
          argv[i][j] = '\0';
          if (!strcmp(argv[i],"n") or !strcmp(argv[i],"num")
              or !strcmp(argv[i],"numvars")) {
            numvars = atoi(&(argv[i][j+1]));
            if (numvars < 3) {
              good_args = false;
              cout << "Invalid number of variables: must be at least 3.\n";
            }
          }
          else if (!strcmp(argv[i],"m") or !strcmp(argv[i],"mod")
                   or !strcmp(argv[i],"modulus"))
          {
            modulus = atoi(&(argv[i][j+1]));
            if (modulus < 2) {
              good_args = false;
              cout << "Invalid modulus; must be at least 2.\n";
            }
          }
          else if (!strcmp(argv[i],"r") or !strcmp(argv[i],"repr")
                   or !strcmp(argv[i], "representation"))
          {
            method = (SPolyCreationFlags )atoi(&(argv[i][j+1]));
            if (
                method <= SPolyCreationFlags::MIN_SPCREATE_FLAG or
                method >= SPolyCreationFlags::MAX_SPCREATE_FLAG
            ) {
              good_args = false;
              cout << "Invalid method; must be at least 1 and at most 3.\n";
            }
          }
          else if (!strcmp(argv[i], "heur") or !strcmp(argv[i],"heuristic")) {
            heuristic = (Dynamic_Heuristic)atoi(&(argv[i][j+1]));
            if (
                heuristic <= Dynamic_Heuristic::MIN_HEURISTIC or
                heuristic >= Dynamic_Heuristic::MAX_HEURISTIC
            ) {
              good_args = false;\
              cout << "Invalid method; must be at least 1 and at most 10.\n";
            }
          }
          else if (!strcmp(argv[i],"strat") or !strcmp(argv[i],"strategy")) {
            char * request = &(argv[i][j+1]);
            if (!strcmp(request, "normal") or !strcmp(request, "norm"))
              strategy = StrategyFlags::NORMAL_STRATEGY;
            else if (!strcmp(request, "sugar") or !strcmp(request, "sug"))
              strategy = StrategyFlags::SUGAR_STRATEGY;
            else if (!strcmp(request, "wsugar") or !strcmp(request, "wsug")) {
              strategy = StrategyFlags::WSUGAR_STRATEGY;
              unsigned n = (homogeneous) ? numvars + 1 : numvars;
            }
            else {
              good_args = false;
              cout << "Strategy must be 'normal' or 'sugar' or 'wsugar'.\n";
            }
          }
          else if (!strcmp(argv[i],"solver")) {
            char * request = &(argv[i][j+1]);
            if (!strcmp(request, "skel") or !strcmp(request, "skeleton"))
              solver = SKELETON_SOLVER;
            else if (!strcmp(request, "glpk") or !strcmp(request, "GLPK"))
              solver = GLPK_SOLVER;
            else if (!strcmp(request, "ppl") or !strcmp(request, "PPL"))
              solver = PPL_SOLVER;
            else {
              good_args = false;
              cout << "Strategy must be 'skeleton' or 'glpk'.\n";
            }
          }
          else {
            cout << "Unrecognized argument.\n"; good_args = false;
          }
        }
      }
    }
  }
  return good_args;
}

void give_help() {
  cout << "Call with options n=<num> m=<mod> r=<repr> [hom] heur=<heur>\n";
  cout << "You *must* specify <num> vars, an integer greater than 2.\n";
  cout << "You can add optional <mod>ulus (please make it prime).\n";
  cout << "The option <hom>ogenize will give you a homogenized ideal.\n";
  cout << "You can also select the <repr>esentation of s-polynomials:\n";
  cout << "\t1) linked lists,\n";
  cout << "\t2) geobuckets, or\n";
  cout << "\t3) double-buffered polynomials.\n";
  cout << "The option <heur>istic allows you to select from the following:\n";
  cout << "\t1) Hilbert heuristic w/ties broken by lex ordering,\n";
  cout << "\t2) Hilbert heuristic w/ties broken by total degree,\n";
  cout << "\t3) total degree w/ties broken by Hilbert heuristic,\n";
  cout << "\t4) graded Hilbert heuristic w/ties broken by lex ordering,\n";
  cout << "\t5) total degree w/ties broken by graded Hilbert heuristic,\n";
  cout << "\t6) aiming for polynomials of smoothest degrees (nearly homogeneous),\n";
  cout << "\t7) aiming for a polynomial with the largest maximal component,\n";
  cout << "\t8) aiming to minimize the number of new critical pairs at min degree,\n";
  cout << "\t9) aiming to minimize the number of new critical pairs at min graded degree,\n";
  cout << "\t10) pseudosignature Hilbert heuristic.\n";
  cout << "So 'test_cyclicn n=6 m=43 r=2' would compute the Groebner basis\n";
  cout << "of the Cyclic-n ideal in 6 variables, modulo 43,";
  cout << "using geobuckets.\n";
  cout << "You can also specify the strategy ('normal', 'sugar', 'wsugar')\n";
  cout << "('wsugar' requires a list of <num> integers, where <num> is as above.";
}
