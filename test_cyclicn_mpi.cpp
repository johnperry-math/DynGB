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

#include <set>
#include <cstring>
#include <iostream>

using std::set;
using std::cout; using std::endl;

#include "system_constants.hpp"

#include "cyclic_n.hpp"
#include "polynomial.hpp"
#include "strategies.hpp"
#include "monomial_ordering.hpp"
#include "particular_orderings.hpp"
#include "algorithm_buchberger_basic.hpp"
#include "algorithm_buchberger_explorer.hpp"
#include "f4_reduction.hpp"

#include "mpi.h"

extern Monomial_Ordering * generic_grevlex_ptr;

// Forward declarations
bool meaningful_arguments(
    int, char **, bool &, int &, int &, SPolyCreationFlags &, bool &,
    Monomial_Ordering **, StrategyFlags &, WT_TYPE **
  );
void give_help();

int main(int argc, char *argv[]) {
  bool homog;
  bool f4 = false;
  int modulus, numvars;
  SPolyCreationFlags method;
  StrategyFlags strategy = NORMAL_STRATEGY;
  WT_TYPE *grading;
  Monomial_Ordering * mord = generic_grevlex_ptr;
  // Multiprocessing
  MPI_Init(nullptr, nullptr);
  double start_time = MPI_Wtime();
  int comm_size, comm_id;
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &comm_id);
  if (not meaningful_arguments(
        argc, argv, homog, modulus, numvars, method, f4,
        &mord, strategy, &grading
      )) {
    if (comm_id == 0) // only need help from the first process
      give_help();
  } else {
    Prime_Field FF = Prime_Field(modulus);
    // set up the basis
    list<Abstract_Polynomial *> B = cyclic_n(numvars, FF, homog, mord);
    vector<Abstract_Polynomial *> F;
    for (auto b : B)
      F.push_back(b);
    // message
    if (comm_id == 0) {
      cout << "Computing a Groebner basis for:\n";
      for (Abstract_Polynomial * f : F)
        cout << '\t' << *f << endl;
    }
    // compute basis
    //set<Constant_Polynomial *, smaller_lm> G = buchberger(F, method);
    list<Constant_Polynomial *> G;
    G = buchberger_explorer(F, method, strategy, grading, comm_id, comm_size);
    cout << comm_id << " returned from GB\n";
    // display basis
    if (comm_id == 0) {
      cout << G.size() << " polynomials in basis:\n";
      /*for (list<Constant_Polynomial *>::const_iterator g = G.begin(); g != G.end(); ++g)
        cout << '\t' << *(*g) << endl;*/
      Polynomial_Ring * R = & (G.front()->base_ring());
      cout << G.size() << " leading monomials:\n";
      for (Constant_Polynomial * g : G) {
        cout << g->leading_monomial() << ' ' << endl;
        delete g;
      }
      for (Abstract_Polynomial * f : F) delete f;
      delete R;
    }
  }
  if (mord != generic_grevlex_ptr)
    delete mord;
  cout << "bye from " << comm_id << "\n";
  MPI_Barrier(MPI_COMM_WORLD);
  if (comm_id == 0)
    cout << "elapsed time: " << MPI_Wtime() - start_time << endl;
  MPI_Finalize();
}

enum order_flags { GENERIC_GREVLEX = 0, GREVLEX, LEX, WGREVLEX };

bool meaningful_arguments(int argc, char *argv[], bool & homogeneous,
                          int & modulus, int & numvars,
                          SPolyCreationFlags & method, bool &f4,
                          Monomial_Ordering ** mord, StrategyFlags & strategy,
                          WT_TYPE ** grading
                         )
{
  modulus = 43;
  method = LINKED_LST;
  homogeneous = false;
  WT_TYPE * weights = nullptr;
  unsigned int order_flag = 0;
  bool good_args = (argc > 1);
  if (good_args) {
    for (int i = 1; good_args and i < argc; ++i) {
      if (!strcmp(argv[i],"hom") or !strcmp(argv[i],"homog")
          or !strcmp(argv[i],"homogeneous"))
        homogeneous = true;
      else if (!strcmp(argv[i],"f4"))
        f4 = true;
      else {
        int j = 0;
        for (/* */; argv[i][j] != '=' and argv[i][j] != '\0'; ++j) { /* */ }
        if (argv[i][j] != '=') {
          good_args = false;
          cout << "Arguments must have form <option>=<value>.\n";
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
            if (method < 1 or method > 3) {
              good_args = false;
              cout << "Invalid method; must be at least 1 and at most 3.\n";
            }
          }
          else if (!strcmp(argv[i],"ord") or !strcmp(argv[i],"order")
                   or !strcmp(argv[i],"ordering"))
          {
            char * request = &(argv[i][j+1]);
            if (!strcmp(request, "generic"))
              order_flag = GENERIC_GREVLEX;
            else if (!strcmp(request, "grevlex"))
              order_flag = GREVLEX;
            else if (!strcmp(request, "lex"))
              order_flag = LEX;
            else if (!strcmp(request, "wgrevlex")) {
              order_flag = WGREVLEX;
              unsigned n = (homogeneous) ? numvars + 1 : numvars;
              weights = new WT_TYPE [n];
              unsigned k = ++i;
              for (/* */; k < i + n and k < argc; ++k)
                weights[k - i] = atoi(argv[k]);
              if (k - i < n)
                good_args = false;
              i = k;
            }
            else {
              good_args = false;
              cout << "Ordering must be 'generic', 'grevlex', "
                   << "'lex', or 'wgrevlex'.\n"
                   << "(Matrix orderings not yet supported via command line.)\n"
                   << "When using 'wgrevlex', follow it with a list of 'n' "
                   << "positive integers.\n";
            }
          }
          else if (!strcmp(argv[i],"strat") or !strcmp(argv[i],"strategy")) {
            char * request = &(argv[i][j+1]);
            if (!strcmp(request, "normal") or !strcmp(request, "norm"))
              strategy = NORMAL_STRATEGY;
            else if (!strcmp(request, "sugar") or !strcmp(request, "sug"))
              strategy = SUGAR_STRATEGY;
            else if (!strcmp(request, "wsugar") or !strcmp(request, "wsug")) {
              strategy = WSUGAR_STRATEGY;
              unsigned n = (homogeneous) ? numvars + 1 : numvars;
              *grading = new WT_TYPE [n];
              unsigned k = ++i;
              for (/* */; k < i + n and k < argc; ++k)
                (*grading)[k - i] = atoi(argv[k]);
              if (k - i < n)
                good_args = false;
              i = k;
            }
            else {
              good_args = false;
              cout << "Strategy must be 'normal' or 'sugar'.";
            }
          }
          else {
            cout << "Unrecognized argument.\n"; good_args = false;
          }
        }
      }
    }
  }
  if (good_args) {
    unsigned n = (homogeneous) ? numvars + 1 : numvars;
    switch (order_flag) {
    case GENERIC_GREVLEX: *mord = generic_grevlex_ptr; break;
    case GREVLEX: *mord = new Grevlex_Ordering(n); break;
    case LEX: *mord = new Lex_Ordering(n); break;
    case WGREVLEX: *mord = new CachedWGrevlex_Ordering(n, weights); break;
    default: *mord = generic_grevlex_ptr;
    }
  }
  return good_args;
}

void give_help() {
  cout << "This is the same as test_cyclicn, only using OpenMPI for multiprocessing.\n";
  cout << "Call with options n=<num> m=<mod> r=<repr>\n";
  cout << "You *must* specify <num> vars, an integer greater than 2.\n";
  cout << "You can add optional <mod>ulus (please make it prime).\n";
  cout << "The option <hom>ogenize will give you a homogenized ideal.\n";
  cout << "You can also select the <repr>esentation of s-polynomials:\n";
  cout << "\t1) linked lists,\n";
  cout << "\t2) geobuckets, or\n";
  cout << "\t3) double-buffered polynomials.\n";
  cout << "So 'test_cyclicn n=6 m=43 r=2' would compute the Groebner basis\n";
  cout << "of the Cyclic-n ideal in 6 variables, modulo 43,";
  cout << "using geobuckets.\n";
  cout << "You can also specify the strategy ('normal', 'sugar', 'wsugar')\n";
  cout << "('wsugar' requires a list of <num> integers, where <num> is as above,";
  cout << "and the term ordering ('generic', 'grevlex', 'lex', 'wgrevlex').\n";
  cout << "('wgrevlex' requires a list of <num> integers, where <num>";
  cout << " is as above.\n";
  cout << "The 'f4' option runs the F4 algorithm.\n";
}
