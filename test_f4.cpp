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
#include "monomial_ordering.hpp"
#include "particular_orderings.hpp"
#include "polynomial_linked_list.hpp"
#include "f4_reduction.hpp"

extern Monomial_Ordering * generic_grevlex_ptr;
extern Grading_Order_Data_Allocator<WT_TYPE> * goda;
extern Grading_Order_Data_Allocator<EXP_TYPE> * moda;
extern Grading_Order_Data_Allocator<Monomial> * monoda;
extern Grading_Order_Data_Allocator<Monomial_Node> * monododa;

// Forward declarations
bool meaningful_arguments(int, char **, bool &, int &, int &);

void give_help();

int main(int argc, char *argv[]) {
  bool homog;
  int modulus, numvars;
  if (not meaningful_arguments(argc, argv, homog, modulus, numvars)) {
    give_help();
  } else {
    int true_numvars = (homog) ? numvars + 1 : numvars;
    Prime_Field FF(modulus);
    // set up the basis
    list<Abstract_Polynomial *> F = cyclic_n(numvars, FF, homog);
    // message
    cout << "Computing a Groebner basis for:\n";
    for (Abstract_Polynomial * f : F)
      cout << '\t' << *f << endl;
    // compute basis
    list<Constant_Polynomial *> G = f4_control(F);
    // display basis
    cout << G.size() << " polynomials in basis:\n";
    /*for (list<Constant_Polynomial *>::const_iterator g = G.begin(); g != G.end(); ++g)
      cout << '\t' << *(*g) << endl;*/
    Polynomial_Ring * R = & (G.front()->base_ring());
    auto mord = static_cast<const CachedWGrevlex_Ordering *>(
        G.front()->monomial_ordering()
    );
    cout << G.size() << " leading monomials:\n";
    for (Constant_Polynomial * g : G) {
      cout << g->leading_monomial() << ", ";
      delete g->strategy();
      delete g;
    }
    cout << endl;
    for (Abstract_Polynomial * f : F) {
      delete f->strategy();
      delete f;
    }
    delete R;
    if (mord != generic_grevlex_ptr) delete mord;
  }
  if (goda != nullptr) delete goda;
  if (moda != nullptr) delete moda;
  if (monoda != nullptr) delete monoda;
  if (monododa != nullptr) delete monododa;
  cout << "bye\n";
}

enum order_flags { GENERIC_GREVLEX = 0, GREVLEX, LEX, WGREVLEX };

bool meaningful_arguments(
    int argc, char *argv[],
    bool & homogeneous, int & modulus, int & numvars
) {
  modulus = 43;
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
  cout << "Call with options n=<num> m=<mod> [hom]\n";
  cout << "You *must* specify <num> vars, an integer greater than 2.\n";
  cout << "You can add optional <mod>ulus (please make it prime).\n";
  cout << "The option <hom>ogenize will give you a homogenized ideal.\n";
  cout << "So 'test_f4 n=6 m=43' would compute the Groebner basis\n";
  cout << "of the Cyclic-n ideal in 6 variables, modulo 43.";
}
