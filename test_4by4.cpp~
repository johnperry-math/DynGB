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
#include <cstdlib>
#include <cstring>
#include <iostream>

using std::set;
using std::cout; using std::endl;

#include "system_constants.hpp"

#include "fields.hpp"
#include "monomial.hpp"
#include "polynomial.hpp"

#include "dynamic_engine.hpp"

#include "algorithm_buchberger_basic.hpp"
#include "algorithm_buchberger_dynamic.hpp"

int main(int argc, char *argv[]) {
  if (argc < 3 or (strcmp(argv[2],"stat") and strcmp(argv[2],"dyn"))) {
    cout << "Need to know method (usually 2), then if dynamic (stat or dyn).\n";
    return 1;
  }
  // obtain method -- don't screw it up b/c we don't check it
  SPolyCreationFlags method = (SPolyCreationFlags )atoi(argv[1]);
  bool static_algorithm = true;
  DynamicSolver solver = GLPK_SOLVER;
  Dynamic_Heuristic heuristic;
  if (!strcmp(argv[2],"dyn")) {
    unsigned heur_opt_location = 4;
    if (
        argc < 4 or
        (strcmp(argv[3],"skel") and strcmp(argv[3],"glpk")
          and strcmp(argv[3],"simpl") and strcmp(argv[3],"ppl"))
    ) {
      cout << "For a dynamic algorithm, we also need the solver:\n";
      cout << "\tskeleton (skel) or simplex (glpk or simpl).\n";
      return 2;
    }
    static_algorithm = false;
    if (!strcmp(argv[3],"skel")) solver = SKELETON_SOLVER;
    else if (!strcmp(argv[3],"ppl")) solver = PPL_SOLVER;
    // obtain heuristic
    if (argc < heur_opt_location + 1 or
        argv[heur_opt_location][0] != 'h' or argv[heur_opt_location][1] != 'e' or
        argv[heur_opt_location][2] != 'u' or argv[heur_opt_location][3] != 'r'
    ) {
      cout << "Need to know heuristic. Use format heur=...\n";
      return 1;
    }
    heuristic = (Dynamic_Heuristic )atoi(&(argv[heur_opt_location][5]));
  }
 // set up the field
  Prime_Field F17 = Prime_Field(32003);
  Polynomial_Ring R(32, F17);
  Prime_Field_Element a = F17.unity();
  // set up our polynomials
  // first poly
  EXP_TYPE e11 [] {0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0}; Monomial t11(32, e11);
  EXP_TYPE e12 [] {0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0}; Monomial t12(32, e12);
  EXP_TYPE e13 [] {0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0}; Monomial t13(32, e13);
  EXP_TYPE e14 [] {0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0}; Monomial t14(32, e14);
  EXP_TYPE e15 [] {0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0}; Monomial t15(32, e15);
  EXP_TYPE e16 [] {0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0}; Monomial t16(32, e16);
  Monomial M1 [] { t11, t12, t13, t14, t15, t16 };
  Prime_Field_Element C1 [] { -a, a, a, a, -a, -a };
  Constant_Polynomial f1(6, R, M1, C1);
  f1.sort_by_order();
  // second poly
  EXP_TYPE e21 [] {0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}; Monomial t21(32, e21);
  EXP_TYPE e22 [] {1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0}; Monomial t22(32, e22);
  EXP_TYPE e23 [] {0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0}; Monomial t23(32, e23);
  EXP_TYPE e24 [] {0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0}; Monomial t24(32, e24);
  EXP_TYPE e25 [] {0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0}; Monomial t25(32, e25);
  EXP_TYPE e26 [] {0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0}; Monomial t26(32, e26);
  EXP_TYPE e27 [] {0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0}; Monomial t27(32, e27);
  EXP_TYPE e28 [] {0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0}; Monomial t28(32, e28);
  Monomial M2 [] { t21, t22, t23, t24, t25, t26, t27, t28 };
  Prime_Field_Element C2 [] { -a, a, -a, a, a, a, -a, -a };
  Constant_Polynomial f2(8, R, M2, C2);
  f2.sort_by_order();
  // third poly
  EXP_TYPE e31 [] {1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0}; Monomial t31(32, e31);
  EXP_TYPE e32 [] {0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0}; Monomial t32(32, e32);
  EXP_TYPE e33 [] {0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0}; Monomial t33(32, e33);
  EXP_TYPE e34 [] {0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}; Monomial t34(32, e34);
  EXP_TYPE e35 [] {0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0}; Monomial t35(32, e35);
  EXP_TYPE e36 [] {0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0}; Monomial t36(32, e36);
  EXP_TYPE e37 [] {0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0}; Monomial t37(32, e37);
  EXP_TYPE e38 [] {0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0}; Monomial t38(32, e38);
  Monomial M3 [] { t31, t32, t33, t34, t35, t36, t37, t38 };
  Prime_Field_Element C3 [] { a, -a, -a, -a, a, a, a, -a };
  Constant_Polynomial f3(8, R, M3, C3);
  f3.sort_by_order();
  // fourth poly
  EXP_TYPE e41 [] {0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0}; Monomial t41(32, e41);
  EXP_TYPE e42 [] {1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0}; Monomial t42(32, e42);
  EXP_TYPE e43 [] {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0}; Monomial t43(32, e43);
  EXP_TYPE e44 [] {0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0}; Monomial t44(32, e44);
  EXP_TYPE e45 [] {0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0}; Monomial t45(32, e45);
  EXP_TYPE e46 [] {0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}; Monomial t46(32, e46);
  EXP_TYPE e47 [] {0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1}; Monomial t47(32, e47);
  EXP_TYPE e48 [] {0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0}; Monomial t48(32, e48);
  Monomial M4 [] { t41, t42, t43, t44, t45, t46, t47, t48 };
  Prime_Field_Element C4 [] { -a, a, -a, a, a, -a, a, -a };
  Constant_Polynomial f4(8, R, M4, C4);
  f4.sort_by_order();
  // fifth poly
  EXP_TYPE e51 [] {1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0}; Monomial t51(32, e51);
  EXP_TYPE e52 [] {0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0}; Monomial t52(32, e52);
  EXP_TYPE e53 [] {0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}; Monomial t53(32, e53);
  EXP_TYPE e54 [] {0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0}; Monomial t54(32, e54);
  EXP_TYPE e55 [] {0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0}; Monomial t55(32, e55);
  EXP_TYPE e56 [] {0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0}; Monomial t56(32, e56);
  EXP_TYPE e57 [] {0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0}; Monomial t57(32, e57);
  EXP_TYPE e58 [] {0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0}; Monomial t58(32, e58);
  Monomial M5 [] { t51, t52, t53, t54, t55, t56, t57, t58 };
  Prime_Field_Element C5 [] { -a, -a, a, -a, a, a, a, -a };
  Constant_Polynomial f5(8, R, M5, C5);
  f5.sort_by_order();
  // sixth poly
  EXP_TYPE e61 [] {0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0}; Monomial t61(32, e61);
  EXP_TYPE e62 [] {0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0}; Monomial t62(32, e62);
  EXP_TYPE e63 [] {0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0}; Monomial t63(32, e63);
  EXP_TYPE e64 [] {0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0}; Monomial t64(32, e64);
  EXP_TYPE e65 [] {0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0}; Monomial t65(32, e65);
  EXP_TYPE e66 [] {0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0}; Monomial t66(32, e66);
  Monomial M6 [] { t61, t62, t63, t64, t65, t66 };
  Prime_Field_Element C6 [] { -a, -a, a, a, a, -a };
  Constant_Polynomial f6(6, R, M6, C6);
  f6.sort_by_order();
  // seventh poly
  EXP_TYPE e71 [] {0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0}; Monomial t71(32, e71);
  EXP_TYPE e72 [] {0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0}; Monomial t72(32, e72);
  EXP_TYPE e73 [] {0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0}; Monomial t73(32, e73);
  EXP_TYPE e74 [] {0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0}; Monomial t74(32, e74);
  EXP_TYPE e75 [] {0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0}; Monomial t75(32, e75);
  EXP_TYPE e76 [] {0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0}; Monomial t76(32, e76);
  EXP_TYPE e77 [] {0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0}; Monomial t77(32, e77);
  EXP_TYPE e78 [] {0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0}; Monomial t78(32, e78);
  Monomial M7 [] { t71, t72, t73, t74, t75, t76, t77, t78 };
  Prime_Field_Element C7 [] { -a, -a, -a, a, a, -a, a, a };
  Constant_Polynomial f7(8, R, M7, C7);
  f7.sort_by_order();
  // eighth poly
  EXP_TYPE e81 [] {0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0}; Monomial t81(32, e81);
  EXP_TYPE e82 [] {0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0}; Monomial t82(32, e82);
  EXP_TYPE e83 [] {0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0}; Monomial t83(32, e83);
  EXP_TYPE e84 [] {0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0}; Monomial t84(32, e84);
  EXP_TYPE e85 [] {0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0}; Monomial t85(32, e85);
  EXP_TYPE e86 [] {0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0}; Monomial t86(32, e86);
  EXP_TYPE e87 [] {0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0}; Monomial t87(32, e87);
  EXP_TYPE e88 [] {0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0}; Monomial t88(32, e88);
  Monomial M8 [] { t81, t82, t83, t84, t85, t86, t87, t88 };
  Prime_Field_Element C8 [] { -a, -a, -a, a, a, -a, a, a };
  Constant_Polynomial f8(8, R, M8, C8);
  f8.sort_by_order();
  // message
  cout << "Computing a Groebner basis for\n\t" << f1 << ",\n\t" << f2
       << "\n\t" << f3 << "\n\t" << f4
       << "\n\t" << f5 << "\n\t" << f6
       << "\n\t" << f7 << "\n\t" << f8
       << endl;
  // compute basis
  list<Abstract_Polynomial *> F;
  F.push_back(&f1); F.push_back(&f2); F.push_back(&f3); F.push_back(&f4);
  F.push_back(&f5); F.push_back(&f6); F.push_back(&f7); F.push_back(&f8);
  list<Constant_Polynomial *> G;
  if (static_algorithm) G = buchberger(F, method, StrategyFlags::SUGAR_STRATEGY);
  else G = buchberger_dynamic(
      F, method, StrategyFlags::SUGAR_STRATEGY, nullptr, heuristic, solver, true
  );
  auto ordering = static_cast<const CachedWGrevlex_Ordering *>(
      G.front()->monomial_ordering()
  );
  cout << "Basis:\n";
  while (not G.empty()) {
    auto g = G.front();
    G.pop_front();
    cout << '\t' << g->leading_monomial() << endl;
    //delete g;
  }
  if (not static_algorithm) {
    delete [] ordering->order_weights();
    delete ordering;
  }
  //check_correctness(G, SUGAR_STRATEGY);
  cout << "bye\n";
}
