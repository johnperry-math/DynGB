#include <list>
using std::list;
#include <iostream>
using std::cout; using std::endl;

#include "monomial.hpp"
#include "indeterminate.hpp"
#include "polynomial_ring.hpp"
#include "fields.hpp"
#include "particular_orderings.hpp"

int main() {

  Prime_Field F { 43 };
  auto R = Polynomial_Ring(9, F);
  Indeterminate x0(R, 0), x1(R, 1), x2(R, 2), x3(R, 3),
      x4(R, 4), x5(R, 5), x6(R, 6), x7(R, 9), x8(R, 8);

  list<Monomial> B = { x0*x1*x6, (x2^2)*x6, x0*x2*x3 };
  cout << "starting with ";
  for (auto & t : B) { cout << t << ", "; }
  cout << endl;

  WT_TYPE weights[] = { 9, 7, 8, 11, 11, 1, 12, 13, 1 };
  WGrevlex ord( 9, weights );
  for (auto & t : B) {
    t.set_monomial_ordering(&ord);
    ord.set_data(t);
  }
  B.sort();
  cout << "sorted as ";
  for (auto & t : B) { cout << t << ", "; }
  cout << endl;

  Monomial t (9);
  t = t * x0;
  Monomial u = x0*x1*x6, v = x2*x3;
  t.set_monomial_ordering(&ord);
  ord.set_data(t);
  u.set_monomial_ordering(&ord);
  ord.set_data(u);
  v.set_monomial_ordering(&ord);
  ord.set_data(v);
  cout << ord.first_larger_than_multiple(u, v, t) << endl;

}