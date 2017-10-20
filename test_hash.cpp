#include <iostream>
using std::cout; using std::endl;

#include <list>
using std::list;

#include "system_constants.hpp"

#include "monomial.hpp"
#include "f4_hash.hpp"

int main() {

  const NVAR_TYPE n=5;

  F4_Hash table(n);
  list<pair<Monomial *, Monomial *> > factors;
  list<Monomial *> products;

  srand(time(NULL));
  for (unsigned i = 0; i < 500; ++i) {
    auto t = new Monomial({
      (unsigned int) rand() % 2, (unsigned int) rand() % 2,
      (unsigned int) rand() % 2, (unsigned int) rand() % 2,
      (unsigned int) rand() % 2
    });
    auto u = new Monomial({
      (unsigned int) rand() % 2, (unsigned int) rand() % 2,
      (unsigned int) rand() % 2, (unsigned int) rand() % 2,
      (unsigned int) rand() % 2
    });
    factors.emplace_back(t, u);
    auto v = new Monomial((*t)*(*u));
    bool found = false;
    auto iter = products.begin();
    for (/* already initialized */; not found and iter != products.end(); /* */) {
      if (v->is_like(**iter))
        found = true;
      else if (*v > **iter) {
        found = true;
        products.insert(iter, v);
      }
      else ++iter;
    }
    if (not found)
      products.push_back(v);
  }

  size_t i = 0;
  for (auto v : products)
    table.add_monomial(v, i++);

  Monomial * M[i] { nullptr };
  i = 0;
  for (auto v : products) {
    cout << i << ": " << *v << endl;
    M[i++] = v;
  }

  cout << table;
  cout << i << " elements\n";

  for (auto P : factors) {
    size_t i = table.lookup_product(*P.first, *P.second);
    cout << *P.first << ", " << *P.second << " (" << i << ") = " << *M[i] << " ? " << (*M[i] == (*P.first * *P.second)) << endl;
  }

}