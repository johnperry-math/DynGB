#ifndef __POLYNOMIAL_HASHED_CPP_
#define __POLYNOMIAL_HASHED_CPP_

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

#include <algorithm>
using std::sort;

#include "polynomial_hashed.hpp"

bool Polynomial_Hashed::sort_indirect(
    pair< unsigned, Prime_Field_Element > & is_first,
    pair< unsigned, Prime_Field_Element > & is_second
) {
  Monomial & t = *M[ is_first.first ], u = *M[ is_second.first ];
  t.set_monomial_ordering(mord);
  u.set_monomial_ordering(mord);
  return mord->first_smaller(u, t);
}

Polynomial_Hashed::Polynomial_Hashed(
    Polynomial_Ring & R,
    vector< Monomial * > & monomials,
    F4_Hash & monomial_hash,
    const vector< pair< unsigned, COEF_TYPE > > & from_row,
    Monomial_Ordering * mord
)
  : Abstract_Polynomial(R, mord), M(monomials), hash(monomial_hash), mord(mord)
{
  terms.reserve(from_row.size());
  auto & F = R.ground_field();
  for (const auto & p : from_row) {
    const auto t = monomials[ p.first ];
    if (p.second == 0)
      cout << "error: trying to save a monomial w/coefficient 0: " << *t << "\n";
    else {
      if (not hash.contains(t)) {
        monomials.emplace_back(new Monomial(*t));
        auto where = monomials.size() - 1;
        hash.add_monomial(monomials[where], where);
      }
      terms.emplace_back(
          hash[*t],
          Prime_Field_Element( p.second , &F )
      );
    }
  }
}

Polynomial_Hashed::Polynomial_Hashed(
    Polynomial_Ring & R,
    vector< Monomial * > & monomials,
    F4_Hash & monomial_hash,
    const vector< pair< unsigned, COEF_TYPE > > & from_row,
    const vector< size_t > & row_monomials,
    Monomial_Ordering * mord
)
  : Abstract_Polynomial(R, mord), M(monomials), hash(monomial_hash), mord(mord)
{
  terms.reserve(from_row.size());
  auto & F = R.ground_field();
//auto d = row_monomials[from_row[0].first]->total_degree();
  for (const auto & p : from_row) {
    const auto t = monomials[ row_monomials[p.first] ];
//if (t->total_degree() != d) cout << "ERROR 1 adding " << *row_monomials[from_row[0].first] << " 's monomial " << *t << endl;
    if (p.second == 0)
      cout << "error: trying to save a monomial w/coefficient 0: " << *t << "\n";
    else {
      if (not hash.contains(t)) {
        monomials.emplace_back(new Monomial(*t));
//cout << "emplacing " << *t << " in " << &hash << " at " << monomials.size() - 1 << endl;
        auto where = monomials.size() - 1;
        hash.add_monomial(monomials[where], where);
      }
//else {
//  cout << "found " << *t << " in " << &hash << endl;
//}
/*if (not monomials[hash[*t]]->is_like(*t)) {
  cout << "hashing error: " << *t << " ( " << t << " ) is not " << *monomials[hash[*t]] << " ( " << monomials[hash[*t]] << " )\n";
  cout << "hash is " << hash[*t] << " ; length is " << monomials.size() << endl;
  cout << "hash reports that it contains " << *t << " ? " << hash.contains(t) << endl;
  hash.vomit(*t);
  cout << "hash reports that it contains " << *monomials[hash[*t]] << " ? " << hash.contains(monomials[hash[*t]])<< "\n";
  hash.vomit(*monomials[hash[*t]]);
}*/
      terms.emplace_back(
          hash[*t],
          Prime_Field_Element( p.second , &F )
      );
/*if (monomials[terms[terms.size() - 1].first]->total_degree() != d) {
   cout << "ERROR 2 adding " << *row_monomials[from_row[0].first] << " 's monomial " << *t << " as " << *monomials[terms[terms.size() - 1].first] << endl;
}*/
    }
  }
}

Polynomial_Hashed::Polynomial_Hashed(
    Abstract_Polynomial & p,
    vector< Monomial * > & M,
    F4_Hash & table,
    Monomial_Ordering * mord
)
  : Abstract_Polynomial(p.base_ring(), p.monomial_ordering()),
    M(M), hash(table), mord(mord)
{
  terms.reserve(p.length());
  Polynomial_Iterator * pi = p.new_iterator();
  while (not pi->fellOff()) {
    size_t location;
    auto & u = pi->currMonomial();
    if (table.contains(&u)) {
      location = hash[u];
    } else {
      Monomial * v = new Monomial(u);
      M.emplace_back(v);
      location = M.size() - 1;
      table.add_monomial(v, location);
    }
    terms.emplace_back( location, pi->currCoeff() );
    pi->moveRight();
  }
  delete pi;
}

Polynomial_Hashed::~Polynomial_Hashed() {}

void Polynomial_Hashed::set_monomial_ordering(
    const Monomial_Ordering * new_ord, bool resort
) {
  if (mord != new_ord) {
    mord = const_cast<Monomial_Ordering *>(new_ord);
    sort_by_order();
  }
}

void Polynomial_Hashed::sort_by_order() {
  sort(
      terms.begin(), terms.end(),
      [this](auto l, auto r) { return sort_indirect(l, r); }
  );
}

Monomial & Polynomial_Hashed::leading_monomial() const {
  auto & result = *M[ terms.begin()->first ];
  result.set_monomial_ordering(mord);
  return result;
}

Prime_Field_Element Polynomial_Hashed::leading_coefficient() const {
  return terms.begin()->second;
}

unsigned Polynomial_Hashed::length() const {
  return terms.size();
}

bool Polynomial_Hashed::is_zero() const {
  return terms.size() == 0;
}

Polynomial_Hashed * Polynomial_Hashed::zero_polynomial() const {
  vector< pair< unsigned, COEF_TYPE> > new_terms;
  return new Polynomial_Hashed(R, M, hash, new_terms, mord);
}

Polynomial_Hashed * Polynomial_Hashed::monomial_multiple(
    const Monomial & t
) const {
  vector< pair< unsigned, COEF_TYPE> > new_terms( terms.size() );
  for (auto p : terms) {
    if (hash.contains_product(t, *M[ p.first ])) {
      new_terms.emplace_back(
          hash.lookup_product(t, *M[ p.first ]), p.second.value()
      );
    } else {
      auto u = new Monomial(t, *M[ p.first ]);
      M.emplace_back(u);
      hash.add_monomial(u, M.size() - 1);
      new_terms.emplace_back(M.size() - 1, p.second.value());
    }
  }
  return new Polynomial_Hashed(R, M, hash, new_terms, mord);
}

Polynomial_Hashed * Polynomial_Hashed::scalar_multiple(
    const Prime_Field_Element & a
) const {
  vector< pair< unsigned, COEF_TYPE > > new_terms( terms.size() );
  auto mod = R.ground_field().modulus();
  auto a0 = a.value();
  for (auto p : terms) {
    new_terms.emplace_back(
        p.first, p.second.value() * a0 % mod
    );
  }
  return new Polynomial_Hashed(R, M, hash, new_terms, mord);
}

Hashed_Polynomial_Iterator * Polynomial_Hashed::new_iterator() const {
  return new Hashed_Polynomial_Iterator(this);
}

Hashed_Polynomial_Iterator * Polynomial_Hashed::begin() const {
  return new Hashed_Polynomial_Iterator(this);
}                                                 

Hashed_Polynomial_Iterator * Polynomial_Hashed::end() const {
  return new Hashed_Polynomial_Iterator(this, true);
}

Hashed_Polynomial_Iterator::Hashed_Polynomial_Iterator(
    const Polynomial_Hashed *p, bool at_end
)
  : p(p), started_at_beginning(not at_end)
{
  p_base = p;
  if (not at_end) i = 0;
  else i = p->terms.size() - 1;
}

Hashed_Polynomial_Iterator::~Hashed_Polynomial_Iterator() {}

void Hashed_Polynomial_Iterator::restart_iteration() {
  if (started_at_beginning) i = 0;
  else i = p->terms.size() - 1;
}

void Hashed_Polynomial_Iterator::moveRight() { ++i; }

void Hashed_Polynomial_Iterator::moveLeft() { --i; }

bool Hashed_Polynomial_Iterator::canMoveRight() const {
  return (i < p->terms.size() - 1);
}

bool Hashed_Polynomial_Iterator::canMoveLeft() const {
  return (i > 0);
}

bool Hashed_Polynomial_Iterator::fellOff() const {
  return not (i >= 0 and i < p->terms.size());
}

const Monomial & Hashed_Polynomial_Iterator::currMonomial() const {
  return *p->M[ p->terms[i].first ] ;
}

const Prime_Field_Element & Hashed_Polynomial_Iterator::currCoeff() const {
  return p->terms[i].second ;
}

#endif
