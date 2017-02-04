#ifndef __ALGORITHM_BUCHBERGER_EXPLORER_CPP_
#define __ALGORITHM_BUCHBERGER_EXPLORER_CPP_

#include <vector>
using std::vector;

#include "hilbert_functions.hpp"

#include "algorithm_buchberger_basic.hpp"
#include "reduction_support.hpp"

#include "algorithm_buchberger_explorer.hpp"

#include "mpi.h"

#define XPLOR_GENERATOR -1
#define PROC_UNASSIGNED -1

/**
  @ingroup GBComputation
  @brief contains information on critical pairs by their index in the basis,
    in addition to the usual information
*/
class Critical_Pair_XPlor : public Critical_Pair_Basic {
public:
  /** @name Construction */
  ///@{
  /** @brief create critical pair (f,0) where f is at index @c i */ 
  Critical_Pair_XPlor(int i, unsigned strategy, Abstract_Polynomial * f)
  : Critical_Pair_Basic(f, strategy)
  {
    pi = i; qi = XPLOR_GENERATOR; proc = PROC_UNASSIGNED;
  }
  /** @brief create critical pair (f,g) where f, g are at indices @c i, @c j */
  Critical_Pair_XPlor(
      int i, int j, unsigned strategy, vector<Abstract_Polynomial *>G
  )
  : Critical_Pair_Basic(G[i], G[j], strategy)
  {
    pi = i; qi = j; proc = PROC_UNASSIGNED;
  }
  /** @brief create critical pair (f,g) where f is at index @c i */
  Critical_Pair_XPlor(
      int i, Abstract_Polynomial *g, unsigned strategy,
      vector<Abstract_Polynomial *>G
  )
  : Critical_Pair_Basic(G[i], g, strategy)
  {
    pi = i; qi = G.size(); proc = PROC_UNASSIGNED;
  }
  ///@}
  /** @name Basic properties */
  ///@{
  /** @brief returns index of first polynomial in pair */
  int first_index() { return pi; }
  /** @brief returns index of second polynomial in pair */
  int second_index() { return qi; }
  /** @brief returns sugar of this pair; use ONLY if with sugar strategy */
  DEG_TYPE sugar() {
    return (static_cast<Pair_Sugar_Data *>(key))->pair_sugar();
  }
  ///@}
  /** @name Multiprocessing data */
  ///@{
  /** @brief record that this pair is assigned to processor @c i */
  void set_processor(int i) { proc = i; }
  /**
    @brief query whether this pair is assigned to a processor, and which
      (nonnegative value indicates assignment, and to which)
  */
  int get_processor() { return proc; }
  ///@}
protected:
  /** @brief first polynomial in the critical pair */
  int pi;
  /** @brief second polynomial in the critical pair */
  int qi;
  /** @brief processor to which this pair has been assigned */
  int proc;
};

/**
  @brief compares the Hilbert functions, as specified by Hilbert polynomial and
    series
  @return  @c true iff the second Hilbert function is measurably better than
    the first.
*/
bool second_HF_smaller(
    Dense_Univariate_Rational_Polynomial * hp1,
    Dense_Univariate_Integer_Polynomial * hn1,
    Dense_Univariate_Rational_Polynomial * hp2,
    Dense_Univariate_Integer_Polynomial * hn2
) {
  // this adapts LessByHilbert in dynamic_engine.cpp
  // in this case, we want the opposite result,
  // so the return statement at the end is negated
  bool result;
  // first check the coefficients of the Hilbert polynomial
  Dense_Univariate_Rational_Polynomial HPdiff(*hp1);
  HPdiff -= *hp2;
  if (not HPdiff.is_zero())
    result = (HPdiff.numerator(HPdiff.degree()) < 0);
  else // use Hilbert series
  {
    Dense_Univariate_Integer_Polynomial * h1 = hn1;
    Dense_Univariate_Integer_Polynomial * h2 = hn2;
    DEG_TYPE i = 0;
    for ( /* already initialized */ ;
          i <= h1->degree() and i <= h2->degree() and (*h1)[i] == (*h2)[i];
          i++)
    { /* taken care of in loop */ }
    if (i > h1->degree())
    {
      if (i > h2->degree())
        // the numerators are equal; second is not measurably better
        result = false;
      else
        result = true;
    }
    else
    {
      if (i > h2->degree()) result = false;
      else result = (*h1)[i] < (*h2)[i];
    }
  }
  //cout << "\tfirst less than second? " << result << endl;
  return not result;
}

/**
  @brief compares the Hilbert function at position @c i (&ldquo;older&rdquo;)
    with the Hilbert function at position @c j (&ldquo;newer&rdquo;).
    Returns @c true iff j&rsquo;s Hilbert function is measurably better than
    i&rsquo;s.
*/
bool newer_HF_smaller(Monomial & t, unsigned i, Monomial & u, unsigned j,
    Dense_Univariate_Rational_Polynomial ** HP,
    Dense_Univariate_Integer_Polynomial ** HFN
) {
  // this adapts LessByHilbert in dynamic_engine.cpp
  // in this case, we want the opposite result,
  // so the return statement at the end is negated
  bool result;
  NVAR_TYPE n = t.num_vars();
  // first check the coefficients of the Hilbert polynomial
  Dense_Univariate_Rational_Polynomial HPdiff(*(HP[i]));
  HPdiff -= *(HP[j]);
  if (not HPdiff.is_zero())
    result = (HPdiff.numerator(HPdiff.degree()) < 0);
  else // use Hilbert series
  {
    Dense_Univariate_Integer_Polynomial * h1 = HFN[i];
    Dense_Univariate_Integer_Polynomial * h2 = HFN[j];
    DEG_TYPE i = 0;
    for ( /* already initialized */ ;
          i <= h1->degree() and i <= h2->degree() and (*h1)[i] == (*h2)[i];
          i++)
    { /* taken care of in loop */ }
    if (i > h1->degree())
    {
      if (i > h2->degree())
      { // the numerators are equal; break tie via lex
        int i = 0;
        while (i < n and t[i] == u[i]) ++i;
        if (i == n) result = false;
        else result = (t.degree(i) > u.degree(i));
      }
      else
        result = true;
    }
    else
    {
      if (i > h2->degree()) result = false;
      else result = (*h1)[i] < (*h2)[i];
    }
  }
  //cout << "\tfirst less than second? " << result << endl;
  return not result;
}

/**
  @brief Implementation of Gebauer-Moeller algorithm, with XPLOR critical pairs.
  Based on description in Becker and Weispfenning (1993).
  @ingroup GBComputation
  @param P list of critical pairs that are not assigned
  @param Pass list of critical pairs that are assigned
  @param Pcancel array of critical pairs that are discovered to be redundant
  @param G current basis
  @param r polynomial to add to basis (and to generate new pairs)
  @param strategy how to sort pairs
*/
void gm_update_explorer(
    list<Critical_Pair_XPlor *> & P,
    list<Critical_Pair_XPlor *> & Pass,
    list<Critical_Pair_XPlor *> * Pcancel,
    vector<Abstract_Polynomial *> & G,
    Abstract_Polynomial * r,
    unsigned strategy
) {
  //cout << "----------------------\n";
  list<Critical_Pair_XPlor *> C;
  //cout << "creating and pruning pairs, starting with:\n";
  //cout << '\t' << P.size() << " pairs\n";
  //cout << '\t' << G.size() << " polynomials\n";
  // critical pairs with new polynomial
  unsigned m = G.size();
  for (unsigned i = 0; i < G.size(); ++i)
    C.push_back(new Critical_Pair_XPlor(i, r, strategy, G));
  //cout << "Created " << C.size() << " pairs with new polynomial\n";
  // apply Buchberger's lcm criterion to new pairs
  list<Critical_Pair_XPlor *> D;
  while (C.size() != 0) {
    Critical_Pair_XPlor * p = C.front();
    C.pop_front();
    if ((p->first()->leading_monomial().is_coprime(
            p->second()->leading_monomial()))
        or (no_triplet(p, C) and no_triplet(p, D))
        )
      D.push_back(p);
    else {
      //cout << "triplet prunes " << *p << endl;
      delete p;
    }
  }
  //cout << "After applying Buchberger's lcm criterion to new pairs, now have "
  //     << D.size() << " new pairs\n";
  // apply Buchberger's gcd criterion
  list<Critical_Pair_XPlor *> E;
  while (D.size() != 0) {
    Critical_Pair_XPlor * p = D.front();
    D.pop_front();
    if (!(p->first()->leading_monomial().is_coprime(
            p->second()->leading_monomial())))
      E.push_back(p);
    else {
      //cout << "gcd prunes " << *p << endl;
      delete p;
    }
  }
  //cout << "After applying Buchberger's gcd criterion to new pairs, now have "
  //     << E.size() << " new pairs\n";
  // apply Buchberger's lcm criterion to old pairs
  list<Critical_Pair_XPlor *> Q;
  while (P.size() != 0) {
    Critical_Pair_XPlor * p = P.front();
    P.pop_front();
    if (!(r->leading_monomial() | p->lcm())
          or lcm_alike(p->first()->leading_monomial(), r->leading_monomial(), p)
          or lcm_alike(p->second()->leading_monomial(), r->leading_monomial(), p)
       )
      Q.push_back(p);
    else {
      //cout << "triplet prunes " << *p << endl;
      delete p;
    }
  }
  //cout << "After applying Buchberger's lcm criterion to new pairs, now have "
  //     << Q.size() << " old pairs\n";
  P = Q;
  list <Critical_Pair_XPlor *> R;
  while (Pass.size() != 0) {
    Critical_Pair_XPlor * p = Pass.front();
    Pass.pop_front();
    if (!(r->leading_monomial() | p->lcm())
          or lcm_alike(p->first()->leading_monomial(), r->leading_monomial(), p)
          or lcm_alike(p->second()->leading_monomial(), r->leading_monomial(), p)
       )
      R.push_back(p);
    else {
      //cout << "\ttriplet prunes " << *p << endl;
      Pcancel[p->get_processor()].push_back(p);
    }
  }
  Pass = R;
  // add new pairs to old pairs
  for (Critical_Pair_XPlor * e : E)
    P.push_back(e);
  /*cout << "All pairs:\n";
  for (list<Critical_Pair_XPlor *>::iterator pi = P.begin(); pi != P.end(); ++pi)
    cout << '\t' << **pi << endl;
  cout << "----------------------\n";*/
}

typedef Critical_Pair_XPlor * Critical_Pair_XPlor_Ptr;

/**
  @brief used to pass inforation on a critical pair from one polynomial to another
*/
typedef struct {
  /** @brief index of first polynomial in the pair */
  int first;
  /** @brief index of second polynomial in the pair */
  int second;
  /** @brief pair&rsquo;s sugar */
  DEG_TYPE sugar;
} Critical_Pair_Communication;

list<Constant_Polynomial *> buchberger_explorer(
    const vector<Abstract_Polynomial *> &F,
    int method,
    unsigned strategy,
    WT_TYPE * strategy_weights,
    const int comm_id,
    const int comm_size
) {
  double r_bcast_time = 0;
  unsigned number_of_spolys = 0;
  vector<Abstract_Polynomial *> G; // basis
  list<Critical_Pair_XPlor *> P; // critical pairs
  list<Critical_Pair_XPlor *> Pass, * Pdel; // critical pair assignments
  list<Critical_Pair_XPlor *> R; // reduced polynomials
  list<Constant_Polynomial *> B; // end result
  Polynomial_Ring & Rx = (*(F.begin()))->base_ring();
  Monomial_Ordering * mord = (*(F.begin()))->monomial_ordering();
  // set up MPI_Datatype
  MPI_Datatype pair_type, pair_basetypes[2];
  MPI_Aint pair_offsets[2], extent, lb;
  int pair_block_counts[2];
  pair_offsets[0] = 0; pair_basetypes[0] = MPI_INT; pair_block_counts[0] = 2;
  MPI_Type_get_extent(MPI_INT, &lb, &extent);
  pair_offsets[1] = 2 * extent; pair_basetypes[1] = MPI_UNSIGNED_LONG_LONG;
  pair_block_counts[1] = 1;
  MPI_Type_create_struct(2, pair_block_counts, pair_offsets, pair_basetypes, &pair_type);
  MPI_Type_commit(&pair_type);
  // set up basis with generators
  vector<Abstract_Polynomial *>::const_iterator Fi = F.begin();
  NVAR_TYPE n = (*Fi)->number_of_variables();
  for (unsigned i = 0; i < F.size(); ++i, ++Fi)
  {
    Abstract_Polynomial * fo = *Fi;
    Constant_Polynomial * f = new Constant_Polynomial(*fo);
    f->set_strategy(new Poly_Sugar_Data(f));
    if (f->strategy() != nullptr) { f->strategy()->at_generation_tasks(); }
    if (comm_id == 0) // only control needs to take care of critical pairs
      P.push_back(new Critical_Pair_XPlor(i, strategy, f));
  }
  // main loop
  bool verbose = false;
  bool very_verbose = false;
  unsigned min_todo = 0;
  if (comm_id == 0) {
    Pdel = new list<Critical_Pair_XPlor *>[comm_size];
    Dense_Univariate_Rational_Polynomial ** HP
        = new Dense_Univariate_Rational_Polynomial *[comm_size];
    Dense_Univariate_Integer_Polynomial ** HFN
        = new Dense_Univariate_Integer_Polynomial *[comm_size];
    Dense_Univariate_Integer_Polynomial ** HSN
        = new Dense_Univariate_Integer_Polynomial *[comm_size];
    min_todo = P.size();
  }
  MPI_Bcast(&min_todo, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
  list<Monomial> T;
  Mutable_Polynomial * s;
  // create, reduce s-polynomials
  Critical_Pair_Communication p_in;
  while (min_todo != 0) {
    /*for (int i = 0; i < comm_size; ++i) {
      if (comm_id == i) {
        cout << comm_id << "'s pairs\n";
        for (Critical_Pair_XPlor * p : R)
          cout << '\t' << comm_id << ' ' << *p << endl;
      }
      MPI_Barrier(MPI_COMM_WORLD);
    }*/
    //double start_time = MPI_Wtime();
    if (comm_id == 0) {
      if (P.size() > 0) 
        sort_pairs_by_strategy(P);
      //cout << "estimate " << min_todo << " steps left\n";
      // first select comm_id pairs for reduction, send to coprocessors, reduce
      int i = 1;
      Critical_Pair_Communication p_new;
      for (/* already initialized */; !P.empty() and i < comm_size; ++i) {
        Critical_Pair_XPlor * p = P.front();
        p_new.first = p->first_index(); p_new.second = p->second_index();
        p_new.sugar = p->sugar();
        //report_front_pair(p, strategy);
        P.pop_front();
        p->set_processor(i);
        Pass.push_back(p);
        //cout << comm_id << " sending "
        //     << p_new.first << ',' << p_new.second << ':' << p_new.sugar
        //     << " to " << i << endl;
        MPI_Send(&p_new, 1, pair_type, i, 0, MPI_COMM_WORLD);
      }
      for (/* already initialized */; P.empty() and i < comm_size; ++i) {
        p_new.first = p_new.second = XPLOR_GENERATOR; p_new.sugar = 0;
        //cout << comm_id << " sending "
        //     << p_new.first << ',' << p_new.second << ':' << p_new.sugar
        //     << " to " << i << endl;
        MPI_Send(&p_new, 1, pair_type, i, 0, MPI_COMM_WORLD);
      }
      if (!P.empty()) {
        Critical_Pair_XPlor * p = P.front();
        P.pop_front();
        //report_front_pair(p, strategy);
        p_in.first = p->first_index(); p_in.second = p->second_index();
        p_in.sugar = p->sugar();
        p->set_processor(comm_id);
        Pass.push_back(p);
        //cout << comm_id << " sending "
        //     << p_in.first << ',' << p_in.second << ':' << p_in.sugar
        //     << " to " << comm_id << endl;
      } else {
        p_in.first = p_in.second = XPLOR_GENERATOR;
      }
    }
    //MPI_Barrier(MPI_COMM_WORLD);
    if (comm_id != 0)
      MPI_Recv(&p_in, 1, pair_type, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    if (p_in.first != XPLOR_GENERATOR) {
      if (p_in.second == XPLOR_GENERATOR)
        //R.push_back(new Critical_Pair_XPlor(p_in.first, p_in.second, p_in.sugar, F));
        R.push_back(new Critical_Pair_XPlor(p_in.first, SUGAR_STRATEGY, F[p_in.first]));
      else
        R.push_back(new Critical_Pair_XPlor(p_in.first, p_in.second, SUGAR_STRATEGY, G));
    }
    //r_bcast_time += MPI_Wtime() - start_time;
    //MPI_Barrier(MPI_COMM_WORLD);
    double start_time = MPI_Wtime();
    for (Critical_Pair_XPlor * p : R) {
      // make s-poly
      if ((s = p->s_polynomial()) == nullptr) {
        s = p->s_polynomial(method, strategy);
        ++number_of_spolys;
      }
      //cout << comm_id << " reducing s-poly "
      //     << p->first_index() << ',' << p->second_index() << endl;
      //double start_time = MPI_Wtime();
      if (!s->is_zero())
        reduce_over_basis<vector<Abstract_Polynomial *>>(&s, G, comm_id);
      //r_bcast_time += MPI_Wtime() - start_time;
      p->set_spoly(s);
    }
    r_bcast_time += MPI_Wtime() - start_time;
    /*for (int i = 0; i < comm_size; ++i) {
      if (comm_id == i)
        cout << comm_id << " finished reduction with " << R.size() << " polynomials\n";
      MPI_Barrier(MPI_COMM_WORLD);
    }*/
    // create, compare Hilbert functions
    //MPI_Barrier(MPI_COMM_WORLD);
    //double start_time = MPI_Wtime();
    unsigned winning_index = R.size() + 1;
    Monomial * wt = nullptr;
    Dense_Univariate_Rational_Polynomial * wHP = nullptr;
    Dense_Univariate_Integer_Polynomial * wHS = nullptr;
    Dense_Univariate_Integer_Polynomial * hn = nullptr;
    Dense_Univariate_Integer_Polynomial * hsn = nullptr;
    Dense_Univariate_Rational_Polynomial * hp = nullptr;
    // loop through critical pairs to find optimal HF
    list<Critical_Pair_XPlor *>::iterator Rbest = R.end();
    //double start_time = MPI_Wtime();
    for (
         list<Critical_Pair_XPlor *>::iterator Ri = R.begin();
         Ri != R.end();
         /* advance manually */
    ) {
      s = (*Ri)->s_polynomial();
      if (s->is_zero()) {
        //cout << '\t' << comm_id << ' ' << (*Ri)->first_index() << ','
        //     << (*Ri)->second_index() << ": reduced to zero\n";
        Critical_Pair_XPlor * p = *Ri;
        if (Ri == R.begin()) {
          R.pop_front();
          Ri = R.begin();
        } else {
          list<Critical_Pair_XPlor *>::iterator Rdel = Ri;
          ++Ri;
          R.erase(Rdel);
        }
        delete p;
      } else {
        T.push_back(s->leading_monomial());
        unsigned n = T.front().num_vars();
        hn = hilbert_numerator_bigatti(T);
        hsn = hilbert_second_numerator(n, hn);
        unsigned d = ideal_dimension(n, hn, hsn);
        hp = hilbert_polynomial(n, d, T, hn, hsn);
        T.pop_back();
        if (wt == nullptr) {
          wt = &(s->leading_monomial()); wHP = hp; wHS = hn; Rbest = Ri;
        } else {
          if (second_HF_smaller(wHP, wHS, hp, hn)) {
            wt = &(s->leading_monomial()); wHP = hp; wHS = hn; Rbest = Ri;
          }
        }
        ++Ri;
      }
    }
    //r_bcast_time = MPI_Wtime() - start_time;
    //MPI_Barrier(MPI_COMM_WORLD);
    // move result of winning reduction to basis; return others to P
    /*cout << comm_id << " about to send Hilbert data\n";
    for (int i = 1; i < comm_size; ++i) {
      if (comm_id == i and wt != nullptr) {
        cout << comm_id << "'s Hilbert data for " << *wt << ":\n";
        cout << '\t' << *wHP << endl;
        cout << '\t' << *wHS << endl;
      }
      MPI_Barrier(MPI_COMM_WORLD);
    }*/
    //MPI_Barrier(MPI_COMM_WORLD);
    //double start_time = MPI_Wtime();
    int64_t size_data = -1;
    if (comm_id != 0) {
      if (wHP != nullptr) {
        size_data = wHP->degree();
        MPI_Send(&size_data, 1, MPI_INT64_T, 0, 0, MPI_COMM_WORLD);
        MPI_Send(wHP->numerator_array(), size_data + 1, MPI_INT64_T, 0, 0, MPI_COMM_WORLD);
        MPI_Send(wHP->denominator_array(), size_data + 1, MPI_UINT64_T, 0, 0, MPI_COMM_WORLD);
        delete wHP;
        size_data = wHS->degree();
        MPI_Send(&size_data, 1, MPI_INT64_T, 0, 0, MPI_COMM_WORLD);
        MPI_Send(wHS->coefficient_array(), size_data + 1, MPI_INT64_T, 0, 0, MPI_COMM_WORLD);
        delete wHS;
      } else {
        size_data = -1;
        MPI_Send(&size_data, 1, MPI_INT64_T, 0, 0, MPI_COMM_WORLD);
      }
    }
    int winners_id = -1;
    //MPI_Barrier(MPI_COMM_WORLD);
    if (comm_id == 0) {
      /*cout << "sorting results\n";
      cout << "My Hilbert data:\n";
      if (wt == nullptr)
        cout << "\tnothing\n";
      else {
        cout << '\t' << *wt << endl;
        cout << '\t' << *wHP << endl;
        cout << '\t' << *wHS << endl;
      }*/
      if (wHP != nullptr) winners_id = 0;
      for (int i = 1; i < comm_size; ++i) {
        int64_t size_data;
        MPI_Recv(&size_data, 1, MPI_INT64_T, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        if (size_data >= 0) {
          int64_t * nums = new int64_t[size_data + 1];
          uint64_t * denoms = new uint64_t[size_data + 1];
          MPI_Recv(nums, size_data + 1, MPI_INT64_T, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
          MPI_Recv(denoms, size_data + 1, MPI_UINT64_T, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
          /*cout << "Received from " << i << ":\n\t";
          for (unsigned k = 0; k < size_data + 1; ++k)
            cout << nums[k] << '/' << denoms[k] << " , ";
          cout << endl;*/
          Dense_Univariate_Rational_Polynomial * hp_in
              = new Dense_Univariate_Rational_Polynomial(size_data, nums, denoms);
          //cout << "Received from " << i << ": " << *hp_in << endl;
          delete [] nums; delete [] denoms;
          MPI_Recv(&size_data, 1, MPI_UINT64_T, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
          int64_t * coefs = new int64_t[size_data + 1];
          MPI_Recv(coefs, size_data + 1, MPI_INT64_T, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
          Dense_Univariate_Integer_Polynomial * hn_in
              = new Dense_Univariate_Integer_Polynomial(size_data, coefs);
          //cout << "Received from " << i << ": " << *hn_in << endl;
          delete [] coefs;
          if (wHP == nullptr) { // process 0 has exhausted its pairs
            winners_id = i; wHP = hp_in; wHS = hn_in;
          }
          else if (second_HF_smaller(wHP, wHS, hp_in, hn_in)) {
            if (wHP != nullptr) delete wHP;
            if (wHS != nullptr) delete wHS;
            winners_id = i; wHP = hp_in; wHS = hn_in;
          } else {
            delete hp_in; delete hn_in;
          }
        }
        /*cout << "Received from " << i << ":\n";
        cout << '\t' << *hp_in << endl;
        cout << '\t' << *hn_in << endl;*/
      }
      if (wHP == nullptr) winners_id = -1;
      /*cout << winners_id << " has the best HF\n";
      if (winners_id != -1) {
        cout << '\t' << *wHP << endl;
        cout << '\t' << *wHS << endl;
      }*/
    }
    MPI_Bcast(&winners_id, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (winners_id >= 0) {
      uint64_t * r_bcast;
      uint64_t r_bcast_size;
      uint64_t r_bcast_sugar;
      Constant_Polynomial * r;
      if (comm_id == winners_id) {
        // broadcast winning polynomial
        s = (*Rbest)->s_polynomial();
        r = new Constant_Polynomial(*s);
        //cout << "selected " << (*Rbest)->first_index()
        //     << ',' << (*Rbest)->second_index() << ": " << r->leading_monomial()
        //     << endl;
        Poly_Sugar_Data * sd = static_cast<Poly_Sugar_Data *>(s->strategy());
        r->set_strategy(sd);
        r_bcast_sugar = sd->poly_sugar();
        s->set_strategy(nullptr);
        delete s;
        Critical_Pair_XPlor * p = *Rbest;
        R.erase(Rbest);
        delete p;
        r_bcast = r->serialized(r_bcast_size);
      }
      //MPI_Barrier(MPI_COMM_WORLD);
      MPI_Bcast(&r_bcast_size, 1, MPI_UINT64_T, winners_id, MPI_COMM_WORLD);
      if (comm_id != winners_id)
        r_bcast = new uint64_t [r_bcast_size];
      //double start_time = MPI_Wtime();
      MPI_Bcast(r_bcast, r_bcast_size, MPI_UINT64_T, winners_id, MPI_COMM_WORLD);
      //r_bcast_time += MPI_Wtime() - start_time;
      MPI_Bcast(&r_bcast_sugar, 1, MPI_UINT64_T, winners_id, MPI_COMM_WORLD);
      // other processors need to create a copy of the new polynomial
      if (comm_id != winners_id) {
        r_bcast_size /= n + 1; // adjust size from # of words to # of terms
        r = new Constant_Polynomial(Rx, mord, r_bcast_size, r_bcast);
        Poly_Sugar_Data * rd = new Poly_Sugar_Data(r);
        rd->force_sugar(r_bcast_sugar);
        r->set_strategy(rd);
      }
      //MPI_Barrier(MPI_COMM_WORLD);
      //r_bcast_time += MPI_Wtime() - start_time;
      delete [] r_bcast;
      //double start_time = MPI_Wtime();
      if (comm_id == 0) {
        //cout << "added " << G.size() << ": " << r->leading_monomial() << endl;
        very_verbose = false;
        if (very_verbose) { cout << "\tin full "; r->println(); }
        very_verbose = false;
        gm_update_explorer(P, Pass, Pdel, G, r, strategy);
        // remove useless pairs that were sent out to processors
        for (int i = 1; i < comm_size; ++i) {
          int outgoing = Pdel[i].size();
          MPI_Send(&outgoing, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
          while (Pdel[i].size() > 0) {
            Critical_Pair_XPlor * p = Pdel[i].front();
            Critical_Pair_Communication p_new;
            Pdel[i].pop_front();
            p_new.first = p->first_index();
            p_new.second = p->second_index();
            MPI_Send(&p_new, 1, pair_type, i, 0, MPI_COMM_WORLD);
            delete p;
          }
        }
        // now delete my own redundant pairs
        while (Pdel[0].size() > 0) {
          Critical_Pair_XPlor * p = Pdel[0].front();
          Pdel[0].pop_front();
          list<Critical_Pair_XPlor *>::iterator Ri = R.begin();
          while (
              Ri != R.end() and
              ((*Ri)->first_index() != p->first_index() or
              (*Ri)->second_index() != p->second_index())
          )
            ++Ri;
          // if the polynomial reduced to zero, it will not be in R
          if (Ri != R.end())
            R.erase(Ri);
          delete p;
        }
      } else {
        int incoming;
        MPI_Recv(&incoming, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        for (int j = 0; j < incoming; ++j) {
          MPI_Recv(&p_in, 1, pair_type, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
          list<Critical_Pair_XPlor *>::iterator Ri = R.begin();
          while (
              Ri != R.end() and
              ((*Ri)->first_index() != p_in.first or
              (*Ri)->second_index() != p_in.second)
          )
            ++Ri;
          // if Ri has reduced to zero, it will not be in R
          if (Ri != R.end()) {
            Critical_Pair_XPlor * p = *Ri;
            R.erase(Ri);
            delete p;
          }
        }
      }
      //r_bcast_time += MPI_Wtime() - start_time;
      //MPI_Barrier(MPI_COMM_WORLD);
      G.push_back(r);
      T.push_back(r->leading_monomial());
    }
    //r_bcast_time += MPI_Wtime() - start_time;
    /*for (int i = 0; i < comm_size; ++i) {
      if (comm_id == i) {
        cout << comm_id << "'s basis\n";
        for (Abstract_Polynomial * g : G)
          cout << '\t' << comm_id << ' ' << *g << endl;
      }
      MPI_Barrier(MPI_COMM_WORLD);
    }*/
    //double start_time = MPI_Wtime();
    if (comm_id == 0) min_todo = P.size();
    unsigned my_todo = R.size();
    unsigned all_todo;
    MPI_Reduce(&my_todo, &all_todo, 1, MPI_UNSIGNED, MPI_MAX, 0, MPI_COMM_WORLD);
    if (comm_id == 0)
      min_todo = min_todo < all_todo ? all_todo : min_todo;
    MPI_Bcast(&min_todo, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    //r_bcast_time += MPI_Wtime() - start_time;
    //cout << comm_id << " understands there to be " << min_todo << " pairs remaining\n";
    /*for (int i = 0; i < comm_size; ++i) {
      if (comm_id == i)
        for (Critical_Pair_XPlor * p : R)
          cout << comm_id << "has " << p->first_index() << ',' << p->second_index() << endl;
      MPI_Barrier(MPI_COMM_WORLD);
    }*/
    //MPI_Barrier(MPI_COMM_WORLD);
  }
  if (comm_id != 0) {
    for (Abstract_Polynomial * g : G)
      delete g;
  }
  MPI_Type_free(&pair_type);
  cout << comm_id << " reduced " << number_of_spolys << endl;
  MPI_Barrier(MPI_COMM_WORLD);
  int total_spolys;
  MPI_Reduce(&number_of_spolys, &total_spolys, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  double max_bcast_time;
  MPI_Reduce(&r_bcast_time, &max_bcast_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  for (unsigned i = 0; i < comm_size; ++i) {
    if (comm_id == i)
      cout << comm_id << " took " << r_bcast_time << " milliseconds on timed segment(s).\n";
    MPI_Barrier(MPI_COMM_WORLD);
  }
  if (comm_id == 0) {
    cout << total_spolys << " s-polynomials computed and reduced\n";
    cout << max_bcast_time << " seconds spent on timed segment(s)\n";
    // cleanup
    cout << G.size() << " polynomials before interreduction\n";
    //check_correctness(G, strategy);
    list<Abstract_Polynomial *> G_final;
    for (Abstract_Polynomial * g : G)
      G_final.push_back(g);
    G_final = reduce_basis(G_final);
    cout << G_final.size() << " polynomials after interreduction\n";
    //set<Constant_Polynomial *, smaller_lm> B;
    for (Abstract_Polynomial * g : G_final) {
      B.push_back(new Constant_Polynomial(*g));
      //if (F.find(g) == F.end()) delete g;
    }
    delete [] Pdel;
  }
  return B;
}

#endif