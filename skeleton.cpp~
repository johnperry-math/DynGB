/**
  \brief implementation of classes for double description method
  \author John Perry
  \version 1.0
  \date October 2014
  @copyright The University of Southern Mississippi
*/

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

#ifndef SKELETON_C
#define SKELETON_C

#include <iostream>
#include <cstdlib>
#include <algorithm>
using std::set_union;

#include "system_constants.hpp"

#include "goda.hpp"
#include "skeleton.hpp"

namespace LP_Solvers {

Edge::Edge(const Ray &first_ray, const Ray &second_ray)
    : first(first_ray), second(second_ray)
{
  if (first < second)
    first.swap(second);
}

Edge::Edge(const Edge &old_edge)
    : first(old_edge.first), second(old_edge.second)
{
  // nothing to do
}

ostream & operator<<(ostream & ostr, const Edge &e)
{
  ostr << "{ " << e.first << " , " << e.second << " }";
  return ostr;
}

bool operator==(const Edge &e1, const Edge &e2)
{ return e1.first == e2.first && e1.second == e2.second; }

bool operator < (const Edge &first, const Edge &second)
{
  bool result = true;
  if (first.first < second.first)
  {
    // do nothing
  }
  else if (second.first < first.first)
    result = false;
  else // first entries equal; look at second
    if (first.second < second.second)
    {
      // do nothing
    }
    else
      result = false;
  return result;
}

Edge & Edge::operator=(const Edge &other)
{
  if (!(*this == other))
  {
    first = other.first;
    second = other.second;
  }
  return *this;
}

void Skeleton::common_initialization(NVAR_TYPE dimension)
{
  //cout << "creating skeleton with dimension " << dimension << endl;
  dim = dimension;
  // initialize the constraints
  CONSTR_TYPE * constr_coords = new CONSTR_TYPE[dim];
  for (NVAR_TYPE i = 0; i < dim; ++i) constr_coords[i] = 0;
  // add the constraint xi >= 0 for each i = 0, ..., dim - 1
  for (NVAR_TYPE i = 0; i < dim; ++i)
  {
    constr_coords[i] = 1;
    constraints.push_back(Constraint(dim, constr_coords));
    constr_coords[i] = 0;
  }
  delete [] constr_coords;
  // initialize the rays, one for each axis
  for (NVAR_TYPE i = 0; i < dim; ++i)
  {
    Ray new_ray(dim, i);
    rays.insert(new_ray);
  }
  // initialize the edges
  // currently, all the rays are adjacent
  for (auto riter = rays.begin(); riter != rays.end(); ++riter)
    for (auto siter = rays.begin(); siter != rays.end(); ++siter)
      if (*riter != *siter)
      {
        Edge new_edge(*riter, *siter);
        edges.insert(new_edge);
      }
}

Skeleton::Skeleton(NVAR_TYPE dimension)
{
  common_initialization(dimension);
}

Skeleton::Skeleton(NVAR_TYPE dimension, const vector<Constraint> &constraints)
        //: Skeleton(dimension)
{
  common_initialization(dimension);
  solve(constraints);
}

Skeleton::Skeleton(const Skeleton &old_skeleton)
        : constraints(old_skeleton.constraints),
          edges(old_skeleton.edges), dim(old_skeleton.dim)
          
{
  rays = old_skeleton.rays;
  // nothing more to do
}

bool Skeleton::copy(const LP_Solver * other) {
  const Skeleton * old_skeleton = dynamic_cast<const Skeleton *>(other);
  if (old_skeleton != nullptr) {
    constraints.clear(); rays.clear(); edges.clear();
    for (const Constraint & c : old_skeleton->constraints)
      constraints.push_back(c);
    for (const Ray & r : old_skeleton->rays)
      rays.emplace(r);
    for (const Edge & e : old_skeleton->edges)
      edges.emplace(e);
    /*constraints = old_skeleton->constraints;
    rays = old_skeleton->rays;
    edges = old_skeleton->edges;*/
    dim = old_skeleton->dim;
  }
  return (old_skeleton != nullptr);
}

Skeleton::~Skeleton()
{
}

bool Skeleton::solve(const Constraint &constraint)
{
  // cout << "processing constraint " << constraint << endl;
  // innocent until proven guilty
  bool consistent = true;
  // sort the rays into the ones above, below, or on the constraint
  set<Ray> rays_above, rays_below, rays_on;
  for (auto riter = rays.begin(); riter != rays.end(); ++riter)
  {
    DOTPROD_TYPE dp = (*riter) * constraint; // overloaded * as dot product :-)
    if (dp > 0)
    {
      rays_above.insert(*riter);
      // cout << *riter << " is above constraint\n";
    }
    else if (dp < 0)
    {
      rays_below.insert(*riter);
      // cout << *riter << " is below constraint\n";
    }
    else
    {
      Ray old_ray = *riter;
      rays_on.insert(old_ray);
      // cout << *riter << " is on constraint\n";
    }
  }
  // cout << rays_above.size() << " rays above; " << rays_below.size() << " rays below; " << rays_on.size() << " rays on\n";
  // check for constitency
  if (rays_above.size() == 0)
  {
    consistent = false;
    //cout << "long test inconsistent\n";
  }
  // proceed only if constraint is consistent, and *not* redundant;
  // redundancy can be checked by making sure
  // that at least one ray is below the constraint
  if (consistent and rays_below.size() != 0)
  {
    set<Edge> edges_above, edges_on;
    for (auto eiter = edges.begin(); eiter != edges.end(); ++eiter)
    {
      Edge e = *eiter;
      Ray u = e.get_first_ray();
      Ray v = e.get_second_ray();
      //cout << "edge " << u << ',' << v << endl;
      // identify the edges that lie above and on this constraint
      if ((u*constraint >= 0) and (v*constraint >= 0))
      {
        // cout << "old edge preserved: " << u << ',' << v << "\n";
        edges_above.insert(e);
      }
    }
    for (auto eiter = edges.begin(); eiter != edges.end(); ++eiter)
    {
      Edge e = *eiter;
      Ray u = e.get_first_ray();
      Ray v = e.get_second_ray();
      DOTPROD_TYPE a = u*constraint;
      DOTPROD_TYPE b = v*constraint;
      // identify edges that pass through the constraint
      // (one ray above, one ray below)
      if (a > 0 and b < 0)
      {
        Ray w = a*v - b*u;
        w.simplify_ray();
        rays_on.insert(w);
        edges_on.insert(Edge(u,w));
        // cout << "new ray (u,v) is " << w << " with constraints " << endl;
      }
      else if (b > 0 and a < 0)
      {
        Ray w = b*u - a*v;
        w.simplify_ray();
        rays_on.insert(w);
        edges_on.insert(Edge(v,w));
        // cout << "new ray (v,u) is " << w << " with constraints " << endl;
      }
    }
    // clear the old rays, add the new ones (above and on the constraint)
    rays.clear();
    for (auto riter = rays_above.begin(); riter != rays_above.end(); ++riter)
    {
      rays.insert(*riter);
    }
    // cout << "inserted rays above; rays on is\n";
    // for (auto riter = rays_on.begin(); riter != rays_on.end(); ++riter) { cout << '\t' << *riter << endl; }
    for (auto riter = rays_on.begin(); riter != rays_on.end(); ++riter)
    {
      // cout << "inserting " << *riter << endl;
      //cout << "return value: " << *(get<0>(rays.insert(*riter)));
      rays.insert(*riter);
      //for (auto siter = rays.begin(); siter != rays.end(); ++siter) { cout << '\t' << *siter << endl; }
    }
    //cout << rays.size() << " rays\n";
    // for (auto riter = rays.begin(); riter != rays.end(); ++riter) { cout << '\t' << *riter << endl; }
    // add the good constraint
    constraints.push_back(constraint);
    // determine new edges
    set<Edge> edges_new = adjacencies_by_graphs(rays_on);
    // combine new edges with old ones that are known to be valid
    edges = union_of_edge_sets(union_of_edge_sets(edges_above, edges_on), edges_new);
    //cout << edges.size() << " edges\n";
    //for (auto eiter = edges.begin(); eiter != edges.end(); ++eiter) { cout << *eiter << ' '; } cout << '\n';
  }
  return consistent;
}

bool Skeleton::solve(const vector<Constraint> &new_constraints)
{
  // innocent until proven guilty
  bool consistent = true;
  //cout << "adding " << new_constraints.size() << "constraints\n";
  //for (const constraint & c : new_constraints)
  //  cout << '\t' << c << endl;
  // process each constraint sequentially
  for (
        auto nciter = new_constraints.begin();
        consistent and nciter != new_constraints.end();
        ++nciter
      )
  {
    // perform short test of consistency first
    consistent = is_consistent(*nciter) and solve(*nciter);
    //if (!consistent)
    //{
    //  cout << "inconsistent\n";
    //  cout << "failed ray: " << *nciter;
    //  cout << "skeleton: \n" << *this;
    //}
  }
  //cout << rays.size() << " corner vectors\n";
  return consistent;
}

int number_of_common_constraints(
    //vector<bool> &a, vector<bool> &b
    bool * a, bool * b, unsigned m
)
{
  int result = 0;
  /*for (auto aiter = a.begin(); aiter != a.end(); ++aiter)
  {
    //cout << "checking " << *aiter << " in other: " << (b.find(*aiter) != b.end()) << endl;
    if (b.find(*aiter) != b.end())
      ++result;
  } */
  for (unsigned i = 0; i < m; ++i)
    if (a[i] and b[i]) ++result;
  return result;
}

//vector<bool> intersections_of_active_constraints(
    //vector<bool> &a, vector<bool> &b
void intersections_of_active_constraints(
    bool * a, bool * b, bool * result, unsigned m
)
{
  // highly unoptimized, but off the top of my head i don't know how to do better
  //vector<bool> result(a.size());
  /*for (auto aiter = a.begin(); aiter != a.end(); ++aiter)
    if (b.find(*aiter) != b.end())
      result.insert(*aiter);*/
  for (unsigned i = 0; i < m; ++i)
    result[i] = (a[i] and b[i]); 
  //return result;
}

bool is_first_subset_of_second(
    //vector<bool> & a, vector<bool> & b
    bool * a, bool * b, unsigned m
)
{
  // highly unoptimized, but off the top of my head i don't know how to do better
  bool result = true;
  for (unsigned i = 0; result and i < m; ++i)
    if (a[i]) result = b[i];
  return result;
}

set<Edge> union_of_edge_sets(const set<Edge> & a, const set<Edge> & b)
{
  // optimized with a hint for the position (riter) of the new element
  set<Edge> result;
  for (const Edge & e : a) result.insert(e);
  for (const Edge & e : b) result.insert(e);
  return result;
}

set<Edge> Skeleton::adjacencies_by_graphs(const set<Ray> & new_rays)
{
  static unsigned long invocations;
  set<Edge> new_edges;
  set<Ray> tested_rays;
  bool *   Zu = new bool [constraints.size()] { false };
  bool *   Zv = new bool [constraints.size()] { false };
  bool * w_active = new bool [constraints.size()] { false };
  bool *      Zuv = new bool [constraints.size()] { false };
  // loop through each new ray, examining active constraints shared with other rays
  for (auto riter = new_rays.begin(); riter != new_rays.end(); ++riter)
  {
    Ray u = *riter;
    tested_rays.insert(u);
    which_constraints_active_at(u, Zu);
    // D's rays have at least dim - 2 active constraints in common with u
    // (see Proposition 3 in Zolotych's paper)
    set<Ray> D;
    for (auto siter = new_rays.begin(); siter != new_rays.end(); ++siter)
      if (*riter != *siter)
      {
        Ray v = *siter;
        which_constraints_active_at(v, Zv);
        //cout << "checking constraints of " << u << " against " << v  << " for " << dim << endl;
        //if (number_of_common_constraints(*Zu, *Zv) >= dim - 2)
        if (number_of_common_constraints(Zu, Zv, constraints.size()) >= dim - 2)
        {
          //cout << "accept " << u << ',' << v << " from active constraints\n";
          D.insert(v);
        } else {
          //cout << "reject " << u << ',' << v << " from active constraints\n";
        }
      }
    // check u with each v in D, making sure their active constraints
    // are not a subset of the active constraints of any w in D
    // (see Proposition 4 (graph test) in Zolotych's paper)
    unsigned ijk = 0;
    for (auto diter = D.begin(); diter != D.end(); ++diter)
    {
      Ray v = *diter;
      if (tested_rays.find(v) == tested_rays.end()) // avoid doubling edges
      {
        which_constraints_active_at(v, Zv);
        // WARNING: I have commented out the following line, because it seems
        // unnecessary: v is in D iff the size of the intersection is at least
        // dim - 2. If there are unexpected bugs, this commenting should be
        // reconsidered.
        // if (intersections_of_active_constraints(Zu, Zv).size() >= dim - 2)
        {
          bool can_be_added = true;
          intersections_of_active_constraints(Zu, Zv, Zuv, constraints.size());
          for (
               auto dditer = D.begin();
               dditer != D.end() and can_be_added;
               ++dditer
              )
          {
            Ray w = *dditer;
            if (!(w == v)) {
              which_constraints_active_at(w, w_active);
              if (is_first_subset_of_second(Zuv, w_active, constraints.size()))
              {
                //cout << "rejecting " << u << ',' << v << " because of " << w << endl;
                can_be_added = false;
              }
            }
          }
          if (can_be_added)
          {
            Edge new_edge(u, v);
            //cout << "edge " << new_edge << " passes all criteria\n";
            new_edges.insert(new_edge);
          }
        }
      }
    }
  }
  delete [] Zu;
  delete [] Zv;
  delete [] w_active;
  delete [] Zuv;
  return new_edges;
}

ostream & operator << (ostream & ostr, const Skeleton &skel)
{
  // header, start constraints
  ostr << "Skeleton defined by constraints" << endl;
  for (
       vector<Constraint>::const_iterator citer=skel.constraints.begin();
       citer != skel.constraints.end();
       ++citer
      )
    ostr << '\t' << *citer << endl;
  // rays
  ostr << "has " << skel.rays.size() << " rays" << endl;
  for (auto riter=skel.rays.begin(); riter != skel.rays.end(); ++riter)
    ostr << '\t' << *riter << endl;
  //edges
  ostr << "connected in " << skel.edges.size() << " edges" << endl;
  for (auto eiter=skel.edges.begin(); eiter != skel.edges.end(); ++eiter)
    ostr << '\t' << *eiter << endl;
  // footer
  ostr << "End of skeleton" << endl;
  return ostr;
}

Skeleton & Skeleton::operator=(const Skeleton & other)
{
  rays.clear();
  edges.clear();
  constraints.clear();
  dim = other.dim;
  for (
       auto siter = other.rays.begin();
       siter != other.rays.end();
       ++siter
      )
    rays.insert(*siter);
  for (
       auto eiter = other.edges.begin();
       eiter != other.edges.end();
       ++eiter
      )
    edges.insert(*eiter);
  for (
        vector<Constraint>::const_iterator citer = other.constraints.begin();
        citer != other.constraints.end();
        ++citer
      )
    constraints.push_back(*citer);
  return *this;
}

}

#endif