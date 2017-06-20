#include <iostream>
using std::cout; using std::endl;

#include "system_constants.hpp"

#include "glpk_solver.hpp"
using LP_Solvers::GLPK_Solver;
using LP_Solvers::Constraint; using LP_Solvers::Ray;

int main() {
  NVAR_TYPE n = 5;
  GLPK_Solver gs(5);
  CONSTR_TYPE cdata[5] { 1, -1, 0, 0, 0 };
  Constraint c1(5, cdata);
  cdata[0] = -1; cdata[1] = 4;
  Constraint c2(5, cdata);
  vector<Constraint> cs;
  cs.push_back(c1); cs.push_back(c2);
  cout << "adding " << c1 << " , " << c2 << endl;
  bool success = gs.solve(cs);
  if (not success)
    cout << "inconsistent system!\n";
  else {
    const set<Ray> & rays = gs.get_rays();
    cout << "consistent system! solutions are:\n";
    for (const Ray & r : rays)
      cout << '\t' << r << endl;
    cout << ray_sum(rays) << endl;
  }
  cs.clear();
  cdata[0] = 1; cdata[1] = 0; cdata[2] = -1;
  Constraint c3(5, cdata);
  cdata[0] = 0; cdata[1] = -1; cdata[2] = 1;
  Constraint c4(5, cdata);
  cdata[0] = 0; cdata[1] = 3; cdata[2] = -1;
  Constraint c5(5, cdata);
  cs.push_back(c3); cs.push_back(c4); cs.push_back(c5);
  cout << "adding " << c3 << " , " << c4 << " , " << c5 << endl;
  success = gs.solve(cs);
  if (not success)
    cout << "inconsistent system!\n";
  else {
    const set<Ray> & rays = gs.get_rays();
    cout << "consistent system! solutions are:\n";
    for (const Ray & r : rays)
      cout << '\t' << r << endl;
    cout << ray_sum(rays) << endl;
  }
  /*cs.clear();
  cdata[0] = -1; cdata[1] = 0; cdata[2] = 0;
  for (CONSTR_TYPE i = 102; i > 1; --i) {
    cdata[3] = i;
    Constraint new_constr(5, cdata);
    bool success = gs.solve(new_constr);
    success = not success;
  }
  const set<Ray> & rays = gs.get_rays();
  cout << "final solution:\n";
  for (const Ray & r : rays)
    cout << '\t' << r << endl;
  cout << ray_sum(rays) << endl;
  cout << endl;*/
  cout << "copying original\n";
  GLPK_Solver gs2(gs);
  cs.clear();
  cdata[0] = -1; cdata[1] = 0; cdata[2] = 0; cdata[3] = 0; cdata[4] = 1;
  Constraint another_constr(5, cdata);
  cout << "adding " << another_constr << endl;
  success = gs2.solve(another_constr);
  cout << "success? " << (success ? 'T' : 'F') << endl;
  const set<Ray> new_rays = gs2.get_rays();
  cout << "new solution:\n";
  for (const Ray & r : new_rays)
    cout << '\t' << r << endl;
  cout << ray_sum(new_rays) << endl;
  cout << endl;
  cout << "now making it inconsistent...\n";
  cs.clear();
  cdata[0] = 1; cdata[1] = -5; cdata[4] = 0;
  Constraint final_constr(5, cdata);
  cout << "adding " << final_constr << endl;
  success = gs2.solve(final_constr);
  cout << "success? " << (success ? 'T' : 'F') << endl;
}