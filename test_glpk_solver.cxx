#include <iostream>
using std::cout; using std::endl;

#include "system_constants.hpp"

#include "glpk_solver.hpp"

int main() {
  NVAR_TYPE n = 5;
  GLPK_Solver gs(5);
  CONSTR_TYPE cdata[5] { 1, -1, 0, 0, 0 };
  constraint c1(5, cdata);
  cdata[0] = -1; cdata[1] = 4;
  constraint c2(5, cdata);
  vector<constraint> cs;
  cs.push_back(c1); cs.push_back(c2);
  cout << "adding " << c1 << " , " << c2 << endl;
  bool success = gs.solve(cs);
  if (not success)
    cout << "inconsistent system!\n";
  else {
    const set<ray> & rays = gs.get_rays();
    cout << "consistent system! solutions are:\n";
    for (const ray & r : rays)
      cout << '\t' << r << endl;
    cout << ray_sum(rays) << endl;
  }
  cs.clear();
  cdata[0] = 1; cdata[1] = 0; cdata[2] = -1;
  constraint c3(5, cdata);
  cdata[0] = 0; cdata[1] = -1; cdata[2] = 1;
  constraint c4(5, cdata);
  cdata[0] = 0; cdata[1] = 3; cdata[2] = -1;
  constraint c5(5, cdata);
  cs.push_back(c3); cs.push_back(c4); cs.push_back(c5);
  cout << "adding " << c3 << " , " << c4 << " , " << c5 << endl;
  success = gs.solve(cs);
  if (not success)
    cout << "inconsistent system!\n";
  else {
    const set<ray> & rays = gs.get_rays();
    cout << "consistent system! solutions are:\n";
    for (const ray & r : rays)
      cout << '\t' << r << endl;
    cout << ray_sum(rays) << endl;
  }
  /*cs.clear();
  cdata[0] = -1; cdata[1] = 0; cdata[2] = 0;
  for (CONSTR_TYPE i = 102; i > 1; --i) {
    cdata[3] = i;
    constraint new_constr(5, cdata);
    bool success = gs.solve(new_constr);
    success = not success;
  }
  const set<ray> & rays = gs.get_rays();
  cout << "final solution:\n";
  for (const ray & r : rays)
    cout << '\t' << r << endl;
  cout << ray_sum(rays) << endl;
  cout << endl;*/
  cout << "copying original\n";
  GLPK_Solver gs2(gs);
  cs.clear();
  cdata[0] = -1; cdata[1] = 0; cdata[2] = 0; cdata[3] = 0; cdata[4] = 1;
  constraint another_constr(5, cdata);
  cout << "adding " << another_constr << endl;
  success = gs2.solve(another_constr);
  cout << "success? " << (success ? 'T' : 'F') << endl;
  const set<ray> new_rays = gs2.get_rays();
  cout << "new solution:\n";
  for (const ray & r : new_rays)
    cout << '\t' << r << endl;
  cout << ray_sum(new_rays) << endl;
  cout << endl;
  cout << "now making it inconsistent...\n";
  cs.clear();
  cdata[0] = 1; cdata[1] = -5; cdata[4] = 0;
  constraint final_constr(5, cdata);
  cout << "adding " << final_constr << endl;
  success = gs2.solve(final_constr);
  cout << "success? " << (success ? 'T' : 'F') << endl;
}