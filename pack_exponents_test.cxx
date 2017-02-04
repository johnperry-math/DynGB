#include <bitset>
#include <cstring>
#include <iostream>

using std::cout; using std::endl; using std::bitset;

int main(int argc, char ** argv) {
  unsigned long * a = new unsigned long[8];
  unsigned long * b = new unsigned long[8];
  unsigned long long * c = new unsigned long long[4] {0} ;
  unsigned long long * d = new unsigned long long[4] {0} ;
  for (unsigned long long i = 0; i < 8; ++i) {
    a[i] = b[i] = i;
    if ((i % 2) == 0)
      c[i/2] = d[i/2] += i;
    else
      c[i/2] = d[i/2] += (i << 32);
  }
  cout << "size of unsigned long " << sizeof(unsigned long) << endl;
  cout << "size of unsigned long long " << sizeof(unsigned long long) << endl;
  cout << endl;
  for (unsigned i = 0; i < 8; ++i)
    cout << bitset<64>(a[i]) << endl << bitset<64>(b[i]) << endl;
  cout << endl;
  for (unsigned i = 0; i < 4; ++i)
    cout << bitset<64>(c[i]) << endl << bitset<64>(d[i]) << endl;
  if (argc > 1) {
    if (!strcmp(argv[1],"1")) {
      for (unsigned i = 0; i < 128000003; ++i)
        for (unsigned j = 0; j < 8; ++j)
          if (a[j] != b[j]) cout << "fault " << j << endl;
    } else {
      for (unsigned i = 0; i < 128000003; ++i)
        for (unsigned j = 0; j < 4; ++j)
          if (c[j] != d[j]) cout << "fault " << j << endl;
    }
  }
}