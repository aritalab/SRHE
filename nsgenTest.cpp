#include <iostream>
#include "nsgen.h"

using namespace std;  

void noise_test (noise& ns, long dim) {
  cout << "noise_test ---------------------------" << endl;
  long test_num = 5;
  ZZ_p::init(ZZ(1024));
  
  cout << "Uniform : p = " << ZZ_p::modulus() << endl;
  for (long c = 0; c < test_num; ++c) {
    vec_ZZ_p e;
    uniform(e, dim);
    cout << "e = " << e << endl;
  }

  cout << "Hamming : hwt = " << ns.hwt << endl;
  for (long c = 0; c < test_num; ++c) {
    vec_ZZ s;
    ns.hamming(s, dim);
    cout << "s = " << s << endl;
  }

  cout << "Gaussian : sigma = " << ns.sigma << endl;
  for (long c = 0; c < test_num; ++c) {  
    vec_ZZ e;
    ns.gaussian(e, dim);
    cout << "e = " << e << endl;
  }
}

int main (int argc, char *argv[]) {
  cout << "# --- nsgenTest ---" << endl;

  long dim = 18;
  noise ns(/*label=*/"noise 3.2", /*sigma=*/3.2, /*hwt=*/4);

  noise_test(ns, dim);
}
