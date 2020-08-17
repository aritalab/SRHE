#include "nsgen.h"

void noise::hamming (vec_ZZ& s, long dim) {
  sample_hwt(s, hwt, dim);
}

void noise::gaussian (vec_ZZ& e, long dim) {
  double f[dim];
  sample_gaussian(f, sigma, dim);

  e.SetLength(dim);
  for (long i = 0; i < dim; ++i)
    e[i] = conv<ZZ>(round(f[i]));  
}

ostream& operator<< (ostream& s, const noise& ns) {
  s << ns.label + ": sigma = " << ns.sigma << ": hwt = " << ns.hwt << endl;
  return s;
}

// Imported from NumbTh.cpp of helib
void sample_gaussian (double& e1, double& e2, const double stdev) {
  static long const bignum = 0xfffffff;
  // THREADS: C++11 guarantees these are initialized only once

  // Uses the Box-Muller method to get two Normal(0,stdev^2) variables
  double r1 = (1+RandomBnd(bignum))/((double)bignum+1);
  double r2 = (1+RandomBnd(bignum))/((double)bignum+1);
  double theta=2*Pi*r1;
  double rr= sqrt(-2.0*log(r2))*stdev;

  assert(rr < 8*stdev); // sanity-check, no more than 8 standard deviations

  // Generate two Gaussians RV's
  e1 = rr*cos(theta);
  e2 = rr*sin(theta);
  // cout << "e1 = " << e1 << ", e2 = " << e2 << endl;
}

// sample gaussian 
void sample_gaussian (double* f, const double stdev, const long dim) {
  for (long j = 0; j < dim; j+=2) {
    double e1, e2;
    sample_gaussian(e1, e2, stdev);
    f[j] = e1;
    if(j+1 < dim) f[j+1] = e2;
  }
}

void uniform (vec_ZZ& e, const ZZ& q, long dim) {
  e.SetLength(dim);
  for (long i = 0; i < dim; ++i) e[i] = RandomBnd(q);
}

void uniform (vec_ZZ_p& e, long dim) {
  e.SetLength(dim);
  for (long i = 0; i < dim; ++i) e[i] = random_ZZ_p();
}

void sample_hwt (vec_ZZ& v, long hwt, long dim) {
  v.SetLength(dim);
  clear(v);
  if (hwt > dim) hwt = dim;

  long i = 0;
  while (i < hwt) {
    long u = RandomBnd(dim);
    if (v[u] == 0) {
      long b = RandomBnd(2) << 1;
      v[u] = b - 1;
      i++;
    }
  }
}




