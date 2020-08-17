#ifndef NSGEN_H
#define NSGEN_H

#include <cassert>
#include <cmath>
#include <NTL/ZZ.h>
#include <NTL/vec_ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/vec_ZZ_p.h>
#include "fft.h"

using namespace std;
using namespace NTL;

class noise {
 public:
  string label;
  double sigma;
  long hwt;

 public:
  noise () {}
  
  noise (string l, double s, long h) {
    label = l;
    sigma = s;
    hwt = h;
  }

  void hamming (vec_ZZ& s, long dim);
  void gaussian (vec_ZZ& e, long dim);
};

ostream& operator<< (ostream& s, const noise& ns);

void uniform (vec_ZZ& e, const ZZ& q, long dim);
void uniform (vec_ZZ_p& e, long dim);

void sample_hwt (vec_ZZ& v, long hwt, long dim);
 
void sample_gaussian (double& e1, double& e2, const double stdev);
void sample_gaussian (double* f, const double stdev, const long dim);

#endif /*NSGEN_H*/
