#include <complex>
#include <iostream>
#include <NTL/ZZ.h>
#include "fft.h"

// n: a power of 2
void fft_forward (complex<double>* b, complex<double>* a, long n) {
  if (n == 1) {
    b[0] = a[0];
    return;
  }

  complex<double> zeta_n = polar(1.0, 2.0 * Pi / (double) n);
  complex<double> zeta = 1.0;

  complex<double> aeven[n/2], aodd[n/2];
  complex<double> beven[n/2], bodd[n/2];
  for (long i = 0; i < n/2; i++) {
    aeven[i] = a[2*i];
    aodd[i] = a[2*i+1];
  }
  fft_forward(beven, aeven, n/2);
  fft_forward(bodd, aodd, n/2);

  for (long k = 0; k < n/2; k++) {
    b[k] = beven[k] + zeta * bodd[k];
    b[k+n/2] = beven[k] - zeta * bodd[k];
    zeta = zeta * zeta_n;
  }
}

// n: a power of 2
void fft_backward_1 (complex<double>* b, complex<double>* a, long n) {
  if (n == 1) {
    b[0] = a[0];
    return;
  }

  complex<double> zeta_n = polar(1.0, -2.0 * Pi / (double) n);
  complex<double> zeta = 1.0;

  complex<double> aeven[n/2], aodd[n/2];
  complex<double> beven[n/2], bodd[n/2];
  for (long i = 0; i < n/2; i++) {
    aeven[i] = a[2*i];
    aodd[i] = a[2*i+1];
  }
  fft_backward_1(beven, aeven, n/2);
  fft_backward_1(bodd, aodd, n/2);

  for (long k = 0; k < n/2; k++) {
    b[k] = beven[k] + zeta * bodd[k];
    b[k+n/2] = beven[k] - zeta * bodd[k];
    zeta = zeta * zeta_n;
  }
}

void fft_backward (complex<double>* b, complex<double>* a, long n) {
  fft_backward_1(b, a, n);
  for (long k = 0; k < n; k++) {
    b[k] /= (double)n;
  }
}

// n: a power of 2
//   c(x) = a(x) * b(x) mod x^n - 1
void fft_convolution (complex<double>* c, complex<double>* a, complex<double>* b, long n) {
  complex<double> ax[n], bx[n], cx[n];

  fft_forward(ax, a, n); 
  fft_forward(bx, b, n); 
  for (long i = 0; i < n; ++i)
    cx[i] = ax[i] * bx[i];
  fft_backward(c, cx, n);
}

void fft_convolution (double* c, double* a, double* b, long n) {
  complex<double> aa[n], bb[n], cc[n];
  for (long i = 0; i < n; ++i) {
    aa[i] = a[i]; bb[i] = b[i];
  }

  fft_convolution(cc, aa, bb, n);
  for (long i = 0; i < n; ++i)
    c[i] = cc[i].real();
}

//--------------------
// c(x) = a(x) * b(x) mod x^dim - 1
void naive_convolution (double* c, double* a, double* b, long dim) {
  for (long k = 0; k < dim; ++k) {
    c[k] = 0;
    for (long i = 0; i <=  k; ++i) {
      double w;
      w = a[i] * b[k - i];
      c[k] += w;
    }
    for (long i = k+1; i < dim; ++i) {
      double w;
      w = a[i] * b[dim + k - i];
      c[k] += w;
    }
  }  
}

// c(x) = a(x) * b(x) mod x^dim - 1
void convolution (double* c, double* a, double* b, long dim) {
  if ((1<<NextPowerOfTwo(dim-1)) == dim) {
    //cout << "dim is powerof2" << endl;
    fft_convolution(c, a, b, dim);
    return;
  }
  
  long dim2 = 1<<NextPowerOfTwo(2*dim-1);
  double aa[dim2], bb[dim2], cc[dim2];

  for (long n = 0; n < dim2; ++n) {
    if (n < dim) aa[n] = a[n];
    else aa[n] = 0;
  }

  for (long n = 0; n < dim2; ++n) {
    if (n == 0) bb[n] = b[n];
    else if (n < dim) bb[n] = b[n];
    else if (n < dim2-dim+1) bb[n] = 0;
    else // n >= dim2-dim+1
      bb[n] = b[n - (dim2-dim)];
  }

  fft_convolution(cc, aa, bb, dim2);
  
  for (long n = 0; n < dim; ++n)
    c[n] = cc[n];
}

