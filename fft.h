#ifndef FFT_H
#define FFT_H

#include <complex>
#include <NTL/ZZ.h>

using namespace std;
using namespace NTL;

// auxiliary constants
static double const Pi=4.0*atan(1.0); // Pi=3.1415..

// n: a power of 2
void fft_forward (complex<double>* b, complex<double>* a, long n);
void fft_backward (complex<double>* b, complex<double>* a, long n);

void fft_convolution (complex<double>* c, complex<double>* a, complex<double>* b, long n);
void fft_convolution (double* c, double* a, double* b, long n);

// n: need not be a power of 2
void naive_convolution (double* c, double* a, double* b, long dim);
void convolution (double* c, double* a, double* b, long dim);
//void convolution (long* c, long* a, long* b, long dim);

#endif /*FFT_H*/
