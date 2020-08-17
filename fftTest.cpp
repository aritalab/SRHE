#include <iostream>
#include <string>
#include <NTL/ZZ.h>

#include "fft.h"

using namespace std;
using namespace NTL;

void vec_print (string label, long* a, long dim) {
  cout << label << "[";
  for (long i = 0; i < dim; i++) cout << a[i] << ",";
  cout << "]" << endl;
}

void vec_print (string label, double* a, long dim) {
  cout << label << "[";
  for (long i = 0; i < dim; i++) cout << a[i] << ",";
  cout << "]" << endl;
}

void vec_print (string label, complex<double>* a, long dim) {
  cout << label << "[";
  for (long i = 0; i < dim; i++) cout << a[i] << ",";
  cout << "]" << endl;
}

void basic_test (long q, long dim) {
  cout << "basic_test ---------------------------" << endl;

  long testNum = 1000;
  
  for (long cnt = 0; cnt < testNum; ++cnt) {  
    complex<double> a[dim];
    for (long i = 0; i < dim; ++i) {
      a[i] = (double) RandomBnd(q);
    }

    //vec_print("a = ", a, dim);

    complex<double> b[dim];
    fft_forward(b, a, dim);
    //vec_print("b = ", b, dim);

    complex<double> c[dim];
    fft_backward(c, b, dim);
    //vec_print("c = ", c, dim);

    for (long i = 0; i < dim; ++i) {
      if (abs(c[i]-a[i]) > 0.001) {
	cout << "Error in basic_test at the " << cnt << "-th test" << endl;
	cout << "a[" << i << "] = " << a[i] << endl;
	cout << "c[" << i << "] = " << c[i] << endl;
	exit(-1);
      }
    }
  }
  cout << "OK: basic_test" << endl;
  
}

void convolution_test (long q, long dim) {
  cout << "convolution_test ---------------------------" << endl;

  long testNum = 1000;
  
  for (long cnt = 0; cnt < testNum; ++cnt) {
    double a[dim], b[dim], c[dim], d[dim];
    
    for (long i = 0; i < dim; ++i) a[i] = (double) RandomBnd(q);
    for (long i = 0; i < dim; ++i) b[i] = (double) RandomBnd(q);

    //vec_print("a = ", a, dim);
    //vec_print("b = ", b, dim);

    naive_convolution(c, a, b, dim);
    convolution(d, a, b, dim);

    //vec_print("c = ", c, dim);
    //vec_print("d = ", d, dim);

    for (long i = 0; i < dim; ++i) {
      if (abs(c[i]-d[i]) > 0.001) {
	cout << "Error in convolution_test at the " << cnt << "-th test" << endl;
	cout << "c[" << i << "] = " << c[i] << endl;
	cout << "d[" << i << "] = " << d[i] << endl;
	cout << "(c-d)[" << i << "] = " << c[i] - d[i] << endl;
	exit(-1);
      }
    }
  }
  cout << "OK: convolution_test" << endl;
}

int main (int argc, char *argv[]) {
  cout << "# --- fftTest ---" << endl;

  long q = 1024-1;
  basic_test(q, 8);
  
  long dim = 1<<NextPowerOfTwo(2*19-1);
  cout << "dim = " << dim << endl;

  convolution_test(q, dim);

  exit(0);
}
