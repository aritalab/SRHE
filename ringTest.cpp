#include <iostream>

#include "ring.h"
#include "params.h"

using namespace std;
using namespace NTL;

void apply_omega_test (param& prm, ZZ& q) {
  cout << "apply_omega_test for " << prm.label << endl;

  ZZ_pContext context;
  context.save();
  ZZ_p::init(q);

  ring& rg = prm.rg;
  ZZ_pX phi = (ZZ_pX(INIT_MONO, rg.mR)-1)/(ZZ_pX(INIT_MONO, 1)-1);
  long g = rg.gR;

  long test_num = 3;
  for (long cnt = 0; cnt < test_num; ++cnt) {
    vec_ZZ_p a, b, w, ax, bx, wx;
    wx.SetLength(g);
    ZZ_pX a1, b1, c, c1;
    
    uniform(a, g);
    //a[0] = ZZ_p(1);
    uniform(b, g);
    //b[0] = ZZ_p(2);

    prm.applyOmega(ax, a);
    //cout << "ax = " << ax << endl;
    prm.applyOmega(bx, b);
    //cout << "bx = " << bx << endl;
    mult(wx, ax, bx);
    //cout << "wx = " << wx << endl;
    prm.applyOmegaInv(w, wx);
    //cout << "w = " << w << endl;
    rg.to_eta(c1, w);

    rg.to_eta(a1, a);
    rg.to_eta(b1, b);
    c = (a1 * b1) % phi;

    if (c != c1) {
      cout << "NG: apply_omega_test" << endl;
      cout << "cnt = " << cnt << endl;
      cout << "c  = (" << c[0] << ", " << c[1] << ",...)" << endl;
      cout << "c1 = (" << c1[0] << ", " << c1[1] << ",...)" << endl;
      exit(-1);
    }
  }
  
  context.restore();
  cout << "OK: apply_omega_test" << endl;
}

void power_decomp_test (param& prm) {
  cout << "--- power_decomp_test for " << prm.label << " --- " << endl;

  long l = prm.l;
  long dim = prm.rg.gR;
  vec_ZZ zero; zero.SetLength(dim);

  long test_num = 3;
  for (long cnt = 0; cnt < test_num; ++cnt) {
    vec_ZZ x = zero;
    for (long i = 0; i < dim; ++i) x[i] = RandomBnd(power(ZZ(prm.rg.pR), l));

    Vec<vec_ZZ> ds;
    power_decomp(ds, l, x);

    for (long i = 0; i < dim; ++i) {
      ZZ z(0);
      for (long j = 0; j < l; ++j) z += ds[j][i] * power(ZZ(prm.rg.pR), j);
      if (x[i] != z) {
	  cout << "NG: power_decomp_test" << endl;
	  cout << "  be: " << x[i] << endl;
	  cout << "  is: " << z << endl;
	  exit(-1);
      }
    }
  }

  cout << "OK: power_decomp_test" << endl;
}

void right_shift_test (param& prm) {
  cout << "--- right_shift_test for " << prm.label << " --- " << endl;

  long l = prm.l;
  long dim = prm.rg.gR;
  vec_ZZ zero; zero.SetLength(dim);

  long test_num = 3;
  for (long cnt = 0; cnt < test_num; ++cnt) {
    vec_ZZ x = zero;
    for (long i = 0; i < dim; ++i) x[i] = RandomBnd(power(ZZ(prm.rg.pR), l));

    long m = 2;
    vec_ZZ z = zero;
    right_shift(z, x, m, l);
    
    for (long i = 0; i < dim; ++i) {
      ZZ w = x[i]/power(ZZ(2),m);
      if (w != z[i]) {
	cout << "NG: right_shift_test" << endl;
	cout << "  be: " << w << endl;
	cout << "  is: " << z[i] << endl;
	exit(-1);
      }
    }
  }

  cout << "OK: right_shift_test" << endl;
}

// Tests if (a x a1) (b x b1) = (a b) x (a1 b1)
// ('x' denotes tensor product)
void kronecker_test (param& prm, ZZ q) {
  cout << "--- kronecker_test for " << prm.label << endl;

  ZZ_p::init(q);

  ring& rg = prm.rg;
  ring& lrg = *(rg.left);
  ring& rrg = *(rg.right);

  long ldim = lrg.gR;
  long rdim = rrg.gR;
  long dim = ldim * rdim;

  vec_ZZ lzero; lzero.SetLength(ldim);
  vec_ZZ_p lzero_p; lzero_p.SetLength(ldim);
  vec_ZZ rzero; rzero.SetLength(rdim);
  vec_ZZ_p rzero_p; rzero_p.SetLength(rdim);
  vec_ZZ zero; zero.SetLength(dim);
  vec_ZZ_p zero_p; zero_p.SetLength(dim);

  long test_cnt = 3;
  for (long c = 0; c < test_cnt; ++c) {
    vec_ZZ a, b;
    uniform(a, q, ldim);
    uniform(b, q, ldim);
    
    vec_ZZ a1, b1;
    uniform(a1, q, rdim);
    uniform(b1, q, rdim);
   
    // A = a x a1, B = b x b1
    vec_ZZ A = zero, B = zero;
    tensor_prod(A, a, a1, q);
    tensor_prod(B, b, b1, q);

    // C = A B
    vec_ZZ_p Ax = zero_p, Bx = zero_p, Cx = zero_p, _C = zero_p;
    rg.applyOmega(Ax, conv<vec_ZZ_p>(A));
    rg.applyOmega(Bx, conv<vec_ZZ_p>(B));
    mult(Cx, Ax, Bx);
    rg.applyOmegaInv(_C, Cx);
    vec_ZZ C = zero;
    C = conv<vec_ZZ>(_C);
    //cout << "C[0]  = " << C[0] << endl;

    // c = a b
    vec_ZZ_p ax = lzero_p, bx = lzero_p, cx = lzero_p, _c = lzero_p;
    lrg.applyOmega(ax, conv<vec_ZZ_p>(a));
    lrg.applyOmega(bx, conv<vec_ZZ_p>(b));
    mult(cx, ax, bx);
    lrg.applyOmegaInv(_c, cx);

    // c1 = a1 b1
    vec_ZZ_p a1x = rzero_p, b1x = rzero_p, c1x = rzero_p, _c1 = rzero_p;
    rrg.applyOmega(a1x, conv<vec_ZZ_p>(a1));
    rrg.applyOmega(b1x, conv<vec_ZZ_p>(b1));
    mult(c1x, a1x, b1x);
    rrg.applyOmegaInv(_c1, c1x);

    vec_ZZ C1 = zero;
    tensor_prod(C1, conv<vec_ZZ>(_c), conv<vec_ZZ>(_c1), q);
    //cout << "C1[0] = " << C1[0] << endl;

    if (C != C1) {
      cout << "NG: kronecker_test" << endl;
      cout << "  be: " << C[0] << endl;
      cout << "  is: " << C1[0] << endl;
      exit(-1);
    }
  }
}

void kronecker_test (param& prm) {
  if (prm.rg.type == prime) return;

  ZZ q = power(ZZ(prm.rg.pR), prm.r);
  kronecker_test(prm, q);

  q = power(ZZ(prm.rg.pR), prm.l);
  kronecker_test(prm, q);

  cout << "OK: kronecker_test" << endl;  
}

void ring_test (param& prm) {
  cout << "--- ring_test for " << prm.label << " --- " << endl;

  if (prm.rg.type == prime) {
    ZZ q = power(ZZ(prm.rg.pR), prm.r);
    apply_omega_test(prm, q);
    
    q = power(ZZ(prm.rg.pR), prm.l);
    apply_omega_test(prm, q);
  }

  power_decomp_test(prm);
  right_shift_test(prm);

  if (prm.rg.type == composite)
    kronecker_test(prm);
}

void do_test (param& prm) {
  cout << prm;
  ring_test(prm);  
}

void do_test (param& prm, param& prm_tmp, param& prm_big) {
  do_test(prm);
  do_test(prm_tmp);
  do_test(prm_big);
}

int main (int argc, char *argv[]) {
  /* Do Test */

  long prec = 2*composed_param_big.r - composed_param_big.l;

  left_param.init(prec);
  left_param_tmp.init(prec);
  left_param_big.init(prec);
  do_test(left_param, left_param_tmp, left_param_big);      
  
  right_param.init(prec);
  right_param_tmp.init(prec);
  right_param_big.init(prec);
  do_test(right_param, right_param_tmp, right_param_big);      

  composed_param.init(prec);
  composed_param_tmp.init(prec);
  composed_param_big.init(prec);
  do_test(composed_param, composed_param_tmp, composed_param_big);      
}
