#include <iostream>
#include <NTL/RR.h>
#include "ring.h"
#include "she.h"
#include "util.h"
#include "params.h"

using namespace std;
using namespace NTL;

void encode_test (param& prm) {
  cout << "--- encode_test for " << prm.label << " ---" << endl;

  long test_num = 5;
  for (long cnt = 0; cnt < test_num; ++cnt) {
    PlainText x(prm);
    x.uniform();

    PlainText z(prm);
    encode_pt(prm, z, x);

    PlainText x1(prm);
    decode_pt(prm, x1, z);

    if (x != x1) {
      cout << "NG: encode_test: (" << cnt << ")" << endl;
      cout << " x = " << x << endl;
      cout << " z = " << z << endl;
      cout << "x1 = " << x1 << endl;
      exit(-1);
    }
  }

  cout << "OK: encode_test" << endl;
}

void dec_test (param& prm, const SKey& s) {
  cout << "--- dec_test for " << prm.label << " ---" << endl;

  for (long i = 0; i < 5; ++i) {
    PlainText m(prm), m1(prm);
    m.uniform();

    CipherText ct(prm);
    encrypt(prm, ct, s, m);

    if (dec_ng(prm, m1, s, ct, m)) {
      cout << "NG: dec" << endl;
      cout << "m  = " << m << endl;
      cout << "m1 = " << m1 << endl;
      exit(-1);
    }
  }

  cout << "OK: dec_test" << endl;
}

void div_by_2_test (param& prm, const SKey& s) {
  cout << "--- div_by_2_test for " << prm.label << " ---" << endl;

  long test_num = 5;
  for (long cnt = 0; cnt < test_num; ++cnt) {
    PlainText m(prm), m1(prm), m2(prm);
    m.uniform();
    for (long i = 0; i < prm.rg.gR; ++i) m[i] = m[i] - (m[i] % 2);
    for (long i = 0; i < prm.rg.gR; ++i) m2[i] = m[i]/2;

    CipherText ct(prm);
    encrypt(prm, ct, s, m);

    CipherText dt(prm);
    div_by_2(dt, ct);

    ZZ v = decrypt(prm, m1, s, dt);
    for (long i = 0; i < dt.dim; ++i) {
      if (((m2[i] - m1[i]) % (ct.t/2)) != 0) {
	cout << "NG: div_by_2" << endl;
	cout << "be: " << VectorCopy(m2.data, 10) << endl;
	cout << "is: " << VectorCopy(m1.data, 10) << endl;
	exit(-1);
      }
    }
  }  

  cout << "OK: div_by_2_test" << endl;
}

void add_test (param& prm, const SKey& s) {
  cout << "--- add_test for " << prm.label << " ---" << endl;

  long test_num = 5;
  for (long cnt = 0; cnt < test_num; ++cnt) {
    PlainText m1(prm), m2(prm), m3(prm), _m(prm);
    m1.uniform();
    m2.uniform();
    add(m3, m1, m2);

    CipherText ct1(prm), ct2(prm), ct3(prm);
    encrypt(prm, ct1, s, m1);
    encrypt(prm, ct2, s, m2);

    add_ct(ct3, ct1, ct2);
    if (dec_ng(prm, _m, s, ct3, m3)) {
      cout << "NG: add" << endl;
      cout << "  be: " << m3 << endl;
      cout << "  is: " << _m << endl;
      exit(-1);
    }
  }

  cout << "OK: add_test" << endl;
}

void plain_add_test (param& prm, const SKey& s) {
  cout << "--- plain_add_test for " << prm.label << " ---" << endl;

  long test_num = 5;
  for (long cnt = 0; cnt < test_num; ++cnt) {
    PlainText m1(prm), m2(prm), m3(prm), _m(prm);
    m1.uniform();
    m2.uniform();
    add(m3, m1, m2);

    CipherText ct2(prm), ct3(prm);
    encrypt(prm, ct2, s, m2);
    
    plain_add_ct(ct3, m1, ct2);
    if (dec_ng(prm, _m, s, ct3, m3)) {
      cout << "NG: plain_add_test" << endl;
      cout << "  be: " << VectorCopy(m3.data, 10) << endl;
      cout << "  is: " << VectorCopy(_m.data, 10)  << endl;
      exit(-1);
    }
  }

  cout << "OK: plain_add_test" << endl;  
}

void plain_mult_test (param& prm, const SKey& s) {
  cout << "--- plain_mult_test for " << prm.label << " ---" << endl;

  long test_num = 5;
  for (long cnt = 0; cnt < test_num; ++cnt) {
    PlainText m1(prm), m2(prm), m3(prm), _m(prm);
    m1.uniform();
    m2.uniform();
    mult(m3, m1, m2);

    CipherText ct2(prm), ct3(prm);
    encrypt(prm, ct2, s, m2);
    
    plain_mult_ct(ct3, m1, ct2);
    if (dec_ng(prm, _m, s, ct3, m3)) {
      cout << "NG: plain_mult_test" << endl;
      cout << "  be: " << m3 << endl;
      cout << "  is: " << _m << endl;
      exit(-1);
    }
  }

  cout << "OK: plain_mult_test" << endl;  
}

void hint_gen_test (param& prm, const SKey& s) {
  cout << "--- hint_gen_test for " << prm.label << " ---" << endl;
  ZZ_p::init(prm.qq);
  ZZ& sq = prm.sR;

  long test_num = 1;
  for (long cnt = 0; cnt < test_num; ++cnt) {
    ZZ w(1);
    for (long j = 0; j < prm.lw; ++j) {
      cout << "-- " << j << " ------" << endl;
      vec_ZZ_p e;
      e.SetLength(prm.rg.gR);

      // e = H[j][0] + H[j][1] s.x
      mult(e, conv<vec_ZZ_p>(s.Hint.H[j][1]), conv<vec_ZZ_p>(s.x));
      add(e, conv<vec_ZZ_p>(s.Hint.H[j][0]), e);
      //prm.applyGammaInv(e, e);
      prm.applyOmegaInv(e, e);

      // f = sR s2 g^T
      vec_ZZ_p f =  conv<ZZ_p>(w) * conv<vec_ZZ_p>(s.s2);
      f = conv<ZZ_p>(sq) * f;

      e = e - f;

      vec_ZZ ee;
      ee.SetLength(prm.rg.gR);
      center_lift(ee, e);
      cout << "|e| = " << inf_norm(ee) << endl;
      
      w *= prm.w;
    }
  }

  //cout << "OK: hint_gen_test" << endl;
}

void hint_gen_test2 (param& prm, const SKey& s) {
  cout << "--- hint_gen_test2 for " << prm.label << " ---" << endl;
  ZZ& q0 = prm.qq;
  ZZ& q = prm.q;
  
  long test_num = 1;
  for (long cnt = 0; cnt < test_num; ++cnt) {
    ZZ w(1);
    for (long j = 0; j < prm.lw; ++j) {
      cout << "-- " << j << " ------" << endl;

      ZZ_p::init(q0);
      
      vec_ZZ_p _e;
      _e.SetLength(prm.rg.gR);

      // e = H[j][0] + H[j][1] s.x
      mult(_e, conv<vec_ZZ_p>(s.Hint.H[j][1]), conv<vec_ZZ_p>(s.x));
      add(_e, conv<vec_ZZ_p>(s.Hint.H[j][0]), _e);
      //prm.applyGammaInv(_e, _e);
      prm.applyOmegaInv(_e, _e);

      vec_ZZ e;
      e.SetLength(prm.rg.gR);
      rescale(e, conv<vec_ZZ>(_e), q0, q);


      ZZ_p::init(q);
      
      // f = s2 g^T
      vec_ZZ_p f =  conv<ZZ_p>(w) * conv<vec_ZZ_p>(s.s2);

      vec_ZZ_p e1 = conv<vec_ZZ_p>(e) - f;
      vec_ZZ ee;
      ee.SetLength(prm.rg.gR);
      center_lift(ee, e1);
      cout << "|e| = " << inf_norm(ee) << endl;
      
      w *= prm.w;
    }
  }

  //cout << "OK: hint_gen_test2" << endl;
}

void decomp_test (param& prm) {
  cout << "--- decomp_test for " << prm.label << " ---" << endl;

  ZZ& q = prm.qq;
  ZZ_p::init(q);
  
  long test_num = 5;
  for (long cnt = 0; cnt < test_num; ++cnt) {
    vec_ZZ c;
    c.SetLength(prm.rg.gR);
    for (long i = 0; i < c.length(); ++i) c[i] = RandomBnd(q);

    Vec<vec_ZZ> d;
    d.SetLength(prm.lw);
    for (long l = 0; l < d.length(); ++l) d[l].SetLength(prm.rg.gR);
    Decomp(prm, d, c);

    vec_ZZ_p y;
    y.SetLength(prm.rg.gR);
    
    ZZ w(1);
    for (long l = 0; l < prm.lw; ++l) {
      vec_ZZ_p dd = conv<ZZ_p>(w) * conv<vec_ZZ_p>(d[l]);
      y += dd;
      w *= prm.w;
    }
    
    if (conv<vec_ZZ>(y) != c) {
      cout << "NG: decomp" << endl;
      cout << "  be: " << c[0] << endl;
      cout << "  is: " << (conv<vec_ZZ>(y))[0] << endl;
      exit(-1);
    }
  }

  cout << "OK: decomp_test" << endl;
}

// Typically, sin = s2x (s2=s^2), sout = sx
// -- H : encryption of sin
// -- d = apply_hint(H, c)  --> d(sout) ~ sin c
void hint_apply_test (param& prm, const hint Hint, const vec_ZZ& sin, const vec_ZZ& sout) {
  cout << "--- hint_apply_test for " << prm.label << " ---" << endl;
  long test_num = 3;
  for (long cnt = 0; cnt < test_num; ++cnt) {
    ZZ& q = prm.q;
    ZZ_p::init(q);

    vec_ZZ a, b, c;
    c.SetLength(prm.rg.gR);
    for (long i = 0; i < c.length(); ++i) c[i] = RandomBnd(q);
    
    apply_hint(prm, a, b, Hint, c, q);

    vec_ZZ_p _a, _b, _e;
    _e.SetLength(prm.rg.gR);
    
    //prm.applyGamma(_a, conv<vec_ZZ_p>(a));
    prm.applyOmega(_a, conv<vec_ZZ_p>(a));
    prm.applyOmega(_b, conv<vec_ZZ_p>(b));
    mult(_e, _b, conv<vec_ZZ_p>(sout));
    _e = _a + _e;

    vec_ZZ_p _c;
    prm.applyOmega(_c, conv<vec_ZZ_p>(c));
    mult(_c, conv<vec_ZZ_p>(sin), _c);

    _c = _c - _e;
    //prm.applyGammaInv(_c, _c);
    prm.applyOmegaInv(_c, _c);

    center_lift(c, _c);
    cout << "inf_norm(noise) = " << inf_norm(c) << endl;
  }
  
  //cout << "OK: hint_apply_test" << endl;
}

void key_switch_test (param& prm, const SKey& s, const SKey& t, const hint& Hint) {
  cout << "--- key_switch_test for " << prm.label << " ---" << endl;
  
  long test_num = 3;
  for (long cnt = 0; cnt < test_num; ++cnt) {
    PlainText m(prm), m1(prm);
    m.uniform();
    CipherText ct(prm);
    encrypt(prm, ct, t, m);

    // ct/t --> dt/s
    CipherText dt(prm);
    key_switch(dt, Hint, ct);
   
    if (dec_ng(prm, m1, s, dt, m)) {
      cout << "NG: dec" << endl;
      cout << "be: " << VectorCopy(m.data, 10) << endl;
      cout << "is: " << VectorCopy(m1.data, 10) << endl;
      exit(-1);
    }
  }
  
  cout << "OK: key_switch_test" << endl;
}

void direct_mult_test (param& prm, const SKey& s) {
  cout << "--- direct_mult_test for " << prm.label << " ---" << endl;
  
  long test_num = 5;
  for (long cnt = 0; cnt < test_num; ++cnt) {
    PlainText m1(prm), m2(prm), m3(prm), _m(prm);
    m1.uniform();
    m2.uniform();
    mult(m3, m1, m2);

    CipherText ct1(prm), ct2(prm);
    encrypt(prm, ct1, s, m1);
    encrypt(prm, ct2, s, m2);

    CipherText2 dt(prm);
    direct_mult_ct(dt, ct1, ct2);

    if (dec_ng2(prm, _m, s, dt, m3)) {
      cout << "NG: direct_mult_test" << endl;
      cout << "  be: " << m3 << endl;
      cout << "  is: " << _m << endl;
      exit(-1);
    }    
  }

  cout << "OK: direct_mult_test" << endl;
}

void mult_test (param& prm, const hint& Hint, const SKey& s) {
  cout << "--- mult_test for " << prm.label << " ---" << endl;
  
  long test_num = 5;
  for (long cnt = 0; cnt < test_num; ++cnt) {
    PlainText m1(prm), m2(prm), m3(prm), _m(prm);
    m1.uniform();
    m2.uniform();
    mult(m3, m1, m2);

    CipherText ct1(prm), ct2(prm);
    encrypt(prm, ct1, s, m1);
    encrypt(prm, ct2, s, m2);

    CipherText dt(prm);
    mult_ct(dt, Hint, ct1, ct2);

    if (dec_ng(prm, _m, s, dt, m3)) {
      cout << "NG: mult_test" << endl;
      cout << "  be: " << VectorCopy(m3.data, 10) << endl;
      cout << "  is: " << VectorCopy(_m.data, 10) << endl;
      exit(-1);
    }    
  }

  cout << "OK: mult_test" << endl;
}

void square_test (param& prm, const hint& Hint, const SKey& s) {
  cout << "--- square_test for " << prm.label << " ---" << endl;
  
  long test_num = 5;
  for (long cnt = 0; cnt < test_num; ++cnt) {
    PlainText m1(prm), m3(prm), _m(prm);
    m1.uniform();
    mult(m3, m1, m1);

    CipherText ct1(prm);
    encrypt(prm, ct1, s, m1);

    CipherText dt(prm);
    //square_ct(dt, Hint, ct1);
    square_ct_debug(dt, Hint, ct1, s);
  
    if (dec_ng(prm, _m, s, dt, m3)) {
      cout << "NG: square_test" << endl;
      cout << "  be: " << VectorCopy(m3.data, 10) << endl;
      cout << "  is: " << VectorCopy(_m.data, 10) << endl;
      exit(-1);
    }    
  }

  cout << "OK: square_test" << endl;
}

void power_test (param& prm, const hint& Hint, const SKey& s, const long n) {
  cout << "--- power_test for " << prm.label << " ---" << endl;
  
  long test_num = 3;
  for (long cnt = 0; cnt < test_num; ++cnt) {
    cout << "[" << cnt << "]-----------------------" << endl;
    
    Vec<PlainText> m;
    for (long i = 0; i < n + 1; ++i) {
      PlainText _m(prm);
      append(m, _m);
    }

    Vec<CipherText> ct;
    for (long i = 0; i < n + 1; ++i) {
      CipherText _ct(prm);
      append(ct, _ct);
    }

    m[0].uniform();
    encrypt(prm, ct[0], s, m[0]);

    PlainText _m(prm);
    if (dec_ng(prm, _m, s, ct[0], m[0])) {
      cout << "NG: power_test" << endl;
      cout << "  be: " << VectorCopy(m[0].data,10) << endl;
      cout << "  is: " << VectorCopy(_m.data,10) << endl;
      exit(-1);
    }

    for (long k = 1; k < n+1; ++k) {
      cout << "-- " << k << " --------------------" << endl;
      mult(m[k], m[k-1], m[k-1]);
      //cout << "m[" << k << "] = " << VectorCopy(m[k].data,10) << endl;

      //square_ct_debug(ct[k], Hint, ct[k-1], s);
      square_ct(ct[k], Hint, ct[k-1]);
      if (dec_ng(prm, _m, s, ct[k], m[k])) {
	cout << "NG: power_test" << endl;
	cout << "  be: " << VectorCopy(m[k].data,10) << endl;
	cout << "  is: " << VectorCopy(_m.data,10) << endl;
	exit(-1);
      }
    }
  }

  cout << "OK: power_test" << endl;
}

void hom_test (param& prm, const hint& Hint, const SKey& s, long L) {
  cout << "--- hom_test for " << prm.label << " --- : level = " << L << endl;
  
  long test_num = 1;
  for (long cnt = 0; cnt < test_num; ++cnt) {
    cout << "[" << cnt << "]-----------------------" << endl;
    
    Vec<Vec<long>> op;
    op.SetLength(L+1);
    for (long j = 0; j < L + 1; ++j) op[j].SetLength(pow(2,L));

    // plain computation
    Vec<Vec<PlainText>> m;
    for (long j = 0; j < L + 1; ++j) {
      Vec<PlainText> _mm;
      for (long i = 0; i < pow(2,L); ++i) {
	PlainText _m(prm);
	append(_mm, _m);
      }
      append(m, _mm);
    }

    for (long i = 0; i < pow(2,L); ++i) m[0][i].uniform();

    for (long j = 1; j < L + 1; ++j) {
      for (long i = 0; i < pow(2,L-j); ++i) {
	op[j][i] = RandomBnd(2);
	if (op[j][i]) add(m[j][i], m[j-1][2*i], m[j-1][2*i+1]);
	else mult(m[j][i], m[j-1][2*i], m[j-1][2*i+1]);
      }
    }

    // homomorphic computation
    Vec<Vec<CipherText>> ct;
    for (long j = 0; j < L + 1; ++j) {
      Vec<CipherText> _cc;
      for (long i = 0; i < pow(2,L); ++i) {
	CipherText _ct(prm);
	append(_cc, _ct);
      }
      append(ct, _cc);
    }

    for (long i = 0; i < pow(2,L); ++i) encrypt(prm, ct[0][i], s, m[0][i]);
    
    for (long j = 1; j < L + 1; ++j) {
      cout << "-- level " << j << endl;
      for (long i = 0; i < pow(2,L-j); ++i) {
	if (op[j][i]) {
	  cout << "ct[" << j-1 << "][" << 2*i << "] + ct[" << j-1 << "][" << 2*i+1 << "]" << endl;
	  add_ct(ct[j][i], ct[j-1][2*i], ct[j-1][2*i+1]);
	}
	else {
	  cout << "ct[" << j-1 << "][" << 2*i << "] * ct[" << j-1 << "][" << 2*i+1 << "]" << endl;
	  mult_ct(ct[j][i], Hint, ct[j-1][2*i], ct[j-1][2*i+1]);
	}

	PlainText _m(prm);
	if (dec_ng(prm, _m, s, ct[j][i], m[j][i])) {
	  cout << "NG: hom_test" << endl;
	  cout << "  be: " << m[j][i] << endl;
	  cout << "  is: " << _m << endl;
	  exit(-1);
	}
      }
    }
  }
  
  cout << "OK: hom_test" << endl;
}

void she_test (param& prm) {
  cout << "--- she_test for " << prm << " ---" << endl;

  cout << "generating a secret key..." << flush;
  SKey s(prm);
  KeyGen(prm, s);
  cout << "done" << endl;

  encode_test(prm);
  dec_test(prm, s);
  
  div_by_2_test(prm, s);
  
  add_test(prm, s);
  plain_add_test(prm, s);

  if (prm.level > 1) {
    plain_mult_test(prm, s);
  
    hint_gen_test(prm, s);
    hint_gen_test2(prm, s);
    decomp_test(prm);

    hint_apply_test(prm, s.Hint, s.s2x, s.x);
  
    SKey t(prm);
    KeyGen(prm, t);

    hint Hint;
    Hint.hint_gen(prm, t.data, s.x);
    hint_apply_test(prm, Hint, t.x, s.x);
    key_switch_test(prm, s, t, Hint);

    direct_mult_test(prm, s);
    mult_test(prm, s.Hint, s);
    square_test(prm, s.Hint, s);

    //hom_test(prm, s.Hint, s, prm.level-1);
    power_test(prm, s.Hint, s, prm.level-1);
  }
}

void check_modulus_switch (param& prm, param& prm_tmp) {
  cout << "check_modulus_switch: " << flush;

  // the original secret key
  SKey s(prm);
  KeyGen(prm, s);

  // the corresponding secret key in prm_tmp
  SKey s_tmp(prm_tmp);
  KeyGen(prm_tmp, s_tmp, s.data);  // sx = s_tmp.x
  
  PlainText m(prm);
  CipherText ct(prm), dt(prm_tmp);

  m.uniform();
  encrypt(prm, ct, s, m);

  modulus_switch(dt, prm_tmp, prm_tmp.q, prm, ct);

  // check the result
  PlainText _m(prm_tmp);
  if (dec_ng(prm_tmp, _m, s_tmp, dt, m)) {
    cout << "NG: check_modulus_switch: " << endl;
    cout << "  be: " << VectorCopy(m.data,10) << endl;
    cout << "  is: " << VectorCopy(_m.data,10) << endl;
    exit(-1); 
  }
  cout << "-- OK" << endl;
}

int main (int argc, char *argv[]) {
  /* Do Test */
  //long prec = 2*left_param_big.r - left_param_big.l;
  //long prec = 2*right_param_big.r - right_param_big.l;
  long prec = 2*composed_param_big.r - composed_param_big.l;
  
  // left_param.init(prec);
  // she_test(left_param);

  // left_param_tmp.init(prec);
  // she_test(left_param_tmp);

  // left_param_big.init(prec);
  // she_test(left_param_big);
  
  // right_param.init(prec);
  // she_test(right_param);

  // right_param_tmp.init(prec);
  // she_test(right_param_tmp);

  // right_param_big.init(prec);
  // she_test(right_param_big);

  composed_param.init(prec);
  she_test(composed_param);

  composed_param_tmp.init(prec);
  she_test(composed_param_tmp);

  composed_param_big.init(prec);
  she_test(composed_param_big);

  check_modulus_switch(composed_param, composed_param_tmp);
}
