#include "ring.h"
#include "she.h"

static long max_digits = 10;

ostream& operator<< (ostream& s, const PlainText& x) {
  s << VectorCopy(x.data, max_digits);
  return s;
}

ostream& operator<< (ostream& s, const CipherText& ct) {
  s << "ct[t = " << ct.t << ", q = " << ct.q 
    << ", a = " << VectorCopy(ct.a,max_digits) << ", b = " << VectorCopy(ct.b,max_digits) << "]";
  return s;
}

bool operator != (const PlainText& x, const PlainText& y) {return x.data != y.data;}

void add (PlainText& x3, const PlainText& x1, const PlainText& x2) {
  ZZ_p::init(x1.p);
  x3.data = conv<vec_ZZ>(conv<vec_ZZ_p>(x1.data) + conv<vec_ZZ_p>(x2.data));
}

void mult (PlainText& x3, const PlainText& x1, const PlainText& x2) {
  ZZ_p::init(x3.p);
  vec_ZZ_p _x3; _x3.SetLength(x3.dim);

  mult(_x3, conv<vec_ZZ_p>(x1.data), conv<vec_ZZ_p>(x2.data));
  x3.data = conv<vec_ZZ>(_x3);
}

void deep_shift (PlainText& m2, string label, const long i, const PlainText& m1) {
  ring& rg = m1.prm.rg;
  VectorCopy(m2.data, m1.data, m1.data.length());
  rg.deep_shift(label, m2.data, i);
}

void copy_ct (CipherText& dt, const CipherText& ct) {
  if (dt.prm.label != ct.prm.label) {
    cout << "Error: copy_ct: param mismatch" << endl;
    exit(-1);
  }
  dt.dim = ct.dim;
  dt.q = ct.q;
  dt.t = ct.t;
  VectorCopy(dt.a, ct.a, ct.dim);
  VectorCopy(dt.b, ct.b, ct.dim);
}

void zero_ct (CipherText& dt, const CipherText& ct) {
  dt.prm = ct.prm;
  dt.dim = ct.dim;
  dt.q = ct.q;
  dt.t = ct.t;
  clear(dt.a);
  clear(dt.b);
}

/***********************************************************************/
void KeyGen (param& prm, SKey& s, const vec_ZZ& sdata) {
  ZZ_p::init(s.p);
  
  VectorCopy(s.data, sdata, sdata.length());
  
  // convert 'sdata' into the gamma-vector
  vec_ZZ_p _s = conv<vec_ZZ_p>(s.data);
  vec_ZZ_p _sx;

  prm.applyOmega(_sx, _s);
  s.x = conv<vec_ZZ>(_sx);

  // compute 's^2' 
  vec_ZZ_p _s2x, _s2;
  _s2x.SetLength(s.data.length()); 

  mult(_s2x, _sx, _sx);
  s.s2x = conv<vec_ZZ>(_s2x);

  prm.applyOmegaInv(_s2, _s2x);
  center_lift(s.s2, _s2);
  cout << "s2 = " << VectorCopy(s.s2, 10) << endl;

  // hint_gen for ciphertext multiplication
  s.Hint.hint_gen(prm, s.s2, s.x);
}

void KeyGen_x (param& prm, SKey& s, const vec_ZZ& sx) {
  ZZ_p::init(s.p);
  
  VectorCopy(s.x, sx, sx.length());
  
  // convert 'sx' into the eta-vector
  vec_ZZ_p _sx = conv<vec_ZZ_p>(s.x);
  vec_ZZ_p _s;

  prm.applyOmegaInv(_s, _sx);
  center_lift(s.data, _s);
  cout << "s = " << VectorCopy(s.data, 10) << endl;

  // compute 's^2' 
  vec_ZZ_p _s2x, _s2;
  _s2x.SetLength(s.data.length()); 

  mult(_s2x, _sx, _sx);
  s.s2x = conv<vec_ZZ>(_s2x);

  //prm.applyGammaInv(_s2, _s2x);
  prm.applyOmegaInv(_s2, _s2x);
  center_lift(s.s2, _s2);
  cout << "s2 = " << VectorCopy(s.s2, 10) << endl;

  // hint_gen for ciphertext multiplication
  s.Hint.hint_gen(prm, s.s2, s.x);
}

void KeyGen (param& prm, SKey& s) {
  // generate a secret key 's'
  ZZ w;
  vec_ZZ sdata; 

  do {
    prm.ns.hamming(sdata, prm.rg.gR);
    w = 0;
    for (long i = 0; i < sdata.length(); ++i) w += sdata[i];
  } while ((w > 3) || (w < -3));

  KeyGen(prm, s, sdata);  
}

/***********************************************************************/
void encode_pt (param& prm, vec_ZZ_p& z, const vec_ZZ_p& x) {
  prm.applyOmegaInv(z, x);
}

void encode_pt (param& prm, PlainText& z, const PlainText& x) {
  ZZ_pContext context;
  context.save();  
  ZZ_p::init(x.p);
  
  vec_ZZ_p _z;
  encode_pt(prm, _z, conv<vec_ZZ_p>(x.data));
  z.data = conv<vec_ZZ>(_z);

  context.restore();
}

void decode_pt (param& prm, vec_ZZ_p& x, const vec_ZZ_p& z) {
  prm.applyOmega(x, z);
}

void decode_pt (param& prm, PlainText& x, const PlainText& z) {
  ZZ_pContext context;
  context.save();  
  ZZ_p::init(z.p);

  vec_ZZ_p _x; 
  decode_pt(prm, _x, conv<vec_ZZ_p>(z.data));
  x.data = conv<vec_ZZ>(_x);

  context.restore();
}

/***********************************************************************/
// returns (a, b) s.t. a + bs = c + t e with some noise e
void RLWE (vec_ZZ& a, vec_ZZ& b, ring& rg, noise& ns, const vec_ZZ& sx, const vec_ZZ& c,
	   const ZZ& t) {
  // t e + c
  vec_ZZ e;
  ns.gaussian(e, rg.gR);

  vec_ZZ_p w;
  w.SetLength(rg.gR);
  
  mul(w, conv<ZZ_p>(t), conv<vec_ZZ_p>(e));
  add(w, conv<vec_ZZ_p>(c), w);
  //rg.applyGamma(w, w);
  rg.applyOmega(w, w);
  
  vec_ZZ_p _a, _b;
  _a.SetLength(rg.gR);

  uniform(_b, rg.gR);
  mult(_a, _b, conv<vec_ZZ_p>(sx));
  sub(_a, w, _a);  // a = (te + c) - bs

  a = conv<vec_ZZ>(_a);
  b = conv<vec_ZZ>(_b);
}

void RLWE (vec_ZZ& a, vec_ZZ& b, param& prm, const vec_ZZ& sx, const vec_ZZ& c, const ZZ& t) {
  RLWE(a, b, prm.rg, prm.ns, sx, c, t);
}

void trivial_RLWE (vec_ZZ& a, vec_ZZ& b, param& prm, const vec_ZZ& c) {
  vec_ZZ_p _a;
  //prm.applyGammaa(_a, conv<vec_ZZ_p>(c));
  prm.applyOmega(_a, conv<vec_ZZ_p>(c));
  a = conv<vec_ZZ>(_a);
  clear(b);
}

void invert_RLWE (vec_ZZ& e, param& prm, const vec_ZZ& sx, const vec_ZZ& a, const vec_ZZ& b) {
  vec_ZZ_p w;
  w.SetLength(prm.rg.gR);

  // w = a + b s mod q
  mult(w, conv<vec_ZZ_p>(b), conv<vec_ZZ_p>(sx));
  add(w, conv<vec_ZZ_p>(a), w);

  //prm.applyGammaInv(w, w);
  prm.applyOmegaInv(w, w);
  e = conv<vec_ZZ>(w);
}

void invert_RLWE2 (vec_ZZ& e, param& prm, const SKey& s,
		   const vec_ZZ& a, const vec_ZZ& b, const vec_ZZ& c) {
  vec_ZZ_p w, v;
  w.SetLength(prm.rg.gR);
  v.SetLength(prm.rg.gR);

  // w = a + b s + c s^2mod q
  mult(w, conv<vec_ZZ_p>(b), conv<vec_ZZ_p>(s.x));
  add(w, conv<vec_ZZ_p>(a), w);
  mult(v, conv<vec_ZZ_p>(c), conv<vec_ZZ_p>(s.s2x));
  add(w, v, w);

  //prm.applyGammaInv(w, w);
  prm.applyOmegaInv(w, w);
  e = conv<vec_ZZ>(w);
}

/***********************************************************************/
void encrypt (param& prm, CipherText& ct, const SKey& s, const PlainText& m) {
  ZZ_p::init(ct.q);
  
  PlainText x(prm);
  encode_pt(prm, x, m);

  ZZ delta = ct.q/ct.t;
  vec_ZZ_p w = conv<ZZ_p>(delta) * conv<vec_ZZ_p>(x.data);
  RLWE(ct.a, ct.b, prm, s.x, conv<vec_ZZ>(w), ZZ(1));
}

void dummy_encrypt (param& prm, CipherText& ct, const PlainText& m) {
  ZZ_p::init(ct.q);
  
  PlainText x(prm);
  encode_pt(prm, x, m);
  
  ZZ delta = ct.q/ct.t;
  vec_ZZ_p w = conv<ZZ_p>(delta) * conv<vec_ZZ_p>(x.data);
  trivial_RLWE(ct.a, ct.b, prm, conv<vec_ZZ>(w));
}

ZZ decrypt (param& prm, PlainText& m, const SKey& s, const CipherText& ct) {
  ZZ_p::init(ct.q);
  vec_ZZ w0, w;
  invert_RLWE(w0, prm, s.x, ct.a, ct.b);
  rescale(w, w0, ct.q, ct.t);

  PlainText x(prm, w);
  //x.p = ct.t;
  decode_pt(prm, m, x);
 
  vec_ZZ_p e;
  e.SetLength(ct.dim);

  ZZ delta = ct.q/ct.t;
  mul(e, conv<ZZ_p>(delta), conv<vec_ZZ_p>(w));
  sub(e, conv<vec_ZZ_p>(w0), e);

  return inf_norm(center_lift(e));
}

ZZ decrypt2 (param& prm, PlainText& m, const SKey& s, CipherText2& ct) {
  ZZ_p::init(ct.q);
  
  vec_ZZ w0, w;
  w.SetLength(ct.dim);
  invert_RLWE2(w0, prm, s, ct.a, ct.b, ct.c);
  rescale(w, w0, ct.q, ct.t);

  PlainText x(prm, w);
  //x.p = ct.t;
  decode_pt(prm, m, x);

  vec_ZZ_p e;
  e.SetLength(ct.dim);

  ZZ delta = ct.q/ct.t;
  mul(e, conv<ZZ_p>(delta), conv<vec_ZZ_p>(w));
  sub(e, conv<vec_ZZ_p>(w0), e);

  return inf_norm(center_lift(e));
}


/***********************************************************************/
// assumes the underlying plaintext of 'ct' is divisible by 2
void div_by_2 (CipherText& dt, const CipherText& ct) {
  VectorCopy(dt.a, ct.a, ct.dim);
  VectorCopy(dt.b, ct.b, ct.dim);
  dt.q = ct.q;
  dt.t = ct.t/2;
}

void add_ct (CipherText& dt, CipherText& ct1, CipherText& ct2) {
  dt.q = ct1.q;
  ZZ_p::init(dt.q);
  
  dt.a = conv<vec_ZZ>(conv<vec_ZZ_p>(ct1.a) + conv<vec_ZZ_p>(ct2.a));
  dt.b = conv<vec_ZZ>(conv<vec_ZZ_p>(ct1.b) + conv<vec_ZZ_p>(ct2.b));
}

void sub_ct (CipherText& dt, CipherText& ct1, CipherText& ct2) {
  dt.q = ct1.q;
  ZZ_p::init(dt.q);
  
  dt.a = conv<vec_ZZ>(conv<vec_ZZ_p>(ct1.a) - conv<vec_ZZ_p>(ct2.a));
  dt.b = conv<vec_ZZ>(conv<vec_ZZ_p>(ct1.b) - conv<vec_ZZ_p>(ct2.b));  
}

void double_ct (CipherText& dt, CipherText& ct) {
  dt.dim = ct.dim;
  dt.q = ct.q;
  dt.t = ct.t;
  
  for (long i = 0; i < ct.dim; i++) {
    dt.a[i] = ct.a[i] * 2;
    if (dt.a[i] >= dt.q) dt.a[i] -= dt.q;
    dt.b[i] = ct.b[i] * 2;
    if (dt.b[i] >= dt.q) dt.b[i] -= dt.q;    
  }
}

/***********************************************************************/
void modulus_switch (CipherText& dt, param& new_prm, ZZ& q1, param& prm, CipherText& ct) {
  ZZ_pContext context;
  context.save();  
  
  ZZ& q = ct.q;
  ZZ_p::init(q);

  vec_ZZ_p wa, wb; 
  prm.applyOmegaInv(wa, conv<vec_ZZ_p>(ct.a));
  prm.applyOmegaInv(wb, conv<vec_ZZ_p>(ct.b));

  vec_ZZ va, vb; 
  rescale(va, conv<vec_ZZ>(wa), q, q1);
  rescale(vb, conv<vec_ZZ>(wb), q, q1);
  
  ZZ_p::init(q1);
  vec_ZZ_p _a, _b;
  new_prm.applyOmega(_a, conv<vec_ZZ_p>(va));
  new_prm.applyOmega(_b, conv<vec_ZZ_p>(vb));

  dt.a = conv<vec_ZZ>(_a);
  dt.b = conv<vec_ZZ>(_b);

  context.restore();
}

/***********************************************************************/
// a key-switching hint for sin under sout
// -- H(sout) = (1, sout) H  ~ sin g^T
// -- sin: eta-vector (typically, s^2)
// -- sout: sigma-vector (typically, s)

void hint::hint_gen (param& prm, const vec_ZZ& sin, const vec_ZZ& sout) {
  ZZ_pContext context;
  context.save();

  ZZ t = prm.t;
  ZZ sR = prm.sR;
  ZZ qqR = prm.qq;
  ZZ wR = prm.w;

  ZZ_p::init(qqR);

  len = prm.lw;
  dim = prm.rg.gR;
  set_len_dim(len, dim);

  ZZ w(1);
  for (long j = 0; j < len; ++j) {
    // v = w^j sin
    vec_ZZ_p v =  conv<ZZ_p>(w) * conv<vec_ZZ_p>(sin);
    v = conv<ZZ_p>(sR) * v;

    RLWE(H[j][0], H[j][1], prm.rg, prm.ns, sout, conv<vec_ZZ>(v), t);
    w *= wR;
  }

  context.restore();
}

void Decomp (param& prm, Vec<vec_ZZ>& d, const vec_ZZ& c) {
  ZZ& w = prm.w;
  vec_ZZ cc(c);

  for (long l = 0; l < prm.lw; ++l) {
    for (long i = 0; i < cc.length(); ++i) {
      d[l][i] = cc[i] % w;
      cc[i] = cc[i] / w;
    }
  } 
}

// Given an encryption H of sin and an element c,
//  returns an `encryption' of c sin
// d = H g^(-1)(c),
// -- H(sout) = (1, sout) H  ~ sin g^T : an encryption of sin
// -- c: eta-vector
// -- d0: delta-vector, d1: eta-vector s.t. d(sout) ~ sin c
void apply_hint (param& prm, vec_ZZ& a, vec_ZZ& b, const hint Hint,
		 const vec_ZZ& c, const ZZ& q) {
  ring& rg = prm.rg;
  ZZ sR = prm.sR;
  ZZ q0 = q * sR;

  ZZ_pContext context;
  context.save();  
  ZZ_p::init(q0);
  
  Vec<vec_ZZ> cc;
  cc.SetLength(prm.lw);
  for (long l = 0; l < cc.length(); ++l) cc[l].SetLength(rg.gR);
  
  Decomp(prm, cc, c);
    
  vec_ZZ_p w, _cc, _a, _b;
  w.SetLength(rg.gR);
  _a.SetLength(rg.gR);
  _b.SetLength(rg.gR);  

  for (long j = 0; j < prm.lw; ++j) {
    prm.applyOmega(_cc, conv<vec_ZZ_p>(cc[j]));

    // a = a + H[j].a cc[j]
    mult(w, _cc, conv<vec_ZZ_p>(Hint.H[j][0]));
    _a = w + _a;

    // b = b + H[j].b cc[j]
    mult(w, _cc, conv<vec_ZZ_p>(Hint.H[j][1]));
    _b = w + _b;
  }

  prm.applyOmegaInv(_a, _a);
  prm.applyOmegaInv(_b, _b);

  rescale(a, conv<vec_ZZ>(_a), q0, q);
  rescale(b, conv<vec_ZZ>(_b), q0, q);

  context.restore();
}

void key_switch (CipherText& dt, const hint& Hint, const CipherText& ct) {
  dt.q = ct.q; 
  ZZ_p::init(dt.q);

  vec_ZZ aa, bb;
  vec_ZZ b;

  dt.prm.applyOmegaInv(b, ct.b);
  apply_hint(dt.prm, aa, bb, Hint, b, dt.q);  // aa + bb s ~ b t

  dt.prm.applyOmega(aa, aa);
  dt.prm.applyOmega(bb, bb);

  dt.a = conv<vec_ZZ>(conv<vec_ZZ_p>(aa) + conv<vec_ZZ_p>(ct.a));
  dt.b = bb;
}

void linearlize (CipherText& ct, const hint& Hint, CipherText2& dt) {
  param& prm = dt.prm;
  ZZ& q = dt.q;
  ZZ_p::init(q);
  
  vec_ZZ_p _a, _b, _c;
  prm.applyOmegaInv(_c, conv<vec_ZZ_p>(dt.c));

  vec_ZZ aa, bb;
  apply_hint(prm, aa, bb, Hint, conv<vec_ZZ>(_c), q);

  ct.q = dt.q;
  ZZ_p::init(q);
  
  prm.applyOmega(_a, conv<vec_ZZ_p>(aa));
  prm.applyOmega(_b, conv<vec_ZZ_p>(bb));

  _a += conv<vec_ZZ_p>(dt.a);
  _b += conv<vec_ZZ_p>(dt.b);

  ct.a = conv<vec_ZZ>(_a);
  ct.b = conv<vec_ZZ>(_b);
}

/***********************************************************************/
void plain_add_ct (CipherText& dt, const PlainText& m, CipherText& ct1) {
  param& prm = ct1.prm;
  ZZ q = ct1.q;
  ZZ_p::init(q);
  
  dt.q = q;

  CipherText ct2(prm);
  dummy_encrypt(prm, ct2, m);  // ct2.b = 0

  dt.a = conv<vec_ZZ>(conv<vec_ZZ_p>(ct1.a) + conv<vec_ZZ_p>(ct2.a));
  dt.b = ct1.b;
}

// requires finer precision
void plain_mult_ct (CipherText& dt, const PlainText& m, CipherText& ct) {
  param& prm = ct.prm;
  ZZ q = ct.q;
  ZZ_p::init(q);

  vec_ZZ a1, b1, a2;
  a1.SetLength(ct.dim);
  b1.SetLength(ct.dim);
  a2.SetLength(ct.dim);
  prm.applyOmegaInv(a1, ct.a);
  prm.applyOmegaInv(b1, ct.b);

  ZZ qq = power(ZZ(prm.rg.pR), 2*prm.r - prm.l);
  ZZ_p::init(qq);

  prm.applyOmega(a1, a1);
  prm.applyOmega(b1, b1);

  ZZ_p::init(q);
  CipherText ct2(prm);
  dummy_encrypt(prm, ct2, m);  // ct2.b = 0
  prm.applyOmegaInv(a2, ct2.a);

  ZZ_p::init(qq);
  prm.applyOmega(a2, a2);

  // a = a2 a1 / delta, 
  vec_ZZ_p w;
  w.SetLength(ct.dim);
  mult(w, conv<vec_ZZ_p>(a1), conv<vec_ZZ_p>(a2));
  prm.applyOmegaInv(w, w);
  rescale(dt.a, conv<vec_ZZ>(w), qq, q);

  ZZ_p::init(q);
  prm.applyOmega(dt.a, dt.a);

  // b = a2 b1 / delta
  ZZ_p::init(qq);
  
  vec_ZZ_p v;
  v.SetLength(ct.dim);
  mult(v, conv<vec_ZZ_p>(b1), conv<vec_ZZ_p>(a2));
  prm.applyOmegaInv(v, v);
  rescale(dt.b, conv<vec_ZZ>(v), qq, q);

  ZZ_p::init(q);
  prm.applyOmega(dt.b, dt.b);
}

/***********************************************************************/
void direct_mult_ct (CipherText2& dt, CipherText& ct1, CipherText& ct2) {
  param& prm = ct1.prm;
  ZZ q = dt.q = ct1.q;
  ZZ t = dt.t = ct1.t;
  ZZ qq = power(ZZ(prm.rg.pR), 2*prm.r - prm.l);
  long dim = ct1.dim;

  // lift the precision of xi-vectors
  vec_ZZ a1, b1, a2, b2;
  a1.SetLength(dim);
  b1.SetLength(dim);
  a2.SetLength(dim);
  b2.SetLength(dim);

  ZZ_p::init(q);
  prm.applyOmegaInv(a1, ct1.a);
  prm.applyOmegaInv(b1, ct1.b);
  prm.applyOmegaInv(a2, ct2.a);
  prm.applyOmegaInv(b2, ct2.b);

  ZZ_p::init(qq);
  prm.applyOmega(a1, a1);
  prm.applyOmega(b1, b1);
  prm.applyOmega(a2, a2);
  prm.applyOmega(b2, b2);
  
  vec_ZZ_p v,w;
  v.SetLength(dim); w.SetLength(dim);

  // a = a1 a2 / delta
  mult(w, conv<vec_ZZ_p>(a1), conv<vec_ZZ_p>(a2));
  prm.applyOmegaInv(w, w);
  rescale(dt.a, conv<vec_ZZ>(w), qq, q);
  ZZ_p::init(q);
  prm.applyOmega(dt.a, dt.a);

  // b = (a1 b2 + a2 b1) / delta
  ZZ_p::init(qq);
  mult(v, conv<vec_ZZ_p>(a1), conv<vec_ZZ_p>(b2));
  mult(w, conv<vec_ZZ_p>(a2), conv<vec_ZZ_p>(b1));
  w = v + w;
  prm.applyOmegaInv(w, w);
  rescale(dt.b, conv<vec_ZZ>(w), qq, q);
  ZZ_p::init(q);
  prm.applyOmega(dt.b, dt.b);

  // c = b1 b2 / delta
  ZZ_p::init(qq);
  mult(w, conv<vec_ZZ_p>(b1), conv<vec_ZZ_p>(b2));
  prm.applyOmegaInv(w, w);
  rescale(dt.c, conv<vec_ZZ>(w), qq, q);
  ZZ_p::init(q);
  prm.applyOmega(dt.c, dt.c);
}

void direct_square_ct (CipherText2& dt, CipherText& ct) {
  param& prm = ct.prm;
  ZZ q = dt.q = ct.q;
  ZZ t = dt.t = ct.t;
  ZZ qq = power(ZZ(prm.rg.pR), 2*prm.r - prm.l);
  long dim = ct.dim;

  // lift the precision of xi-vectors
  vec_ZZ a, b;
  a.SetLength(dim);
  b.SetLength(dim);

  ZZ_p::init(q);
  prm.applyOmegaInv(a, ct.a);
  prm.applyOmegaInv(b, ct.b);

  ZZ_p::init(qq);
  prm.applyOmega(a, a);
  prm.applyOmega(b, b);
  
  vec_ZZ_p v,w;
  v.SetLength(dim); w.SetLength(dim);

  // a = a^2 / delta
  sqr(w, conv<vec_ZZ_p>(a));
  prm.applyOmegaInv(w, w);
  rescale(dt.a, conv<vec_ZZ>(w), qq, q);
  ZZ_p::init(q);
  prm.applyOmega(dt.a, dt.a);

  // b = 2 a b / delta
  ZZ_p::init(qq);
  mult(w, conv<vec_ZZ_p>(a), conv<vec_ZZ_p>(b));
  w = w + w;
  prm.applyOmegaInv(w, w);
  rescale(dt.b, conv<vec_ZZ>(w), qq, q);
  ZZ_p::init(q);
  prm.applyOmega(dt.b, dt.b);

  // c = b^2 / delta
  ZZ_p::init(qq);
  sqr(w, conv<vec_ZZ_p>(b));
  prm.applyOmegaInv(w, w);
  rescale(dt.c, conv<vec_ZZ>(w), qq, q);
  ZZ_p::init(q);
  prm.applyOmega(dt.c, dt.c);  
}

void mult_ct (CipherText& ct, const hint& Hint, CipherText& ct1, CipherText& ct2) {
  CipherText2 dt(ct1.prm);
  direct_mult_ct(dt, ct1, ct2);

  linearlize(ct, Hint, dt);
  ct.level = min(ct1.level, ct2.level) - 1;
}

void square_ct (CipherText& ct, const hint& Hint, CipherText& ct1) {
  CipherText2 dt(ct1.prm);
  
  direct_square_ct(dt, ct1); 
  linearlize(ct, Hint, dt);
  ct.level = ct1.level - 1;
}

void square_ct_debug (CipherText& ct, const hint& Hint, CipherText& ct1, const SKey& s) {
  CipherText2 dt(ct1.prm);
  
  direct_square_ct(dt, ct1);
  {
    CipherText2 dt2(dt);
    PlainText _m(dt2.prm);
    ZZ v = decrypt2(dt2.prm, _m, s, dt2);
    cout << "-- direct: log2(noise) = " << log2(conv<double>(v)) << endl;
  }    
  
  linearlize(ct, Hint, dt);
  {
    CipherText ct2(ct);
    PlainText _m(ct2.prm);
    ZZ v = decrypt(ct2.prm, _m, s, ct2);
    cout << "-- linearlize: log2(noise) = " << log2(conv<double>(v)) << endl;
  }

  ct.level = ct1.level - 1;
}

