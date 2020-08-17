#ifndef SHE_H
#define SHE_H

#include <NTL/vec_ZZ_p.h>

#include "nsgen.h"
#include "ring.h"

using namespace std;
using namespace NTL;

class PlainText {
 public:
  param& prm;
  ZZ p;  // modulus
  long dim;   // dimension
  vec_ZZ data;  // plaintext data

 public:
  PlainText (param& _prm) : prm(_prm) {
    p = power(ZZ(prm.rg.pR), prm.l);
    dim = prm.rg.gR;
    data.SetLength(dim);
  }

  PlainText (param& _prm, vec_ZZ& _data) : prm(_prm) {
    p = power(ZZ(prm.rg.pR), prm.l);
    dim = prm.rg.gR;
    data.SetLength(dim);
    for (long i = 0; i < dim; ++i) data[i] = _data[i] % p;
  }

  PlainText& operator = (const PlainText& m) {
    if (this == &m) return *this;
    prm = m.prm;
    p = m.p;
    dim = m.dim;

    data.SetLength(dim);
    VectorCopy(data, m.data, dim);

    return *this;
  }
  
  ZZ& operator [] (long i) {return data[i];}

  void uniform () {
    for (long i = 0; i < dim; ++i) data[i] = RandomBnd(p);
  }
};

ostream& operator<< (ostream& s, const PlainText& x);
bool operator != (const PlainText& x, const PlainText& y);

void add (PlainText& x3, const PlainText& x1, const PlainText& x2);
void mult (PlainText& x3, const PlainText& x1, const PlainText& x2);
void deep_shift (PlainText& m2, string label, const long i, const PlainText& m1);

void encode_pt (param& prm, PlainText& z, const PlainText& x);
void decode_pt (param& prm, PlainText& x, const PlainText& z);

class hint {
 public:
  long len;  // length
  long dim;  // dimension
  Vec<Vec<vec_ZZ>> H;

 public:
  hint () {};

  void set_len_dim (long _len, long _dim) {
    len = _len;
    dim = _dim;
    H.SetLength(len);
    for (long i = 0; i < len; ++i) {
      H[i].SetLength(2);
      H[i][0].SetLength(dim);
      H[i][1].SetLength(dim);
    }  
  }

  void hint_gen (param& prm, const vec_ZZ& sin, const vec_ZZ& sout);
  
  hint& operator = (const hint& ht) {
    if (this == &ht) return *this;
    set_len_dim(ht.len, ht.dim);
    for (long i = 0; i < len; ++i) {
      for (long j = 0; j < dim; ++i) {
      	H[i][0][j] = ht.H[i][0][j];
	      H[i][1][j] = ht.H[i][1][j];
      }
    }
    return *this;
  }
};

class SKey {
 public:
  ZZ p;  // modulus
  long dim;   // dimension

  vec_ZZ data;  // key data in eta vector
  vec_ZZ x;  // key data in xi vector
  vec_ZZ s2;  // s^2 in eta vector
  vec_ZZ s2x;  // s^2 in eta vector

  hint Hint;
  
 public:
  SKey (param& prm) {
    p = power(ZZ(prm.rg.pR), prm.r + prm.s);
    dim = prm.rg.gR;
    data.SetLength(dim);
    x.SetLength(dim);
    s2.SetLength(dim);
    s2x.SetLength(dim);
  }  
};

void KeyGen (param& prm, SKey& s, const vec_ZZ& sdata);
void KeyGen_x (param& prm, SKey& s, const vec_ZZ& sx);
void KeyGen (param& prm, SKey& s);

class CipherText {
 public:
  param& prm;
  long level;

  long dim;   // dimension
  ZZ q;  // ciphertext modulus
  ZZ t;  // plaintext modulus

  vec_ZZ a;
  vec_ZZ b;

 public:
  CipherText (param& _prm) : prm(_prm) {
    level = prm.level;
    dim = prm.rg.gR;
    q = power(ZZ(prm.rg.pR), prm.r);
    t = power(ZZ(prm.rg.pR), prm.l);
    a.SetLength(dim);
    b.SetLength(dim);
  };

  // copy constructor
  CipherText (const CipherText& dt) : prm(dt.prm) {
    level = dt.level;
    dim = dt.dim;
    q = dt.q;
    t = dt.t;
    
    a.SetLength(dim);
    VectorCopy(a, dt.a, dim);
    b.SetLength(dim);
    VectorCopy(b, dt.b, dim);
  }

  // assignment
  CipherText& operator = (const CipherText& dt) {
    if (this == &dt) return *this;
    level = dt.level;
    prm = dt.prm;
    dim = dt.dim;
    q = dt.q;
    t = dt.t;
    
    a.SetLength(dim);
    VectorCopy(a, dt.a, dim);
    b.SetLength(dim);
    VectorCopy(b, dt.b, dim);

    return *this;
  }
};

class CipherText2 {
 public:
  param& prm;

  long dim;   // dimension
  ZZ q;  // ciphertext modulus
  ZZ t;  // plaintext modulus

  vec_ZZ a;
  vec_ZZ b;
  vec_ZZ c;

 public:
  CipherText2 (param& _prm) : prm(_prm) {
    dim = prm.rg.gR;
    q = power(ZZ(prm.rg.pR), prm.r);
    t = power(ZZ(prm.rg.pR), prm.l);

    a.SetLength(dim);
    b.SetLength(dim);
    c.SetLength(dim);
  };

  // copy constructor
  CipherText2 (const CipherText2& dt) : prm(dt.prm) {
    dim = dt.dim;
    q = dt.q;
    t = dt.t;

    a.SetLength(dim);
    VectorCopy(a, dt.a, dim);
    b.SetLength(dim);
    VectorCopy(b, dt.b, dim);
    c.SetLength(dim);
    VectorCopy(c, dt.c, dim);
  }
};

ostream& operator<< (ostream& s, const CipherText& ct);
void copy_ct (CipherText& dt, const CipherText& ct);
void zero_ct (CipherText& dt, const CipherText& ct);

/***********************************************************************/
void encrypt (param& prm, CipherText& ct, const SKey& s, const PlainText& m);
void dummy_encrypt (param& prm, CipherText& ct, const PlainText& m);
ZZ decrypt (param& prm, PlainText& m, const SKey& s, const CipherText& ct);
ZZ msd_decrypt (param& prm, PlainText& m, const SKey& s, const CipherText& ct);
ZZ decrypt2 (param& prm, PlainText& m, const SKey& s, CipherText2& ct);
void div_by_2 (CipherText& dt, const CipherText& ct);
void add_ct (CipherText& dt, CipherText& ct1, CipherText& ct2);
void sub_ct (CipherText& dt, CipherText& ct1, CipherText& ct2);
void double_ct (CipherText& dt, CipherText& ct);

/***********************************************************************/
void modulus_switch (CipherText& dt, param& new_prm, ZZ& q1, param& prm, CipherText& ct);

/***********************************************************************/
void hint_gen (param& prm, Vec<Vec<vec_ZZ>>& H, const vec_ZZ& sin, const vec_ZZ& sout);
void Decomp (param& prm, Vec<vec_ZZ>& d, const vec_ZZ& c);
void apply_hint (param& prm, vec_ZZ& a, vec_ZZ& b, const hint Hint,
		 const vec_ZZ& c, const ZZ& q);
void key_switch (CipherText& dt, const hint& Hint, const CipherText& ct);

/***********************************************************************/
void plain_add_ct (CipherText& dt, const PlainText& m, CipherText& ct1);
void plain_mult_ct (CipherText& dt, const PlainText& m, CipherText& ct1);
void direct_mult_ct (CipherText2& dt, CipherText& ct1, CipherText& ct2);
void direct_square_ct (CipherText2& dt, CipherText& ct);
void mult_ct (CipherText& dt, const hint& Hint, CipherText& ct1, CipherText& ct2);
void square_ct (CipherText& ct, const hint& Hint, CipherText& ct1);
void square_ct_debug (CipherText& ct, const hint& Hint, CipherText& ct1, const SKey& s);
  
#endif /*SHE_H*/
