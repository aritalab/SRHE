#ifndef RING_H
#define RING_H

#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <NTL/ZZ_pXFactoring.h>
#include "nsgen.h"

using namespace std;
using namespace NTL;

enum ring_type {
  prime = 0,
  composite = 1,
};

enum text_type {
  plain = 0,
  cipher = 1,
};

class ring {
 public:
  bool  initialized;
  
  string label;
  ring_type type;
  long mR;
  long pR;
  long dR;
  long gR;
  long tR;

  // resolution of one
  ZZX res; 

  // matrices regarding with canonical embeddings
  // used only in prime case
  vec_ZZ omega;
  vec_ZZ omega_inv;
  
  long depth;  // the number of prime rings
  Vec<ring> subrings;

  ring* left;
  ring* right;  

 public:
  ring (string l, long m, long p, long d, long g, long t) {
    initialized = false;
    label = l; mR = m; pR = p; dR = d; gR = g; tR = t;
    type = prime; // mR must be prime
    left = 0;
    right = 0;
    }

  ring (string l, ring* lt, ring* rt) {
    initialized = false;
    label = l;
    mR = lt->mR * rt->mR;
    pR = lt->pR;
    dR = -1;
    gR = lt->gR * rt->gR;
    tR = -1;
    type = composite;
    left = lt;
    right = rt;
  }

  ring& operator = (const ring& rg) {
    initialized = rg.initialized;
    
    label = rg.label;
    type = rg.type;
    mR = rg.mR;
    pR = rg.pR;
    dR = rg.dR;
    gR = rg.gR;
    tR = rg.tR;

    res = rg.res;
       
    omega = rg.omega;
    omega_inv = rg.omega_inv;
    
    left = rg.left;
    right = rg.right;
    
    return *this;
  }

  long GZ (long i) { // i be in range(dR)
    return PowerMod(pR, i, mR);
  }

  long GofZ (long i) { // i be in range(gR)
    return PowerMod(tR, i, mR);
  }

  void to_eta (ZZ_pX& w, vec_ZZ_p& a);
  void to_delta (ZZ_pX& w, vec_ZZ_p& a);
  ZZX eta (long j);
  ZZ_pX delta (long j);
  ZZX apply_conj (ZZX& a, long i);
  ZZX trace (ZZX& w);
  ZZ scalar_quo (ZZX& a, ZZX& b, ZZ& p);
  
  void resolution_of_one (ZZX& res, const long m, const ZZ& p);
  void lift_resolution (ZZX& res, const long m, const ZZ& r, const long l);
  
  void canonical_resolution_of_one (ZZX& res, long prec);

  void init (long prec, noise& ns);
  void init_resolution (long prec);

  void init_omega (ZZ q);
  void init_E (noise& ns);

  void comp_omega (vec_ZZ& o, ZZ& q);
  void comp_omega_inv (vec_ZZ& o, const vec_ZZ& _omega, const ZZ& q);

  void comp_gamma (vec_ZZ& o, const vec_ZZ& _omega, const ZZ& q);
  void comp_gamma_inv (vec_ZZ& o, const vec_ZZ& _omega, const ZZ& q);
  
  void applyOmega (vec_ZZ_p& b, const vec_ZZ_p& a);
  void applyOmega (vec_ZZ& b, const vec_ZZ& a);
  void applyOmegaInv (vec_ZZ_p& b, const vec_ZZ_p& a);
  void applyOmegaInv (vec_ZZ& b, const vec_ZZ& a);

  long get_subrings (Vec<ring>& subrings);
  long get_subring_index (string label);

  void deep_shift (string _label, vec_ZZ& v, long i);
  void deep_shift_right (string _label, vec_ZZ& v, long i);

  vec_ZZ index_convert (long j);

  void embed_0 (vec_ZZ& y, vec_ZZ& x);
  vec_ZZ embed_0 (vec_ZZ& x);

  void embed_1 (vec_ZZ& y, vec_ZZ& x);
  vec_ZZ embed_1 (vec_ZZ& x);

  void embed_0 (vec_ZZ& y, string subring_label, vec_ZZ& x);
  vec_ZZ embed_0 (string subring_label, vec_ZZ& x);
  void embed_1 (vec_ZZ& y, string subring_label, vec_ZZ& x);
  vec_ZZ embed_1 (string subring_label, vec_ZZ& x);
  void embed_11 (vec_ZZ& y, string subring_label, vec_ZZ& x);
  vec_ZZ embed_11 (string subring_label, vec_ZZ& x);

  void trace (vec_ZZ& y, vec_ZZ& x, ZZ& q);
  void trace (vec_ZZ& y, string subring_label, vec_ZZ& x, ZZ& q);

  void deep_omega (vec_ZZ& extended_omega, string _label, long i);
  void deep_omega_inv (vec_ZZ& extended_omega_inv, string _label, long i);
};

ZZ compose_integer (const ZZ& r1, const ZZ& q1, const ZZ& r2, const ZZ& q2);

void mult (vec_ZZ_p& x, const vec_ZZ_p& a, const vec_ZZ_p& b);
void sqr (vec_ZZ_p& x, const vec_ZZ_p& a);
void add_scalar (vec_ZZ_p& x, const vec_ZZ_p& a, const ZZ_p& b);
void add_scalar (vec_ZZ_p& x, const ZZ_p& a, const vec_ZZ_p& b);
void zzx_convolution (vec_ZZ_p& c, const vec_ZZ_p& a, const vec_ZZ_p& b);

ZZ inf_norm (const vec_ZZ& a);
void center_lift (vec_ZZ& b, const vec_ZZ_p& a);
void center_lift (vec_ZZ& b, const vec_ZZ& a, ZZ& p);
vec_ZZ center_lift (const vec_ZZ_p& a);
vec_ZZ center_lift (const vec_ZZ& a, ZZ& p);
void rescale (vec_ZZ &b, const vec_ZZ &a, const ZZ& q, const ZZ& q1);
void lsd_rescale (vec_ZZ &b, const vec_ZZ &a, const ZZ& P, const ZZ& t);
  
class param {
 public:
  string label;
  ring& rg;
  long level;
  long l;  // plaintext modulus = pR^l
  long r;  // ciphertext modulus = pR^r
  long s;  // special modulus = pR^s
  long lw;  // word length
  ZZ t;
  ZZ qq;
  ZZ q;
  ZZ sR;
  ZZ w;  
  noise ns;
  
 public:
  param (string _label, ring& _rg, long _level, long _l, long _r, long _s, long _lw, noise& _ns)
    : label(_label), rg(_rg), level(_level), l(_l), r(_r), s(_s), lw(_lw), ns(_ns) {
    t = power(ZZ(rg.pR), l);
    qq = power(ZZ(rg.pR), r + s);
    q = power(ZZ(rg.pR), r);
    sR = power(ZZ(rg.pR), s);
    w = power(ZZ(rg.pR), ceil(((double)r + (double)s)/((double)lw)));  /*@@*/
  }

  param& operator = (const param& prm) {
    if (this == &prm) return *this;
    label = prm.label;
    rg = prm.rg;
    level = prm.level;
    l = prm.l;
    r = prm.r;
    s = prm.s;
    lw = prm.lw;
    t = prm.t;
    qq = prm.qq;
    q = prm.q;
    sR = prm.sR;    
    w = prm.w;    
    ns = prm.ns;
    return *this;
  }
  
  void init (long prec);

  long get_subrings (Vec<ring>& subrings);

  void applyOmega (vec_ZZ_p& b, const vec_ZZ_p& a) {
    rg.applyOmega(b, a);
  }
  void applyOmega (vec_ZZ& b, const vec_ZZ& a) {
    rg.applyOmega(b, a);
  }
  void applyOmegaInv (vec_ZZ_p& b, const vec_ZZ_p& a) {
    rg.applyOmegaInv(b, a);
  }
  void applyOmegaInv (vec_ZZ& b, const vec_ZZ& a) {
    rg.applyOmegaInv(b, a);
  }

  void gaussian (vec_ZZ& e) {    
    ns.gaussian(e, rg.gR);
  }
  
};

ostream& operator<< (ostream& s, const ring& rg);
ostream& operator<< (ostream& s, const param& prm);

void mat (Vec<vec_ZZ_p>& aa, const vec_ZZ_p& a, const long n, const long m);
void mat (Vec<vec_ZZ>& aa, const vec_ZZ& a, const long n, const long m);
void transpose (Vec<vec_ZZ_p>& b, const Vec<vec_ZZ_p>& a);
void transpose (Vec<vec_ZZ>& b, const Vec<vec_ZZ>& a);
void vec (vec_ZZ_p& b, const Vec<vec_ZZ_p>& a);
void transpose_vec (vec_ZZ_p& b, const Vec<vec_ZZ_p>& a);

void applyMatrix (vec_ZZ_p& b, const vec_ZZ_p& m, const vec_ZZ_p& a);
void applyMatrix (vec_ZZ_p& b, const vec_ZZ_p& m_l, const vec_ZZ_p& m_r, const vec_ZZ_p& a);

/***********************************************************************/
/* returns the power decomp of x mod 2^n  */
void power_decomp (Vec<vec_ZZ>& z, long n, const vec_ZZ& x);
void right_shift (vec_ZZ& z, const vec_ZZ& x, long m, long n);
void right_shift (vec_ZZ_p& z, const vec_ZZ_p& x, long m, long n);

/***********************************************************************/
/* returns the tensor product of a and a1 */
void diagonal (vec_ZZ& MD, const long k, const vec_ZZ& M);
void diagonal (vec_ZZ_p& MD, const long k, const vec_ZZ_p& M);
void tensor_prod (vec_ZZ& A, const vec_ZZ& a, const vec_ZZ& a1, const ZZ& q);
void shift (vec_ZZ& b, const  long i, const vec_ZZ& a);
void shift (vec_ZZ_p& b, const  long i, const vec_ZZ_p& a);
void shift (vec_ZZ_p& a);
void shift_right (vec_ZZ& b, const long i, const vec_ZZ& a);
void shift_right (vec_ZZ_p& b, const long i, const vec_ZZ_p& a);
void right_shift (vec_ZZ_p& a);

/***********************************************************************/
void init_qR_list (long *qR_list, long m, long min_2_exp, long min_prime_size, long chain_length);

#endif /*RING_H*/
