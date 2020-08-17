#include <NTL/ZZ.h>
#include <NTL/ZZ_pXFactoring.h>
#include "nsgen.h"
#include "ring.h"

void print_poly (ZZX& f) {
  long d = deg(f);
  if (d == -1) {
    cout << "0" << endl;
  }
  else {
    for (long i = d; i >= 0; --i) {
      if (f[i] != 0)
	cout << f[i] << "*" << "X^" << i << " + ";
    }
    cout << endl;
  }
}

/***********************************************************************/
ZZX ring::apply_conj (ZZX& a, long i) {
  long d = deg(a);
  if (d == -1) return ZZX(0);

  ZZX phi = (ZZX(INIT_MONO, mR)-1)/(ZZX(INIT_MONO, 1)-1);
  ZZX w = ZZX(0);
  
  for (long j = 0; j <= d; j++) {
    if (a[j] != 0) {
      long t = (j*i) % mR;
      ZZX z = ZZX(INIT_MONO, t);
      rem(z, z, phi);
      z *= ZZX(a[j]);
      w += z;
    }
  }
  
  return w;
}

ZZX ring::trace (ZZX& w) {
  ZZX r = ZZX(0);

  for (long a = 0; a < dR; ++a) {
    ZZX v = apply_conj(w, GZ(a));
    r += v;
  }
  return r;
}

ZZX ring::eta (long j) {
  ZZX z = ZZX(INIT_MONO, 1);
  ZZX phi = (ZZX(INIT_MONO, mR)-1)/(ZZX(INIT_MONO, 1)-1);

  long t = GofZ(j);
  ZZX w = ZZX(INIT_MONO, t);  // w = X^t
  rem(w, w, phi);

  return trace(w);
}

ZZ_pX ring::delta (long j) {
  return conv<ZZ_pX>(eta(j) - dR);
}

ZZ ring::scalar_quo (ZZX& a, ZZX& b, ZZ& p) {
  long d = max(deg(a), deg(b));
  a.SetLength(d+1);
  b.SetLength(d+1);

  // ZZ_pContext context;
  // context.save();
  ZZ_p::init(p);

  ZZ_pX a1, b1;
  a1.SetLength(d+1);
  a1 = conv<ZZ_pX>(a);
  b1.SetLength(d+1);
  b1 = conv<ZZ_pX>(b);
  
  long i;
  ZZ w;
  for (i = 0; i <= d; ++i) {
    if (GCD(b[i]%p, p) == 1) {
      w = (a[i] * InvMod(b[i]%p, p)) % p;
      break;
    }
  }

  // ZZ_pX c1 = conv<ZZ_pX>(w) * b1;
  // if (c1 != a1) {
  //   cout << "Error: unexpected behaviour in scalar_quo" << endl;
  //   cout << "i = " << i << ", p = " << p << endl;
  //   // cout << "a1 = " << a1 << endl;
  //   // cout << "b1 = " << b1 << endl;
  //   // cout << "c1 = " << c1 << endl;
  //   exit(-1);
  // }

  //context.restore();
  return w;
}

/***********************************************************************/
ZZ compose_integer (const ZZ& r1, const ZZ& q1, const ZZ& r2, const ZZ& q2) {
  ZZ a, b, c;
  XGCD(a, b, c, q1, q2);  // 1 = a = b*q1 + c*q2
  ZZ w1 = c * q2;
  ZZ w2 = b * q1;
  return (r1*w1 + r2*w2) % (q1*q2);
}

void compose_vector (vec_ZZ& omega, vec_ZZ& omega1, ZZ& q1, vec_ZZ& omega2, ZZ& q2) {
  for (long i = 0; i < omega1.length(); i++)
    omega[i] = compose_integer(omega1[i], q1, omega2[i], q2);
}

void compose (ZZX& res, ring& rg, ZZX& res0, ZZ& q0, ZZX& res1, ZZ& q1) {
  vec_ZZ v0 = VectorCopy(res0, rg.mR-1);
  vec_ZZ v1 = VectorCopy(res1, rg.mR-1);
  vec_ZZ v; v.SetLength(rg.mR-1);
  compose_vector(v, v0, q0, v1, q1);
  res.SetLength(rg.mR-1);
  for (long i = 0; i < rg.mR-1; i++) SetCoeff(res, i, v[i]);
}

// assumes that 'm' is prime
// modulus of the resolution is p^r
void ring::resolution_of_one (ZZX& res, const long m, const ZZ& p) {
  ZZ_pContext context;
  context.save();
  ZZ_p::init(p);
  
  ZZ_pX phi = (ZZ_pX(INIT_MONO, m)-1)/(ZZ_pX(INIT_MONO, 1)-1);

  cout << "." << flush;

  vec_ZZ_pX factors;
  //SFBerlekamp(factors, phi);
  if (p==ZZ(2))
    SFCanZass(factors, phi);
  else
    RootEDF(factors, phi);

  ZZ_pX G = ZZ_pX(1);
  for (long i = 1; i < factors.length(); ++i) G *= factors[i];
  ZZ_pX a, b, c;
  XGCD(a, b, c, factors[0], G); // 1 = a = b * factors[0] + c * G
  ZZ_pX res1 = c*G;
  
  // ZZ_pX res2 = (res1*res1) % phi;
  // if (res1 != res2) {
  //   cout << "Error in resolution_of_one" << endl;
  //   cout << "res   = " << res1 << endl;
  //   cout << "res2  = " << res2 << endl;
  //   exit(-1);
  // }
  
  res = conv<ZZX>(res1);
  context.restore();
}

void ring::lift_resolution (ZZX& res, const long m, const ZZ& r, const long l) {
  ZZ_p::init(power(r, l));
  ZZ_pX phi = (ZZ_pX(INIT_MONO, m)-1)/(ZZ_pX(INIT_MONO, 1)-1);

  ZZ_pX x = conv<ZZ_pX>(-res);
  for (long j = 0; j < l-1; ++j) {
    if ((j % 100) == 0) cout << "." << flush;
    x = (power(x+1, conv<long>(r)) - 1) % phi;
  }
  res = conv<ZZX>(-x);
}

void ring::canonical_resolution_of_one (ZZX& res, long prec) {
  resolution_of_one(res, mR, ZZ(pR));
  lift_resolution(res, mR, ZZ(pR), prec);    
}

void ring::init_resolution (long prec) {
  if (type == prime) {
    canonical_resolution_of_one(res, prec);
  }
  else { // ring is composite
    if (left->initialized == false)
      left->init_resolution(prec);
    if (right->initialized == false)
      right->init_resolution(prec);
  }
}

void ring::init (long prec, noise& ns) {
  if (initialized == true) return;

  cout << "Initializing ring " << label << flush;
  
  init_resolution(prec);
  
  cout << "computing omega..." << flush;
  init_omega(power(ZZ(2), prec));

  depth = get_subrings(subrings);

  initialized = true;
  cout << "done" << endl;
}

/***********************************************************************/
// Transpose a (m*n)-dim. vector 'a'
// into a m x n matrix (m-dim vector of n-dim vector)'w'
void mat (Vec<double*>& aa, const double* a, const long n, const long m) {
  for (long j = 0; j < m; ++j)
    for (long i = 0; i < n; ++i)
      aa[j][i] = a[i + n*j];
}

void transpose (Vec<double*>& b, Vec<double*>& a) {
  long m = a.length();
  long n = b.length();

  for (long i = 0; i < m; ++i)
    for (long j = 0; j < n; ++j)
      b[j][i] = a[i][j];  
}

void vec (double* b, const Vec<double*>& a, long n) {
  long m = a.length();

  for (long i = 0; i < m; ++i)
    for (long j = 0; j < n; ++j)
      b[i*n+j] = a[i][j];
}

/***************************************************************************/
// (Omega)_{i,j} = (rho_{t_i}(eta_j))
void ring::comp_omega (vec_ZZ& o, ZZ& q) {
  ZZX phi = (ZZX(INIT_MONO, mR)-1)/(ZZX(INIT_MONO, 1)-1);

  ZZ_p::init(q);
  
  o.SetLength(gR);
  for (long i = 0; i < gR; ++i) {
    if ((i % 100) == 0) cout << "." << flush;
    ZZ_pX e = conv<ZZ_pX>(eta(i));
    ZZ_pX a = (e * conv<ZZ_pX>(res)) % conv<ZZ_pX>(phi);
    ZZX aa = conv<ZZX>(a);
    o[i] = scalar_quo(aa, res, q);
  }
}

void ring::comp_omega_inv (vec_ZZ& o, const vec_ZZ& _omega, const ZZ& q) {
  o.SetLength(gR);
  
  for (long i = 0; i < gR; ++i) {
    ZZ w, conj_eta;
    if ((dR % 2) == 1) conj_eta = _omega[(i+gR/2) % gR];
    else conj_eta = _omega[i];
    w = conj_eta - dR;
    w = (w * InvMod(ZZ(mR)%q, q)) % q;
    o[i] = w;
  }  
}

// (Gamma)_{i,j} = (rho_{t_i}(eta_j - d))
void ring::comp_gamma (vec_ZZ& o, const vec_ZZ& _omega, const ZZ& q) {
  o.SetLength(gR);

  for (long i = 0; i < gR; ++i) {
    o[i] = (omega[i] - dR) % q;
  }  
}

void ring::comp_gamma_inv (vec_ZZ& o, const vec_ZZ& _omega, const ZZ& q) {
  o.SetLength(gR);
  
  for (long i = 0; i < gR; ++i) {
    ZZ conj_eta;
    if ((dR % 2) == 1) conj_eta = _omega[(i+gR/2) % gR];
    else conj_eta = _omega[i];
    o[i] = (conj_eta * InvMod(ZZ(mR)%q, q)) % q;
  }  
}

void ring::init_omega (ZZ q) {
  if (type == prime) {
    comp_omega(omega, q);
    comp_omega_inv(omega_inv, omega, q);
  }
  else { // ring type is composite
    if (left->initialized == false)
      left->init_omega(q);
    if (right->initialized == false)
      right->init_omega(q);
  }
}

void param::init (long prec) {
  rg.init(prec, ns);  
}

