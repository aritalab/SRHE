#include <cassert>
#include <NTL/RR.h>
#include "ring.h"

/**********************************************************************/
ostream& operator<< (ostream& s, const ring& rg) {
  s << rg.label << ": type = " << rg.type << ", mR = " << rg.mR << ", pR = " << rg.pR
    << ", gR = " << rg.gR << ", tR = " << rg.tR << ", depth = " << rg.depth << endl;
  
  return s;
}

ostream& operator<< (ostream& s, const param& prm) {
  s << endl << "<" << prm.label << ">" << endl;
  s << "l = " << prm.l << ", ";
  s << "r = " << prm.r << ", ";
  s << "s = " << prm.s << ", ";
  s << "lw = " << prm.lw << endl;
  s << prm.rg << prm.ns;
  
  return s;
}

/**********************************************************************/
void ring::to_eta (ZZ_pX& w, vec_ZZ_p& a) {
  if (type != prime) {
    cout << "Error: to_eta" << endl;
    cout << "- not yet implemented for " << type << "-type" << endl;
  }

  for (long i = 0; i < gR; ++i) w += a[i] * conv<ZZ_pX>(eta(i));
}

void ring::to_delta (ZZ_pX& w, vec_ZZ_p& a) {
  if (type != prime) {
    cout << "Error: to_eta" << endl;
    cout << "- not yet implemented for " << type << "-type" << endl;
  }

  for (long i = 0; i < gR; ++i) w += a[i] * delta(i);
}

void mult (vec_ZZ_p& x, const vec_ZZ_p& a, const vec_ZZ_p& b) {
  for (long i = 0; i < a.length(); ++i) x[i] = a[i] * b[i];
}

void sqr (vec_ZZ_p& x, const vec_ZZ_p& a) {
  for (long i = 0; i < a.length(); ++i) x[i] = sqr(a[i]);
}

void add_scalar (vec_ZZ_p& x, const vec_ZZ_p& a, const ZZ_p& b) {
  for (long i = 0; i < a.length(); ++i) x[i] = a[i] + b;
}

void add_scalar (vec_ZZ_p& x, const ZZ_p& a, const vec_ZZ_p& b) {
  for (long i = 0; i < b.length(); ++i) x[i] = a + b[i];
}

void zzx_convolution (vec_ZZ_p& c, const vec_ZZ_p& a, const vec_ZZ_p& b) {
  ZZ_pX aX, bX;
  long g = a.length();

  aX.SetLength(g);
  bX.SetLength(g);
  c.SetLength(g);

  aX[0] = a[0];
  for (long i = 1; i < g; ++i) aX[i] = a[g - i];
  for (long i = 0; i < g; ++i) bX[i] = b[i];  
  
  ZZ_pX w = aX * bX;
  w.SetLength(2*g-1);
  
  // c = w mod m, m = x^gR - 1
  for (long i = 0; i < g-1; i++) c[i] = w[i] + w[i + g];
  c[g-1] = w[g-1];
}

// Transpose a (m*n)-dim. vector 'a'
// into a m x n matrix (m-dim vector of n-dim vector)'w'
void mat (Vec<vec_ZZ_p>& aa, const vec_ZZ_p& a, const long n, const long m) {
  aa.SetLength(m);
  for (long j = 0; j < m; ++j) aa[j].SetLength(n);
  
  for (long j = 0; j < m; ++j)
    for (long i = 0; i < n; ++i)
      aa[j][i] = a[i + n*j];
}

// returns m-dim vector of n-dim vectors
void mat (Vec<vec_ZZ>& aa, const vec_ZZ& a, const long n, const long m) {
  aa.SetLength(m);
  for (long j = 0; j < m; ++j) aa[j].SetLength(n);
  
  for (long j = 0; j < m; ++j)
    for (long i = 0; i < n; ++i)
      aa[j][i] = a[i + n*j];
}

void transpose (Vec<vec_ZZ_p>& b, const Vec<vec_ZZ_p>& a) {
  long m = a.length();
  long n = a[0].length();

  b.SetLength(n);
  for (long j = 0; j < n; ++j) b[j].SetLength(m);

  for (long j = 0; j < n; ++j)
    for (long i = 0; i < m; ++i)
      b[j][i] = a[i][j];
}

void transpose (Vec<vec_ZZ>& b, const Vec<vec_ZZ>& a) {
  long m = a.length();
  long n = a[0].length();

  b.SetLength(n);
  for (long j = 0; j < n; ++j) b[j].SetLength(m);

  for (long i = 0; i < m; ++i)
    for (long j = 0; j < n; ++j)
      b[j][i] = a[i][j];
}


void vec (vec_ZZ_p& b, const Vec<vec_ZZ_p>& a) {
  long m = a.length();
  long n = a[0].length();

  b.SetLength(m*n);

  for (long i = 0; i < m; ++i)
    for (long j = 0; j < n; ++j)
      b[i*n+j] = a[i][j];
}

void transpose_vec (vec_ZZ_p& b, const Vec<vec_ZZ_p>& a) {
  long m = a.length();
  long n = a[0].length();

  b.SetLength(m*n);

  for (long i = 0; i < n; ++i)
    for (long j = 0; j < m; ++j)
      b[i*m+j] = a[j][i];
}

void transpose_vec (vec_ZZ& b, const Vec<vec_ZZ>& a) {
  long m = a.length();
  long n = a[0].length();

  b.SetLength(m*n);

  for (long i = 0; i < n; ++i)
    for (long j = 0; j < m; ++j)
      b[i*m+j] = a[j][i];
}

void diagonal (vec_ZZ& MD, const long k, const vec_ZZ& M) {
  long dim = M.length();
  MD.SetLength(dim);

  for (long i = 0; i < dim; ++i) MD[i] = M[(i+(i+k)) % dim];
}

void diagonal (vec_ZZ_p& MD, const long k, const vec_ZZ_p& M) {
  long dim = M.length();
  MD.SetLength(dim);

  for (long i = 0; i < dim; ++i) MD[i] = M[(i+(i+k)) % dim];
}

void naiveApplyMatrix (vec_ZZ_p& b, const vec_ZZ_p& m, const vec_ZZ_p& a) {
  long dim = m.length();
  vec_ZZ_p S; S.SetLength(dim);

  vec_ZZ_p aa; //aa.SetLength(dim);
  VectorCopy(aa, a, dim);

  for (long i=0; i < dim; ++i) {
    vec_ZZ_p md;
    diagonal(md, i, m);
    
    vec_ZZ_p w; w.SetLength(dim);
    mult(w, md, aa);

    shift(aa);
    S += w;
  }

  b = S;
}

void applyMatrix (vec_ZZ_p& b, const vec_ZZ_p& m, const vec_ZZ_p& a) {
  zzx_convolution(b, a, m);
  //naiveApplyMatrix(b, m, a);  // slow, for debug only
}

void ring::applyOmega (vec_ZZ_p& b, const vec_ZZ_p& a) {
  if (type == prime) {
    vec_ZZ_p w = conv<vec_ZZ_p>(omega);
    applyMatrix(b, w, a);
  }
  else {  // type == composite
    long gl = left->gR;
    long gr = right->gR;
    
    Vec<vec_ZZ_p> w;
    mat(w, a, gr, gl);
    for (long i = 0; i < gl; ++i) right->applyOmega(w[i], w[i]);

    Vec<vec_ZZ_p> wt;
    transpose(wt, w);
    for (long j = 0; j < gr; ++j) left->applyOmega(wt[j], wt[j]);
    
    transpose_vec(b, wt);
  }
}

void ring::applyOmega (vec_ZZ& b, const vec_ZZ& a) {
  vec_ZZ_p _b;
  applyOmega(_b, conv<vec_ZZ_p>(a));
  b = conv<vec_ZZ>(_b);
}

void ring::applyOmegaInv (vec_ZZ_p& b, const vec_ZZ_p& a) {
  if (type == prime) {
    vec_ZZ_p w = conv<vec_ZZ_p>(omega_inv);
    applyMatrix(b, w, a);
  }
  else {  // type == composite
    long gl = left->gR;
    long gr = right->gR;
    
    Vec<vec_ZZ_p> w;
    mat(w, a, gr, gl);
    for (long i = 0; i < gl; ++i) right->applyOmegaInv(w[i], w[i]);

    Vec<vec_ZZ_p> wt;
    transpose(wt, w);
    for (long j = 0; j < gr; ++j) left->applyOmegaInv(wt[j], wt[j]);
    
    transpose_vec(b, wt);
  }
}

void ring::applyOmegaInv (vec_ZZ& b, const vec_ZZ& a) {
  vec_ZZ_p _b;
  applyOmegaInv(_b, conv<vec_ZZ_p>(a));
  b = conv<vec_ZZ>(_b);
}

/**********************************************************************/
ZZ inf_norm (const vec_ZZ& a) {
  ZZ w(0);
  for (long i = 0; i < a.length(); ++i) 
    if (w < abs(a[i])) w = abs(a[i]);
  return w;
}

void center_lift (vec_ZZ& b, const vec_ZZ_p& a) {
  b.SetLength(a.length());
  ZZ q = ZZ_p::modulus();
  for (long i=0; i < a.length(); ++i) {
    if (conv<ZZ>(a[i]) >= q/2) b[i] = conv<ZZ>(a[i]) - q;
    else b[i] = conv<ZZ>(a[i]);
  }
}

void center_lift (vec_ZZ& b, const vec_ZZ& a, ZZ& p) {
  ZZ_pContext context;
  context.save();  

  ZZ_p::init(p);
  center_lift(b, conv<vec_ZZ_p>(a));

  context.restore();
}

vec_ZZ center_lift (const vec_ZZ_p& a) {
  vec_ZZ b;
  center_lift(b, a);
  return b;
}

vec_ZZ center_lift (const vec_ZZ& a, ZZ& p) {
  vec_ZZ b;
  center_lift(b, a, p);
  return b;
}

// q -> q1
void rescale (vec_ZZ &b, const vec_ZZ &a, const ZZ& q, const ZZ& q1) {
  b.SetLength(a.length());
  for (long i = 0; i < a.length(); ++i) {
    b[i] = (a[i] * q1 + q/2) / q;
    if (b[i] >= q1) b[i] -= q1;
    //b[i] = b[i] % q1;  /*@@*/
  }
}

void lsd_rescale (vec_ZZ &b, const vec_ZZ &a, const ZZ& P, const ZZ& t) {
  long dim = a.length();
  b.SetLength(dim);
  vec_ZZ delta;
  delta.SetLength(dim);

  ZZ _g, u, v;
  XGCD(_g, u, v, P, t);  // 1 = w = u*P + v*t
  ZZ w1 = v*t; // = 1 mod P, = 0 mod t
  
  for (long i = 0; i < dim; ++i) {
    delta[i] = (((- a[i]) % P) * w1) % (P*t);
    b[i] = (a[i] + delta[i]) / P;
  }  
}

/***********************************************************************/
/* returns the power decomp of x mod 2^n  */
void power_decomp (Vec<vec_ZZ>& z, long n, const vec_ZZ& x) {
  long dim = x.length();
  vec_ZZ zero; zero.SetLength(dim);

  z.SetLength(n);
  for (long j = 0; j < n; ++j) z[j] = zero;

  Vec<Vec<vec_ZZ>> y;
  y.SetLength(n);
  for (long i = 0; i < n; ++i) {
    y[i].SetLength(n-i);
    for (long j = 0; j < n-i; ++j) y[i][j] = zero;
  }

  for (long i = 0; i < n; ++i) {
    // y[i][0] = (x - sum_{j=0}^{i-1} 2^j * y[j][i-j]) / 2^i
    y[i][0] = x;
    for (long j = 0; j < i; ++j) {
      for (long k = 0; k < dim; ++k) {
	y[i][0][k] = (y[i][0][k] - y[j][i-j][k]) % power(ZZ(2), n);
	y[i][0][k] /= ZZ(2);
      }
    }

    // y[i][j] = y[i][j-1]^2  (j = 1..(n-1-i))
    for (long j = 1; j < n-i; ++j) {
      for (long k = 0; k < dim; ++k) y[i][j][k] = (y[i][j-1][k] * y[i][j-1][k]) % power(ZZ(2), n);
    }
  }

  for (long i = 0; i < n; ++i) {
    z[i] = y[i][n-1-i];
    for (long k = 0; k < dim; ++k) z[i][k] = z[i][k] % power(ZZ(2), n-i);
  }  
}

// returns (x >> m) mod 2^(n-m) given (x mod 2^n)
void right_shift (vec_ZZ& z, const vec_ZZ& x, long m, long n) {
  long dim = x.length();
  vec_ZZ zero; zero.SetLength(dim);

  long l = n - m;

  Vec<vec_ZZ> y; y.SetLength(n);
  for (long j = 0; j < n; ++j) y[j] = zero;
  power_decomp(y, n, x);

  // Drop the least m bits (LowerClear)
  for (long i = 0; i < l; ++i) y[i] = y[m+i];

  // z = y[0] + 2*y[1] + ... + 2^(l-1)*y[l-1]
  z = zero;
  for (long i = 0; i < l; ++i) {
    // z = (z + (2^i)*y[i]) % 2^l
    vec_ZZ w = zero;
    for (long j = 0; j < dim; ++j) {
      w[j] = (power(ZZ(2),i) * y[i][j]) % power(ZZ(2),l);
      z[j] = (z[j] + w[j]) % power(ZZ(2),l);
    }
  }
}

void right_shift (vec_ZZ_p& z, const vec_ZZ_p& x, long m, long n) {
  vec_ZZ _z; _z.SetLength(x.length());
  right_shift(_z, conv<vec_ZZ>(x), m, n);
  z = conv<vec_ZZ_p>(_z);
}

// returns the (i,j)-th component of the cyclic matrix with the first row of 'm'
ZZ_p to_mat (const vec_ZZ_p& m, long i, long j) {
  long dim = m.length();
  assert(0 <= i && i < dim && 0 <= j && j < dim);
  return m[(i + j) % dim];
}

/***********************************************************************/

/* returns the tensor product of a and a1 */
void tensor_prod (vec_ZZ& A, const vec_ZZ& a, const vec_ZZ& a1, const ZZ& q) {
  long gout = a.length();
  long gin = a1.length();

  for (long i = 0; i < gout; ++i) {
    for (long j = 0; j < gin; ++j) {
      A[i*gin + j] = (a[i] * a1[j]) % q;
    }
  }
}

// left shift of vectors
void shift (vec_ZZ& b, const long i, const vec_ZZ& a) {
  long dim = a.length();
  vec_ZZ w;
  w.SetLength(dim);
  for (long j = 0; j < dim; ++j) w[j] = a[(j + i) % dim];
  b = w;
}

void shift (vec_ZZ_p& b, const long i, const vec_ZZ_p& a) {
  long dim = a.length();
  vec_ZZ_p w;
  w.SetLength(dim);
  for (long j = 0; j < dim; ++j) w[j] = a[(j + i) % dim];
  b = w;
}

vec_ZZ shift (const long i, const vec_ZZ& a) {
  vec_ZZ w;
  shift(w, i, a);
  return w;
}

// one-position left shift
void shift (vec_ZZ_p& a) {
  long dim = a.length();
  ZZ_p w = a[0];
  for (long i = 1; i < dim; ++i) a[i-1] = a[i];
  a[dim-1] = w;
}

void shift_right (vec_ZZ& b, const long i, const vec_ZZ& a) {
  shift(b, a.length()-i, a);
}

void shift_right (vec_ZZ_p& b, const long i, const vec_ZZ_p& a) {
  shift(b, a.length()-i, a);
}

// one-position right shift
void right_shift (vec_ZZ_p& a) {
  long dim = a.length();
  ZZ_p w = a[dim-1];
  for (long i = dim-1; i > 0; --i) a[i] = a[i-1];
  a[0] = w;
}

void ring::deep_shift (string _label, vec_ZZ& v, long i) {
  if ((type == prime) && (label == _label)) {
    shift(v, (i<0)?(i+gR):i, v);
  }
  else if ((type == prime) && (label != _label)) {
    return;
  }
  else {  // composite ring
    long gl = left->gR;
    long gr = right->gR;
    
    Vec<vec_ZZ> w;
    mat(w, v, gr, gl);
    for (long j = 0; j < gl; ++j) right->deep_shift(_label, w[j], i);

    Vec<vec_ZZ> wt;
    transpose(wt, w);
    for (long j = 0; j < gr; ++j) left->deep_shift(_label, wt[j], i);

    transpose_vec(v, wt);
  }
}

void ring::deep_shift_right (string _label, vec_ZZ& v, long i) {
  deep_shift(_label, v, -i);
}

/***********************************************************************/
void get_subrings (ring& rg, Vec<ring>& subrings) {
  if (rg.type == prime) {
    append(subrings, rg);
  }
  else {  // type == composite
    get_subrings(*(rg.right), subrings);
    get_subrings(*(rg.left), subrings);
  }  
}

long ring::get_subrings (Vec<ring>& subrings) {
  ::get_subrings(*this, subrings);
  return subrings.length();
}

long param::get_subrings (Vec<ring>& subrings) {
  return rg.get_subrings(subrings);
}

long ring::get_subring_index (string label) {
  Vec<ring> subrings;
  long depth = get_subrings(subrings);

    for (long i = 0; i < depth; ++i) {
      if (subrings[i].label == label) return i;
    }

    return -1;
}

vec_ZZ ring::index_convert (long j) {
  vec_ZZ tensor_index;
  tensor_index.SetLength(depth);

  long jj = j;
  
  for (long i = 0; i < depth; ++i) {
    long a = jj % subrings[i].gR;
    tensor_index[i] = a;
    jj = (jj - a) / subrings[i].gR;
  }

  return tensor_index;
}

// y = (x, x, ..., x) where x is assumed to be an element of the right subring
void ring::embed_1 (vec_ZZ& y, vec_ZZ& x) {
  y.SetLength(gR);
  for (long i = 0; i < left->gR; ++i) {
    for (long j = 0; j < right->gR; ++j) {
      y[i*right->gR + j] = x[j];
    }
  }
}

vec_ZZ ring::embed_1 (vec_ZZ& x) {
  vec_ZZ y;
  embed_1(y, x);
  return y;
}

// y = (x, 0, 0, ..., 0) where x is assumed to be an element of the right subring
void ring::embed_0 (vec_ZZ& y, vec_ZZ& x) {
  y.SetLength(gR);
  for (long j = 0; j < right->gR; ++j) {
    y[j] = x[j];
  }
}

vec_ZZ ring::embed_0 (vec_ZZ& x) {
  vec_ZZ y;
  embed_0(y, x);
  return y;
}

// y = (x, 0, 0, ..., 0) where x is assumed to be an element of the prime subring with label subring_label
void ring::embed_0 (vec_ZZ& y, string subring_label, vec_ZZ& x) {
  y.SetLength(gR);
  long d = get_subring_index(subring_label);

  for (long j = 0; j < gR; ++j) {
    vec_ZZ tensor_index = index_convert(j);
    long skip = 0;
    for (long k = 0; k < tensor_index.length(); ++k) {
      if (k != d && tensor_index[k] != 0) {
        skip = 1;
        break;
      }
    }
    if (!skip) y[j] = x[conv<long>(tensor_index[d])];
  }
}

vec_ZZ ring::embed_0 (string subring_label, vec_ZZ& x) {
  vec_ZZ y;
  embed_0(y, subring_label, x);
  return y;
}

// y = (x, x, ..., x) where x is assumed to be an element of the prime subring with label subring_label
void ring::embed_1 (vec_ZZ& y, string subring_label, vec_ZZ& x) {
  y.SetLength(gR);
  long d = get_subring_index(subring_label);

  for (long j = 0; j < gR; ++j) {
    vec_ZZ tensor_index = index_convert(j);
    y[j] = x[conv<long>(tensor_index[d])];
  }
}

vec_ZZ ring::embed_1 (string subring_label, vec_ZZ& x) {
  vec_ZZ y;
  embed_1(y, subring_label, x);
  return y;
}

// y = (-x, -x, ..., -x) where x is assumed to be an element of the prime subring with label subring_label
void ring::embed_11 (vec_ZZ& y, string subring_label, vec_ZZ& x) {
  y.SetLength(gR);
  long d = get_subring_index(subring_label);

  for (long j = 0; j < gR; ++j) {
    vec_ZZ tensor_index = index_convert(j);
    y[j] = -x[conv<long>(tensor_index[d])];
  }
}

vec_ZZ ring::embed_11 (string subring_label, vec_ZZ& x) {
  vec_ZZ y;
  embed_11(y, subring_label, x);
  return y;
}

// trace to the right ring
void ring::trace (vec_ZZ& y, vec_ZZ& x, ZZ& q) {
  y.SetLength(right->gR);
  for (long j = 0; j < x.length(); ++j) {
    long i = j % right->gR;
    y[i] = (y[i] + x[j]) % q;
  }
}

// trace to the subring of label 'subring_label'
// -- the subring must be prime
void ring::trace (vec_ZZ& y, string subring_label, vec_ZZ& x, ZZ& q) {
  //y.SetLength(gR);
  long d = get_subring_index(subring_label);
  y.SetLength(subrings[d].gR);

  for (long j = 0; j < x.length(); ++j) {
    vec_ZZ tensor_index = index_convert(j);
    long i = conv<long>(tensor_index[d]);
    y[i] = (y[i] + x[j]) % q;
  }
}

// returns the i-th diagonal vector of the extended Omega matrix,
// corresponding to the subring with label '_label'
void ring::deep_omega (vec_ZZ& extended_omega, string _label, long i) {
  Vec<ring> subrings;
  get_subrings(subrings);
  long d = get_subring_index(_label);
  
  vec_ZZ omg;
  diagonal(omg, i, subrings[d].omega);
  
  //embed_0(extended_omega, _label, omg); /*@@*/
  embed_1(extended_omega, _label, omg); /*@@*/
}

// returns the i-th diagonal vector of the extended Omega_inv matrix,
// corresponding to the subring with label '_label'
void ring::deep_omega_inv (vec_ZZ& extended_omega_inv, string _label, long i) {
  Vec<ring> subrings;
  get_subrings(subrings);
  long d = get_subring_index(_label);
  
  vec_ZZ omg_inv;
  diagonal(omg_inv, i, subrings[d].omega_inv);

  //embed_0(extended_omega_inv, _label, omg_inv); /*@@*/
  embed_1(extended_omega_inv, _label, omg_inv); /*@@*/
}

