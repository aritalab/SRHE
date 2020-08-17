#include <sys/time.h>
#include <ctime>
#include "util.h"

/* Returns the amount of milliseconds elapsed since the UNIX epoch. Works on both
 * windows and linux. */

uint64 GetTimeMs64 ()
{
 struct timeval tv;

 gettimeofday(&tv, NULL);

 uint64 ret = tv.tv_usec;
 /* Convert from micro seconds (10^-6) to milliseconds (10^-3) */
 ret /= 1000;

 /* Adds the seconds (10^0) after converting them to milliseconds (10^-3) */
 ret += (tv.tv_sec * 1000);

 return ret;
}

//----------------------------------------------------------------------
void debug_dec (param& prm, PlainText& m, const SKey& s, const CipherText& ct) {
  CipherText dt(prm);
  copy_ct(dt, ct);
  
  ZZ v = decrypt(prm, m, s, dt);

  cout << "log2(noise) = " << log(conv<RR>(v))/log(conv<RR>(ZZ(2)))
       << " with ratio = " << log(conv<RR>(v)/conv<RR>(ct.q))/log(conv<RR>(ZZ(2))) << endl;
}

bool dec_ng (param& prm, PlainText& m1,
	     const SKey& s, const CipherText& ct, const PlainText& m) {
  CipherText _ct(ct);
  ZZ v = decrypt(prm, m1, s, _ct);
  if (v < 1) v = 1;

  //cout << "v = " << v << ", ";
  cout << "log2(noise) = " << log(conv<RR>(v))/log(conv<RR>(ZZ(2))) << " with ratio = "
       << log(conv<RR>(v)/conv<RR>(_ct.q))/log(conv<RR>(ZZ(2))) << endl;
  return (m1 != m);
}

bool dec_ng2 (param& prm, PlainText& m1,
	     const SKey& s, CipherText2& ct, const PlainText& m) {
  CipherText2 _ct(ct);
  ZZ v = decrypt2(prm, m1, s, _ct);
  if (v < 1) v = 1;

  cout << "log2(noise) = " << log(conv<RR>(v))/log(conv<RR>(ZZ(2))) << " with ratio = "
       << log(conv<RR>(v)/conv<RR>(_ct.q))/log(conv<RR>(ZZ(2))) << endl;
  return (m1 != m);
}
