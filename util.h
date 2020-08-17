#ifndef UTIL_H
#define UTIL_H

#include <NTL/vec_ZZ_p.h>
#include <NTL/RR.h>

#include "nsgen.h"
#include "ring.h"
#include "she.h"

using namespace std;
using namespace NTL;

/* Remove if already defined */
typedef long long int64; typedef unsigned long long uint64;

uint64 GetTimeMs64 ();

void debug_dec (param& prm, PlainText& m, const SKey& s, const CipherText& ct);

bool dec_ng (param& prm, PlainText& m1,
	     const SKey& s, const CipherText& ct, const PlainText& m);

bool dec_ng2 (param& prm, PlainText& m1,
	     const SKey& s, CipherText2& ct, const PlainText& m);

#endif /*UTIL_H*/
