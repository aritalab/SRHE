/* Parameter Setting */

ring left_ring(/*label=*/"Ring-127", /*mR=*/127, /*pR=*/2, /*dR=*/7, /*gR=*/18, /*tR=*/3);
ring right_ring(/*label=*/"Ring-1801", /*mR=*/1801, /*pR=*/2, /*dR=*/25, /*gR=*/72, /*tR=*/11);

// composite rings
ring composed_ring(/*label=*/"Ring-127-1801", &left_ring, &right_ring);

noise ns_32(/*label=*/"noise 3.2", /*sigma=*/3.2, /*hwt=*/128);
noise ns_64(/*label=*/"noise 6.4", /*sigma=*/6.4, /*hwt=*/128);

param left_param(/*label=*/"left_param", /*rg=*/left_ring, /*level=*/5,
		  /*l=*/4,/*r=*/12*4+10, /*s=*/11, /*lw=*/6, /*ns=*/ns_32);
param left_param_tmp(/*label=*/"left_param_tmp", /*rg=*/left_ring,  /*level=*/1,
		      /*l=*/4,/*r=*/9, /*s=*/1, /*lw=*/4, /*ns=*/ns_32);
param left_param_big(/*label=*/"left_param_big", /*rg=*/left_ring,  /*level=*/45,
		      /*l=*/9,/*r=*/17*44+10, /*s=*/17*5, /*lw=*/10, /*ns=*/ns_32);

param right_param(/*label=*/"right_param", /*rg=*/right_ring, /*level=*/5,
		  /*l=*/4,/*r=*/24*4+10, /*s=*/12, /*lw=*/6, /*ns=*/ns_32);
param right_param_tmp(/*label=*/"right_param_tmp", /*rg=*/right_ring,  /*level=*/1,
		      /*l=*/4,/*r=*/17, /*s=*/1, /*lw=*/4, /*ns=*/ns_32);
param right_param_big(/*label=*/"right_param_big", /*rg=*/right_ring,  /*level=*/45,
		      /*l=*/17,/*r=*/37*44+28, /*s=*/37*4, /*lw=*/10, /*ns=*/ns_32);

param composed_param(/*label=*/"composed_param", /*rg=*/composed_ring, /*level=*/5,
		  /*l=*/4,/*r=*/28*4+18, /*s=*/12, /*lw=*/6, /*ns=*/ns_32);
param composed_param_tmp(/*label=*/"composed_param_tmp", /*rg=*/composed_ring,  /*level=*/1,
		      /*l=*/4,/*r=*/21, /*s=*/1, /*lw=*/4, /*ns=*/ns_32);
param composed_param_big(/*label=*/"composed_param_big", /*rg=*/composed_ring,  /*level=*/45,
		      /*l=*/22,/*r=*/46*44+30, /*s=*/45*4, /*lw=*/21, /*ns=*/ns_32);
