#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "bio++.H"
#include "libmeryl.H"


u32bit *
collectCounts(char *name, u32bit base) {
  merylStreamReader   *A = new merylStreamReader(name);
  u32bit              *C = new u32bit [4*4*4];
  char                 S[32];
  u32bit               code = 0;

  for (u32bit i=0; i<4*4*4; i++)
    C[i] = 0;

  while (A->nextMer()) {
    A->theFMer().merToString(S);

    code   = 0;
    code  |= letterToBits[S[base]];
    code <<= 2;
    code  |= letterToBits[S[base+1]];
    code <<= 2;
    code  |= letterToBits[S[base+2]];

    C[code] += A->theCount();

    kMer R = A->theFMer();
    R.reverseComplement();
    R.merToString(S);

    code   = 0;
    code  |= letterToBits[S[base]];
    code <<= 2;
    code  |= letterToBits[S[base+1]];
    code <<= 2;
    code  |= letterToBits[S[base+2]];

    C[code] += A->theCount();
  }

  delete A;

  return(C);
}


void
showBias(u32bit base=5) {
  u32bit  *A = collectCounts("CNPT3", base);
  u32bit  *B = collectCounts("25.errorless", base);
  u32bit  *C = collectCounts("25.errorless.simulated", base);

  for (u32bit i=0; i<4*4*4; i++) {
    double   bp = 0.0;
    double   cp = 0.0;

    if (A[i] > 0) {
      bp = (double)B[i] / (double)A[i];
      cp = (double)C[i] / (double)A[i];
    }

    fprintf(stdout, "%c%c%c "u32bitFMTW(3)" A "u32bitFMTW(6)" B "u32bitFMTW(6)" %.5f C "u32bitFMTW(6)" %.5f\n",
            bitsToLetter[(i >> 4) & 0x00000003],
            bitsToLetter[(i >> 2) & 0x00000003],
            bitsToLetter[(i >> 0) & 0x00000003],
            i,
            A[i],
            B[i],
            bp,
            C[i],
            cp);
  }
}


double
computeRMSD(u32bit base) {
  u32bit  *A = collectCounts("CNPT3", base);
  u32bit  *B = collectCounts("25.errorless", base);
  u32bit  *C = collectCounts("25.errorless.simulated", base);

  double rmsd = 0;

  for (u32bit i=0; i<4*4*4; i++) {
    double   bp = 0.0;
    double   cp = 0.0;

    if (A[i] > 0) {
      bp = (double)B[i] / (double)A[i];
      cp = (double)C[i] / (double)A[i];
    }

    rmsd += (bp - cp) * (bp - cp);
  }

  rmsd /= 4*4*4;

  return(sqrt(rmsd));
}





int
main(int argc, char **argv) {

  showBias(5);

  //for (u32bit i=0; i<23; i++)
  //  fprintf(stdout, "rmsd "u32bitFMTW(2)" %f\n", i, computeRMSD(i));
}
