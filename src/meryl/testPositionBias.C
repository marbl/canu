#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "bio++.H"
#include "libmeryl.H"


uint32 *
collectCounts(char *name, uint32 base) {
  merylStreamReader   *A = new merylStreamReader(name);
  uint32              *C = new uint32 [4*4*4];
  char                 S[32];
  uint32               code = 0;

  for (uint32 i=0; i<4*4*4; i++)
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
showBias(uint32 base=5) {
  uint32  *A = collectCounts("CNPT3", base);
  uint32  *B = collectCounts("25.errorless", base);
  uint32  *C = collectCounts("25.errorless.simulated", base);

  for (uint32 i=0; i<4*4*4; i++) {
    double   bp = 0.0;
    double   cp = 0.0;

    if (A[i] > 0) {
      bp = (double)B[i] / (double)A[i];
      cp = (double)C[i] / (double)A[i];
    }

    fprintf(stdout, "%c%c%c "uint32FMTW(3)" A "uint32FMTW(6)" B "uint32FMTW(6)" %.5f C "uint32FMTW(6)" %.5f\n",
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
computeRMSD(uint32 base) {
  uint32  *A = collectCounts("CNPT3", base);
  uint32  *B = collectCounts("25.errorless", base);
  uint32  *C = collectCounts("25.errorless.simulated", base);

  double rmsd = 0;

  for (uint32 i=0; i<4*4*4; i++) {
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

  //for (uint32 i=0; i<23; i++)
  //  fprintf(stdout, "rmsd "uint32FMTW(2)" %f\n", i, computeRMSD(i));
}
