#include <stdio.h>
#include <stdlib.h>


#include "../libutil/util++.H"
#include "hashIterator.H"


int
main(int argc, char **argv) {

  u64bit   merorig    = 0;

  u64bit   hashorig   = 0;
  u64bit   hasherrors = 0;

  u32bit   skipped = 0;
  u32bit   usable  = 0;

  //  1,000,000 hashes with three errors in 156u

  for (hashorig = 0; hashorig < 1000000; hashorig++) {
    hashIterator  iter(26, 3);

    while (iter.next()) {

      //  40 == 01000000 in binary
      for (u64bit errors=0; errors < 0x00000040; errors++) {

        //  Build the hash with errors
        hasherrors  = hashorig;
        hasherrors &= ~(iter.error1() | iter.error2() | iter.error3());

        hasherrors |= (errors & 0x00000003) << (iter.shift1() - 0);
        hasherrors |= (errors & 0x0000000c) << (iter.shift2() - 2);
        hasherrors |= (errors & 0x00000030) << (iter.shift3() - 4);

        //  If any of those errors aren't errors, we don't want to use
        //  this hash, as it represents a hash with fewer errors that we
        //  expected.  The problem is that other error patterns would do
        //  the same work; for example, if the error mask is 011001, and
        //  the last 1 isn't an error, the mask 011010 does the same
        //  work all over again.

        if (((hasherrors & iter.error1()) == (hashorig & iter.error1())) ||
            ((hasherrors & iter.error2()) == (hashorig & iter.error2())) ||
            ((hasherrors & iter.error3()) == (hashorig & iter.error3()))) {
          //  Skip it.
          skipped++;
          //fprintf(stdout, "hash:"u64bitHEX" errs:"u64bitHEX" e1:"u64bitHEX" e2:"u64bitHEX" e3:"u64bitHEX" SKIP\n",
          //        hasherrors, errors, iter.error1(), iter.error2(), iter.error3());
        } else {
          usable++;
          //fprintf(stdout, "hash:"u64bitHEX" errs:"u64bitHEX" e1:"u64bitHEX" e2:"u64bitHEX" e3:"u64bitHEX"\n",
          //        hasherrors, errors, iter.error1(), iter.error2(), iter.error3());
        }
      }
    }
  }

  fprintf(stderr, "skipped "u32bitFMT" usable "u32bitFMT"\n", skipped, usable);

  return(0);
}
