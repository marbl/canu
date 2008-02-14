#include "positionDB.H"
#include "bio++.H"

#include "hashIterator.H"

bool
positionDB::getMismatch(u64bit   mer,
                        u64bit   numMismatches,
                        u64bit*& posn,
                        u64bit&  posnMax,
                        u64bit&  posnLen) {

  posnLen = 0;

  if (posnMax == 0) {
    posnMax = 16384;
    try {
      posn    = new u64bit [posnMax];
    } catch (...) {
      fprintf(stderr, "positionDB::get()-- Can't allocate space for initial positions, requested "u64bitFMT" u64bit's.\n", posnMax);
      abort();
    }
  }

  u64bit  hashorig   = HASH(mer);
  u64bit  hasherrors = 0;

  hashIterator  iter(_tableSizeInBits, numMismatches);

  //  All errors are seen when we iterate to this number.
  //
  u64bit errorsMax = u64bitONE << (2 * numMismatches);

  while (iter.next()) {
    for (u64bit errors=0; errors < errorsMax; errors++) {

      //  Shift the errors to the proper spot in the mer.  We assume
      //  that three ops are faster than a conditional, and this
      //  should be easy for the compiler to pipeline.  But, what do I
      //  know...
      //
      u64bit e1 = ((errors >> 0) & 0x00000003) << (iter.shift1());
      u64bit e2 = ((errors >> 2) & 0x00000003) << (iter.shift2());
      u64bit e3 = ((errors >> 4) & 0x00000003) << (iter.shift3());
      u64bit e4 = ((errors >> 6) & 0x00000003) << (iter.shift4());

      //  If any of those errors aren't errors, we don't want to use
      //  this hash, as it represents a hash with fewer errors that we
      //  expected.  The problem is that other error patterns would do
      //  the same work; for example, if the error mask is 011001, and
      //  the last 1 isn't an error, the mask 011010 does the same
      //  work all over again.
      //
      bool   s1 = ((iter.error1() != 0) && (e1 == (hashorig & iter.error1())));
      bool   s2 = ((iter.error2() != 0) && (e2 == (hashorig & iter.error2())));
      bool   s3 = ((iter.error3() != 0) && (e3 == (hashorig & iter.error3())));
      bool   s4 = ((iter.error4() != 0) && (e4 == (hashorig & iter.error4())));

      if (s1 || s2 || s3 || s4) {
        //  Skip it.

        //fprintf(stdout, "hash:"u64bitHEX" errs:"u64bitHEX" e1:"u64bitHEX" e2:"u64bitHEX" e3:"u64bitHEX" SKIP\n",
        //        hasherrors, errors, iter.error1(), iter.error2(), iter.error3());
      } else {

        //fprintf(stdout, "hash:"u64bitHEX" errs:"u64bitHEX" e1:"u64bitHEX" e2:"u64bitHEX" e3:"u64bitHEX"\n",
        //        hasherrors, errors, iter.error1(), iter.error2(), iter.error3());

        //  Build the hash with errors.
        //
        hasherrors  = hashorig;
        hasherrors &= ~(iter.error1() | iter.error2() | iter.error3() | iter.error4());
        hasherrors |= e1;
        hasherrors |= e2;
        hasherrors |= e3;
        hasherrors |= e4;

        u64bit st = getDecodedValue(_hashTable, hasherrors * _hashWidth,              _hashWidth);
        u64bit ed = getDecodedValue(_hashTable, hasherrors * _hashWidth + _hashWidth, _hashWidth);

        for (u64bit i=st, J=st * _wFin; i<ed; i++, J += _wFin) {

          //  Rebuild the mer from the hash and its check code.
          //
          u64bit  chk = REBUILD(hasherrors, getDecodedValue(_buckets, J, _chckWidth));

          //  Compare the rebuilt mer and the original mer -- if there
          //  are exactly N errors, it's a hit!  (if there are fewer
          //  than N, we'll find it when we look for N-1 errors).

          u64bit  diffs = chk ^ mer;
          u64bit  d1    = diffs & u64bitNUMBER(0x5555555555555555);
          u64bit  d2    = diffs & u64bitNUMBER(0xaaaaaaaaaaaaaaaa);
          u64bit  count = 0;

          if (countNumberOfSetBits64(d1 | (d2 >> 1)) == numMismatches)
            loadPositions(J, posn, posnMax, posnLen, count);
        }
      }
    }
  }

  //  True if we find matches.
  return(posnLen > 0);
}
