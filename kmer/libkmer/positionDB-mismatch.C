#include "positionDB.H"
#include "bio++.H"

#include "hashIterator.H"

#define DEBUG_MISMATCH

//  Returns hits with _AT_MOST_ numMismatches mistakes.
bool
positionDB::getUpToNMismatches(u64bit   mer,
                               u64bit   numMismatches,
                               u64bit*& posn,
                               u64bit&  posnMax,
                               u64bit&  posnLen) {

  posnLen = 0;

  fprintf(stderr, "HASH("u64bitHEX") = "u64bitHEX"\n", mer, HASH(mer));

  u64bit count = 0;
  if (getExact(mer, posn, posnMax, posnLen, count)) {
    fprintf(stderr, "mer:"u64bitHEX" -- EXACT HIT\n", mer);
  }

  for (u32bit i=1; i<=numMismatches; i++)
    getExactlyNMismatches(mer, i, posn, posnMax, posnLen);

  return(posnLen > 0);
}


//  ARGH!  So, one error in the mer can show up in the hash in up to
//  three places.
//
//  If our mer is 25 bp long, and we toss out the first letter because
//  it's biased by the letter before it (the adapter), we have a 24
//  mer.  The posDB wants to give us 22 tblBits, leaving shifts of 12
//  and 24 bits.
//
//     [-----------1-------------]
//           [-----------1-------------]
//                 [------------1------------]
//
//  The two shifted mers (second and third lines) both contribute the
//  error to the hash.  The first line barely misses contributing.
//
//  Can three ever occur?

//  Returns hits with _EXACTLY_ numMismatches mistakes.
bool
positionDB::getExactlyNMismatches(u64bit   mer,
                                  u64bit   numMismatches,
                                  u64bit*& posn,
                                  u64bit&  posnMax,
                                  u64bit&  posnLen) {

  if (posnMax == 0) {
    posnMax = 16384;
    try {
      posn    = new u64bit [posnMax];
    } catch (...) {
      fprintf(stderr, "positionDB::get()-- Can't allocate space for initial positions, requested "u64bitFMT" u64bit's.\n", posnMax);
      abort();
    }
  }

  u64bit        hashorig   = HASH(mer);
  hashIterator  iter(_tableSizeInBits, numMismatches);

  //  The 'errors' loop iterates over all possible errors.  All errors
  //  are seen when we iterate to this number.  The 'errors' codes
  //  what we want to make the hash be at some position;
  //  errors=aabbccdd -- the first error is encoded in dd, etc.
  //
  u64bit errorsMax = u64bitONE << (2 * numMismatches);

  //  The hashIterator 'iter' gives us all placements of N errors.  In
  //  those positions, we insert the decoded 'errors'.
  //
  while (iter.next()) {
    for (u64bit errors=0; errors < errorsMax; errors++) {

      //  Shift the errors to the proper spot in the mer.  We assume
      //  that three ops are faster than a conditional, and this
      //  should be easy for the compiler to pipeline.  But, what do I
      //  know...
      //
      //  e? is now the thing we want to make the hash be.
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
      //  iter.error?() is the mask for that error.  If the error we
      //  want to put in there is the same as what's already there,
      //  it's not an error, and we should skip it.
      //
      bool   s1 = ((iter.error1() != 0) && (e1 == (hashorig & iter.error1())));
      bool   s2 = ((iter.error2() != 0) && (e2 == (hashorig & iter.error2())));
      bool   s3 = ((iter.error3() != 0) && (e3 == (hashorig & iter.error3())));
      bool   s4 = ((iter.error4() != 0) && (e4 == (hashorig & iter.error4())));

      //  Ah!  No, stupid.  Just because the doesn't change doesn't
      //  mean the mer doesn't change.  Two errors in the same
      //  'column' will keep the hash the same.
      //
      //  Or maybe not.
      s1 = s2 = s3 = s4 = false;


      if (s1 || s2 || s3 || s4) {
        //  Skip it.
#ifdef DEBUG_MISMATCH
        fprintf(stdout, "hash:------------------ errs:"u64bitHEX" e1:"u64bitHEX" e2:"u64bitHEX" e3:"u64bitHEX" SKIP\n",
                errors, iter.error1(), iter.error2(), iter.error3());
#endif
      } else {

        //  Build the hash with errors.
        //
        u64bit hasherrors = hashorig;

        hasherrors &= ~(iter.error1() | iter.error2() | iter.error3() | iter.error4());
        hasherrors |= e1;
        hasherrors |= e2;
        hasherrors |= e3;
        hasherrors |= e4;

        u64bit st = getDecodedValue(_hashTable, hasherrors * _hashWidth,              _hashWidth);
        u64bit ed = getDecodedValue(_hashTable, hasherrors * _hashWidth + _hashWidth, _hashWidth);

#ifdef DEBUG_MISMATCH
        fprintf(stdout, "hash:"u64bitHEX" errs:"u64bitHEX" e1:"u64bitHEX" e2:"u64bitHEX" e3:"u64bitHEX" st="u64bitFMT" ed="u64bitFMT"%s\n",
                hasherrors, errors, iter.error1(), iter.error2(), iter.error3(), st, ed, (st == ed) ? " NOTHING" : "");
#endif

        for (u64bit i=st, J=st * _wFin; i<ed; i++, J += _wFin) {

          //  Rebuild the mer from the hash and its check code.
          //
          //  Compare the rebuilt mer and the original mer -- if there
          //  are exactly N errors, it's a hit!  (if there are fewer
          //  than N, we'll find it when we look for N-1 errors).

          //  Before rebuilding, compute diffs on the chckBits only --
          //  if things are wildly different (the usual case) we'll
          //  get enough difference here to abort.

          u64bit  chck  = getDecodedValue(_buckets, J, _chckWidth);
          u64bit  chk   = chck;
          u64bit  diffs = chk ^ (mer & _mask2);
          u64bit  d1    = diffs & u64bitNUMBER(0x5555555555555555);
          u64bit  d2    = diffs & u64bitNUMBER(0xaaaaaaaaaaaaaaaa);

          fprintf(stderr, "mer "u64bitHEX" -- chk "u64bitHEX" (partial)\n", mer & _mask2, chk);

          if (countNumberOfSetBits64(d1 | (d2 >> 1)) <= numMismatches) {
            chk   = REBUILD(hasherrors, chck);
            diffs = chk ^ mer;
            d1    = diffs & u64bitNUMBER(0x5555555555555555);
            d2    = diffs & u64bitNUMBER(0xaaaaaaaaaaaaaaaa);

            fprintf(stderr, "mer "u64bitHEX" -- chk "u64bitHEX" (full)\n", mer, chk);

            if (countNumberOfSetBits64(d1 | (d2 >> 1)) == numMismatches) {
              fprintf(stderr, "mer:"u64bitHEX" -- "u64bitFMT" MISMATCHES\n", mer, numMismatches);
              u64bit c = 0;
              loadPositions(J, posn, posnMax, posnLen, c);
            }
          }
        }
      }
    }
  }

  //  True if we find matches.
  return(posnLen > 0);
}
