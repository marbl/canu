#include <stdio.h>
#include <stdlib.h>

#include "util++.H"
#include "bio++.H"
#include "libmeryl.H"
#include "existDB.H"

#include "seqCache.H"
#include "seqStream.H"
#include "merStream.H"

//  1) Use meryl to find the list of mers in common between HMISSING and NCBI.
//  2) Read mers from that meryl database into an existDB.
//  3) Stream each sequence from HMISSING.  Mask out any mer in HMISSING.

#define DEBUG

int
main(int argc, char **argv) {
  const char   *merName = "/project/huref4/assembly-mapping/missing/missing0/HMISSING-and-B35LC";
  const char   *seqName = "/project/huref4/assembly-mapping/missing/missing0/HMISSING-ge64.uid.fasta";
  u32bit        merSize = 28;

  fprintf(stderr, "Build existDB.\n");

  existDB *exist = 0L;
  if (0) {
    exist = new existDB(merName, merSize, existDBnoFlags, 0, ~u32bitZERO);
    exist->saveState("/project/huref4/assembly-mapping/missing/missing0/HMISSING-and-B35LC.existDB");
  } else {
    exist = new existDB("/project/huref4/assembly-mapping/missing/missing0/HMISSING-and-B35LC.existDB");
  }

  seqCache     *F     = new seqCache(seqName);
  seqInCore    *S     = 0L;

  u32bit   maskLen = 1048576;
  bool    *mask    = new bool [maskLen];
  char    *maskSeq = new char [maskLen];
  char    *maskBit = new char [maskLen];

  fprintf(stderr, "Begin.\n");

  while ((S = F->getSequenceInCore()) != 0L) {
    //fprintf(stderr, "iid="u32bitFMT" len="u32bitFMT" %s\n", S->getIID(), S->sequenceLength(), S->header());

    if (maskLen <= S->sequenceLength() + 1024) {
      maskLen = S->sequenceLength() + S->sequenceLength() + 1024;

      delete [] mask;
      delete [] maskSeq;
      delete [] maskBit;

      mask    = new bool [maskLen];
      maskSeq = new char [maskLen];
      maskBit = new char [maskLen];
    }

    for (u32bit i=0; i<S->sequenceLength(); i++)
      mask[i] = false;

    //  Build the initial masking
    //
    merStream    *MS = new merStream(new kMerBuilder(merSize),
                                     new seqStream(S->sequence(), S->sequenceLength()),
                                     true, true);
    while (MS->nextMer())
      if (exist->exists(MS->theFMer()) || exist->exists(MS->theRMer()))
        mask[MS->thePositionInSequence()] = true;
    delete    MS;


#ifdef PRINT_BITS
    for (u32bit i=0; i<S->sequenceLength(); i++) {
      if (mask[i])
        maskBit[i] = '1';
      else
        maskBit[i] = '0';
    }
    maskBit[S->sequenceLength()] = 0;
    fprintf(stdout, "%s\n",     maskBit);
#endif


    //  Clean up the masking.
    //
    //  If there are no spurious mer matchs, a real mismatch (and also an insert) should look like:
    //
    //    AAAAAAAAAAAAAAAAAAAAAAAAAXAAAAAAAA.....
    //    1111111111111111111110000011111111
    //                   -----
    //                    -----
    //                     -----
    //                      -----
    //                       -----
    //                        -----
    //                              -----
    //                               -----
    //  The mers cover all letters, except for the mismatch.
    //
    //
    //  A break in contiguous sequence (or a deletion) is very similar, except
    //  mers cover all letters, they just don't span the gap.
    //    AAAAAAAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCC....
    //    111111111111111111111000011111111
    //                   -----
    //                    -----
    //                     -----
    //                      -----
    //                       -----
    //                        -----
    //                             -----
    //                              -----
    //
    //
    //  Thus, A block of fewer than merSize 0's indicates that we have
    //  spurious mer hits, and don't find any contiguous sequence
    //  match.
    //
    //  This block marks the zeros and ms-1 bases on either side as questionable
    //  sequence.
    //
    //  I'm not sure if the blocks outside are questionable -- we do
    //  find the mer, and overlapping mers on the other side.  So
    //  possibly the unmapped sequence is just:
    //    11111111111111100110010001111111111111111111
    //                   ----------
    //
    //
    //  Alg is now to look for blocks of less than merSize zeros,
    //  and....do something.  Certainly, chain blocks together if
    //  there is no intervening block of merSize 1's in between.
    //  Possibly extend the block merSize-1 bases to either side.
    //
    //
    for (u32bit i=0; i<S->sequenceLength(); i++) {

      //  Now we're potentially at the start of a small block.  The
      //  goal is to see if this is an isolated block of spurious
      //  match, and if so, merge it into the surrounding blocks of
      //  mismatch.
      //
      //  If it's a small block of 1's, we extend the block into the
      //  next small blocks, stopping at the start of a large block.
      //  e.g.: (the --'s are the end of the small blocks):
      //
      //   ...0000000--1111110000011111--000000000000000000000...
      //   ...0000000--1110001100110010--111111111111111111111...
      //   ...1111111--0000111100100000--111111111111111111111...
      //   ...1111111--0000111100100011--000000000000000000000...
      //

      bool    isIsolated = true;
      u32bit  blockEnd   = i;
      u32bit  bs         = 0;
      u32bit  nextStart  = i;

      while (isIsolated) {
        u32bit j  = blockEnd;

        while ((j < S->sequenceLength()) &&
               (mask[j] == mask[blockEnd]))
          j++;

        bs = j - blockEnd;

        //  Remember how far we got.
        nextStart += j;

        //  If a big block terminate.
        //
        if (bs >= merSize)
          isIsolated = false;
        else
          blockEnd += bs;

        //  If now at the end of the sequence, terminate.
        //
        if (blockEnd >= S->sequenceLength())
          isIsolated = false;
      }


      //  But if we stop because we find a big block of zeros, then
      //  the previous block of ones is still valid.  Our big block of
      //  zeros indicates a true mismatch.
      //
      //  .....11111000110000111100000000000000000000000000111
      //

      //  If we did find a small block
      if (i != blockEnd) {

        //  And that block doesn't run up to the end
        if (blockEnd < S->sequenceLength()) {

          //  And the last thing in that block is 1's, move back to the previous block of zeros.
          //
          //  EXCEPT, that this big block of zeros might REALLY be
          //  unmatched sequence, not just an isolated error.  So we
          //  don't reverse if the big block of zeros is big.
          //
          if (bs < merSize + 6)
            while ((blockEnd > 0) && (mask[blockEnd-1] == true))
              blockEnd--;
        }

        fprintf(stdout, "mask "u32bitFMT" "u32bitFMT"\n", i, blockEnd);

        //  Supposedly, our block of unmasked sequence is now from i to
        //  blockEnd.  Flip all that stuff to unmasked.
        //
        for (u32bit j=i; j<blockEnd; j++)
          mask[j] = false;
      }


      //  And move up to the next starting position.
      i = nextStart - 1;
    }


#ifdef PRINT_BITS
    for (u32bit i=0; i<S->sequenceLength(); i++) {
      if (mask[i])
        maskBit[i] = '1';
      else
        maskBit[i] = '0';
    }
    maskBit[S->sequenceLength()] = 0;
    fprintf(stdout, "%s\n",     maskBit);
#endif



    //  Convert the mer-markings to base-markings
    //
    u32bit isMasking = 0;
    for (u32bit i=0; i<S->sequenceLength(); i++) {
      if (mask[i])
        isMasking = merSize;

      if (isMasking > 0) {
        mask[i] = true;
        isMasking--;
      }
    }


#ifdef PRINT_BITS
    for (u32bit i=0; i<S->sequenceLength(); i++) {
      if (mask[i])
        maskBit[i] = '1';
      else
        maskBit[i] = '0';
    }
    maskBit[S->sequenceLength()] = 0;
    fprintf(stdout, "%s\n",     maskBit);
#endif



    char *seq          = S->sequence();

    for (u32bit i=0; i<S->sequenceLength(); i++) {
      if (mask[i]) {
        maskSeq[i] = 'n';
      } else {
        maskSeq[i] = seq[i];
      }
    }

    maskSeq[S->sequenceLength()] = 0;

    fprintf(stdout, "%s\n%s\n", S->header(), maskSeq);

    delete S;
  }

  delete [] maskSeq;
  delete [] mask;
}
