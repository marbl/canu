#include <stdio.h>
#include <stdlib.h>

#include "bio++.H"

//  Splits a sequence into itty-bitty pieces.
//
//  By default, splits into non-overlapping pieces of length L.
//  Pieces will not start with nor end with N, but may have embedded N's.
//
//  All pieces will be at least L long.  Most pieces will be exactly L
//  long.  All pieces will be less than 2L long.
//
//  If a piece has more than (currently) 50 N's, it will be broken --
//  the first piece and last piece will be saved, and the middle (with
//  the N's) will be discarded.


//#define DEBUG


int
main(int argc, char **argv) {

  if (argc != 5) {
    fprintf(stderr, "usage: %s <L> <input-file> <out-file> <log-file>\n", argv[0]);
    exit(1);
  }

  u32bit   desiredLength = strtou32bit(argv[1], 0L);
  char    *seqFile       = argv[2];

  FastAWrapper         *F = new FastAWrapper(seqFile);
  FastASequenceInCore  *B;

  u32bit       fragmentIndex = 0;

  FILE        *O = fopen(argv[3], "w");
  FILE        *L = fopen(argv[4], "w");

  while (F->eof() == false) {
    B = F->getSequence();

#ifdef DEBUG
    fprintf(stderr, "working on %s\n", B->header());
#endif

    char     *seq = (char *)B->sequence();

    u32bit    pos = 0;
    u32bit    max = 4 * (B->sequenceLength() / desiredLength + 1);
    u32bit   *sta = new u32bit [max];
    u32bit   *end = new u32bit [max];


    //  step 1: Build a list of regions to output.  we now just make a
    //  list of all length L segments, regardless of N content.  We'll
    //  fix them later.

    for (u32bit s = 0; s < B->sequenceLength() ;) {
      sta[pos]  = s;
      s        += desiredLength;

      if (s > B->sequenceLength())
        s = B->sequenceLength();

      end[pos]  = s;
      pos++;

      if (pos >= max) {
        fprintf(stderr, "pos >= max at line %d\n", __LINE__);
        exit(1);
      }
    }

#ifdef DEBUG
    fprintf(stderr, "step1: created %d regions\n", pos);
#endif

    if (pos > 0) {

    //  step 2: examine each region, trimming N's on the ends

    for (u32bit p=0; p<pos; p++) {
      while ((seq[sta[p]] == 'N') && (sta[p] < end[p]))
        sta[p]++;

      while ((seq[end[p]] == 'N') && (sta[p] < end[p]))
        end[p]--;
    }

    //  step 3: remove empty things

    u32bit p=0;
    u32bit q=0;

    for (p=0, q=0; p < pos; p++) {
      if (sta[p] != end[p]) {
        if (p != q) {
          sta[q] = sta[p];
          end[q] = end[p];
        }
        q++;
      }
    }
    pos = q;

#ifdef DEBUG
    fprintf(stderr, "step3: removed empty regions, left with %d regions\n", pos);
#endif

    if (pos > 0) {

    //  step 4: split things with "significant" amount N these should
    //  be things with contained blocks, or low-quality things

    u32bit    pos2 = 0;
    u32bit   *sta2 = new u32bit [max];
    u32bit   *end2 = new u32bit [max];

    for (p=0; p<pos; p++) {

      u32bit  numN = 0;
      for (u32bit s=sta[p]; s<end[p]; s++) {
        if (seq[s] == 'N')
          numN++;
      }

      //  If there are more N than allowed, see if we can split it.
      //
      //  (u32bit)(0.1 * (end[p] - sta[p]))
      //
      if (numN > 50) {
        u32bit s = sta[p];
        u32bit e = end[p] - 1;

        while (seq[s] != 'N')
          s++;

        while (seq[--e] != 'N')
          ;
        e++;

        //  s is on the first N after the start
        //  e is on the first non-N after the last N


        //  We want to save (0...sta[p]) and (e...end[p]) Probably
        //  don't want to keep the stuff in the middle, but that's not
        //  clear.

        sta2[pos2] = sta[p];
        end2[pos2] = s;
        pos2++;

        sta2[pos2] = e;
        end2[pos2] = end[p];
        pos2++;

#if 0
        fprintf(L, "discard %9d to %9d -- %7d N's (%8.5f%%)\n",
                s, e, numN, 100.0 * (double)numN / (double)(e-s));
        fflush(L);
#endif

      } else {
        sta2[pos2] = sta[p];
        end2[pos2] = end[p];
        pos2++;
      }

      if (pos2 >= max) {
        fprintf(stderr, "pos >= max at line %d\n", __LINE__);
        exit(1);
      }

    }

#ifdef DEBUG
    fprintf(stderr, "step4: split/removed significant N regions, left with %d regions\n", pos2);
#endif

    //  step 5: Merge any adjacent regions that are less than 2*L

    u32bit merged = 0;

    for (p=0, q=1; q<pos2;) {
      if ((end2[p] == sta2[q]) &&
          (end2[q] - sta2[p] < 2*desiredLength)) {
        // merge q into p
        //
        merged++;
        end2[p] = end2[q];
        q++;
      } else {
        // save q on the next position
        //
        p++;
        sta2[p] = sta2[q];
        end2[p] = end2[q];
        q++;
      }

#if 0
      fprintf(stdout, "p=%5d [%6d-%6d]  q=%5d [%6d-%6d]\n",
              p, sta2[p], end2[p], q, sta2[q], end2[q]);
#endif
    }

    delete [] sta;
    delete [] end;

    sta = sta2;
    end = end2;
    pos = p + 1;


#ifdef DEBUG
    fprintf(stderr, "step5: merged %d pairs of regions, left with %d regions\n", merged, pos);
#endif

    //  step 6: write regions

    char *buf = new char [desiredLength + desiredLength + 1];

    for (p=0; p<pos; p++) {
      memcpy(buf, seq + sta[p], sizeof(char) * (end[p] - sta[p]));
      buf[ end[p] - sta[p] ] = 0;
      fprintf(O, "%s sta="u32bitFMT" end="u32bitFMT"\n%s\n",
              B->header(), sta[p], end[p], buf);
      fprintf(L, u32bitFMT" : "u32bitFMT"["u32bitFMT"-"u32bitFMT"]\n",
              fragmentIndex++,
              B->getIID(),
              sta[p],
              end[p]);
    }

    delete [] buf;

    }
    }

    delete [] sta;
    delete [] end;

    delete B;
  }

  fclose(L);
  fclose(O);
}
