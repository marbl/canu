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


void
usage(char *name) {
  fprintf(stderr, "usage: %s [-overlap len] -length len -input X.fasta -output Y.fasta -log T.log\n",
          name);
  exit(1);
}


int
main(int argc, char **argv) {
  u32bit                desiredLength = 0;
  u32bit                overlapLength = 0;
  bool                  beVerbose = false;
  seqFile              *F = 0L;
  seqInCore            *B = 0L;
  FILE                 *O = 0L;
  FILE                 *L = 0L;

  u32bit                fragmentIndex = 0;

  int arg=1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-length") == 0) {
      desiredLength = strtou32bit(argv[++arg], 0L);
    } else if (strcmp(argv[arg], "-overlap") == 0) {
      overlapLength = strtou32bit(argv[++arg], 0L);
    } else if (strcmp(argv[arg], "-input") == 0) {
      F = openSeqFile(argv[++arg]);
    } else if (strcmp(argv[arg], "-output") == 0) {
      errno = 0;
      O = fopen(argv[++arg], "w");
      if (errno)
        fprintf(stderr, "ERROR:  Can't open output file '%s': %s\n", argv[arg], strerror(errno)), exit(1);
    } else if (strcmp(argv[arg], "-log") == 0) {
      errno = 0;
      L = fopen(argv[++arg], "w");
      if (errno)
        fprintf(stderr, "ERROR:  Can't open log file '%s': %s\n", argv[arg], strerror(errno)), exit(1);
    } else if (strcmp(argv[arg], "-verbose") == 0) {
      beVerbose = true;
    } else {
      usage(argv[arg]);
    }

    arg++;
  }

  if ((F == 0L) || (O == 0L) || (L == 0L))
    usage(argv[0]);

  B = F->getSequenceInCore();
  while (B) {
    if (beVerbose)
      fprintf(stderr, "working on %s\n", B->header());

    char     *seq = (char *)B->sequence();

    u32bit    pos = 0;
    u32bit    max = 16384;
    u32bit   *sta = new u32bit [max];
    u32bit   *end = new u32bit [max];

    //  step 1: build a list of regions to output.  Scan the sequence,
    //  making a new region if we see a significant chunk of N, or if we
    //  hit the desiredLength.
    //
    u32bit s = 0;
    u32bit e = 0;
    while (s < B->sequenceLength()) {

      //  Skip any N at the start
      while ((seq[s] == 'n') || (seq[s] == 'N') && (s < B->sequenceLength()))
        s++;

      //  Construct the preliminary block.
      //
      e = s + desiredLength;
      if (e > B->sequenceLength())
        e = B->sequenceLength();

      fprintf(stderr, "got block1 "u32bitFMT" - "u32bitFMT"\n", s, e);

      //  Scan from s to e, looking for significant N.  If we find it,
      //  reset e and stop.
      //
      u32bit  numN = 0;
      for (u32bit i=s; i<e; i++) {
        if ((seq[i] == 'n') || (seq[i] == 'N')) {
          numN++;
        } else {
          numN = 0;
        }
        if (numN >= 50) {
          e = i;
          break;
        }
      }

      fprintf(stderr, "got block2 "u32bitFMT" - "u32bitFMT"\n", s, e);

      //  Back up e until we hit the first non-N
      if ((s < e) && ((seq[e] == 'n') || (seq[e] == 'N'))) {
        while ((s <= e) && ((seq[e] == 'n') || (seq[e] == 'N')))
          e--;
        e++;
      }

      fprintf(stderr, "got block3 "u32bitFMT" - "u32bitFMT"\n", s, e);

      //  Add this region
      //
      if (s > e) {
        fprintf(stderr, "ERROR!  s>e!  "u32bitFMT" "u32bitFMT"\n", s, e);
      }
      if (s != e) {
        fprintf(stderr, "ADD ["u32bitFMTW(3)"] "u32bitFMTW(9)" "u32bitFMTW(9)" length "u32bitFMTW(9)"\n", pos, s, e, e-s);
        sta[pos] = s;
        end[pos] = e;
        pos++;
        if (pos >= max) {
          fprintf(stderr, "ERROR!  max exceeded!\n");
        }
      }

      s = e;
    }


    //  If we're supposed to be overlapping, fiddle with the begin position to make it so.
    //
    if (overlapLength > 0) {
      for (u32bit p=1; p<pos; p++) {
        if (end[p-1] == sta[p]) {
          sta[p] -= overlapLength;
          fprintf(stderr, "ADJ ["u32bitFMTW(3)"] "u32bitFMTW(9)" "u32bitFMTW(9)" length "u32bitFMTW(9)"\n",
                  p, sta[p], end[p], end[p] - sta[p]);
        }
      }
    }


    if (beVerbose)
      fprintf(stderr, "created %d regions\n", pos+1);


    for (u32bit p=0; p<pos; p++) {

#if 1
      fprintf(O, "%s begin "u32bitFMT" end "u32bitFMT" length "u32bitFMT"\n",
              B->header(), sta[p], end[p], end[p] - sta[p]);
      fwrite(seq+sta[p], sizeof(char), end[p] - sta[p], O);
      fprintf(O, "\n");
#endif

      fprintf(L, u32bitFMT" : "u32bitFMT"["u32bitFMT"-"u32bitFMT"]\n",
              fragmentIndex++,
              B->getIID(),
              sta[p],
              end[p]);
    }

    delete [] sta;
    delete [] end;

    delete B;
    B = F->getSequenceInCore();
  }

  fclose(L);
  fclose(O);
}
