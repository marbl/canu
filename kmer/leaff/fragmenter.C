#include <stdio.h>
#include <stdlib.h>

#include "bio++.H"
#include "seqCache.H"

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
  uint32                desiredLength = 0;
  uint32                overlapLength = 0;
  bool                  beVerbose = false;
  seqCache             *F   = 0L;
  seqInCore            *B   = 0L;
  uint32                Bid = 0;
  FILE                 *O   = 0L;
  FILE                 *L   = 0L;

  uint32                fragmentIndex = 0;

  int arg=1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-length") == 0) {
      desiredLength = strtouint32(argv[++arg], 0L);
    } else if (strcmp(argv[arg], "-overlap") == 0) {
      overlapLength = strtouint32(argv[++arg], 0L);
    } else if (strcmp(argv[arg], "-input") == 0) {
      F = new seqCache(argv[++arg]);
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

  B = F->getSequenceInCore(Bid);
  while (B) {
    if (beVerbose)
      fprintf(stderr, "working on %s\n", B->header());

    char     *seq = (char *)B->sequence();

    uint32    pos = 0;
    uint32    max = 16384;
    uint32   *sta = new uint32 [max];
    uint32   *end = new uint32 [max];

    //  step 1: build a list of regions to output.  Scan the sequence,
    //  making a new region if we see a significant chunk of N, or if we
    //  hit the desiredLength.
    //
    uint32 s = 0;
    uint32 e = 0;
    while (s < B->sequenceLength()) {

      //  Skip any N at the start
      while ((seq[s] == 'n') || (seq[s] == 'N') && (s < B->sequenceLength()))
        s++;

      //  Construct the preliminary block.
      //
      e = s + desiredLength;
      if (e > B->sequenceLength())
        e = B->sequenceLength();

      fprintf(stderr, "got block1 "uint32FMT" - "uint32FMT"\n", s, e);

      //  Scan from s to e, looking for significant N.  If we find it,
      //  reset e and stop.
      //
      uint32  numN = 0;
      for (uint32 i=s; i<e; i++) {
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

      fprintf(stderr, "got block2 "uint32FMT" - "uint32FMT"\n", s, e);

      //  Back up e until we hit the first non-N
      if ((s < e) && ((seq[e] == 'n') || (seq[e] == 'N'))) {
        while ((s <= e) && ((seq[e] == 'n') || (seq[e] == 'N')))
          e--;
        e++;
      }

      fprintf(stderr, "got block3 "uint32FMT" - "uint32FMT"\n", s, e);

      //  Add this region
      //
      if (s > e) {
        fprintf(stderr, "ERROR!  s>e!  "uint32FMT" "uint32FMT"\n", s, e);
      }
      if (s != e) {
        fprintf(stderr, "ADD ["uint32FMTW(3)"] "uint32FMTW(9)" "uint32FMTW(9)" length "uint32FMTW(9)"\n", pos, s, e, e-s);
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
      for (uint32 p=1; p<pos; p++) {
        if (end[p-1] == sta[p]) {
          sta[p] -= overlapLength;
          fprintf(stderr, "ADJ ["uint32FMTW(3)"] "uint32FMTW(9)" "uint32FMTW(9)" length "uint32FMTW(9)"\n",
                  p, sta[p], end[p], end[p] - sta[p]);
        }
      }
    }


    if (beVerbose)
      fprintf(stderr, "created %d regions\n", pos+1);


    for (uint32 p=0; p<pos; p++) {

#if 1
      fprintf(O, "%s begin "uint32FMT" end "uint32FMT" length "uint32FMT"\n",
              B->header(), sta[p], end[p], end[p] - sta[p]);
      fwrite(seq+sta[p], sizeof(char), end[p] - sta[p], O);
      fprintf(O, "\n");
#endif

      fprintf(L, uint32FMT" : "uint32FMT"["uint32FMT"-"uint32FMT"]\n",
              fragmentIndex++,
              B->getIID(),
              sta[p],
              end[p]);
    }

    delete [] sta;
    delete [] end;

    delete B;
    B = F->getSequenceInCore(++Bid);
  }

  fclose(L);
  fclose(O);
}
