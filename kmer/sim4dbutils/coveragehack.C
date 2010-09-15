#include <stdio.h>
#include <stdlib.h>
#include "util++.H"
#include "bio++.H"
#include "sim4.H"

//  g++ -o coveragehack coveragehack.C -I../libutil -I../libbio -I../libsim4 -L../libutil -L../libbio -L../libsim4 -lsim4 -lbio -lutil

//  Flag that tells which side of the alignment our contaminated assembly is on.
//    1 (R) -- if atac     the contaminant is on the left, the assembly is on the right
//
u32bit  orientation = 1;

//
//  WARNING!  This is stale code.  It does not compile.  The fasta interface has changed.
//

void
readATAC(intervalList **coverage, char *path) {
  char             line[1024] = {0};
  splitToWords     S(line);

  errno = 0;
  FILE *F = fopen(path, "r");
  if (errno)
    fprintf(stderr, "Failed to open '%s': %s\n", path, strerror(errno)), exit(1);

  while (!feof(F)) {
    fgets(line, 1024, F);

    if ((line[0] == 'M') && (line[2] == 'u')) {
      S.split(line);

      u32bit  taglength = 0;
      while (S[8][taglength] != ':')
        taglength++;
      u32bit  idx = atoi(S[8] + taglength + 1);
      u32bit  beg = atoi(S[9]);
      u32bit  len = atoi(S[10]);

      if (orientation == 2) {
        while (S[4][taglength] != ':')
          taglength++;
        idx = atoi(S[4] + taglength + 1);
        beg = atoi(S[5]);
        len = atoi(S[6]);
      }

      if (coverage[idx] == 0L)
        coverage[idx] = new intervalList();

      coverage[idx]->add(beg, len);
    }
  }

  fclose(F);
}


void
readSIM4(intervalList **coverage, int which, char *path) {

  errno = 0;
  FILE *F = fopen(path, "r");
  if (errno)
    fprintf(stderr, "Failed to open '%s': %s\n", path, strerror(errno)), exit(1);

  while (!feof(F)) {
    sim4polish *p = new sim4polish(F);

    if (p) {

      switch (which) {
        case 1:
          //  The query are contaminant reads, the genomic is the assembly
          if ((p->_percentIdentity >= 94) && (p->_querySeqIdentity >= 80)) {
            u32bit  idx = p->_genID;

            if (coverage[idx] == 0L)
              coverage[idx] = new intervalList();

            coverage[idx]->add(p->_exons[0]._genFrom,
                               p->_exons[0]._genTo - p->_exons[0]._genFrom + 1);
          }
          break;
        case 2:
          //  The query are assembly scaffolds, the genomic is the contaminant assembly (one or a few contigs)
          //
          u32bit  idx = p->_estID;

          if (coverage[idx] == 0L)
            coverage[idx] = new intervalList();

          if (p->_matchOrientation == SIM4_MATCH_FORWARD) {
            coverage[idx]->add(p->_exons[0]._estFrom,
                               p->_exons[0]._estTo - p->_exons[0]._estFrom + 1);
          } else {
            coverage[idx]->add(p->_estLen - p->_exons[0]._estTo + 1,
                               p->_exons[0]._estTo - p->_exons[0]._estFrom + 1);
          }
          break;
      }

      delete p;
    }
  }
  fclose(F);
}


#define MAXSCAFFOLD 200000

int
main(int argc, char **argv) {
  intervalList    **coverage = new intervalList* [MAXSCAFFOLD];
  intervalList    **gaps     = new intervalList* [MAXSCAFFOLD];
  FastAWrapper     *W        = 0L;
  u32bit            minCov   = 80;

  bool              includeGapsAsContamination = true;

  for (u32bit i=0; i<MAXSCAFFOLD; i++)
    coverage[i] = 0L;

  int arg=1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-a") == 0) {
      readATAC(coverage, argv[++arg]);
    } else if (strcmp(argv[arg], "-s") == 0) {
      //  SNAPPER of query fragments onto the contaminant.
      readSIM4(coverage, 1, argv[++arg]);
    } else if (strcmp(argv[arg], "-S") == 0) {
      //  SNAPPER of query contaminant onto scaffolds
      readSIM4(coverage, 2, argv[++arg]);
    } else if (strcmp(argv[arg], "-f") == 0) {
      W = new seqCache(argv[++arg]);
    } else if (strcmp(argv[arg], "-c") == 0) {
      minCov = strtou32bit(argv[++arg], 0L);
    } else if (strcmp(argv[arg], "-r") == 0) {
      orientation = 1;
    } else if (strcmp(argv[arg], "-l") == 0) {
      orientation = 2;
    } else if (strcmp(argv[arg], "-g") == 0) {
      includeGapsAsContamination = false;
    } else {
      fprintf(stderr, "%s: unknown arg %s\n", argv[0], argv[arg]);
    }
    arg++;
  }

  if (W == 0L) {
    fprintf(stderr, "usage: %s [-a atacmapping] [-s sim4db] [-g] -f seq.fasta\n", argv[0]);
    fprintf(stderr, "  -g   don't count gaps in scaffolds as contamination.\n");
    exit(1);
  }

  u32bit  sumOfLengths = 0;
  u32bit  sequences   = 0;

  for (u32bit i=0; i<MAXSCAFFOLD; i++) {
    if (coverage[i]) {

      W->find(i);
      FastASequenceInCore  *S = W->getSequence();

      intervalList          gaps;

      //  Compute how much of the scaffold is gap.

      u32bit  gapBeg = W->sequenceLength(i);
      char   *seq    = S->sequence();

      for (u32bit beg=0, len=W->sequenceLength(i); beg<len; beg++) {
        if ((seq[beg] == 'N') || (seq[beg] == 'n')) {
          if (gapBeg > beg)
            gapBeg = beg;
        } else {
          if (gapBeg < beg) {
            gaps.add(gapBeg, beg-gapBeg);
            gapBeg = W->sequenceLength(i);
          }
        }
      }

      //  Geez!  I suppose we could have just directly counted ACGT above!

      gaps.merge();
      coverage[i]->merge();

      u32bit   coveredLength = coverage[i]->sumOfLengths();
      u32bit   gapLength     = gaps.sumOfLengths();
      u32bit   totalLength   = W->sequenceLength(i) - gapLength;

      if (100 * coveredLength > minCov * totalLength) {

        sumOfLengths += coveredLength;
        sequences++;

        double cov = 100.0 * coveredLength / (double)totalLength;

        fprintf(stderr, "sequence ["u32bitFMT"] %s covered "u32bitFMT" out of "u32bitFMT" (%7.3f)\n",
                i,
                S->header(),
                coveredLength,
                totalLength,
                cov);

        delete S;
      }

      //  Dump a special scaffold
      if (i == 4796) {
        for (u32bit z=0; z<coverage[i]->numberOfIntervals(); z++) {
          fprintf(stderr, "interval[%3d] %6d - %6d\n", z, coverage[i]->lo(z), coverage[i]->hi(z));
        }

      }

    }
  }

  fprintf(stderr, "Found "u32bitFMT" bases in "u32bitFMT" scaffolds.\n", sumOfLengths, sequences);
}

