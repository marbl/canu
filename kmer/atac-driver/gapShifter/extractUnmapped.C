// This file is part of A2Amapper.
// Copyright (c) 2006 J. Craig Venter Institute
// Author: Brian Walenz
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received (LICENSE.txt) a copy of the GNU General Public 
// License along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "bio++.H"
#include "atac.H"


//  Reads a set of matches and outputs two sequence files containing sequence that
//  is not matched.


void
writeGaplessSequence(FILE                *output,
                     FastASequenceInCore *S,
                     u32bit               beg,
                     u32bit               end,
                     u32bit               extend,
                     atacMatch           *l,
                     atacMatch           *r) {
  char  *s = S->sequence();


  //  Skip any N's starting where we are currently
  //
  while ((beg < end) &&
         (toUpper[s[beg]] == 'N'))
    beg++;

  while ((beg < end) &&
         (toUpper[s[end-1]] == 'N'))
    end--;

  if (beg >= end)
    return;

  //  Extend the ends up to 'extend' positions, as long as we don't
  //  hit a gap.
  //
  for (u32bit x=0; ((x < extend) &&
                    (beg > 0) &&
                    (toUpper[s[beg-1]] != 'N')); x++)
    beg--;

  for (u32bit x=0; ((x < extend) &&
                    (end < S->sequenceLength()) &&
                    (toUpper[s[end]] != 'N')); x++)
    end++;

  //  Just make sure we're still in bounds!
  if (end > S->sequenceLength())
    end = S->sequenceLength();

  //  Over the whole sequence
  //
  while (beg < end) {

    //  Skip any N's starting where we are currently
    //
    while ((beg < end) &&
           (toUpper[s[beg]] == 'N'))
      beg++;
        
    //  Move our current up to here
    u32bit cur = beg;

    //  If we're at the end of the sequence, this block doesn't
    //  exist; it's solid N.
    //
    if (beg < end) {

      //  Move cur up to the next N
      //
      while ((cur < end) &&
             (toUpper[s[cur]] != 'N'))
        cur++;

      //  And output whatever this block is
      //
      char  *lmuid = "none";
      char  *lpuid = "none";
      char  *rmuid = "none";
      char  *rpuid = "none";

      if (l) {
        lmuid = l->matchuid;
        lpuid = l->parentuid;
      }
      if (r) {
        rmuid = r->matchuid;
        rpuid = r->parentuid;
      }

      fprintf(output, "%s extracted from iid "u32bitFMT" pos "u32bitFMT" "u32bitFMT" between match %s(%s) and %s(%s)\n",
              S->header(), S->getIID(), beg, cur,
              lmuid, lpuid, rmuid, rpuid);

      fwrite(S->sequence() + beg, sizeof(char), cur-beg, output);
      fprintf(output, "\n");

    }

    //  Move to the next block.
    beg = cur;
  }
}





//  COPIED from libatac/matchList.C, but need to dereference the atacMatch again.
//
static
int
sort1_(const void *a, const void *b) {
  const atacMatch *A = *(const atacMatch **)a;
  const atacMatch *B = *(const atacMatch **)b;

  if (A->iid1 < B->iid1)  return(-1);
  if (A->iid1 > B->iid1)  return(1);
  if (A->pos1 < B->pos1)  return(-1);
  if (A->pos1 > B->pos1)  return(1);
  if (A->len1 > B->len1)  return(-1);
  if (A->len1 < B->len1)  return(1);
  if (A->iid2 < B->iid2)  return(-1);
  if (A->iid2 > B->iid2)  return(1);
  if (A->pos2 < B->pos2)  return(-1);
  if (A->pos2 > B->pos2)  return(1);
  if (A->len2 > B->len2)  return(-1);
  if (A->len2 < B->len2)  return(1);
  return(0);
}

static
int
sort2_(const void *a, const void *b) {
  const atacMatch *A = *(const atacMatch **)a;
  const atacMatch *B = *(const atacMatch **)b;

  if (A->iid2 < B->iid2)  return(-1);
  if (A->iid2 > B->iid2)  return(1);
  if (A->pos2 < B->pos2)  return(-1);
  if (A->pos2 > B->pos2)  return(1);
  if (A->len2 > B->len2)  return(-1);
  if (A->len2 < B->len2)  return(1);
  if (A->iid1 < B->iid1)  return(-1);
  if (A->iid1 > B->iid1)  return(1);
  if (A->pos1 < B->pos1)  return(-1);
  if (A->pos1 > B->pos1)  return(1);
  if (A->len1 > B->len1)  return(-1);
  if (A->len1 < B->len1)  return(1);

  return(0);
}



//  New method, uses an intervalList to find the unmapped regions for
//  each sequence.
//
class  extractMatchList {

public:
  extractMatchList() {
    matchesLen = 0;
    matchesMax = 16;
    matches    = new atacMatch * [matchesMax];
  };
  ~extractMatchList() {
    delete [] matches;
  };

  atacMatch *operator[](u32bit idx) {
    return(matches[idx]);
  };

  u32bit len(void) {
    return(matchesLen);
  };

  void   add(atacMatch *m) {
    if (matchesLen >= matchesMax) {
      matchesMax *= 2;
      atacMatch **M = new atacMatch * [matchesMax];
      memcpy(M, matches, sizeof(atacMatch *) * matchesLen);
      delete [] matches;
      matches = M;
    }
    matches[matchesLen++] = m;
  };

  void   sort1(void) {
    qsort(matches, matchesLen, sizeof(atacMatch*), sort1_);
  };
  void   sort2(void) {
    qsort(matches, matchesLen, sizeof(atacMatch*), sort2_);
  };

private:
  atacMatch  **matches;
  u32bit       matchesLen;
  u32bit       matchesMax;
};




void
extractUnmapped(FastACache *A, FastACache *B,
                FILE *Aoutput, FILE *Boutput,
                u32bit extend, atacMatchList &ML) {
  u32bit   numSeqsA = ML.fastaA()->getNumberOfSequences();
  u32bit   numSeqsB = ML.fastaB()->getNumberOfSequences();

  extractMatchList  *coveredA = new extractMatchList [numSeqsA];
  extractMatchList  *coveredB = new extractMatchList [numSeqsB];

  //  Populate the intervals with the mapping
  //
  for (u32bit x=0; x<ML.numMatches(); x++) {
    atacMatch *m = ML[x];

    coveredA[m->iid1].add(m);
    coveredB[m->iid2].add(m);
  }

  //  Sort the intervals, manually invert the interval -- remembering
  //  what matches are where.
  //
  for (u32bit seq=0; seq<numSeqsA; seq++) {
    coveredA[seq].sort1();

    //ML.fastaA()->find(seq);
    //FastASequenceInCore  *S = ML.fastaA()->getSequence();

    FastASequenceInCore *S = A->getSequence(seq);

    if (coveredA[seq].len() == 0) {
      //  Hey!  This sequence has NO matches at all!
      //
      writeGaplessSequence(Aoutput,
                           S,
                           0,
                           ML.fastaA()->sequenceLength(seq),
                           extend,
                           0L, 0L);
    } else {
      if (0 < coveredA[seq][0]->pos1) {
        writeGaplessSequence(Aoutput,
                             S,
                             0,
                             coveredA[seq][0]->pos1,
                           extend,
                             0L, coveredA[seq][0]);
      }

      for (u32bit i=1; i<coveredA[seq].len(); i++) {
        if (coveredA[seq][i-1]->pos1 + coveredA[seq][i-1]->len1 < coveredA[seq][i]->pos1) {
          writeGaplessSequence(Aoutput,
                               S,
                               coveredA[seq][i-1]->pos1 + coveredA[seq][i-1]->len1,
                               coveredA[seq][i]->pos1,
                           extend,
                               coveredA[seq][i-1], coveredA[seq][i]);
        }
      }

      u32bit last = coveredA[seq].len()-1;
      if (coveredA[seq][last]->pos1) {
        writeGaplessSequence(Aoutput,
                             S,
                             coveredA[seq][last]->pos1 + coveredA[seq][last]->len1,
                             ML.fastaA()->sequenceLength(seq),
                           extend,
                             coveredA[seq][0], 0L);
      }
    }
  }



  //  DUPLICATION OF THE ABOVE!  (Replace 1 with 2, A with B)


  //  Sort the intervals, manually invert the interval -- remembering
  //  what matches are where.
  //
  for (u32bit seq=0; seq<numSeqsB; seq++) {
    coveredB[seq].sort2();

    //ML.fastaB()->find(seq);
    //FastASequenceInCore  *S = ML.fastaB()->getSequence();

    FastASequenceInCore *S = B->getSequence(seq);

    if (coveredB[seq].len() == 0) {
      //  Hey!  This sequence has NO matches at all!
      //
      writeGaplessSequence(Boutput,
                           S,
                           0,
                           ML.fastaB()->sequenceLength(seq),
                           extend,
                           0L, 0L);
    } else {
      if (0 < coveredB[seq][0]->pos2) {
        writeGaplessSequence(Boutput,
                             S,
                             0,
                             coveredB[seq][0]->pos2,
                           extend,
                             0L, coveredB[seq][0]);
      }

      for (u32bit i=1; i<coveredB[seq].len(); i++) {
        if (coveredB[seq][i-1]->pos2 + coveredB[seq][i-1]->len2 < coveredB[seq][i]->pos2) {
          writeGaplessSequence(Boutput,
                               S,
                               coveredB[seq][i-1]->pos2 + coveredB[seq][i-1]->len2,
                               coveredB[seq][i]->pos2,
                           extend,
                               coveredB[seq][i-1], coveredB[seq][i]);
        }
      }

      u32bit last = coveredB[seq].len()-1;
      if (coveredB[seq][last]->pos2) {
        writeGaplessSequence(Boutput,
                             S,
                             coveredB[seq][last]->pos2 + coveredB[seq][last]->len2,
                             ML.fastaB()->sequenceLength(seq),
                           extend,
                             coveredB[seq][0], 0L);
      }
    }
  }









}


void
extractUnmappedRuns(FastACache *A, FastACache *B,
                    FILE *ARoutput, FILE *BRoutput,
                    u32bit extend, atacMatchList &ML) {
  FastASequenceInCore  *S1 = 0L;
  FastASequenceInCore  *S2 = 0L;

  atacMatchOrder  MO(ML);
  MO.sortA();

  for (u32bit i=1; i<MO.numMatches(); i++) {
    atacMatch *l = MO[i-1];
    atacMatch *r = MO[i];

    if (l->iid1 != r->iid1)
      continue;
    if (l->iid2 != r->iid2)
      continue;

    //  Extract from (l->pos1 + l->len1) to (r->pos1), if it's longer than 20bp

    bool  lengthOK = true;
    if (l->pos1 + l->len1 + 20 >= r->pos1)
      lengthOK = false;
    if ((l->fwd2 == true) && (l->pos2 + l->len2 + 20 >= r->pos2))
      lengthOK = false;
    if ((l->fwd2 == false) && (r->pos2 + r->len2 + 20 >= l->pos2))
      lengthOK = false;

    //  Extract if our two matches are in the same run.
    //
    if ((lengthOK) &&
        (strcmp(l->parentuid, r->parentuid) == 0)) {

#if 0
      if (l->iid1 != S1->getIID()) {
        delete S1;
        W1->find(l->iid1);
        S1 = W1->getSequence();
      }

      if (l->iid2 != S2->getIID()) {
        delete S2;
        W2->find(l->iid2);
        S2 = W2->getSequence();
      }
#else
      S1 = A->getSequence(l->iid1);
      S2 = B->getSequence(l->iid2);
#endif

      writeGaplessSequence(ARoutput,
                           S1,
                           l->pos1 + l->len1,
                           r->pos1,
                           extend,
                           l, r);

      //  Need to deal with reverse matches here!  In run matches
      //  should be the same orientation, but we'll still check.
      //      
      if (l->fwd2 != r->fwd2) {
        fprintf(stderr, "WOAH!  Matches of different orientation in a run?!?\n");
        exit(1);
      }

      if (l->fwd2) {
        writeGaplessSequence(BRoutput,
                             S2,
                             l->pos2 + l->len2,
                             r->pos2,
                           extend,
                             l, r);
      } else {
        writeGaplessSequence(BRoutput,
                             S2,
                             r->pos2 + r->len2,
                             l->pos2,
                           extend,
                             l, r);
      }
    }
  }
}








void
usage(char *name) {
  fprintf(stderr, "usage: %s [-OP output.fasta] [-t trfile] -m matches\n", name);
  fprintf(stderr, "   OP\n");
  fprintf(stderr, "   -a        extract all unmapped sequence in A\n");
  fprintf(stderr, "   -b        extract all unmapped sequence in B\n");
  fprintf(stderr, "   -ar       extract within run unmapped sequence in A\n");
  fprintf(stderr, "   -br       extract within run unmapped sequence in B\n");
  fprintf(stderr, "             BOTH -ar and -br need to be specified!\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "   -t        mask out tandem repeats listed in trfile\n");
}

FILE *
openOutputFile(char *name) {
  errno = 0;
  FILE *R = fopen(name, "w");
  if (errno)
    fprintf(stderr, "Failed to open '%s': %s\n", name, strerror(errno)), exit(1);
  return(R);
}

int
main(int argc, char *argv[]) {
  char         *matchesFile = 0L;
  FILE         *Aoutput = 0L;
  FILE         *Boutput = 0L;
  FILE         *ARoutput = 0L;
  FILE         *BRoutput = 0L;
  u32bit        extend = 0;
  char         *trFile = 0L;

  int arg=1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-m") == 0) {
      matchesFile = argv[++arg];
    } else if (strcmp(argv[arg], "-a") == 0) {
      Aoutput = openOutputFile(argv[++arg]);
    } else if (strcmp(argv[arg], "-b") == 0) {
      Boutput = openOutputFile(argv[++arg]);
    } else if (strcmp(argv[arg], "-ar") == 0) {
      ARoutput = openOutputFile(argv[++arg]);
    } else if (strcmp(argv[arg], "-br") == 0) {
      BRoutput = openOutputFile(argv[++arg]);
    } else if (strcmp(argv[arg], "-e") == 0) {
      extend = strtou32bit(argv[++arg], 0L);
    } else if (strcmp(argv[arg], "-t") == 0) {
      trFile = argv[++arg];
    } else {
      usage(argv[0]);
      exit(1);
    }
    arg++;
  }

  if (matchesFile == 0L)
    usage(argv[0]), exit(1);

  atacMatchList  ML(matchesFile, 'm');

  //  Build caches for both sequences, then modify that sequence to
  //  mask out tandem repeats.
  //
  FastACache  *A = new FastACache(ML.assemblyFileA(), 0, true, true);
  FastACache  *B = new FastACache(ML.assemblyFileB(), 0, true, true);

  if (trFile) {
    errno =0;
    FILE         *F = fopen(trFile, "r");
    if (errno)
      fprintf(stderr, "Error opening '%s': %s\n", trFile, strerror(errno));

    char          L[1024] = { 0 };
    splitToWords  W(L);

    fprintf(stderr, "Masking repeats in '%s'\n", trFile);

    u32bit  statidx = 0;
    u32bit  stats[2] = { 0 };

    while (!feof(F)) {
      fgets(L, 1024, F);
      W.split(L);

      char   source = W[0][0];
      u32bit iid    = strtou32bit(W[1], 0L);
      u32bit pos    = strtou32bit(W[2], 0L);
      u32bit len    = strtou32bit(W[3], 0L);
      bool   fwd    = (W[4][0] != '-');

      FastASequenceInCore  *S = 0L;
      char                 *s = 0L;

      if (source == 'B') {
        S = A->getSequence(iid);
        s = A->getSequence(iid)->sequence();
        statidx = 0;
      } else if (source == 'H') {
        S = B->getSequence(iid);
        s = B->getSequence(iid)->sequence();
        statidx = 1;
      } else {
        fprintf(stderr, "Unknown source '%c'\n", source);
        exit(1);
      }

      //fprintf(stderr, "Masking %c "u32bitFMTW(8)" from "u32bitFMTW(9)" to "u32bitFMTW(9)" on strand %c\r",
      //        source, iid, pos, pos+len, (fwd) ? 'f' : 'r');

      if (fwd) {
        s += pos;
      } else {
        s += S->sequenceLength() - pos - len;
      }

      for (u32bit i=0; i<len; i++) {
        if (toUpper[s[i]] != 'N')
          stats[statidx]++;
        s[i] = 'N';
      }
    }
    fclose(F);

    fprintf(stderr, "Done masking.  "u32bitFMT" in A, "u32bitFMT" in B.\n", stats[0], stats[1]);
  }


  if (Aoutput && Boutput) {
    extractUnmapped(A, B, Aoutput, Boutput, extend, ML);
    fclose(Aoutput);
    fclose(Boutput);
  }

  if (ARoutput && BRoutput) {
    extractUnmappedRuns(A, B, ARoutput, BRoutput, extend, ML);
    fclose(ARoutput);
    fclose(BRoutput);
  }

  return(0);
}
