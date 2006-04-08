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
                     atacMatch           *l,
                     atacMatch           *r) {
  char  *s = S->sequence();

  //  Over the whole sequence
  //
  while (beg < end) {

    //  Skip any N's starting where we are currently
    //
    while ((beg < end) &&
           ((s[beg] == 'N') || (s[beg] == 'n')))
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
             ((s[cur] != 'N') && (s[cur] != 'n')))
        cur++;

      //  And output whatever this block is
      //
      //writeSequence(output, S->header(), S->getIID(), beg, cur,
      //              S->sequence() + beg);

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



void
extractUnmapped1(FILE *Aoutput, atacMatchList &ML) {
  FastASequenceInCore  *S = 0L;
  FastAWrapper         *W = ML._seq1;

  ML.sort1();

  W->find(ML[0]->iid1);
  S = W->getSequence();
  for (u32bit i=1; i<ML.numMatches(); i++) {
    atacMatch *l = ML[i-1];
    atacMatch *r = ML[i];

    if (l->iid1 != r->iid1)
      continue;

    if (l->iid1 != S->getIID()) {
      delete S;
      W->find(l->iid1);
      S = W->getSequence();
    }

    //  Extract from (l->pos1 + l->len1) to (r->pos1), if it's longer than 20bp
    //
    if (l->pos1 + l->len1 + 20 < r->pos1)
      writeGaplessSequence(Aoutput,
                           S,
                           l->pos1 + l->len1,
                           r->pos1,
                           l, r);
  }
}


void
extractUnmapped2(FILE *Boutput, atacMatchList &ML) {
  FastASequenceInCore  *S = 0L;
  FastAWrapper         *W = ML._seq2;

  ML.sort2();

  W->find(ML[0]->iid2);
  S = W->getSequence();
  for (u32bit i=1; i<ML.numMatches(); i++) {
    atacMatch *l = ML[i-1];
    atacMatch *r = ML[i];

    if (l->iid2 != r->iid2)
      continue;

    if (l->iid2 != S->getIID()) {
      delete S;
      W->find(l->iid2);
      S = W->getSequence();
    }

    //  Extract from (l->pos2 + l->len2) to (r->pos2), if it's longer than 20bp
    //
    if (l->pos2 + l->len2 + 20 < r->pos2)
      writeGaplessSequence(Boutput,
                           S,
                           l->pos2 + l->len2,
                           r->pos2,
                           l, r);
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
extractUnmapped(FILE *Aoutput, FILE *Boutput, atacMatchList &ML) {
  u32bit   numSeqsA = ML._seq1->getNumberOfSequences();
  u32bit   numSeqsB = ML._seq2->getNumberOfSequences();

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

    ML._seq1->find(seq);
    FastASequenceInCore  *S = ML._seq1->getSequence();

    if (coveredA[seq].len() == 0) {
#if 1
      writeGaplessSequence(Aoutput,
                           S,
                           0,
                           ML._seq1->sequenceLength(seq),
                           0L, 0L);
#endif
    } else {
      if (0 < coveredA[seq][0]->pos1) {
        writeGaplessSequence(Aoutput,
                             S,
                             0,
                             coveredA[seq][0]->pos1,
                             0L, coveredA[seq][0]);
      }

      for (u32bit i=1; i<coveredA[seq].len(); i++) {
        if (coveredA[seq][i-1]->pos1 + coveredA[seq][i-1]->len1 < coveredA[seq][i]->pos1) {
          writeGaplessSequence(Aoutput,
                               S,
                               coveredA[seq][i-1]->pos1 + coveredA[seq][i-1]->len1,
                               coveredA[seq][i]->pos1,
                               coveredA[seq][i-1], coveredA[seq][i]);
        }
      }

      u32bit last = coveredA[seq].len()-1;
      if (coveredA[seq][last]->pos1) {
        writeGaplessSequence(Aoutput,
                             S,
                             coveredA[seq][last]->pos1 + coveredA[seq][last]->len1,
                             ML._seq1->sequenceLength(seq),
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

    ML._seq2->find(seq);
    FastASequenceInCore  *S = ML._seq2->getSequence();

    if (coveredB[seq].len() == 0) {
#if 1
      writeGaplessSequence(Boutput,
                           S,
                           0,
                           ML._seq2->sequenceLength(seq),
                           0L, 0L);
#endif
    } else {
      if (0 < coveredB[seq][0]->pos2) {
        writeGaplessSequence(Boutput,
                             S,
                             0,
                             coveredB[seq][0]->pos2,
                             0L, coveredB[seq][0]);
      }

      for (u32bit i=1; i<coveredB[seq].len(); i++) {
        if (coveredB[seq][i-1]->pos2 + coveredB[seq][i-1]->len2 < coveredB[seq][i]->pos2) {
          writeGaplessSequence(Boutput,
                               S,
                               coveredB[seq][i-1]->pos2 + coveredB[seq][i-1]->len2,
                               coveredB[seq][i]->pos2,
                               coveredB[seq][i-1], coveredB[seq][i]);
        }
      }

      u32bit last = coveredB[seq].len()-1;
      if (coveredB[seq][last]->pos2) {
        writeGaplessSequence(Boutput,
                             S,
                             coveredB[seq][last]->pos2 + coveredB[seq][last]->len2,
                             ML._seq2->sequenceLength(seq),
                             coveredB[seq][0], 0L);
      }
    }
  }









}


void
extractUnmappedRuns(FILE *ARoutput, FILE *BRoutput, atacMatchList &ML) {
  FastASequenceInCore  *S1 = 0L;
  FastAWrapper         *W1 = ML._seq1;
  FastASequenceInCore  *S2 = 0L;
  FastAWrapper         *W2 = ML._seq2;

  ML.sort1();

  W1->find(ML[0]->iid1);
  S1 = W1->getSequence();
  W2->find(ML[0]->iid2);
  S2 = W2->getSequence();

  for (u32bit i=1; i<ML.numMatches(); i++) {
    atacMatch *l = ML[i-1];
    atacMatch *r = ML[i];

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

      writeGaplessSequence(ARoutput,
                           S1,
                           l->pos1 + l->len1,
                           r->pos1,
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
                             l, r);
      } else {
        writeGaplessSequence(BRoutput,
                             S2,
                             r->pos2 + r->len2,
                             l->pos2,
                             l, r);
      }
    }
  }
}








void
usage(char *name) {
  fprintf(stderr, "usage: %s [-OP output.fasta] -m matches\n", name);
  fprintf(stderr, "   OP\n");
  fprintf(stderr, "   -a        extract all unmapped sequence in A\n");
  fprintf(stderr, "   -b        extract all unmapped sequence in B\n");
  fprintf(stderr, "   -ar       extract within run unmapped sequence in A\n");
  fprintf(stderr, "   -br       extract within run unmapped sequence in B\n");
  fprintf(stderr, "             BOTH -ar and -br need to be specified!\n");
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
    } else {
      usage(argv[0]);
      exit(1);
    }
    arg++;
  }

  if (matchesFile == 0L)
    usage(argv[0]), exit(1);

  atacMatchList  ML(matchesFile, 'm', false);

#if 0
  if (Aoutput) {
    extractUnmapped1(Aoutput, ML);
    fclose(Aoutput);
  }

  if (Boutput) {
    extractUnmapped2(Boutput, ML);
    fclose(Boutput);
  }
#else
  if (Aoutput && Boutput) {
    extractUnmapped(Aoutput, Boutput, ML);
    fclose(Aoutput);
    fclose(Boutput);
  }
#endif

  if (ARoutput && BRoutput) {
    extractUnmappedRuns(ARoutput, BRoutput, ML);
    fclose(ARoutput);
    fclose(BRoutput);
  }

  return(0);
}
