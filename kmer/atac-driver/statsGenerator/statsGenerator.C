// This file is part of A2Amapper.
// Copyright (c) 2005, 2006 J. Craig Venter Institute
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

//  Compute some simple statistics on a set of matches

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include "atac.H"
#include "util++.H"
#include "bio++.H"
#include "seqCache.H"

bool  noHistogramPlots = true;


//  Sort uint32 backwards
int
uint32compare(const void *a, const void *b) {
  const uint32 A = *((const uint32 *)a);
  const uint32 B = *((const uint32 *)b);
  if (A < B) return(1);
  if (A > B) return(-1);
  return(0);
}


class histogram {
public:
  histogram(uint64 blockSize, uint64 maxSize) {
    _b = blockSize;
    _m = maxSize;
    _l = 0;
    _h = new uint32 [maxSize / blockSize + 1];
    for (uint32 i=0; i<maxSize / blockSize + 1; i++)
      _h[i] = 0;
    _eLen = 0;
    _eMax = 10240;
    _e = new uint32 [_eMax];
  };
  ~histogram() {
    delete [] _h;
    delete [] _e;
  };

  void       add(uint64 x) {
    if (_eLen >= _eMax) {
      _eMax *= 2;
      uint32 *e = new uint32 [_eMax];
      memcpy(e, _e, sizeof(uint32) * _eLen);
      delete [] _e;
      _e = e;
    }
    _e[_eLen++] = x;

    if (x > _m)
      _l++;
    else
      _h[x/_b]++;
  };

  void       show(char const *label) {
    double  average = 0;
    double  stddev  = 0;

    for (uint32 i=0; i<_eLen; i++)
      average += _e[i];
    average /= _eLen;

    for (uint32 i=0; i<_eLen; i++)
      stddev += (_e[i] - average) * (_e[i] - average);
    stddev = sqrt(stddev / _eLen);

    fprintf(stdout, "histogram %s "uint32FMT" items %8.3f average %8.3f std.dev.\n",
            label, _eLen, average, stddev);
  };

  void       dump(char const *prefix, char const *label) {
    if (noHistogramPlots)
      return;
    char   filename[1024];
    sprintf(filename, "%s.%s.histogramdat", prefix, label);
    FILE *out = fopen(filename, "w");
    for (uint64 i=0; i<_m / _b; i++)
      fprintf(out, uint64FMT" "uint32FMT"\n", i * _b, _h[i]);
    fprintf(out, ">"uint64FMT" "uint32FMT"\n", _m, _l);
    fclose(out);
  }

  void       plot(char const *prefix, char const *label) {

    if (noHistogramPlots)
      return;

    //  Find max's of the data
    uint64  maxx = 0;
    uint64  maxy = 0;

    for (uint64 i=0; i<_m / _b; i++) {
      if (_h[i] > 0)
        maxx = i * _b;
      if (maxy < _h[i])
        maxy = _h[i];
    }

    if ((maxx == 0) || (maxy == 0))
      return;

    char   filename[1024];
    sprintf(filename, "%s.%s.histogram.gnuplot", prefix, label);
    FILE *out = fopen(filename, "w");
    fprintf(out, "set terminal postscript color\n");
    fprintf(out, "set output \"%s.%s.histogram.ps\"\n", prefix, label);
    fprintf(out, "set xlabel \"length bp\"\n");
    fprintf(out, "set ylabel \"number of matches\"\n");
    fprintf(out, "plot [0:"uint64FMT"][0:"uint64FMT"] \"%s.%s.histogramdat\" using 2 with lines\n",
            maxx, maxy, prefix, label);
    fprintf(out, "set output \"%s.%s.histogram.closeup.ps\"\n", prefix, label);
    fprintf(out, "plot [0:"uint64FMT"][0:"uint64FMT"] \"%s.%s.histogramdat\" using 2 with lines\n",
            maxx/10, maxy, prefix, label);
    fprintf(out, "quit\n");
    fclose(out);
    sprintf(filename, "gnuplot < %s.%s.histogram.gnuplot", prefix, label);
    if (system(filename))
      fprintf(stderr, "Failed to execute '%s'\n", filename);
  };


private:
  uint64    _b;  //  blockSize
  uint64    _m;  //  maximum element size
  uint32    _l;  //  number of things bigger than _m
  uint32   *_h;  //  the histogram
  uint32    _eMax;
  uint32    _eLen;
  uint32   *_e;  //  the elements -- for computing the stats;
};



//  Compute the total gapped and ungapped length of the input
//  sequences.  Uses atacMatchList only to access the underlying fasta
//  sequences.
//
void
totalLength(atacFile &AF, seqCache *A, seqCache *B) {
  uint64   length1 = 0;
  uint64   length2 = 0;

  for (uint32 i=0; i<A->getNumberOfSequences(); i++)
    length1 += A->getSequenceLength(i);
  for (uint32 i=0; i<B->getNumberOfSequences(); i++)
    length2 += B->getSequenceLength(i);

  fprintf(stdout, "totalLength    %s "uint64FMT"  %s "uint64FMT" # all letters, including N\n",
          AF.labelA(), length1,
          AF.labelB(), length2);

  length1 = 0;
  length2 = 0;
  for (uint32 i=0; i<A->getNumberOfSequences(); i++) {
    seqInCore   *S = A->getSequenceInCore(i);
    char        *s = S->sequence();
    for (uint32 j=0; j<S->sequenceLength(); j++)
      if (letterToBits[s[j]] != 0xff)
        length1++;
  }
  for (uint32 i=0; i<B->getNumberOfSequences(); i++) {
    seqInCore   *S = B->getSequenceInCore(i);
    char        *s = S->sequence();
    for (uint32 j=0; j<S->sequenceLength(); j++)
      if (letterToBits[s[j]] != 0xff)
        length2++;
  }

  fprintf(stdout, "totalLength    %s "uint64FMT"  %s "uint64FMT" # ACGT only\n",
          AF.labelA(), length1,
          AF.labelB(), length2);
}




uint64
tandemRepeatACGTLength(intervalList &il,
                       uint64       *offset,
                       seqCache     *A) {

  //  s -- the sequence
  //  i -- the interval list index

  il.sort();  //  Both should be done already.
  il.merge();
  uint64 length = 0;
  uint64 unknown[256] = {0};
  for (uint32 i=0, s=0; i<il.numberOfIntervals(); i++) {
    while ((offset[s + 1]) <= il.lo(i)) 
      s++;

    char *S = A->getSequenceInCore(s)->sequence();

    uint64 lo = il.lo(i) - offset[s];
    uint64 hi = il.hi(i) - offset[s];

    for (uint64 j=lo; j < hi; j++)
      if (letterToBits[S[j]] != 0xff)
        length++;
      else
        unknown[S[j]]++;
  }

  //fprintf(stderr, "tandemRepeatACGTLength: "uint64FMT"\n", length);
  //for (uint32 i=0; i<256; i++)
  //  if (unknown[i] > 0)
  //    fprintf(stderr, "tandemRepeatACGTLength["uint32FMT"] = "uint64FMT" (%c)\n", i, unknown[i], i);

  return(length);
}


uint64 *
buildOffset(seqCache *F) {
  uint64  *offset = new uint64 [F->getNumberOfSequences() + 1];
  offset[0] = 1000000;
  for (uint32 i=0; i<F->getNumberOfSequences(); i++)
    offset[i+1] = offset[i] + F->getSequenceLength(i) + 1;
  return(offset);
}


void
tandemRepeatStats(atacFileStream   &featuresA,
                  atacFileStream   &featuresB,
                  atacFile         &AF,
                  seqCache         *A,
                  seqCache         *B) {
  intervalList  ifa, ifb;
  intervalList  ima, imb;
  intervalList  mma, mmb;

  atacMatchList    &matches = *AF.matches();

  uint64  *offset1 = buildOffset(A);
  uint64  *offset2 = buildOffset(B);

  //  ifa, ifb are intervalLists, storing the intervals labeled as
  //  tandem repeats.  They are using the offset[] to encode the
  //  entire sequence as one consecutive string.
  //
  atacFeature  *f = 0L;
  while ((f = featuresA.nextFeature("tr")) != 0L)
    ifa.add(offset1[f->iid] + f->pos, f->len);
  while ((f = featuresB.nextFeature("tr")) != 0L)
    ifb.add(offset2[f->iid] + f->pos, f->len);


  //  ima, imb, like if?, encode the matches in one string.
  //
  for (uint32 m=0; m<matches.numberOfMatches(); m++)
    ima.add(offset1[matches[m]->iid1] + (uint64)matches[m]->pos1, (uint64)matches[m]->len1);
  for (uint32 m=0; m<matches.numberOfMatches(); m++)
    imb.add(offset2[matches[m]->iid2] + (uint64)matches[m]->pos2, (uint64)matches[m]->len2);


  fprintf(stdout, "\nTANDEM REPEATS in %s\n", AF.labelA());
  fprintf(stdout, "numberOfItems "uint64FMT"\n", (uint64)ifa.numberOfIntervals());
  fprintf(stdout, "totalLength   "uint64FMT" # sum of lengths of all features\n", ifa.sumOfLengths());
  ifa.merge();
  fprintf(stdout, "numberOfItems "uint64FMT" # after merging overlapping regions\n", (uint64)ifa.numberOfIntervals());
  fprintf(stdout, "coveredLength "uint64FMT" # sequence covered by a feature, including N\n", ifa.sumOfLengths());
  fprintf(stdout, "coveredLength "uint64FMT" # sequence covered by a feature, ACGT only\n", tandemRepeatACGTLength(ifa, offset1, A));
  mma.intersect(ifa, ima);
  fprintf(stdout, "numberOfItems "uint64FMT" # after merging overlapping regions, only in matches\n", (uint64)mma.numberOfIntervals());
  fprintf(stdout, "inMatches     "uint64FMT" # sequence covered by a feature and in a match, including N\n", mma.sumOfLengths());
  fprintf(stdout, "inMatches     "uint64FMT" # sequence covered by a feature and in a match, ACGT only\n", tandemRepeatACGTLength(mma, offset1, A));


  fprintf(stdout, "\nTANDEM REPEATS in %s\n", AF.labelB());
  fprintf(stdout, "numberOfItems "uint64FMT"\n", (uint64)ifb.numberOfIntervals());
  fprintf(stdout, "totalLength   "uint64FMT" # sum of lengths of all features\n", ifb.sumOfLengths());
  ifb.merge();
  fprintf(stdout, "numberOfItems "uint64FMT" # after merging overlapping regions\n", (uint64)ifb.numberOfIntervals());
  fprintf(stdout, "coveredLength "uint64FMT" # sequence covered by a feature, including N\n", ifb.sumOfLengths());
  fprintf(stdout, "coveredLength "uint64FMT" # sequence covered by a feature, ACGT only\n", tandemRepeatACGTLength(ifb, offset2, B));
  mmb.intersect(ifb, imb);
  fprintf(stdout, "numberOfItems "uint64FMT" # after merging overlapping regions, only in matches\n", (uint64)mmb.numberOfIntervals());
  fprintf(stdout, "inMatches     "uint64FMT" # sequence covered by a feature and in a match, including N\n", mmb.sumOfLengths());
  fprintf(stdout, "inMatches     "uint64FMT" # sequence covered by a feature and in a match, ACGT only\n", tandemRepeatACGTLength(mmb, offset2, B));

  delete [] offset1;
  delete [] offset2;
}



void
mappedLengths(atacFile &AF, atacMatchList &matches, seqCache *A, seqCache *B, char *prefix) {
  histogram  h1(100, 1000000);
  histogram  h2(100, 1000000);

  //  For the coverage to work correctly, we need to either have one
  //  intervalList per input sequence, or build a table of the chained
  //  sequence positions.
  //
  uint64  *offset1 = buildOffset(AF.fastaA());
  uint64  *offset2 = buildOffset(AF.fastaB());

  intervalList  intervalA;
  intervalList  intervalB;

  for (uint32 m=0; m<matches.numberOfMatches(); m++) {
    intervalA.add(offset1[matches[m]->iid1] + (uint64)matches[m]->pos1, (uint64)matches[m]->len1);
    intervalB.add(offset2[matches[m]->iid2] + (uint64)matches[m]->pos2, (uint64)matches[m]->len2);

    h1.add(matches[m]->len1);
    h2.add(matches[m]->len2);
  }

  fprintf(stdout, "numberOfItems "uint64FMT"\n", (uint64)matches.numberOfMatches());

  fprintf(stdout, "matchLength   %s "uint64FMT"  %s "uint64FMT" # Sum of lengths of sequence in matches\n",
          AF.labelA(), (uint64)intervalA.sumOfLengths(),
          AF.labelB(), (uint64)intervalB.sumOfLengths());

  h1.show("AmatchLength");
  h2.show("BmatchLength");
  h1.dump(prefix, "AmatchLength");    h1.plot(prefix, "AmatchLength");
  h2.dump(prefix, "BmatchLength");    h2.plot(prefix, "BmatchLength");

  intervalA.merge();
  intervalB.merge();

  fprintf(stdout, "coveredLength  %s "uint64FMT"  %s "uint64FMT" # sequence covered by a match, including N\n",
          AF.labelA(), (uint64)intervalA.sumOfLengths(),
          AF.labelB(), (uint64)intervalB.sumOfLengths());

  fprintf(stdout, "coveredLength  %s "uint64FMT"  %s "uint64FMT" # sequence covered by a match, ACGT only (new)\n",
          AF.labelA(), tandemRepeatACGTLength(intervalA, offset1, A),
          AF.labelB(), tandemRepeatACGTLength(intervalB, offset2, B));

  delete [] offset1;
  delete [] offset2;
}



//  Generate an Nx plot
void
NxOfMapped(atacFile       &AF,
           atacMatchList  &matches,
           uint64          genomeSize,
           char           *prefix) {

  uint32   *n50             = new uint32 [matches.numberOfMatches()];

  for (uint32 i=0; i<matches.numberOfMatches(); i++)
    n50[i] = matches[i]->len1;

  //  Compute the total length of the sequence
  uint64 totalLength = 0;
  switch (genomeSize) {
    case 0:
      for (uint32 i=0; i<AF.fastaA()->getNumberOfSequences(); i++)
        totalLength += AF.fastaA()->getSequenceLength(i);
      break;
    case 1:
      for (uint32 i=0; i<AF.fastaB()->getNumberOfSequences(); i++)
        totalLength += AF.fastaB()->getSequenceLength(i);
      break;
    default:
      totalLength = genomeSize;
      break;
  }

  //  Sort the n50 list of lengths
  qsort(n50, matches.numberOfMatches(), sizeof(uint32), uint32compare);

  //  It's slow and obvious and, yes, there is a better way.  Dump the
  //  Nx plot as it's being generated.
  //
  char   filename[1024];
  sprintf(filename, "%s.Nx", prefix);
  FILE *out = fopen(filename, "w");

  for (uint64 n=1; n<100; n++) {
    uint64  limit = totalLength / 100 * n;
    uint64  iter  = 0;
    uint64  sum   = 0;

    while ((sum < limit) && (iter < matches.numberOfMatches()))
      sum += n50[iter++];

    fprintf(out, uint64FMT" "uint32FMT"\n", n, n50[iter-1]);
  }

  fclose(out);

  //  Now plot it.
  //
  if (noHistogramPlots == false) {
    sprintf(filename, "%s.Nx.gnuplot", prefix);
    out = fopen(filename, "w");
    fprintf(out, "set terminal postscript color\n");
    fprintf(out, "set output \"%s.Nx.ps\"\n", prefix);
    fprintf(out, "set xlabel \"N\"\n");
    fprintf(out, "set ylabel \"match length\"\n");
    fprintf(out, "plot \"%s.Nx\" using 2 with lines\n", prefix);
    fclose(out);
    sprintf(filename, "gnuplot < %s.Nx.gnuplot", prefix);
    if (system(filename))
      fprintf(stderr, "Failed to execute '%s'\n", filename);
  }

  delete [] n50;
}


//  Computes the percentage of each chromosome (assumes chromosomes are A)
//  that is mapped, with and without N's.
//
void
MappedByChromosome(atacFile      &AF,
                   atacMatchList &matches,
                   seqCache      *A,
                   seqCache      *B,
                   char          *prefix) {

  uint32         maxIID1 = A->getNumberOfSequences();
  intervalList  *il1full;
  intervalList  *il1acgt;
  histogram    **hist1full;
  histogram    **hist1acgt;

  if (A->getNumberOfSequences() > 24) {
    fprintf(stderr, "WARNING: too many sequences to be chromosomes, only using the first 24.\n");
    maxIID1 = 24;
  }

  //  We could cache this when we compute the totalLength() above
  uint64   *nonNlength = new uint64 [maxIID1+1];
  for (uint32 i=0; i<maxIID1; i++) {
    seqInCore   *S = A->getSequenceInCore(i);
    char        *s = S->sequence();
    nonNlength[i] = 0;
    for (uint32 j=0; j<S->sequenceLength(); j++)
      if (letterToBits[s[j]] != 0xff)
        nonNlength[i]++;
  }

  il1full = new intervalList [maxIID1 + 1];
  il1acgt = new intervalList [maxIID1 + 1];

  hist1full = new histogram * [maxIID1 + 1];
  hist1acgt = new histogram * [maxIID1 + 1];

  for (uint32 i=0; i<maxIID1; i++) {
    hist1full[i] = new histogram(100, 1000000);
    hist1acgt[i] = new histogram(100, 1000000);
  }

  for (uint32 m=0; m<matches.numberOfMatches(); m++) {
    if (matches[m]->iid1 < maxIID1) {
      il1full[matches[m]->iid1].add(matches[m]->pos1, matches[m]->len1);
      hist1full[matches[m]->iid1]->add(matches[m]->len1);

      seqInCore   *Sa = A->getSequenceInCore(matches[m]->iid1);
      char        *sa = Sa->sequence() + matches[m]->pos1;

      uint32               length = 0;

      for (uint32 j=0; j<matches[m]->len1; j++) {
        bool invalid = (letterToBits[sa[j]] == 0xff);

        if (!invalid)
          length++;

        if (length && invalid) {        //  Last time we were ACGT, this time not.
          il1acgt[matches[m]->iid1].add(matches[m]->pos1 + j - length, length);
          hist1acgt[matches[m]->iid1]->add(length);
          length = 0;
        }
      }
      if (length) {
        il1acgt[matches[m]->iid1].add(matches[m]->pos1 + matches[m]->len1 - length, length);
        hist1acgt[matches[m]->iid1]->add(length);
      }
    }
  }

  for (uint32 c=0; c<maxIID1; c++) {
    fprintf(stdout, "chrCoveredLength["uint32FMTW(2)"]  %s "uint64FMT" "uint64FMT" %6.2f%%   "uint64FMT" "uint64FMT" %6.2f%% # seqCov, totalSeq for both ALL and ACGTonly\n",
            c, AF.labelA(),
            il1full[c].sumOfLengths(), (uint64)A->getSequenceLength(c), 100.0 * il1full[c].sumOfLengths() / A->getSequenceLength(c),
            il1acgt[c].sumOfLengths(), nonNlength[c], 100.0 * il1acgt[c].sumOfLengths() / nonNlength[c]);
  }

  for (uint32 c=0; c<maxIID1; c++) {
    char  label[1024];

    sprintf(label, "chr"uint32FMTW(02)"full", c);
    hist1full[c]->dump(prefix, label);
    hist1full[c]->plot(prefix, label);

    sprintf(label, "chr"uint32FMTW(02)"acgt", c);
    hist1acgt[c]->dump(prefix, label);
    hist1acgt[c]->plot(prefix, label);
  }

  delete [] il1full;
  delete [] il1acgt;
  for (uint32 i=0; i<maxIID1; i++) {
    delete hist1full[i];
    delete hist1acgt[i];
  }
  delete [] hist1full;
  delete [] hist1acgt;
  delete [] nonNlength;
}






void
statsInACGT(seqInCore       *S,
            uint32           beg,
            uint32           len,
            intervalList    *IL,
            histogram       *HI) {
  char     *s = S->sequence() + beg;
  uint32    length = 0;

  for (uint32 j=0; j<len; j++) {
    bool invalid = (letterToBits[s[j]] == 0xff);

    if (!invalid)
      length++;

    if (length && invalid) {        //  Last time we were ACGT, this time not.
      if (IL) IL->add(beg + j - length, length);
      if (HI) HI->add(length);
      length = 0;
    }
  }
  if (length) {
    if (IL)  IL->add(beg + len - length, length);
    if (HI)  HI->add(length);
  }
}




//  Computes the amount of ACGT in runs that is unmapped
//
void
unmappedInRuns(atacFile &AF, seqCache *A, seqCache *B, char *prefix) {

  atacMatchList &matches = *AF.matches();

  //  We must sort by the location and not the parentID; when we
  //  stream through, we check that the pair of matches are in the
  //  same parent.
  //
  atacMatchOrder  MO(matches);
  MO.sortA();

  intervalList  il1full, il2full;
  intervalList  il1acgt, il2acgt;

  histogram     hist1full(100, 1000000), hist2full(100, 1000000);
  histogram     hist1acgt(100, 1000000), hist2acgt(100, 1000000);

  for (uint32 i=1; i<MO.numberOfMatches(); i++) {
    if (strcmp(MO[i-1]->parentuid, MO[i]->parentuid) == 0) {
      uint32  l1, r1, l2, r2;

      if (MO[i]->fwd2 == 1) {
        l1 = MO[i-1]->pos1 + MO[i-1]->len1;
        r1 = MO[i]->pos1;
        l2 = MO[i-1]->pos2 + MO[i-1]->len2;
        r2 = MO[i]->pos2;
      } else {
        l1 = MO[i-1]->pos1 + MO[i-1]->len1;
        r1 = MO[i]->pos1;
        l2 = MO[i]->pos2 + MO[i]->len2;
        r2 = MO[i-1]->pos2;
      }

      il1full.add(l1, r1-l1);
      il2full.add(l2, r2-l2);

      hist1full.add(r1-l1);
      hist2full.add(r2-l2);

      statsInACGT(A->getSequenceInCore(MO[i]->iid1),
                  l1,
                  r1-l1,
                  &il1acgt,
                  &hist1acgt);
      statsInACGT(B->getSequenceInCore(MO[i]->iid2),
                  l2,
                  r2-l2,
                  &il2acgt,
                  &hist2acgt);
    }
  }

  //  Dump the stats

  fprintf(stdout, "runMissingFull  %s "uint64FMT"  %s "uint64FMT" # sequence in run, not covered, including N\n",
          AF.labelA(), (uint64)il1full.sumOfLengths(),
          AF.labelB(), (uint64)il2full.sumOfLengths());
  fprintf(stdout, "runMissingFull  %s "uint64FMT"  %s "uint64FMT" # sequence in run, not covered, ACGT only\n",
          AF.labelA(), (uint64)il1acgt.sumOfLengths(),
          AF.labelB(), (uint64)il2acgt.sumOfLengths());

  hist1full.dump(prefix, "ARunMissingFull");
  hist1full.plot(prefix, "ARunMissingFull");

  hist2full.dump(prefix, "BRunMissingFull");
  hist2full.plot(prefix, "BRunMissingFull");

  hist1acgt.dump(prefix, "ARunMissingACGT");
  hist1acgt.plot(prefix, "ARunMissingACGT");

  hist2acgt.dump(prefix, "BRunMissingACGT");
  hist2acgt.plot(prefix, "BRunMissingACGT");
}



int
main(int argc, char **argv) {
  uint64  genomeSize         = 0;
  char   *atacFileName       = 0L;
  char   *prefix             = 0L;
  char   *trFile1  = 0L;
  char   *trFile2  = 0L;
  char    prefixFull[1024];
  bool    error              = false;

  int arg=1;
  while (arg < argc) {
    if (strcmp(argv[arg], "-g") == 0) {
      ++arg;
      if        (argv[arg][0] == 'A') {
        genomeSize = 0;
      } else if (argv[arg][0] == 'B') {
        genomeSize = 1;
      } else {
        genomeSize = strtouint64(argv[arg], 0L);
      }
    } else if (strcmp(argv[arg], "-a") == 0) {
      atacFileName = argv[++arg];
    } else if (strcmp(argv[arg], "-p") == 0) {
      prefix = argv[++arg];
    } else if (strcmp(argv[arg], "-ta") == 0) {
      trFile1 = argv[++arg];
    } else if (strcmp(argv[arg], "-tb") == 0) {
      trFile2 = argv[++arg];
    } else {
      error = true;
    }
    arg++;
  }

  if (!atacFileName || !prefix || error) {
    fprintf(stderr, "usage: %s -a <file.atac> -p <outprefix> [-ta trfile] [-tb trfile] [-g {A | B | g}]\n", argv[0]);
    fprintf(stderr, "  -a          read input from 'file.atac'\n");
    fprintf(stderr, "  -p          write stats to files prefixed with 'outprefix'\n");
    fprintf(stderr, "  -g          use a genome size of g for the Nx computation, defaults to\n");
    fprintf(stderr, "              the length of the A sequence.  Or use the actual length\n");
    fprintf(stderr, "              of sequence A or B.\n");
    fprintf(stderr, "  -ta         read tandem repeats for A from trfile\n");
    fprintf(stderr, "  -tb         read tandem repeats for B from trfile\n");
    exit(1);
  }
  
  atacFile           AF(atacFileName);
  atacMatchList     &matches = *AF.matches();
  atacMatchList     &runs    = *AF.runs();
  atacMatchList     &clumps  = *AF.clumps();

  //  We end up using sequences a lot here, so just bite it and load them in a cache.
  //
  seqCache  *A = new seqCache(AF.assemblyFileA(), 0, true);
  seqCache  *B = new seqCache(AF.assemblyFileB(), 0, true);

  A->loadAllSequences();
  B->loadAllSequences();

  fprintf(stdout, "\nSEQUENCE\n");
  totalLength(AF, A, B);

  if (trFile1 && trFile2) {
    atacFileStream     tr1(trFile1);
    atacFileStream     tr2(trFile2);
    tandemRepeatStats(tr1, tr2, AF, A, B);
  }

  //  XXX unmappedInRuns only works on runs, and if we have clumps in
  //  the input it fails.
  //
  if ((runs.numberOfMatches() > 0) && (clumps.numberOfMatches() == 0)) {
    fprintf(stdout, "\nMATCHES IN RUNS\n");
    unmappedInRuns(AF, A, B, prefix);
  }

  if (matches.numberOfMatches() > 0) {
    fprintf(stdout, "\nMATCHES\n");
    sprintf(prefixFull, "%s-matches", prefix);
    mappedLengths(AF, matches, A, B, prefixFull);
    NxOfMapped(AF, matches, genomeSize, prefixFull);
    MappedByChromosome(AF, matches, A, B, prefixFull);
  }

  if (runs.numberOfMatches() > 0) {
    fprintf(stdout, "\nRUNS\n");
    sprintf(prefixFull, "%s-runs", prefix);
    mappedLengths(AF, runs, A, B, prefixFull);
    NxOfMapped(AF, runs, genomeSize, prefixFull);
    MappedByChromosome(AF, runs, A, B, prefixFull);
  }

  if (clumps.numberOfMatches() > 0) {
    fprintf(stdout, "\nCLUMPS\n");
    sprintf(prefixFull, "%s-clumps", prefix);
    mappedLengths(AF, clumps, A, B, prefixFull);
    NxOfMapped(AF, clumps, genomeSize, prefixFull);
    MappedByChromosome(AF, clumps, A, B, prefixFull);
  }

  delete A;
  delete B;

  return(0);
}
