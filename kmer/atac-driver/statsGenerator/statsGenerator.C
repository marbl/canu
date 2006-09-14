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

#include "bio++.H"
#include "util++.H"
#include "atac.H"


//  Sort u32bit backwards
int
u32bitcompare(const void *a, const void *b) {
  const u32bit A = *((const u32bit *)a);
  const u32bit B = *((const u32bit *)b);
  if (A < B) return(1);
  if (A > B) return(-1);
  return(0);
}


class histogram {
public:
  histogram(u64bit blockSize, u64bit maxSize) {
    _b = blockSize;
    _m = maxSize;
    _l = 0;
    _h = new u32bit [maxSize / blockSize + 1];
    for (u32bit i=0; i<maxSize / blockSize + 1; i++)
      _h[i] = 0;
    _eLen = 0;
    _eMax = 100000;
    _e = new u32bit [_eMax];
  };
  ~histogram() {
    delete [] _h;
  };

  void       add(u64bit x) {
    if (_eLen >= _eMax) {
      _eMax *= 2;
      u32bit *e = new u32bit [_eMax];
      memcpy(e, _e, sizeof(u32bit) * _eLen);
      delete [] _e;
      _e = e;
    }
    _e[_eLen++] = x;

    if (x > _m)
      _l++;
    else
      _h[x/_b]++;
  };

  void       show(char *label) {
    double  average = 0;
    double  stddev  = 0;

    for (u32bit i=0; i<_eLen; i++)
      average += _e[i];
    average /= _eLen;

    for (u32bit i=0; i<_eLen; i++)
      stddev += (_e[i] - average) * (_e[i] - average);
    stddev = sqrt(stddev / _eLen);

    fprintf(stdout, "histogram %s "u32bitFMT" items %8.3f average %8.3f std.dev.\n",
            label, _eLen, average, stddev);
  };

  void       dump(char *prefix, char *label) {
    char   filename[1024];
    sprintf(filename, "%s.%s.histogramdat", prefix, label);
    FILE *out = fopen(filename, "w");
    for (u64bit i=0; i<_m / _b; i++)
      fprintf(out, u64bitFMT" "u32bitFMT"\n", i * _b, _h[i]);
    fprintf(out, ">"u64bitFMT" "u32bitFMT"\n", _m, _l);
    fclose(out);
  }

  void       plot(char *prefix, char *label) {

    //  Find max's of the data
    u64bit  maxx = 0;
    u64bit  maxy = 0;

    for (u64bit i=0; i<_m / _b; i++) {
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
    fprintf(out, "plot [0:"u64bitFMT"][0:"u64bitFMT"] \"%s.%s.histogramdat\" using 2 with lines\n",
            maxx, maxy, prefix, label);
    fprintf(out, "set output \"%s.%s.histogram.closeup.ps\"\n", prefix, label);
    fprintf(out, "plot [0:"u64bitFMT"][0:"u64bitFMT"] \"%s.%s.histogramdat\" using 2 with lines\n",
            maxx/10, maxy, prefix, label);
    fprintf(out, "quit\n");
    fclose(out);
    sprintf(filename, "gnuplot < %s.%s.histogram.gnuplot", prefix, label);
    if (system(filename))
      fprintf(stderr, "Failed to execute '%s'\n", filename);
  };


private:
  u64bit    _b;  //  blockSize
  u64bit    _m;  //  maximum element size
  u32bit    _l;  //  number of things bigger than _m
  u32bit   *_h;  //  the histogram
  u32bit    _eMax;
  u32bit    _eLen;
  u32bit   *_e;  //  the elements -- for computing the stats;
};



//  Compute the total gapped and ungapped length of the input
//  sequences.  Uses atacMatchList only to access the underlying fasta
//  sequences.
//
void
totalLength(atacMatchList &matches, FastACache *A, FastACache *B) {
  u64bit   length1 = 0;
  u64bit   length2 = 0;

  for (u32bit i=0; i<matches.fastaA()->getNumberOfSequences(); i++)
    length1 += matches.fastaA()->sequenceLength(i);
  for (u32bit i=0; i<matches.fastaB()->getNumberOfSequences(); i++)
    length2 += matches.fastaB()->sequenceLength(i);

  fprintf(stdout, "totalLength    %s "u64bitFMT"  %s "u64bitFMT" # all letters, including N\n",
          matches.labelA(), length1,
          matches.labelB(), length2);

  length1 = 0;
  length2 = 0;
  for (u32bit i=0; i<matches.fastaA()->getNumberOfSequences(); i++) {
    FastASequenceInCore *S = A->getSequence(i);
    char                *s = S->sequence();
    for (u32bit j=0; j<S->sequenceLength(); j++)
      if (validSymbol[s[j]])
        length1++;
  }
  for (u32bit i=0; i<matches.fastaB()->getNumberOfSequences(); i++) {
    FastASequenceInCore *S = B->getSequence(i);
    char                *s = S->sequence();
    for (u32bit j=0; j<S->sequenceLength(); j++)
      if (validSymbol[s[j]])
        length2++;
  }

  fprintf(stdout, "totalLength    %s "u64bitFMT"  %s "u64bitFMT" # ACGT only\n",
          matches.labelA(), length1,
          matches.labelB(), length2);
}




u64bit
tandemRepeatACGTLength(intervalList &mma,
                       u64bit       *offset1,
                       FastACache   *A) {

  mma.sort();  //  Both should be done already.
  mma.merge();
  u64bit length1 = 0;
  for (u32bit i=0, s=0; i<mma.numberOfIntervals(); i++) {
    while ((offset1[s + 1]) < mma.lo(i)) 
      s++;
    char *S = A->getSequence(s)->sequence();
    for (u64bit j=mma.lo(i)-offset1[s]; j < mma.hi(i)-offset1[s]; j++)
      if (validSymbol[S[j]])
        length1++;
  }

  return(length1);
}


void
tandemRepeatStats(atacFeatureList &AF, atacFeatureList &BF,
                  atacMatchList &matches,
                  FastACache *A, FastACache *B) {
  intervalList  ifa, ifb;
  intervalList  ima, imb;
  intervalList  mma, mmb;

  //  XXX  Stolen from mappedLengths() below
  //
  u64bit  *offset1 = new u64bit [matches.fastaA()->getNumberOfSequences() + 1];
  u64bit  *offset2 = new u64bit [matches.fastaB()->getNumberOfSequences() + 1];
  offset1[0] = 1000000;
  for (u32bit i=1; i<=matches.fastaA()->getNumberOfSequences(); i++)
    offset1[i] = offset1[i-1] + matches.fastaA()->sequenceLength(i-1) + 1;
  offset2[0] = 1000000;
  for (u32bit i=1; i<=matches.fastaB()->getNumberOfSequences(); i++)
    offset2[i] = offset2[i-1] + matches.fastaB()->sequenceLength(i-1) + 1;

  for (u32bit i=0; i<AF.numberOfFeatures(); i++)
    ifa.add(offset1[AF[i]->iid] + AF[i]->pos, AF[i]->len);
  for (u32bit i=0; i<BF.numberOfFeatures(); i++)
    ifb.add(offset2[BF[i]->iid] + BF[i]->pos, BF[i]->len);

  for (u32bit m=0; m<matches.numberOfMatches(); m++)
    ima.add(offset1[matches[m]->iid1] + (u64bit)matches[m]->pos1, (u64bit)matches[m]->len1);
  for (u32bit m=0; m<matches.numberOfMatches(); m++)
    imb.add(offset2[matches[m]->iid2] + (u64bit)matches[m]->pos2, (u64bit)matches[m]->len2);


  fprintf(stdout, "\nTANDEM REPEATS in %s\n", AF.label());
  fprintf(stdout, "numberOfItems "u64bitFMT"\n", (u64bit)AF.numberOfFeatures());
  fprintf(stdout, "totalLength   "u64bitFMT" # sum of lengths of all features\n", ifa.sumOfLengths());
  ifa.merge();
  fprintf(stdout, "numberOfItems "u64bitFMT" # after merging overlapping regions\n", (u64bit)ifa.numberOfIntervals());
  fprintf(stdout, "coveredLength "u64bitFMT" # sequence covered by a feature, including N\n", ifa.sumOfLengths());
  fprintf(stdout, "coveredLength "u64bitFMT" # sequence covered by a feature, ACGT only\n",
          tandemRepeatACGTLength(ifa, offset1, A));
  mma.intersect(ifa, ima);
  fprintf(stdout, "numberOfItems "u64bitFMT" # after merging overlapping regions, only in matches\n", (u64bit)mma.numberOfIntervals());
  fprintf(stdout, "inMatches     "u64bitFMT" # sequence covered by a feature and in a match, including N\n", mma.sumOfLengths());
  fprintf(stdout, "inMatches     "u64bitFMT" # sequence covered by a feature and in a match, ACGT only\n",
          tandemRepeatACGTLength(mma, offset1, A));


  fprintf(stdout, "\nTANDEM REPEATS in %s\n", BF.label());
  fprintf(stdout, "numberOfItems "u64bitFMT"\n", (u64bit)BF.numberOfFeatures());
  fprintf(stdout, "totalLength   "u64bitFMT" # sum of lengths of all features\n", ifb.sumOfLengths());
  ifb.merge();
  fprintf(stdout, "numberOfItems "u64bitFMT" # after merging overlapping regions\n", (u64bit)ifb.numberOfIntervals());
  fprintf(stdout, "coveredLength "u64bitFMT" # sequence covered by a feature, including N\n", ifb.sumOfLengths());
  fprintf(stdout, "coveredLength "u64bitFMT" # sequence covered by a feature, ACGT only\n",
          tandemRepeatACGTLength(ifb, offset2, B));
  mmb.intersect(ifb, imb);
  fprintf(stdout, "numberOfItems "u64bitFMT" # after merging overlapping regions, only in matches\n", (u64bit)mmb.numberOfIntervals());
  fprintf(stdout, "inMatches     "u64bitFMT" # sequence covered by a feature and in a match, including N\n", mmb.sumOfLengths());
  fprintf(stdout, "inMatches     "u64bitFMT" # sequence covered by a feature and in a match, ACGT only\n",
          tandemRepeatACGTLength(mmb, offset2, B));


  delete [] offset1;
  delete [] offset2;
}



void
mappedLengths(atacMatchList &matches, char *prefix) {
  histogram  h1(100, 1000000);
  histogram  h2(100, 1000000);

  //  For the coverage to work correctly, we need to either have one
  //  intervalList per input sequence, or build a table of the chained
  //  sequence positions.
  //
  u64bit  *offset1 = new u64bit [matches.fastaA()->getNumberOfSequences() + 1];
  u64bit  *offset2 = new u64bit [matches.fastaB()->getNumberOfSequences() + 1];
  offset1[0] = 1000000;
  for (u32bit i=1; i<=matches.fastaA()->getNumberOfSequences(); i++)
    offset1[i] = offset1[i-1] + matches.fastaA()->sequenceLength(i-1) + 1;
  offset2[0] = 1000000;
  for (u32bit i=1; i<=matches.fastaB()->getNumberOfSequences(); i++)
    offset2[i] = offset2[i-1] + matches.fastaB()->sequenceLength(i-1) + 1;

  intervalList  intervalA;
  intervalList  intervalB;

  for (u32bit m=0; m<matches.numberOfMatches(); m++) {
    intervalA.add(offset1[matches[m]->iid1] + (u64bit)matches[m]->pos1, (u64bit)matches[m]->len1);
    intervalB.add(offset2[matches[m]->iid2] + (u64bit)matches[m]->pos2, (u64bit)matches[m]->len2);

    h1.add(matches[m]->len1);
    h2.add(matches[m]->len2);
  }

  delete [] offset1;
  delete [] offset2;

  fprintf(stdout, "numberOfItems "u64bitFMT"\n", (u64bit)matches.numberOfMatches());

  fprintf(stdout, "matchLength   %s "u64bitFMT"  %s "u64bitFMT" # Sum of lengths of sequence in matches\n",
          matches.labelA(), (u64bit)intervalA.sumOfLengths(),
          matches.labelB(), (u64bit)intervalB.sumOfLengths());

  h1.show("AmatchLength");
  h2.show("BmatchLength");
  h1.dump(prefix, "AmatchLength");    h1.plot(prefix, "AmatchLength");
  h2.dump(prefix, "BmatchLength");    h2.plot(prefix, "BmatchLength");

  intervalA.merge();
  intervalB.merge();

  fprintf(stdout, "coveredLength  %s "u64bitFMT"  %s "u64bitFMT" # sequence covered by a match, including N\n",
          matches.labelA(), (u64bit)intervalA.sumOfLengths(),
          matches.labelB(), (u64bit)intervalB.sumOfLengths());
}


//  Returns the amount of N covered by a match (only possible for runs, we hope!)
//
void
mappedNs(atacMatchList &matches, FastACache *A, FastACache *B, char *prefix) {
  histogram  h1(100, 1000000);
  histogram  h2(100, 1000000);

  u64bit   length  = 0;
  u64bit   length1 = 0;
  u64bit   length2 = 0;
  u64bit   length1n = 0;
  u64bit   length2n = 0;

  for (u32bit m=0; m<matches.numberOfMatches(); m++) {
    FastASequenceInCore *Sa = A->getSequence(matches[m]->iid1);
    FastASequenceInCore *Sb = B->getSequence(matches[m]->iid2);

    char                *sa = Sa->sequence() + matches[m]->pos1;
    char                *sb = Sb->sequence() + matches[m]->pos2;

    length = 0;
    for (u32bit j=0; j<matches[m]->len1; j++) {
      bool valid = (validSymbol[sa[j]] != 0);
      if (valid) {
        length1++;
      } else {
        length1n++;
        length++;
      }
      if (length && valid) {        //  Last time we were N, this time not.
        h1.add(length);
        length = 0;
      }
    }
    if (length)
      h1.add(length);

    length = 0;
    for (u32bit j=0; j<matches[m]->len2; j++) {
      bool valid = (validSymbol[sb[j]] != 0);
      if (valid) {
        length2++;
      } else {
        length2n++;
        length++;
      }
      if (length && valid) {        //  Last time we were N, this time not.
        h2.add(length);
        length = 0;
      }
    }
    if (length)
      h2.add(length);
  }

  fprintf(stdout, "coveredLength  %s "u64bitFMT"  %s "u64bitFMT" # sequence covered by a match, ACGT only\n",
          matches.labelA(), length1,
          matches.labelB(), length2);
  fprintf(stdout, "coveredLength  %s "u64bitFMT"  %s "u64bitFMT" # sequence covered by a match, non ACGT\n",
          matches.labelA(), length1n,
          matches.labelB(), length2n);

  h1.show("AcoveredN");
  h2.show("BcoveredN");
  h1.dump(prefix, "AcoveredN");  h1.plot(prefix, "AcoveredN");
  h2.dump(prefix, "BcoveredN");  h2.plot(prefix, "BcoveredN");
}



//  Generate an Nx plot
void
NxOfMapped(atacMatchList &matches,
           u64bit  genomeSize,
           char *prefix) {

  u32bit   *n50             = new u32bit [matches.numberOfMatches()];

  for (u32bit i=0; i<matches.numberOfMatches(); i++)
    n50[i] = matches[i]->len1;

  //  Compute the total length of the sequence
  u64bit totalLength = 0;
  switch (genomeSize) {
    case 0:
      for (u32bit i=0; i<matches.fastaA()->getNumberOfSequences(); i++)
        totalLength += matches.fastaA()->sequenceLength(i);
      break;
    case 1:
      for (u32bit i=0; i<matches.fastaB()->getNumberOfSequences(); i++)
        totalLength += matches.fastaB()->sequenceLength(i);
      break;
    default:
      totalLength = genomeSize;
      break;
  }

  //  Sort the n50 list of lengths
  qsort(n50, matches.numberOfMatches(), sizeof(u32bit), u32bitcompare);

  //  It's slow and obvious and, yes, there is a better way.  Dump the
  //  Nx plot as it's being generated.
  //
  char   filename[1024];
  sprintf(filename, "%s.Nx", prefix);
  FILE *out = fopen(filename, "w");

  for (u64bit n=1; n<100; n++) {
    u64bit  limit = totalLength / 100 * n;
    u64bit  iter  = 0;
    u64bit  sum   = 0;

    while ((sum < limit) && (iter < matches.numberOfMatches()))
      sum += n50[iter++];

    fprintf(out, u64bitFMT" "u32bitFMT"\n", n, n50[iter-1]);
  }

  fclose(out);

  //  Now plot it.
  //
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

  delete [] n50;
}


//  Computes the percentage of each chromosome (assumes chromosomes are A)
//  that is mapped, with and without N's.
//
void
MappedByChromosome(atacMatchList &matches, FastACache *A, FastACache *B, char *prefix) {
  u32bit         maxIID1 = matches.fastaA()->getNumberOfSequences();
  intervalList  *il1full;
  intervalList  *il1acgt;
  histogram    **hist1full;
  histogram    **hist1acgt;

  if (matches.fastaA()->getNumberOfSequences() > 24) {
    fprintf(stderr, "WARNING: too many sequences to be chromosomes, only using the first 24.\n");
    maxIID1 = 24;
  }

  //  We could cache this when we compute the totalLength() above
  u64bit   *nonNlength = new u64bit [maxIID1+1];
  for (u32bit i=0; i<maxIID1; i++) {
    FastASequenceInCore *S = A->getSequence(i);
    char                *s = S->sequence();
    nonNlength[i] = 0;
    for (u32bit j=0; j<S->sequenceLength(); j++)
      if (validSymbol[s[j]])
        nonNlength[i]++;
  }

  il1full = new intervalList [maxIID1 + 1];
  il1acgt = new intervalList [maxIID1 + 1];

  hist1full = new histogram * [maxIID1 + 1];
  hist1acgt = new histogram * [maxIID1 + 1];

  for (u32bit i=0; i<maxIID1; i++) {
    hist1full[i] = new histogram(100, 1000000);
    hist1acgt[i] = new histogram(100, 1000000);
  }

  for (u32bit m=0; m<matches.numberOfMatches(); m++) {
    if (matches[m]->iid1 < maxIID1) {
      il1full[matches[m]->iid1].add(matches[m]->pos1, matches[m]->len1);
      hist1full[matches[m]->iid1]->add(matches[m]->len1);

      FastASequenceInCore *Sa = A->getSequence(matches[m]->iid1);
      char                *sa = Sa->sequence() + matches[m]->pos1;

      u32bit               length = 0;

      for (u32bit j=0; j<matches[m]->len1; j++) {
        bool invalid = (validSymbol[sa[j]] == 0);

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

  for (u32bit c=0; c<maxIID1; c++) {
    fprintf(stdout, "chrCoveredLength["u32bitFMTW(2)"]  %s "u64bitFMT" "u64bitFMT" %6.2f%%   "u64bitFMT" "u64bitFMT" %6.2f%% # seqCov, totalSeq for both ALL and ACGTonly\n",
            c, matches.labelA(),
            il1full[c].sumOfLengths(), (u64bit)matches.fastaA()->sequenceLength(c), 100.0 * il1full[c].sumOfLengths() / matches.fastaA()->sequenceLength(c),
            il1acgt[c].sumOfLengths(), nonNlength[c], 100.0 * il1acgt[c].sumOfLengths() / nonNlength[c]);
  }

  for (u32bit c=0; c<maxIID1; c++) {
    char  label[1024];

    sprintf(label, "chr%02dfull", c);
    hist1full[c]->dump(prefix, label);
    hist1full[c]->plot(prefix, label);

    sprintf(label, "chr%02dacgt", c);
    hist1acgt[c]->dump(prefix, label);
    hist1acgt[c]->plot(prefix, label);
  }

  delete [] il1full;
  delete [] il1acgt;
  for (u32bit i=0; i<maxIID1; i++) {
    delete hist1full[i];
    delete hist1acgt[i];
  }
  delete [] hist1full;
  delete [] hist1acgt;
  delete [] nonNlength;
}






void
statsInACGT(FastASequenceInCore  *S,
            u32bit                beg,
            u32bit                len,
            intervalList         *IL,
            histogram            *HI) {
  char     *s = S->sequence() + beg;
  u32bit    length = 0;

  for (u32bit j=0; j<len; j++) {
    bool invalid = (validSymbol[s[j]] == 0);

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
unmappedInRuns(atacMatchList &matches, FastACache *A, FastACache *B, char *prefix) {

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

  for (u32bit i=1; i<MO.numberOfMatches(); i++) {
    if (strcmp(MO[i-1]->parentuid, MO[i]->parentuid) == 0) {
      u32bit  l1, r1, l2, r2;

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

      //  Crimeny!  I really should put this in a function....  Lessee...it needs the
      //  sequence, the begin and length, the il and the histogram.

      statsInACGT(A->getSequence(MO[i]->iid1),
                  l1,
                  r1-l1,
                  &il1acgt,
                  &hist1acgt);
      statsInACGT(B->getSequence(MO[i]->iid2),
                  l2,
                  r2-l2,
                  &il2acgt,
                  &hist2acgt);
    }
  }

  //  Dump the stats

  fprintf(stdout, "runMissingFull  %s "u64bitFMT"  %s "u64bitFMT" # sequence in run, not covered, including N\n",
          matches.labelA(), (u64bit)il1full.sumOfLengths(),
          matches.labelB(), (u64bit)il2full.sumOfLengths());
  fprintf(stdout, "runMissingFull  %s "u64bitFMT"  %s "u64bitFMT" # sequence in run, not covered, ACGT only\n",
          matches.labelA(), (u64bit)il1acgt.sumOfLengths(),
          matches.labelB(), (u64bit)il2acgt.sumOfLengths());


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
  u64bit  genomeSize         = 0;
  char   *atacFile           = 0L;
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
        genomeSize = strtou64bit(argv[arg], 0L);
      }
    } else if (strcmp(argv[arg], "-a") == 0) {
      atacFile = argv[++arg];
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

  if (!atacFile || !prefix || error) {
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

  atacMatchList      matches(atacFile, 'm');
  atacMatchList      runs(atacFile,    'r');

  //  We end up using sequences a lot here, so just bite it and load them in a cache.
  //
  FastACache  *A = new FastACache(matches.assemblyFileA(), 0, true, true);
  FastACache  *B = new FastACache(matches.assemblyFileB(), 0, true, true);

  fprintf(stdout, "\nSEQUENCE\n");
  totalLength(matches, A, B);

  if (trFile1 && trFile2) {
    atacFeatureList    tr1(trFile1);
    atacFeatureList    tr2(trFile2);

    fprintf(stdout, "\nTANDEM REPEATS\n");
    tandemRepeatStats(tr1, tr2, matches, A, B);
  }

  fprintf(stdout, "\nMATCHES IN RUNS\n");
  unmappedInRuns(matches, A, B, prefix);

  fprintf(stdout, "\nMATCHES\n");
  sprintf(prefixFull, "%s-matches", prefix);
  mappedLengths(matches, prefixFull);
  mappedNs(matches, A, B, prefixFull);
  NxOfMapped(matches, genomeSize, prefixFull);
  MappedByChromosome(matches, A, B, prefixFull);

  fprintf(stdout, "\nRUNS\n");
  sprintf(prefixFull, "%s-runs", prefix);
  mappedLengths(runs, prefixFull);
  mappedNs(runs, A, B, prefixFull);
  NxOfMapped(runs, genomeSize, prefixFull);
  MappedByChromosome(runs, A, B, prefixFull);
  
  delete A;
  delete B;

  return(0);
}
