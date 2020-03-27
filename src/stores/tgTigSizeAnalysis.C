
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' r4587 (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' r1994 (http://kmer.sourceforge.net)
 *
 *  Except as indicated otherwise, this is a 'United States Government Work',
 *  and is released in the public domain.
 *
 *  File 'README.licenses' in the root directory of this distribution
 *  contains full conditions and disclaimers.
 */

#include "tgTigSizeAnalysis.H"

#include <math.h>

#include <map>
#include <set>
#include <vector>
#include <algorithm>

using namespace std;


tgTigSizeAnalysis::tgTigSizeAnalysis(uint64 genomeSize_) {
  genomeSize = genomeSize_;
}

tgTigSizeAnalysis::~tgTigSizeAnalysis() {
}

void
tgTigSizeAnalysis::evaluateTig(tgTig *tig) {
  uint32  length = tig->length();

  if (tig->_suggestRepeat)
    lenSuggestRepeat.push_back(length);

  if (tig->_suggestCircular)
    lenSuggestCircular.push_back(length);

  switch (tig->_class) {
    case tgTig_unassembled:
      lenUnassembled.push_back(length);
      break;
    case tgTig_contig:
      if (tig->_suggestBubble)
         lenBubble.push_back(length);
      else
         lenContig.push_back(length);
      break;
    default:
      break;
  }
}

void
tgTigSizeAnalysis::finalize(void) {

  sort(lenSuggestRepeat.begin(),   lenSuggestRepeat.end(),   greater<uint32>());
  sort(lenSuggestCircular.begin(), lenSuggestCircular.end(), greater<uint32>());

  sort(lenUnassembled.begin(), lenUnassembled.end(), greater<uint32>());
  sort(lenBubble.begin(),      lenBubble.end(),      greater<uint32>());
  sort(lenContig.begin(),      lenContig.end(),      greater<uint32>());
}

void
tgTigSizeAnalysis::printSummary(FILE *out, char const *description, vector<uint32> &data) {
  uint64  cnt = data.size();
  uint64  sum = 0;
  uint64  tot = 0;
  uint64  nnn = 10;
  uint64  siz = 0;

  //  Duplicates AS_BAT_Instrumentation.C reportN50().

  if (cnt == 0)
    return;

  for (uint64 i=0; i<cnt; i++)
    tot += data[i];

  if (genomeSize > 0)
    siz = genomeSize;
  else
    siz = tot;

  for (uint64 i=0; i<cnt; i++) {
    sum += data[i];

    while (siz * nnn / 100 < sum) {
      fprintf(out, "%s ng%-3" F_U64P " %10" F_U32P " bp   lg%-3" F_U64P " %6" F_U64P "   sum %10" F_U64P " bp\n",
              description,
              nnn, data[i],
              nnn, i+1,
              sum);

      nnn += 10;
    }
  }

  fprintf(out, "%s sum %10" F_U64P " (genomeSize " F_U64 ")\n", description, tot, genomeSize);
  fprintf(out, "%s num %10" F_U64P "\n", description, cnt);
  fprintf(out, "%s ave %10" F_U64P "\n", description, tot / cnt);
}


void
tgTigSizeAnalysis::printSummary(FILE *out) {
  printSummary(out, "lenSuggestRepeat",   lenSuggestRepeat);
  printSummary(out, "lenSuggestCircular", lenSuggestCircular);

  printSummary(out, "lenUnassembled",     lenUnassembled);
  printSummary(out, "lenBubble",          lenBubble);
  printSummary(out, "lenContig",          lenContig);
}
