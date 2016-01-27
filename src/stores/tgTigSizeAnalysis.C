
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' (http://kmer.sourceforge.net)
 *  both originally distributed by Applera Corporation under the GNU General
 *  Public License, version 2.
 *
 *  Canu branched from Celera Assembler at its revision 4587.
 *  Canu branched from the kmer project at its revision 1994.
 *
 *  This file is derived from:
 *
 *    src/AS_CNS/MultiAlignSizeAnalysis.C
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2012-MAR-26 to 2013-OCT-24
 *      are Copyright 2012-2013 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz on 2014-DEC-22
 *      are Copyright 2014 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2015-OCT-30
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
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
tgTigSizeAnalysis::evaluateTig(tgTig *tig, bool useGapped) {

  //  Try to get the ungapped length.
  //  But revert to the gapped length if that doesn't exist.  This should
  //  only occur for pre-consensus unitigs.

  uint32  length = tig->length(useGapped);

  if (tig->_suggestRepeat)
    lenSuggestRepeat.push_back(length);

  if (tig->_suggestCircular)
    lenSuggestCircular.push_back(length);

  switch (tig->_class) {
    case tgTig_unassembled:   lenUnassembled.push_back(length);  break;
    case tgTig_bubble:        lenBubble.push_back(length);       break;
    case tgTig_contig:        lenContig.push_back(length);       break;
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
tgTigSizeAnalysis::printSummary(FILE *out, char *description, vector<uint32> &data) {
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
      fprintf(out, "%s ng%-3"F_U64P" %10"F_U32P" bp   lg%-3"F_U64P" %6"F_U64P"   sum %10"F_U64P" bp\n",
              description,
              nnn, data[i],
              nnn, i+1,
              sum);

      nnn += 10;
    }
  }

  fprintf(out, "%s sum %10"F_U64P" (genomeSize "F_U64")\n", description, tot, genomeSize);
  fprintf(out, "%s num %10"F_U64P"\n", description, cnt);
  fprintf(out, "%s ave %10"F_U64P"\n", description, tot / cnt);
}


void
tgTigSizeAnalysis::printSummary(FILE *out) {
  printSummary(out, "lenSuggestRepeat",   lenSuggestRepeat);
  printSummary(out, "lenSuggestCircular", lenSuggestCircular);

  printSummary(out, "lenUnassembled",     lenUnassembled);
  printSummary(out, "lenBubble",          lenBubble);
  printSummary(out, "lenContig",          lenContig);
}
