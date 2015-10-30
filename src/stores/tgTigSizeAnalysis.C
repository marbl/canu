
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
 *    Brian P. Walenz beginning on 2014-DEC-22
 *      are Copyright 2014 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

static const char *rcsid = "$Id$";

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

  if (tig->_suggestUnique)
    lenSuggestUnique.push_back(length);

  if (tig->_suggestCircular)
    lenSuggestCircular.push_back(length);

  if (tig->_suggestHaploid)
    lenSuggestHaploid.push_back(length);

  if (tig->numberOfChildren() == 1)
    lenSingleton.push_back(length);
  else
    lenAssembled.push_back(length);
}

void
tgTigSizeAnalysis::finalize(void) {

  sort(lenSuggestRepeat.rbegin(),   lenSuggestRepeat.rend());
  sort(lenSuggestUnique.rbegin(),   lenSuggestUnique.rend());
  sort(lenSuggestCircular.rbegin(), lenSuggestCircular.rend());
  sort(lenSuggestHaploid.rbegin(),  lenSuggestHaploid.rend());

  sort(lenSingleton.rbegin(), lenSingleton.rend());
  sort(lenAssembled.rbegin(), lenAssembled.rend());
}

void
tgTigSizeAnalysis::printSummary(FILE *out, char *description, vector<uint32> &data) {
  uint64  cnt = data.size();
  uint64  sum = 0;
  uint64  tot = 0;
  uint64  nnn = 10;
  uint64  siz = 0;

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
      fprintf(out, "%s n%2"F_U64P" siz %10"F_U32P" sum %10"F_U64P" idx %10"F_U64P"\n",
              description, nnn, data[i], sum, i);

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
  printSummary(out, "lenSuggestUnique",   lenSuggestUnique);
  printSummary(out, "lenSuggestCircular", lenSuggestCircular);
  printSummary(out, "lenSuggestHaploid",  lenSuggestHaploid);

  printSummary(out, "lenSingleton",  lenSingleton);
  printSummary(out, "lenAssembled",  lenAssembled);
}
