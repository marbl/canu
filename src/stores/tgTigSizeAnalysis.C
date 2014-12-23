
/**************************************************************************
 * Copyright (C) 2011, J Craig Venter Institute. All rights reserved.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received (LICENSE.txt) a copy of the GNU General Public
 * License along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *************************************************************************/

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
tgTigSizeAnalysis::evaluateTig(tgTig *tig) {

  //  Try to get the ungapped length.
  //  But revert to the gapped length if that doesn't exist.  This should
  //  only occur for pre-consensus unitigs.

  uint32 length = tig->ungappedLength();

  if (length == 0)
    length = tig->layoutLength();


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
