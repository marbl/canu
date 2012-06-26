
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

static const char *rcsid = "$Id: MultiAlignSizeAnalysis.C,v 1.2 2012-06-26 14:08:30 brianwalenz Exp $";

#include "MultiAlignSizeAnalysis.H"

#include <math.h>

#include <map>
#include <set>
#include <vector>
#include <algorithm>

using namespace std;

//  stored in the store, not the ma
//
//  unitig_status
//  AS_UNASSIGNED
//  AS_UNIQUE
//  AS_SEP
//  AS_NOTREZ
//
//  contig_status
//  AS_UNPLACED
//  AS_PLACED
//
//  stored in the u_list of each contig
//
//  unitig_type
//  AS_UNIQUE_UNITIG
//  AS_SINGLE_UNITIG
//  AS_ROCK_UNITIG
//  AS_STONE_UNITIG
//  AS_PEBBLE_UNITIG


sizeAnalysis::sizeAnalysis(uint64 genomeSize_) {
  genomeSize = genomeSize_;
}

sizeAnalysis::~sizeAnalysis() {
}

void
sizeAnalysis::evaluateTig(MultiAlignT *ma, bool isUnitig) {

  //  Try to get the ungapped length.
  uint32 length = GetMultiAlignUngappedLength(ma);

  if (length == 0)
    //  But revert to the gapped length if that doesn't exist.  This should
    //  only occur for pre-consensus unitigs.
    length = GetMultiAlignLength(ma);

  if      ((isUnitig == true) && (ma->data.unitig_status == AS_UNASSIGNED))
    utgLenUnassigned.push_back(length);

  else if ((isUnitig == true) && (ma->data.unitig_status == AS_UNIQUE))
    utgLenUnique.push_back(length);

  else if ((isUnitig == true) && (ma->data.unitig_status == AS_SEP))
    utgLenSep.push_back(length);

  else if ((isUnitig == true) && (ma->data.unitig_status == AS_NOTREZ))
    utgLenNotRez.push_back(length);

  else if ((isUnitig == false) && (ma->data.contig_status == AS_UNPLACED))
    ctgLenUnplaced.push_back(length);

  else if ((isUnitig == false) && (ma->data.contig_status == AS_PLACED))
    ctgLenPlaced.push_back(length);

  if (ma->data.num_frags == 1)
    tigLenSingleton.push_back(length);
  else
    tigLenAssembled.push_back(length);
}

void
sizeAnalysis::finalize(void) {
  sort(utgLenUnassigned.rbegin(), utgLenUnassigned.rend());
  sort(utgLenUnique.rbegin(),     utgLenUnique.rend());
  sort(utgLenSep.rbegin(),        utgLenSep.rend());
  sort(utgLenNotRez.rbegin(),     utgLenNotRez.rend());

  sort(ctgLenPlaced.rbegin(),     ctgLenPlaced.rend());
  sort(ctgLenUnplaced.rbegin(),   ctgLenUnplaced.rend());

  sort(tigLenSingleton.rbegin(),  tigLenSingleton.rend());
  sort(tigLenAssembled.rbegin(),  tigLenAssembled.rend());

}

void
sizeAnalysis::printSummary(FILE *out, char *description, vector<uint32> &data) {
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

    if (siz * nnn / 100 < sum) {
      fprintf(out, "%s n%2"F_U64P" siz %10"F_U32P" sum %10"F_U64P" idx %10"F_U64P"\n",
              description, nnn, data[i], sum, i);
      
      nnn += 10;
    }
  }

  fprintf(out, "%s sum %10"F_U64P" (genomeSize "F_U64")\n",
          description, tot, genomeSize);
  fprintf(out, "%s num %10"F_U64P"\n",
          description, cnt);
  fprintf(out, "%s ave %10"F_U64P"\n",
          description, tot / cnt);
}


void
sizeAnalysis::printSummary(FILE *out) {
  printSummary(out, "utgLenUnassigned", utgLenUnassigned);
  printSummary(out, "utgLenUnique",     utgLenUnique);
  printSummary(out, "utgLenSep",        utgLenSep);
  printSummary(out, "utgLenNotRez",     utgLenNotRez);

  printSummary(out, "ctgLenPlaced",     ctgLenPlaced);
  printSummary(out, "ctgLenUnplaced",   ctgLenUnplaced);

  printSummary(out, "tigLenSingleton",  tigLenSingleton);
  printSummary(out, "tigLenAssembled",  tigLenAssembled);
}
