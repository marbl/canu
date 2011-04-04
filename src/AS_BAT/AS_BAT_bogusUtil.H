
/**************************************************************************
 * Copyright (C) 2010, J Craig Venter Institute. All rights reserved.
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

#ifndef INCLUDE_AS_BAT_BOGUSUTIL
#define INCLUDE_AS_BAT_BOGUSUTIL

static const char *rcsid_INCLUDE_AS_BAT_BOGUSUTIL = "$Id: AS_BAT_bogusUtil.H,v 1.7 2011-04-04 14:24:00 brianwalenz Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>

#include <string>
#include <vector>
#include <map>
#include <algorithm>

using namespace std;

#include "AS_global.h"
#include "AS_UTL_splitToWords.H"
#include "AS_UTL_intervalList.H"


class genomeAlignment {
public:
  genomeAlignment() {
    frgIID = frgBgn = frgEnd = 0;
    genIID = 0;
    genBgn = genEnd = 0;
    identity = 0.0;
    isDeleted = isReverse = isSpanned = isRepeat = false;
  };

  int32   frgIID;
  int32   frgBgn;
  int32   frgEnd;

  int32   genIID;      //  Position in the actual sequence
  int32   genBgn;
  int32   genEnd;

  int32   chnBgn;      //  Position in the chained sequences
  int32   chnEnd;

  double  identity;    //  Percent identity of the alignment

  bool    isDeleted;   //  Used by bogusness

  bool    isReverse;
  bool    isSpanned;
  bool    isRepeat;
};


class referenceSequence {
public:
  referenceSequence(int32 cb, int32 ce, int32 rl, char *cn) {
    rschnBgn = cb;
    rschnEnd = ce;
    rsrefLen = rl;
    assert(strlen(cn) < 256);
    strcpy(rsrefName, cn);
  };

  int32  rschnBgn;
  int32  rschnEnd;
  int32  rsrefLen;
  char   rsrefName[256];
};


bool byFragmentID(const genomeAlignment &A, const genomeAlignment &B);
bool byGenomePosition(const genomeAlignment &A, const genomeAlignment &B);


void addAlignment(vector<genomeAlignment>   &genome,
                  int32  frgIID,
                  int32  frgBgn, int32  frgEnd, bool  isReverse,
                  int32  chnBgn, int32  chnEnd,
                  double identity,
                  int32  genIID,
                  int32  genBgn, int32  genEnd);

void loadNucmer(char                       *nucmerName,
                vector<genomeAlignment>    &genome,
                map<string, int32>         &IIDmap,
                vector<string>             &IIDname,
                vector<referenceSequence>  &refList,
                map<string,uint32>         &refMap);

void loadSnapper(char                       *snapperName,
                 vector<genomeAlignment>    &genome,
                 map<string, int32>         &IIDmap,
                 vector<string>             &IIDname,
                 vector<referenceSequence>  &refList,
                 map<string,uint32>         &refMap);

void
loadReferenceSequence(char                       *refName,
                      vector<referenceSequence>  &refList,
                      map<string,uint32>         &refMap);

#endif  //  INCLUDE_AS_BAT_BOGUSUTIL
