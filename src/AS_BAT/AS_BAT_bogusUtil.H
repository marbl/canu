
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

static const char *rcsid_INCLUDE_AS_BAT_BOGUSUTIL = "$Id: AS_BAT_bogusUtil.H,v 1.4 2010-12-17 09:57:12 brianwalenz Exp $";

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
    frgIID = frgLen = frgBgn = frgEnd = 0;
    genBgn = genEnd = 0;
    isDeleted = isReverse = isSpanned = isRepeat = false;
  }

  int32   frgIID;
  int32   frgLen;  //  UNUSED
  int32   frgBgn;
  int32   frgEnd;

  int32   genBgn;
  int32   genEnd;

  bool    isDeleted;  //  Used by bogusness

  bool    isReverse;
  bool    isSpanned;
  bool    isRepeat;
};


bool byFragmentID(const genomeAlignment &A, const genomeAlignment &B);
bool byGenomePosition(const genomeAlignment &A, const genomeAlignment &B);


void addAlignment(vector<genomeAlignment>   &genome,
                  int32  frgIID, int32  frgLen,
                  int32  frgBgn, int32  frgEnd, bool  isReverse,
                  int32  genBgn, int32  genEnd);

void loadNucmer(char                      *nucmerName,
                vector<genomeAlignment>   &genome,
                map<string, int32>        &IIDmap,
                vector<string>            &IIDname,
                vector<uint32>            &IIDcount);

void loadSnapper(char                      *snapperName,
                 vector<genomeAlignment>   &genome,
                 map<string, int32>        &IIDmap,
                 vector<string>            &IIDname,
                 vector<uint32>            &IIDcount);

#endif  //  INCLUDE_AS_BAT_BOGUSUTIL
