
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

static const char *rcsid_INCLUDE_AS_BAT_BOGUSUTIL = "$Id: AS_BAT_bogusUtil.H,v 1.1 2010-11-24 01:03:31 brianwalenz Exp $";

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

class longestAlignment {
public:
  longestAlignment() {
    bgn = end = len = num = 0;
    rptBgn = rptEnd = frgLen = 0;
  };

  int32   bgn;  //  Begin coord on the fragment for this region
  int32   end;  //  End coord on the fragment for this region
  int32   len;  //  Length of the match
  int32   num;  //  Number of matches on this region

  int32   rptBgn;
  int32   rptEnd;

  int32   frgLen;
};

class genomeAlignment {
public:
  genomeAlignment() {
    frgIID = frgBgn = frgEnd = 0;
    genBgn = genEnd = 0;
    isReverse = isSpanned = isRepeat = false;
  }
  int32   frgIID;
  int32   frgBgn;
  int32   frgEnd;

  int32   genBgn;
  int32   genEnd;

  bool    isReverse;
  bool    isSpanned;
  bool    isRepeat;
};


bool byFragmentID(const genomeAlignment &A, const genomeAlignment &B);
bool byGenomePosition(const genomeAlignment &A, const genomeAlignment &B);


void addAlignment(vector<genomeAlignment>   &genome,
                  int32  frgIID,
                  int32  frgBgn, int32  frgEnd, bool  isReverse,
                  int32  genBgn, int32  genEnd);

void loadNucmer(char                      *nucmerName,
                vector<longestAlignment>  &longest,
                vector<genomeAlignment>   &genome,
                map<string, int32>        &IIDmap,
                vector<string>            &IIDname,
                int32                     &IIDnext,
                FILE                      *outputFile);

void loadSnapper(char                      *snapperName,
                 vector<longestAlignment>  &longest,
                 vector<genomeAlignment>   &genome,
                 map<string, int32>        &IIDmap,
                 vector<string>            &IIDname,
                 int32                     &IIDnext,
                 FILE                      *outputFile);

#endif  //  INCLUDE_AS_BAT_BOGUSUTIL
