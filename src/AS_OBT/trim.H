
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 2005-2007, J. Craig Venter Institute.
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <time.h>
#include <errno.h>

#include "AS_global.H"
#include "AS_PER_gkpStore.H"

void  doTrim(gkFragment *fr,
             double      minQuality,
             uint32     &left,
             uint32     &right);


//  A simple initialized array -- performs a quality letter -> quality
//  value translation.
//
class qualityLookup {
public:
  qualityLookup() {
    for (uint32 i=0; i<='0'; i++)
      q[i] = 1.0;

    for (uint32 i='1'; i<255; i++)
      q[i] = 1 / pow(10, (i - '0') / 10.0);
  };
  ~qualityLookup() {
  };

  double lookupChar(char x)      { return(q[(uint32)x]);  };
  double lookupNumber(uint32 x)  { return(q['0' + x]);    };

private:
  double   q[255];
};


extern
qualityLookup  qual;

