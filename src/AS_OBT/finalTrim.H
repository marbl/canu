
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 2005-2011, J. Craig Venter Institute.
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

#ifndef FINALTRIM_H
#define FINALTRIM_H

static const char *rcsid_FINALTRIM_H = "$Id$";

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <assert.h>

#include "trim.H"

#include "AS_global.H"
#include "AS_PER_gkpStore.H"
#include "AS_OVS_overlapStore.H"
#include "AS_UTL_intervalList.H"


#define OBT_MODE_WIGGLE      (5)

#define OBT_CQ_LENGTH        (100)
#define OBT_CQO_LENGTH       (200)
#define OBT_CQO_OVERLAP      (100)
#define OBT_CQ_SHORT         (5)

#define OBT_QLT_CLOSE_5      (10)  // 5,6  use 5'mode, use 5'min>1
#define OBT_QLT_FAR_5        (50)  // 11   use min5'
#define OBT_QLT_MODE3        (150) // 9    use 3'mode
#define OBT_QLT_CLOSE_MAXM3  (30)  // 14   use max>1 close to max
#define OBT_QLT_CLOSE_MAX3   (100) // 12   use max3'


bool
withoutOverlaps(OVSoverlap  *ovl,
                uint32       ovlLen,
                gkFragment  &fr,
                uint32       ibgn,
                uint32       iend,
                uint32      &fbgn,
                uint32      &fend,
                char        *logMsg,
                uint32       errorRate,
                double       errorLimit,
                bool         qvTrimAllowed);

bool
largestCovered(OVSoverlap  *ovl,
               uint32       ovlLen,
               gkFragment  &fr,
               uint32       ibgn,
               uint32       iend,
               uint32      &fbgn,
               uint32      &fend,
               char        *logMsg,
               uint32       errorRate,
               double       errorLimit,
               bool         qvTrimAllowed,
               uint32       minOverlap,
               uint32       minCoverage);

bool
bestEdge(OVSoverlap  *ovl,
         uint32       ovlLen,
         gkFragment  &fr,
         uint32       ibgn,
         uint32       iend,
         uint32      &fbgn,
         uint32      &fend,
         char        *logMsg,
         uint32       errorRate,
         double       errorLimit,
         bool         qvTrimAllowed,
         uint32       minOverlap,
         uint32       minCoverage);

bool
evidenceBased(OVSoverlap  *ovl,
              uint32       ovlLen,
              gkFragment  &fr,
              uint32       ibgn,
              uint32       iend,
              uint32      &fbgn,
              uint32      &fend,
              char        *logMsg,
              uint32       errorRate,
              double       errorLimit,
              bool         qvTrimAllowed);

#endif  //  FINALTRIM_H
