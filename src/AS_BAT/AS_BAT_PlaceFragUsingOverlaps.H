
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 1999-2004, The Venter Institute. All rights reserved.
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

#ifndef INCLUDE_AS_BAT_PLACEFRAGUSINGOVERLAPS
#define INCLUDE_AS_BAT_PLACEFRAGUSINGOVERLAPS

static const char *rcsid_INCLUDE_AS_BAT_PLACEFRAGUSINGOVERLAPS = "$Id: AS_BAT_PlaceFragUsingOverlaps.H,v 1.6 2011-04-04 14:25:31 brianwalenz Exp $";

#include "AS_BAT_OverlapCache.H"

class overlapPlacement {
public:
  overlapPlacement() {
    frgID                 = 0;
    ovlID                 = 0;
    
    tigID                 = 0;
    position.bgn          = 0;
    position.end          = 0;

    verified.bgn          = 0;
    verified.end          = 0;

    bgnStdDev             = 0.0;
    endStdDev             = 0.0;

    clusterID             = 0;

    fCoverage             = 0.0;

    errors                = 0.0;
    aligned               = 0;
    covered.bgn           = 0;
    covered.end           = 0;

    frag5p                = FragmentEnd();
    frag3p                = FragmentEnd();
  };
  ~overlapPlacement() {
  };

public:
  uint32            frgID;      //  ID of the fragment this position is for.
  uint32            ovlID;      //  ID if the fragment we used to place this frag.

  uint32            tigID;      //  Unitig ID of this placement
  SeqInterval       position;   //  Unitig position of this placement

  SeqInterval       verified;   //  Unitig position of this placement, verified by overlaps

  double            bgnStdDev;  //  Standard deviation of position.bgn
  double            endStdDev;  //  Standard deviation of position.end

  int32             clusterID;

  double            fCoverage;  //  Coverage of the fragment

  double            errors;     //  number of errors in alignments
  uint32            aligned;    //  number of bases in alignments
  SeqInterval       covered;    //  Position of the overlap on the read

  FragmentEnd       frag5p;     //  First unitig fragment that supports this placement
  FragmentEnd       frag3p;     //  Last unitig fragment that supports this placement
};


//  Sort by:
//    tigID
//    orientation
//    position
//
//  The extra ovlID, aligned and errors are only to break ties.  Without them, byCluster() returns
//  ties that change the results of RETURN_ALL_PLACEMENTS in AS_BAT_MergeSplitJoin.C
//
inline
bool
overlapPlacement_byLocation(const overlapPlacement &A, const overlapPlacement &B) {
  if (A.tigID != B.tigID)
    return(A.tigID < B.tigID);
  if (isReverse(A.position) != isReverse(B.position))
    return(isReverse(A.position) < isReverse(B.position));
  if (A.position != B.position)
    return(A.position < B.position);
  if (A.ovlID != B.ovlID)
    return(A.ovlID < B.ovlID);
  if (A.aligned != B.aligned)
    return(A.aligned < B.aligned);
  return(A.errors < B.errors);
}


inline
bool
overlapPlacement_byCluster(const overlapPlacement &A, const overlapPlacement &B) {
  if (A.clusterID != B.clusterID)
    return(A.clusterID < B.clusterID);
  return(overlapPlacement_byLocation(A, B));
}


bool
placeFragUsingOverlaps(UnitigVector             &unitigs,
                       Unitig                   *target,
                       uint32                    fid,
                       vector<overlapPlacement> &placements);

#endif  //  INCLUDE_AS_BAT_PLACEFRAGUSINGOVERLAPS
