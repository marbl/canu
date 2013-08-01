
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

static const char *rcsid_INCLUDE_AS_BAT_PLACEFRAGUSINGOVERLAPS = "$Id$";

#include "AS_BAT_OverlapCache.H"

class overlapPlacement {
public:
  overlapPlacement() {
    frgID                 = 0;
    refID                 = 0;
    
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
  uint32            refID;      //  ID if the fragment we used to place this frag (the reference).

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


//  Sort by:  tigID, orientation, position
//
//  This sort is used to cluster the reads into overlapping regions.  We don't care
//  about ties.
inline
bool
overlapPlacement_byLocation(const overlapPlacement &A, const overlapPlacement &B) {
  if (A.tigID != B.tigID)
    return(A.tigID < B.tigID);
  if (isReverse(A.position) != isReverse(B.position))
    return(isReverse(A.position) < isReverse(B.position));
  return(A.position < B.position);
}


//  Sort by:
//    cluster
//
//  This sort is used to group reads by cluster.  We don't care about ties, but they
//  can change the results if the input overlaps change.
inline
bool
overlapPlacement_byCluster(const overlapPlacement &A, const overlapPlacement &B) {
  return(A.clusterID < B.clusterID);
}


bool
placeFragUsingOverlaps(UnitigVector             &unitigs,
                       Unitig                   *target,
                       AS_IID                    fid,
                       vector<overlapPlacement> &placements);

void
placeFragInBestLocation(UnitigVector   &unitigs,
                        AS_IID          fid);


#endif  //  INCLUDE_AS_BAT_PLACEFRAGUSINGOVERLAPS
