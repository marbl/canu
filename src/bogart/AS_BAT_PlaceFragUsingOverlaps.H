
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
 *    src/AS_BAT/AS_BAT_PlaceFragUsingOverlaps.H
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2010-DEC-06 to 2013-AUG-01
 *      are Copyright 2010-2013 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2014-DEC-19 to 2015-JUN-03
 *      are Copyright 2014-2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2016-JAN-11
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#ifndef INCLUDE_AS_BAT_PLACEFRAGUSINGOVERLAPS
#define INCLUDE_AS_BAT_PLACEFRAGUSINGOVERLAPS

#include "AS_BAT_OverlapCache.H"
#include "AS_BAT_BestOverlapGraph.H"  //  For FragmentEnd
#include "AS_BAT_Unitig.H"            //  For SeqInterval
#include "AS_BAT_UnitigVector.H"


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
  if (A.position.isReverse() != B.position.isReverse())
    return(A.position.isReverse() < B.position.isReverse());
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
                       double                    erate,
                       Unitig                   *target,
                       uint32                    fid,
                       vector<overlapPlacement> &placements);


#endif  //  INCLUDE_AS_BAT_PLACEFRAGUSINGOVERLAPS
