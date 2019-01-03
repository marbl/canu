
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
 *    src/bogart/AS_BAT_PlaceFragUsingOverlaps.H
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

#ifndef INCLUDE_AS_BAT_PLACEREADUSINGOVERLAPS
#define INCLUDE_AS_BAT_PLACEREADUSINGOVERLAPS

#include "AS_BAT_OverlapCache.H"
#include "AS_BAT_BestOverlapGraph.H"  //  For ReadEnd
#include "AS_BAT_Unitig.H"            //  For SeqInterval
#include "AS_BAT_TigVector.H"


class overlapPlacement {
public:
  overlapPlacement(uint32 fi=0) {
    frgID         = fi;

    refID         = 0;
    tigID         = 0;

    clusterID     = 0;

    position      = SeqInterval();
    verified      = SeqInterval();
    covered       = SeqInterval();

    fCoverage     = 0.0;

    errors        = 0.0;
    aligned       = 0;

    tigFidx       = UINT32_MAX;
    tigLidx       = 0;
  };

  overlapPlacement(uint32 fid, overlapPlacement &op) {
    frgID          = fid;

    refID          = UINT32_MAX;       //  Not valid in the output overlapPlacement.
    tigID          = op.tigID;

    clusterID      = op.clusterID;     //  Useless to track forward.

    position.bgn   = 0;
    position.end   = 0;

    verified.bgn   = 0;
    verified.end   = 0;

    covered.bgn    = op.covered.bgn;
    covered.end    = op.covered.end;

    fCoverage      = 0.0;

    errors         = 0.0;
    aligned        = 0;

    tigFidx        = UINT32_MAX;
    tigLidx        = UINT32_MAX;
  };

  ~overlapPlacement() {
  };


public:
  uint32            frgID;      //  Read ID of the read this position is for.

  uint32            refID;      //  Read ID of the overlapping read were placed with.
  uint32            tigID;      //  Unitig ID of this placement

  int32             clusterID;

  SeqInterval       position;   //  Unitig position of this placement
  SeqInterval       verified;   //  Unitig position of this placement, verified by overlaps
  SeqInterval       covered;    //  Position of the overlap on the read

  double            fCoverage;  //  Coverage of the read

  double            errors;     //  number of errors in alignments
  uint32            aligned;    //  number of bases in alignments

  uint32            tigFidx;    //  First unitig read that supports this placement
  uint32            tigLidx;    //  Last  unitig read that supports this placement
};


//  Sort by:  clusterID, tigID, orientation, position
//
//  This sort is used to cluster the reads into overlapping regions.  We don't care
//  about ties.
//
//  clusterID is UINT32_MAX if the placement should be ignored.
//
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


const uint32  placeRead_all        = 0x00;   //  Return all alignments
const uint32  placeRead_fullMatch  = 0x01;   //  Return only alignments for the whole read
const uint32  placeRead_noExtend   = 0x02;   //  Return only alignments contained in the tig

bool
placeReadUsingOverlaps(TigVector                &tigs,
                       Unitig                   *target,
                       uint32                    fid,
                       vector<overlapPlacement> &placements,
                       uint32                    flags = placeRead_all);


#endif  //  INCLUDE_AS_BAT_PLACEREADUSINGOVERLAPS
