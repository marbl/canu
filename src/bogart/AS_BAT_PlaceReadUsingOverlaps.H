
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' r4587 (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' r1994 (http://kmer.sourceforge.net)
 *
 *  Except as indicated otherwise, this is a 'United States Government Work',
 *  and is released in the public domain.
 *
 *  File 'README.licenses' in the root directory of this distribution
 *  contains full conditions and disclaimers.
 */

#ifndef INCLUDE_AS_BAT_PLACEREADUSINGOVERLAPS
#define INCLUDE_AS_BAT_PLACEREADUSINGOVERLAPS

#include "AS_BAT_Unitig.H"            //  For SeqInterval
#include "AS_BAT_TigVector.H"


class overlapPlacement {
public:
  overlapPlacement() {
  };
  overlapPlacement(uint32 f, uint32 t, SeqInterval c) :
    frgID(f),
    tigID(t),
    covered(c) {
  };

  uint32        frgID     = 0;            //  Read ID of the read this position is for.

  uint32        refID     = 0;            //  Read ID of the overlapping read were placed with.
  uint32        tigID     = 0;            //  Unitig ID of this placement

  int32         clusterID = UINT32_MAX;   //  Used for grouping individual placements into mergable order.
                                              //  UINT32_MAX indicates we don't want this placement in a cluster.

  SeqInterval   position;                 //  Unitig position of this placement
  SeqInterval   verified;                 //  Unitig position of this placement, verified by overlaps
  SeqInterval   covered;                  //  Position of the overlap on the read

  double        fCoverage = 0.0;          //  Coverage of the read

  double        errors    = 0.0;          //  number of errors in alignments
  uint32        aligned   = 0;            //  number of bases in alignments

  double        erate(void) const   { return((double)errors / aligned); };

  uint32        tigFidx   = UINT32_MAX;   //  First unitig read that supports this placement
  uint32        tigLidx   = UINT32_MAX;   //  Last  unitig read that supports this placement
};



const uint32  placeRead_all        = 0x00;   //  Return all alignments
const uint32  placeRead_fullMatch  = 0x01;   //  Return only alignments for the whole read
const uint32  placeRead_noExtend   = 0x02;   //  Return only alignments contained in the tig

bool
placeReadUsingOverlaps(TigVector                &tigs,
                       Unitig                   *target,
                       uint32                    fid,
                       vector<overlapPlacement> &placements,
                       uint32                    flags = placeRead_all,
                       double                    errorLimit = 1.0);


#endif  //  INCLUDE_AS_BAT_PLACEREADUSINGOVERLAPS
