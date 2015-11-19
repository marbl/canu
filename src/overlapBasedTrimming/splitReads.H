
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
 *  Modifications by:
 *
 *    Brian P. Walenz from 2015-APR-20 to 2015-JUN-16
 *      are Copyright 2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#ifndef SPLIT_READS_H
#define SPLIT_READS_H

#include "AS_global.H"

#include "gkStore.H"
#include "ovStore.H"

#include "adjustOverlaps.H"
#include "clearRangeFile.H"

#include "intervalList.H"

#include "AS_UTL_decodeRange.H"

#include <vector>
#include <map>

using namespace std;


//  But then use only overlaps larger than this for some of the more questionable overlaps.
#define MIN_INTERVAL_OVERLAP  60

//  Same orient overlaps closer than this are evidence of PacBio subreads.
#define SUBREAD_LOOP_MAX_SIZE  500
#define SUBREAD_LOOP_EXT_SIZE  2000

//  Original
//#define SUBREAD_LOOP_MAX_SIZE  100
//#define SUBREAD_LOOP_EXT_SIZE  0


//  WITH_REPORT_FULL will ALL overlap evidence.
//  REPORT_OVERLAPS  will print the incoming overlaps in the log.
//
#undef WITH_REPORT_FULL
#undef REPORT_OVERLAPS

#undef DEBUG_ISLINKER
#undef DEBUG_INTERVAL
#undef DEBUG_FILTER



//  An overlap, adjusted for trimming done already.
//
class adjOverlap {
public:
  uint32   a_iid;
  uint32   b_iid;

  uint32   flipped;           //  The b read is 3' to 5'

  uint32   aclrbgn, aclrend;  //  Clear range of the two reads
  uint32   bclrbgn, bclrend;  //  B read coords are adjusted for any flip

  uint32   aovlbgn, aovlend;  //  Overlap extent, adjusted for the clear range
  uint32   bovlbgn, bovlend;

  uint32   aOVLbgn, aOVLend;  //  Original unadjusted overlap
  uint32   bOVLbgn, bOVLend;
};


const uint32 badType_nothing = 0;
const uint32 badType_5spur   = 1;
const uint32 badType_3spur   = 2;
const uint32 badType_chimera = 3;
const uint32 badType_subread = 4;


//  A region that is bad.
//
class badRegion {
public:
  badRegion() {
    id   = 0;
    type = 0;
    bgn  = 0;
    end  = 0;
  };

  badRegion(uint32 id_, uint32 type_, uint32 bgn_, uint32 end_) {
    id   = id_;
    type = type_;
    bgn  = bgn_;
    end  = end_;
  };

  const char  *typeName(void) const {
    const char *N[5] = { "nothing", "5'spur", "3'spur", "chimera", "subread" };

    return(N[type]);
  }


public:
  uint32    id;    //  Read ID for this bad region
  uint32    type;  //  Type of bad region - 5'spur, 3'spur, chimera, subread
  uint32    bgn;   //
  uint32    end;   //
};






class workUnit {
public:
  workUnit(uint32 id_=0, uint32 iniBgn_=0, uint32 iniEnd_=0) {
    clear(id_, iniBgn_, iniEnd_);

    logMsg[0] = 0;

    blist.clear();

    adjLen = 0;
    adjMax = 0;
    adj    = NULL;
  };
  ~workUnit() {
    delete [] adj;
  };


  void          clear(uint32 id_, uint32 iniBgn_, uint32 iniEnd_) {
    id          = id_;
    iniBgn      = iniBgn_;
    iniEnd      = iniEnd_;

    isOK        = true;
    isBad       = false;

    isSpur5     = false;
    isSpur3     = false;
    isChimera   = false;
    isSubread   = false;


    clrBgn      = iniBgn_;
    clrEnd      = iniEnd_;

    logMsg[0]   = 0;

    blist.clear();

    adjLen      = 0;
    //adjMax      = 0;      //  Do NOT reset the allocated space.
    //adj         = NULL;
  };

  void          addAndFilterOverlaps(gkStore *gkp,
                                     clearRangeFile *finClr,
                                     double errorRate,
                                     ovOverlap *ovl, uint32 ovlLen);

public:
  uint32        id;         //  Read ID
  uint32        iniBgn;     //  The input clear range
  uint32        iniEnd;

  //  Results

  bool          isOK;       //  Read is acceptable, might have been fixed as below
  bool          isBad;      //  Read couldn't be fixed

  bool          isSpur5;    //  Spur trimmed on 5' end
  bool          isSpur3;    //  Spur trimmed on 3' end
  bool          isChimera;  //  Read split because of suspected chimeric region
  bool          isSubread;  //  Read split because of suspected subreads

  uint32        clrBgn;     //  The final clear range
  uint32        clrEnd;

  char          logMsg[1024];

  //  Work space

  vector<badRegion>  blist;

  //  Overlaps

  uint32        adjLen;
  uint32        adjMax;
  adjOverlap   *adj;
};







void
detectSpur(gkStore               *gkp,
           workUnit              *w,
           FILE                  *subreadFile,
           bool                   subreadFileVerbose);

void
detectChimer(gkStore               *gkp,
             workUnit              *w,
             FILE                  *subreadFile,
             bool                   subreadFileVerbose);

void
detectSubReads(gkStore               *gkp,
               workUnit              *w,
               FILE                  *subreadFile,
               bool                   subreadFileVerbose);

void
trimBadInterval(gkStore               *gkp,
                workUnit              *w,
                uint32                 minReadLength,
                FILE                  *subreadFile,
                bool                   subreadFileVerbose);



#endif  //  SPLIT_READS_H
