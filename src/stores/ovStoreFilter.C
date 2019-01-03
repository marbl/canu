
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
 *    src/AS_OVS/AS_OVS_overlapStore.C
 *    src/AS_OVS/AS_OVS_overlapStore.c
 *    src/stores/ovStore.C
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2007-MAR-08 to 2013-AUG-01
 *      are Copyright 2007-2013 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Sergey Koren on 2007-MAY-08
 *      are Copyright 2007 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Sergey Koren from 2011-JUN-02 to 2011-JUN-03
 *      are Copyright 2011 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Gregory Sims from 2012-FEB-01 to 2012-FEB-14
 *      are Copyright 2012 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2014-DEC-09 to 2015-AUG-14
 *      are Copyright 2014-2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2015-OCT-12
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *    Sergey Koren beginning on 2015-DEC-15
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "ovStore.H"


#define OBT_FAR5PRIME        (29)
#define OBT_MIN_LENGTH       (75)




ovStoreFilter::ovStoreFilter(sqStore *seq_, double maxErate_, bool beVerbose_) {
  seq             = seq_;
  maxID           = seq->sqStore_getNumReads();
  maxEvalue       = AS_OVS_encodeEvalue(maxErate_);

  beVerbose       = beVerbose_;

  resetCounters();

  skipReadOBT     = new char [maxID + 1];
  skipReadDUP     = new char [maxID + 1];

  uint32  numSkipOBT = 0;
  uint32  numSkipDUP = 0;

  for (uint64 iid=0; iid<=maxID; iid++) {
    uint32     Lid = seq->sqStore_getRead(iid)->sqRead_libraryID();
    sqLibrary *L   = seq->sqStore_getLibrary(Lid);

    skipReadOBT[iid] = false;
    skipReadDUP[iid] = false;

    if ((L->sqLibrary_removeDuplicateReads()     == false) &&
        (L->sqLibrary_finalTrim()                == SQ_FINALTRIM_NONE) &&
        (L->sqLibrary_removeSpurReads()          == false) &&
        (L->sqLibrary_removeChimericReads()      == false) &&
        (L->sqLibrary_checkForSubReads()         == false)) {
      numSkipOBT++;
      skipReadOBT[iid] = true;
    }

    if (L->sqLibrary_removeDuplicateReads() == false) {
      numSkipDUP++;
      skipReadDUP[iid] = true;
    }
  }

  if (numSkipOBT > 0)
    fprintf(stderr, "-  Marked " F_U32 " reads to skip trimming.\n", numSkipOBT);

  if (numSkipDUP > 0)
    fprintf(stderr, "-  Marked " F_U32 " reads to skip deduplication.\n", numSkipDUP);
}



ovStoreFilter::~ovStoreFilter() {
  delete [] skipReadOBT;
  delete [] skipReadDUP;
}





//  Are the 5' end points very different?  If the overlap is flipped, then, yes, they are.
static
bool
isOverlapDifferent(ovOverlap &ol) {
  bool   isDiff = true;

  if (ol.flipped() == false) {
    if (ol.a_bgn() > ol.b_bgn())
      isDiff = ((ol.a_bgn() - ol.b_bgn()) > OBT_FAR5PRIME) ? (true) : (false);
    else
      isDiff = ((ol.b_bgn() - ol.a_bgn()) > OBT_FAR5PRIME) ? (true) : (false);
  }

  return(isDiff);
}


//  Is the overlap long?
static
bool
isOverlapLong(ovOverlap &ol) {
  int32 ab    = ol.a_bgn();
  int32 ae    = ol.a_end();
  int32 bb    = ol.b_bgn();
  int32 be    = ol.b_end();

  int32 Alength = ae - ab;
  int32 Blength = be - bb;

  if (be < bb)
    Blength = bb - be;

  return(((Alength > OBT_MIN_LENGTH) && (Blength > OBT_MIN_LENGTH)) ? (true) : (false));
}



void
ovStoreFilter::filterOverlap(ovOverlap       &foverlap,
                             ovOverlap       &roverlap) {

  //  GREATLY annoy the poor user that asked for 'overly verbose' mode.

  if (beVerbose) {
    char ovlstr[256];

    fprintf(stderr, "%s\n", foverlap.toString(ovlstr, ovOverlapAsUnaligned, false));
  }

  //  Quick sanity check on IIDs.

  if ((foverlap.a_iid == 0) ||
      (foverlap.b_iid == 0) ||
      (foverlap.a_iid > maxID) ||
      (foverlap.b_iid > maxID)) {
    char ovlstr[256];

    fprintf(stderr, "Overlap has IDs out of range (maxID " F_U32 "), possibly corrupt input data.\n", maxID);
    fprintf(stderr, "  coords -- %s\n", foverlap.toString(ovlstr, ovOverlapAsCoords, false));
    fprintf(stderr, "  hangs  -- %s\n", foverlap.toString(ovlstr, ovOverlapAsHangs, false));
    exit(1);
  }

  //  Make the reverse overlap (important, AFTER resetting the erate-based 'for' flags).

  roverlap.swapIDs(foverlap);

  //  Ignore high error overlaps

  if ((foverlap.evalue() > maxEvalue)) {
    foverlap.dat.ovl.forUTG = false;
    foverlap.dat.ovl.forOBT = false;
    foverlap.dat.ovl.forDUP = false;

    roverlap.dat.ovl.forUTG = false;
    roverlap.dat.ovl.forOBT = false;
    roverlap.dat.ovl.forDUP = false;

    skipERATE++;
    skipERATE++;
  }

  //  Ignore opposite oriented overlaps
#ifdef IGNORE_FLIPPED_OVERLAPS
  if ((foverlap.flipped() == true)) {
    foverlap.dat.ovl.forUTG = false;
    foverlap.dat.ovl.forOBT = false;
    foverlap.dat.ovl.forDUP = false;

    roverlap.dat.ovl.forUTG = false;
    roverlap.dat.ovl.forOBT = false;
    roverlap.dat.ovl.forDUP = false;

    skipFLIPPED++;
    skipFLIPPED++;
  }
#endif

  //  Don't OBT if not requested.

  if ((foverlap.dat.ovl.forOBT == false) && (skipReadOBT[foverlap.a_iid] == true)) {
    foverlap.dat.ovl.forOBT = false;
    skipOBT++;
  }

  if ((roverlap.dat.ovl.forOBT == false) && (skipReadOBT[roverlap.a_iid] == true)) {
    roverlap.dat.ovl.forOBT = false;
    skipOBT++;
  }

  //  If either overlap is good for either obt or dup, compute if it is different and long.  These
  //  are the same for both foverlap and roverlap.

  bool  isDiff = isOverlapDifferent(foverlap);
  bool  isLong = isOverlapLong(foverlap);

  //  Remove the bad-for-OBT overlaps.

  if ((isDiff == false) && (foverlap.dat.ovl.forOBT == true)) {
    foverlap.dat.ovl.forOBT = false;
    skipOBTbad++;
  }

  if ((isDiff == false) && (roverlap.dat.ovl.forOBT == true)) {
    roverlap.dat.ovl.forOBT = false;
    skipOBTbad++;
  }

  //  Remove the too-short-for-OBT overlaps.

  if ((isLong == false) && (foverlap.dat.ovl.forOBT == true)) {
    foverlap.dat.ovl.forOBT = false;
    skipOBTshort++;
  }

  if ((isLong == false) && (roverlap.dat.ovl.forOBT == true)) {
    roverlap.dat.ovl.forOBT = false;
    skipOBTshort++;
  }

  //  Don't dedupe if not requested.

  if ((foverlap.dat.ovl.forDUP == true) && (skipReadDUP[foverlap.a_iid] == true)) {
    foverlap.dat.ovl.forDUP = false;
    skipDUP++;
  }

  if ((roverlap.dat.ovl.forDUP == true) && (skipReadDUP[roverlap.b_iid] == true)) {
    roverlap.dat.ovl.forDUP = false;
    skipDUP++;
  }

  //  Remove the bad-for-DUP overlaps.

#if 0
  //  Nah, do this in dedupe, since parameters can change.
  if ((isDiff == true) && (foverlap.dat.ovl.forDUP == true)) {
    foverlap.dat.ovl.forDUP = false;
    skipDUPdiff++;
  }

  if ((isDiff == true) && (roverlap.dat.ovl.forDUP == true)) {
    roverlap.dat.ovl.forDUP = false;
    skipDUPdiff++;
  }
#endif

  //  Can't have duplicates between libraries.

  if (((foverlap.dat.ovl.forDUP == true) ||
       (roverlap.dat.ovl.forDUP == true)) &&
      (seq->sqStore_getRead(foverlap.a_iid)->sqRead_libraryID() != seq->sqStore_getRead(foverlap.b_iid)->sqRead_libraryID())) {

    if ((foverlap.dat.ovl.forDUP == true)) {
      foverlap.dat.ovl.forDUP = false;
      skipDUPlib++;
    }

    if ((roverlap.dat.ovl.forDUP == true)) {
      roverlap.dat.ovl.forDUP = false;
      skipDUPlib++;
    }
  }

  //  All done with the filtering, record some counts.

  if (foverlap.dat.ovl.forUTG == true)  saveUTG++;
  if (foverlap.dat.ovl.forOBT == true)  saveOBT++;
  if (foverlap.dat.ovl.forDUP == true)  saveDUP++;

  if (roverlap.dat.ovl.forUTG == true)  saveUTG++;
  if (roverlap.dat.ovl.forOBT == true)  saveOBT++;
  if (roverlap.dat.ovl.forDUP == true)  saveDUP++;
}



void
ovStoreFilter::resetCounters(void) {
  saveUTG         = 0;
  saveOBT         = 0;
  saveDUP         = 0;

  skipERATE       = 0;

  skipFLIPPED     = 0;

  skipOBT         = 0;
  skipOBTbad      = 0;
  skipOBTshort    = 0;

  skipDUP         = 0;
  skipDUPdiff     = 0;
  skipDUPlib      = 0;
}
