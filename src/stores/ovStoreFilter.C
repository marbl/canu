
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

#include "ovStore.H"



ovStoreFilter::ovStoreFilter(sqStore *seq_, double maxErate_) {
  seq             = seq_;
  maxID           = seq->sqStore_lastReadID();
  maxEvalue       = AS_OVS_encodeEvalue(maxErate_);

  resetCounters();

  skipReadOBT     = new char [maxID + 1];

  for (uint64 iid=0; iid<=maxID; iid++)
    skipReadOBT[iid] = false;

  uint32  numSkipOBT = 0;

  //  Query the library to decide if this overlap should be filtered out
  //  because it won't be used by the algorithms processing this data.
  //
  //  Example: if these are overlaps for trimming, but we're not trimming
  //  this read, drop all it's overlaps.
  //
  //  But since this isn't implemented yet, emit all overlaps all the time.
  //
#if 0
  for (uint64 iid=0; iid<=maxID; iid++) {
    sqLibrary *L = seq->sqStore_getLibraryForRead(iid);

    if ((L->sqLibrary_finalTrim()                == SQ_FINALTRIM_NONE) &&
        (L->sqLibrary_removeSpurReads()          == false) &&
        (L->sqLibrary_removeChimericReads()      == false) &&
        (L->sqLibrary_checkForSubReads()         == false)) {
      numSkipOBT++;
      skipReadOBT[iid] = true;
    }
  }
#endif

  if (numSkipOBT > 0)
    fprintf(stderr, "-  Marked " F_U32 " reads to skip trimming.\n", numSkipOBT);
}



ovStoreFilter::~ovStoreFilter() {
  delete [] skipReadOBT;
}



void
ovStoreFilter::filterOverlap(ovOverlap       &foverlap,
                             ovOverlap       &roverlap) {

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

  //  Do NOT initialize the 'for' flags; overlapper is already flagging
  //  overlaps as bad for unitigging if they aren't dovetail.

  //foverlap.dat.ovl.forUTG = true;
  //foverlap.dat.ovl.forOBT = true;
  //foverlap.dat.ovl.forDUP = true;

  //  Make the reverse overlap.

  roverlap.swapIDs(foverlap);

  //  Ignore high error overlaps.

  if ((foverlap.evalue() > maxEvalue)) {
    foverlap.dat.ovl.forUTG = false;
    foverlap.dat.ovl.forOBT = false;
    foverlap.dat.ovl.forDUP = false;
    skipERATE++;

    roverlap.dat.ovl.forUTG = false;
    roverlap.dat.ovl.forOBT = false;
    roverlap.dat.ovl.forDUP = false;
    skipERATE++;
  }

  //  Ignore opposite oriented overlaps
#ifdef IGNORE_FLIPPED_OVERLAPS
  if ((foverlap.flipped() == true)) {
    foverlap.dat.ovl.forUTG = false;
    foverlap.dat.ovl.forOBT = false;
    foverlap.dat.ovl.forDUP = false;
    skipFLIPPED++;

    roverlap.dat.ovl.forUTG = false;
    roverlap.dat.ovl.forOBT = false;
    roverlap.dat.ovl.forDUP = false;
    skipFLIPPED++;
  }
#endif

  //  Don't OBT if not requested.  We allow non-OBT-able reads to be used as
  //  evidence, but if the read isn't eligible for trimming, we can discard
  //  all of its overlaps.

  if ((foverlap.dat.ovl.forOBT == false) && (skipReadOBT[foverlap.a_iid] == true)) {
    foverlap.dat.ovl.forOBT = false;
    skipOBT++;
  }

  if ((roverlap.dat.ovl.forOBT == false) && (skipReadOBT[roverlap.a_iid] == true)) {
    roverlap.dat.ovl.forOBT = false;
    skipOBT++;
  }

  //  All done with the filtering, record some counts.

  if (foverlap.dat.ovl.forUTG == true)  saveUTG++;
  if (foverlap.dat.ovl.forOBT == true)  saveOBT++;
  //if (foverlap.dat.ovl.forDUP == true)  not-used;

  if (roverlap.dat.ovl.forUTG == true)  saveUTG++;
  if (roverlap.dat.ovl.forOBT == true)  saveOBT++;
  //if (roverlap.dat.ovl.forDUP == true)  not-used;
}



void
ovStoreFilter::resetCounters(void) {
  saveUTG         = 0;
  saveOBT         = 0;

  skipERATE       = 0;
  skipFLIPPED     = 0;
  skipOBT         = 0;
}
