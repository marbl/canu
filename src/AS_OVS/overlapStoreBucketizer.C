
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 2007, J. Craig Venter Institute. All rights reserved.
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

const char *mainid = "$Id: overlapStoreBucketizer.C,v 1.1 2012-04-02 10:58:04 brianwalenz Exp $";

#include "AS_PER_gkpStore.h"

#include "overlapStore.h"

#include "AS_OVS_overlap.h"
#include "AS_OVS_overlapFile.h"
#include "AS_OVS_overlapStore.h"
#include "AS_OBT_acceptableOverlap.h"

#include <ctype.h>
#include <unistd.h>  //  sysconf()

#include <vector>
#include <algorithm>

using namespace std;

#define WITH_GZIP 1

static
void
writeToFile(OVSoverlap          *overlap,
            BinaryOverlapFile  **dumpFile,
            uint32               dumpFileMax,
            uint64              *dumpLength,
            uint32               iidPerBucket,
            char                *ovlName,
            uint32               jobIndex) {

  uint32 df = overlap->a_iid / iidPerBucket + 1;

  if (df > dumpFileMax) {
    char   olapstring[256];
    
    fprintf(stderr, "\n");
    fprintf(stderr, "Too many bucket files when adding overlap:\n");
    fprintf(stderr, "  %s\n", AS_OVS_toString(olapstring, *overlap));
    fprintf(stderr, "\n");
    fprintf(stderr, "bucket       = "F_U32"\n", df);
    fprintf(stderr, "iidPerBucket = "F_U32"\n", iidPerBucket);
    fprintf(stderr, "dumpFileMax  = "F_U32"\n", dumpFileMax);
    fprintf(stderr, "\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "This might be a corrupt input file, or maybe you simply need to supply more\n");
    fprintf(stderr, "memory with the runCA option ovlStoreMemory.\n");
    fprintf(stderr, "\n");
    exit(1);
  }

  if (dumpFile[df] == NULL) {
    char name[FILENAME_MAX];

    sprintf(name, "%s/unsorted%04d/tmp.sort.%03d%S", ovlName, jobIndex, df, (WITH_GZIP) ? ".gz" : "");
    dumpFile[df]   = AS_OVS_createBinaryOverlapFile(name, FALSE);
    dumpLength[df] = 0;
  }

  AS_OVS_writeOverlap(dumpFile[df], overlap);
  dumpLength[df]++;
}



static
void
markLoad(gkStore *gkp, uint32 maxIID, char *&skipFragment, uint32 *&iidToLib) {
  gkStream    *gks = new gkStream(gkp, 0, 0, GKFRAGMENT_INF);
  gkFragment   fr;

  fprintf(stderr, "Reading gatekeeper to build a map from fragment ID to library ID.\n");

  skipFragment = new char [maxIID];
  iidToLib     = new uint32 [maxIID];

  memset(skipFragment, 0, sizeof(char)   * maxIID);
  memset(iidToLib,     0, sizeof(uint32) * maxIID);

  while (gks->next(&fr))
    iidToLib[fr.gkFragment_getReadIID()] = fr.gkFragment_getLibraryIID();

  delete gks;
}



static
void
markOBT(gkStore *gkp, uint32 maxIID, char *skipFragment, uint32 *iidToLib) {
  uint64  numMarked = 0;

  if (skipFragment == NULL)
    return;

  fprintf(stderr, "Marking fragments to skip overlap based trimming.\n");

  for (uint64 iid=0; iid<maxIID; iid++) {
    gkLibrary *L = gkp->gkStore_getLibrary(iidToLib[iid]);

    if (L == NULL)
      continue;

    if ((L->doRemoveDuplicateReads     == false) &&
        (L->doTrim_finalLargestCovered == false) &&
        (L->doTrim_finalEvidenceBased  == false) &&
        (L->doRemoveSpurReads          == false) &&
        (L->doRemoveChimericReads      == false)) {
      numMarked++;
      skipFragment[iid] = true;
    }
  }

  fprintf(stderr, "Marked "F_U64" fragments.\n", numMarked);
}


static
void
markDUP(gkStore *gkp, uint32 maxIID, char *skipFragment, uint32 *iidToLib) {
  uint64  numMarked = 0;

  if (skipFragment == NULL)
    return;

  fprintf(stderr, "Marking fragments to skip deduplication.\n");

  for (uint64 iid=0; iid<maxIID; iid++) {
    gkLibrary *L = gkp->gkStore_getLibrary(iidToLib[iid]);

    if (L == NULL)
      continue;

    if (L->doRemoveDuplicateReads == false) {
      numMarked++;
      skipFragment[iid] = true;
    }
  }

  fprintf(stderr, "Marked "F_U64" fragments.\n", numMarked);
}



int
main(int argc, char **argv) {
  char           *ovlName      = NULL;
  char           *gkpName      = NULL;
  uint32          fileLimit    = 512;

  Ovl_Skip_Type_t ovlSkipOpt   = PLC_NONE;
  uint32          doFilterOBT  = 0;

  uint32          jobIndex     = 0;

  double          maxErrorRate = 1.0;
  uint64          maxError     = AS_OVS_encodeQuality(maxErrorRate);

  char           *ovlInput     = NULL;

  argc = AS_configure(argc, argv);

  int err=0;
  int arg=1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-o") == 0) {
      ovlName = argv[++arg];

    } else if (strcmp(argv[arg], "-g") == 0) {
      gkpName = argv[++arg];

    } else if (strcmp(argv[arg], "-F") == 0) {
      fileLimit = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-plc") == 0) {
      //  Former -i option
      //  PLC_NONE, PLC_ALL, PLC_INTERNAL
      ovlSkipOpt = PLC_NONE;

    } else if (strcmp(argv[arg], "-obt") == 0) {
      doFilterOBT = 1;

    } else if (strcmp(argv[arg], "-dup") == 0) {
      doFilterOBT = 2;

    } else if (strcmp(argv[arg], "-job") == 0) {
      jobIndex = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-i") == 0) {
      ovlInput = argv[++arg];

    } else if (strcmp(argv[arg], "-e") == 0) {
      maxError = atof(argv[++arg]);
      maxError = AS_OVS_encodeQuality(maxErrorRate);

    } else {
      fprintf(stderr, "ERROR: unknown option '%s'\n", argv[arg]);
    }

    arg++;
  }
  if (ovlName == NULL)
    err++;
  if (gkpName == NULL)
    err++;
  if (ovlInput == NULL)
    err++;
  if (jobIndex == 0)
    err++;
  if (err) {
    exit(1);
  }


  {
    char name[FILENAME_MAX];

    if (AS_UTL_fileExists(ovlName, TRUE, FALSE) == false)
      AS_UTL_mkdir(ovlName);

    sprintf(name, "%s/unsorted%04d", ovlName, jobIndex);
    if (AS_UTL_fileExists(name, TRUE, FALSE) == false)
      AS_UTL_mkdir(name);
  }


  gkStore *gkp         = new gkStore(gkpName, FALSE, FALSE);

  uint64  maxIID       = gkp->gkStore_getNumFragments() + 1;
  uint64  iidPerBucket = (uint64)ceil((double)maxIID / (double)fileLimit);

  uint32                   dumpFileMax  = sysconf(_SC_OPEN_MAX) + 1;
  BinaryOverlapFile      **dumpFile     = new BinaryOverlapFile * [dumpFileMax];
  uint64                  *dumpLength   = new uint64              [dumpFileMax];

  memset(dumpFile,   0, sizeof(BinaryOverlapFile *) * dumpFileMax);
  memset(dumpLength, 0, sizeof(uint64)              * dumpFileMax);

  if (maxIID / iidPerBucket + 1 > dumpFileMax - 16) {
    fprintf(stderr, "ERROR:\n");
    fprintf(stderr, "ERROR:  Operating system limit of %d open files.  The current -F setting\n", dumpFileMax);
    fprintf(stderr, "ERROR:  will need to create "F_U64" files to construct the store.\n", maxIID / iidPerBucket + 1);
    exit(1);
  }

  //  Read the gkStore to determine which fragments we care about.
  //
  //  If doFilterOBT == 0, we care about all overlaps (we're not processing for OBT).
  //
  //  If doFilterOBT == 1, then we care about overlaps where either fragment is in a doNotOBT == 0
  //  library.
  //
  //  If doFilterOBT == 2, then we care about overlaps where both fragments are in the same
  //  library, and that library is marked doRemoveDuplicateReads == 1

  char    *skipFragment = NULL;
  uint32  *iidToLib     = NULL;

  uint64   skipERATE       = 0;
  uint64   skipOBT1LQ      = 0;
  uint64   skipOBT2HQ      = 0;
  uint64   skipOBT2LIB     = 0;
  uint64   skipOBT2NODEDUP = 0;

  if (doFilterOBT != 0)
    markLoad(gkp, maxIID, skipFragment, iidToLib);

  if (doFilterOBT == 1)
    markOBT(gkp, maxIID, skipFragment, iidToLib);

  if (doFilterOBT == 2)
    markDUP(gkp, maxIID, skipFragment, iidToLib);
  
  BinaryOverlapFile  *inputFile;
  OVSoverlap          fovrlap;
  OVSoverlap          rovrlap;

  int                 df;

  fprintf(stderr, "Bucketizing %s\n", ovlInput);

  inputFile = AS_OVS_openBinaryOverlapFile(ovlInput, FALSE);

  while (AS_OVS_readOverlap(inputFile, &fovrlap)) {

    //  Quick sanity check on IIDs.
    if ((fovrlap.a_iid == 0) ||
        (fovrlap.b_iid == 0) ||
        (fovrlap.a_iid >= maxIID) ||
        (fovrlap.b_iid >= maxIID)) {
      char ovlstr[256];

      fprintf(stderr, "Overlap has IDs out of range (maxIID "F_U64"), possibly corrupt input data.\n", maxIID);
      fprintf(stderr, "  %s\n", AS_OVS_toString(ovlstr, fovrlap));
      exit(1);
    }

    //  Ignore high error overlaps
    if (fovrlap.dat.ovl.orig_erate > maxError ) {
      skipERATE++;
      continue;
    }

    //  If filtering for OBT, skip the crap.
    if ((doFilterOBT == 1) && (AS_OBT_acceptableOverlap(fovrlap) == 0)) {
      skipOBT1LQ++;
      continue;
    }

    //  If filtering for OBT, skip overlaps that we're never going to use.
    //  (for now, we allow everything through -- these are used for just about everything)

    //  If filtering for OBTs dedup, skip the good
    if ((doFilterOBT == 2) && (AS_OBT_acceptableOverlap(fovrlap) == 1)) {
      skipOBT2HQ++;
      continue;
    }

    //  If filtering for OBTs dedup, skip things we don't dedup, and overlaps between libraries.
    if ((doFilterOBT == 2) && (iidToLib[fovrlap.a_iid] != iidToLib[fovrlap.b_iid])) {
      skipOBT2LIB++;
      continue;
    }

    if ((doFilterOBT == 2) && (skipFragment[fovrlap.a_iid])) {
      skipOBT2NODEDUP++;
      continue;
    }

    if (doFilterOBT == 0) {
      int firstIgnore  = (gkp->gkStore_getFRGtoPLC(fovrlap.a_iid) != 0 ? TRUE : FALSE);
      int secondIgnore = (gkp->gkStore_getFRGtoPLC(fovrlap.b_iid) != 0 ? TRUE : FALSE);
         
      // option means don't ignore them at all
      if (ovlSkipOpt == PLC_NONE) {
      }
      // option means don't overlap them at all
      else if (ovlSkipOpt == PLC_ALL && ((firstIgnore == TRUE || secondIgnore == TRUE))) {
        continue;
      }
      // option means let them overlap other reads but not each other
      else if (ovlSkipOpt == PLC_INTERNAL && ((firstIgnore == TRUE && secondIgnore == TRUE))) {
        continue;
      }
    }


    writeToFile(&fovrlap, dumpFile, dumpFileMax, dumpLength, iidPerBucket, ovlName, jobIndex);

    //  flip the overlap -- copy all the dat, then fix whatever
    //  needs to change for the flip.

    switch (fovrlap.dat.ovl.type) {
	
      case AS_OVS_TYPE_OVL:
        // This inverts the overlap.
        rovrlap.a_iid = fovrlap.b_iid;
        rovrlap.b_iid = fovrlap.a_iid;
        rovrlap.dat   = fovrlap.dat;
        if (fovrlap.dat.ovl.flipped) {
          rovrlap.dat.ovl.a_hang = fovrlap.dat.ovl.b_hang;
          rovrlap.dat.ovl.b_hang = fovrlap.dat.ovl.a_hang;
        } else {
          rovrlap.dat.ovl.a_hang = -fovrlap.dat.ovl.a_hang;
          rovrlap.dat.ovl.b_hang = -fovrlap.dat.ovl.b_hang;
        }

        writeToFile(&rovrlap, dumpFile, dumpFileMax, dumpLength, iidPerBucket, ovlName, jobIndex);
        break;
      case AS_OVS_TYPE_OBT:
        rovrlap.a_iid = fovrlap.b_iid;
        rovrlap.b_iid = fovrlap.a_iid;
        rovrlap.dat   = fovrlap.dat;
        if (fovrlap.dat.obt.fwd) {
          rovrlap.dat.obt.a_beg    = fovrlap.dat.obt.b_beg;
          rovrlap.dat.obt.a_end    = (fovrlap.dat.obt.b_end_hi << 9) | fovrlap.dat.obt.b_end_lo;
          rovrlap.dat.obt.b_beg    = fovrlap.dat.obt.a_beg;
          rovrlap.dat.obt.b_end_hi = fovrlap.dat.obt.a_end >> 9;
          rovrlap.dat.obt.b_end_lo = fovrlap.dat.obt.a_end & 0x1ff;
        } else {
          rovrlap.dat.obt.a_beg    = (fovrlap.dat.obt.b_end_hi << 9) | fovrlap.dat.obt.b_end_lo;
          rovrlap.dat.obt.a_end    = fovrlap.dat.obt.b_beg;
          rovrlap.dat.obt.b_beg    = fovrlap.dat.obt.a_end;
          rovrlap.dat.obt.b_end_hi = fovrlap.dat.obt.a_beg >> 9;
          rovrlap.dat.obt.b_end_lo = fovrlap.dat.obt.a_beg & 0x1ff;
        }

        writeToFile(&rovrlap, dumpFile, dumpFileMax, dumpLength, iidPerBucket, ovlName, jobIndex);
        break;
      case AS_OVS_TYPE_MER:
        //  Not needed; MER outputs both overlaps
        break;
      default:
        assert(0);
        break;
    }
  }

  AS_OVS_closeBinaryOverlapFile(inputFile);

  for (uint32 i=0; i<dumpFileMax; i++)
    AS_OVS_closeBinaryOverlapFile(dumpFile[i]);

  {
    char name[FILENAME_MAX];

    sprintf(name, "%s/unsorted%04d/bucketSizes", ovlName, jobIndex);

    FILE *F = fopen(name, "w");
    if (errno)
      fprintf(stderr, "ERROR:  Failed to open %s: %s\n", name, strerror(errno)), exit(1);

    AS_UTL_safeWrite(F, dumpLength, "dumpLength", sizeof(uint64), dumpFileMax);

    fclose(F);
  }

  fprintf(stderr, "overlaps skipped:\n");
  fprintf(stderr, "%16"F_U64P" ERR - low quality, more than %.2f fraction error\n", skipERATE, maxErrorRate);
  fprintf(stderr, "%16"F_U64P" OBT - low quality\n", skipOBT1LQ);
  fprintf(stderr, "%16"F_U64P" DUP - non-duplicate overlap\n", skipOBT2HQ);
  fprintf(stderr, "%16"F_U64P" DUP - different library\n", skipOBT2LIB);
  fprintf(stderr, "%16"F_U64P" DUP - dedup not requested\n", skipOBT2NODEDUP);

  delete [] skipFragment;  skipFragment = NULL;
  delete [] iidToLib;      iidToLib     = NULL;
}
