
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

static const char *rcsid = "$Id: overlapStore_build.c,v 1.37 2011-12-29 09:26:03 brianwalenz Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>
#include <fcntl.h>
#include <assert.h>
#include <time.h>
#include <sys/stat.h>

#include "AS_global.h"
#include "AS_UTL_fileIO.h"
#include "AS_UTL_qsort_mt.h"
#include "AS_OBT_acceptableOverlap.h"
#include "AS_OVS_overlap.h"
#include "AS_OVS_overlapFile.h"
#include "AS_OVS_overlapStore.h"

#include "overlapStore.h"

#include <math.h>

#include <vector>

using namespace std;

static
uint64
computeIIDperBucket(uint32 fileLimit, uint64 memoryLimit, uint32 maxIID, uint32 fileListLen, char **fileList) {
  uint64  numOverlaps = 0;

  if (fileLimit > 0) {
    uint64  iidPerBucket = (uint64)ceil((double)maxIID / (double)fileLimit);

    fprintf(stderr, "Explicit bucket count supplied, memory sizing disabled.  I'll put "F_U64" IIDs into each of "F_U32" buckets.\n",
            iidPerBucket, fileLimit);
    return(iidPerBucket);
  }

  if (fileList[0][0] == '-') {
    fileLimit = sysconf(_SC_OPEN_MAX) - 16;
    uint64  iidPerBucket = (uint64)ceil((double)maxIID / (double)fileLimit);

    fprintf(stderr, "Reading overlaps from stdin, memory sizing disabled.  I'll put "F_U64" IIDs into each of "F_U32" buckets.\n",
            iidPerBucket, fileLimit);
    return(iidPerBucket);
  }

  fprintf(stderr, "Scanning overlap files to count the number of overlaps.\n");

  for (uint32 i=0; i<fileListLen; i++) {
    uint64  no = AS_UTL_sizeOfFile(fileList[i]);
    if (no == 0)
      fprintf(stderr, "WARNING:  No overlaps found (or file not found) in '%s'.\n", fileList[i]);

    numOverlaps += 2 * no / sizeof(OVSoverlap);
  }

  fprintf(stderr, "Found %.3f million overlaps.\n", numOverlaps / 1000000.0);
  assert(numOverlaps > 0);

  //  Why the +1 below?  Consider the case when the number of overlaps is less than the number of
  //  fragments.  This value is used to figure out how many IIDs we can fit into a single bucket,
  //  and making it too large means we'll get maybe one more bucket and the buckets will be smaller.
  //  Yeah, we probably could have just used ceil.
  //
  double  overlapsPerBucket   = (double)memoryLimit / (double)sizeof(OVSoverlap);
  double  overlapsPerIID      = (double)numOverlaps / (double)maxIID;

  uint64  iidPerBucket        = (uint64)(overlapsPerBucket / overlapsPerIID) + 1;

  fileLimit = maxIID / iidPerBucket + 1;

  fprintf(stderr, "Memory limit "F_U64"MB supplied.  I'll put "F_U64" IIDs (%.2f million overlaps) into each of "F_U32" buckets.\n",
          memoryLimit / (uint64)1048576,
          iidPerBucket,
          overlapsPerBucket / 1000000.0,
          fileLimit);

  return(iidPerBucket);
}




static
void
markLoad(OverlapStore *storeFile, uint32 maxIID, char *&skipFragment, uint32 *&iidToLib) {
  gkStream    *gks = new gkStream(storeFile->gkp, 0, 0, GKFRAGMENT_INF);
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
markOBT(OverlapStore *storeFile, uint32 maxIID, char *skipFragment, uint32 *iidToLib) {
  uint64  numMarked = 0;

  if (skipFragment == NULL)
    return;

  fprintf(stderr, "Marking fragments to skip overlap based trimming.\n");

  for (uint64 iid=0; iid<maxIID; iid++) {
    gkLibrary *L = storeFile->gkp->gkStore_getLibrary(iidToLib[iid]);

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
markDUP(OverlapStore *storeFile, uint32 maxIID, char *skipFragment, uint32 *iidToLib) {
  uint64  numMarked = 0;

  if (skipFragment == NULL)
    return;

  fprintf(stderr, "Marking fragments to skip deduplication.\n");

  for (uint64 iid=0; iid<maxIID; iid++) {
    gkLibrary *L = storeFile->gkp->gkStore_getLibrary(iidToLib[iid]);

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
OVSoverlap_sort(const void *a, const void *b) {
  OVSoverlap const *A = (OVSoverlap const *)a;
  OVSoverlap const *B = (OVSoverlap const *)b;
  if (A->a_iid   < B->a_iid)    return(-1);
  if (A->a_iid   > B->a_iid)    return(1);
  if (A->b_iid   < B->b_iid)    return(-1);
  if (A->b_iid   > B->b_iid)    return(1);
  if (A->dat.dat < B->dat.dat)  return(-1);
  if (A->dat.dat > B->dat.dat)  return(1);
  return(0);
}


static
void
writeToDumpFile(OVSoverlap          *overlap,
                BinaryOverlapFile  **dumpFile,
                uint32               dumpFileMax,
                uint64              *dumpLength,
                uint32               iidPerBucket,
                char                *storeName) {

  uint32 df = overlap->a_iid / iidPerBucket;

  if (df >= dumpFileMax) {
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
    sprintf(name, "%s/tmp.sort.%03d", storeName, df);
    dumpFile[df]   = AS_OVS_createBinaryOverlapFile(name, FALSE);
    dumpLength[df] = 0;
  }

  AS_OVS_writeOverlap(dumpFile[df], overlap);
  dumpLength[df]++;
}





void
buildStore(char *storeName, 
           char *gkpName, 
           uint64 memoryLimit, 
           uint32 fileLimit,
           uint32 nThreads, 
           uint32 doFilterOBT, 
           uint32 fileListLen, 
           char **fileList, 
           Ovl_Skip_Type_t ovlSkipOpt) {

  if (gkpName == NULL) {
    fprintf(stderr, "overlapStore: The '-g gkpName' parameter is required.\n");
    exit(1);
  }

  //  We create the store early, allowing it to fail if it already
  //  exists, or just cannot be created.
  //
  OverlapStore    *storeFile = AS_OVS_createOverlapStore(storeName, TRUE);

  storeFile->gkp = new gkStore(gkpName, FALSE, FALSE);

  uint64  maxIID              = storeFile->gkp->gkStore_getNumFragments() + 1;


  //  Decide on some sizes.  We need to decide on how many IID's to
  //  put in each bucket.  Except for running out of file descriptors
  //  (an OS limit), there isn't much of a penalty for having lots of
  //  buckets -- our BinaryOverlapFile buffers writes, and, in fact,
  //  we could open/close the file each time if things get too bad.
  //
  //  The 2x multiplier isn't really true -- MER overlaps don't need
  //  to be flipped, and so mer overlaps count the true number.
  //  Maybe.
  //
  uint64                   iidPerBucket = computeIIDperBucket(fileLimit, memoryLimit, maxIID, fileListLen, fileList);

  uint32                   dumpFileMax  = sysconf(_SC_OPEN_MAX) - 16;
  BinaryOverlapFile      **dumpFile     = (BinaryOverlapFile **)safe_calloc(sizeof(BinaryOverlapFile *), dumpFileMax);
  uint64                  *dumpLength   = (uint64 *)safe_calloc(sizeof(uint64), dumpFileMax);

  if (maxIID / iidPerBucket + 1 >= dumpFileMax) {
    fprintf(stderr, "ERROR:\n");
    fprintf(stderr, "ERROR:  Operating system limit of %d open files.  The current -M and -F settings\n", dumpFileMax);
    fprintf(stderr, "ERROR:  will need to create "F_U64" files to construct the store.\n", maxIID / iidPerBucket);
    fprintf(stderr, "ERROR:  Increase runCA option ovlStoreMemory.\n");
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

  uint64   skipOBT1LQ      = 0;
  uint64   skipOBT2HQ      = 0;
  uint64   skipOBT2LIB     = 0;
  uint64   skipOBT2NODEDUP = 0;

  if (doFilterOBT != 0)
    markLoad(storeFile, maxIID, skipFragment, iidToLib);

  if (doFilterOBT == 1)
    markOBT(storeFile, maxIID, skipFragment, iidToLib);

  if (doFilterOBT == 2)
    markDUP(storeFile, maxIID, skipFragment, iidToLib);


  

  for (uint32 i=0; i<fileListLen; i++) {
    BinaryOverlapFile  *inputFile;
    OVSoverlap          fovrlap;
    OVSoverlap          rovrlap;
    int                 df;

    fprintf(stderr, "bucketizing %s\n", fileList[i]);

    inputFile = AS_OVS_openBinaryOverlapFile(fileList[i], FALSE);

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
         int firstIgnore = (storeFile->gkp->gkStore_getFRGtoPLC(fovrlap.a_iid) != 0 ? TRUE : FALSE);
         int secondIgnore = (storeFile->gkp->gkStore_getFRGtoPLC(fovrlap.b_iid) != 0 ? TRUE : FALSE);
         
         // option means don't ignore them at all
         if (ovlSkipOpt == NONE) {
         }
         // option means don't overlap them at all
         else if (ovlSkipOpt == ALL && ((firstIgnore == TRUE || secondIgnore == TRUE))) {
            continue;
         }
         // option means let them overlap other reads but not each other
         else if (ovlSkipOpt == INTERNAL && ((firstIgnore == TRUE && secondIgnore == TRUE))) {
            continue;
         }
      }

      writeToDumpFile(&fovrlap, dumpFile, dumpFileMax, dumpLength, iidPerBucket, storeName);

      //  flip the overlap -- copy all the dat, then fix whatever
      //  needs to change for the flip.

      switch (fovrlap.dat.ovl.type) {
        case AS_OVS_TYPE_OVL:
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

          writeToDumpFile(&rovrlap, dumpFile, dumpFileMax, dumpLength, iidPerBucket, storeName);
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

          writeToDumpFile(&rovrlap, dumpFile, dumpFileMax, dumpLength, iidPerBucket, storeName);
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
  }

  for (uint32 i=0; i<dumpFileMax; i++)
    AS_OVS_closeBinaryOverlapFile(dumpFile[i]);

  fprintf(stderr, "bucketizing DONE!\n");

  fprintf(stderr, "overlaps skipped:\n");
  fprintf(stderr, "%16"F_U64P" OBT - low quality\n", skipOBT1LQ);
  fprintf(stderr, "%16"F_U64P" DUP - non-duplicate overlap\n", skipOBT2HQ);
  fprintf(stderr, "%16"F_U64P" DUP - different library\n", skipOBT2LIB);
  fprintf(stderr, "%16"F_U64P" DUP - dedup not requested\n", skipOBT2NODEDUP);

  delete [] skipFragment;  skipFragment = NULL;
  delete [] iidToLib;      iidToLib     = NULL;

  //
  //  Read each bucket, sort it, and dump it to the store
  //

  uint64 dumpLengthMax = 0;
  for (uint32 i=0; i<dumpFileMax; i++)
    if (dumpLengthMax < dumpLength[i])
      dumpLengthMax = dumpLength[i];

  OVSoverlap         *overlapsort = NULL;
  overlapsort = (OVSoverlap *)safe_malloc(sizeof(OVSoverlap) * dumpLengthMax);

  time_t  beginTime = time(NULL);

  for (uint32 i=0; i<dumpFileMax; i++) {
    char                name[FILENAME_MAX];
    BinaryOverlapFile  *bof = NULL;

    if (dumpLength[i] == 0)
      continue;

    //  We're vastly more efficient if we skip the AS_OVS interface
    //  and just suck in the whole file directly....BUT....we can't do
    //  that because the AS_OVS interface is rearranging the data to
    //  make sure the store is cross-platform compatible.

    sprintf(name, "%s/tmp.sort.%03d", storeName, i);
    fprintf(stderr, "reading %s (%ld)\n", name, time(NULL) - beginTime);

    bof = AS_OVS_openBinaryOverlapFile(name, FALSE);

    uint64 numOvl = 0;
    while (AS_OVS_readOverlap(bof, overlapsort + numOvl))
      numOvl++;

    AS_OVS_closeBinaryOverlapFile(bof);

    assert(numOvl == dumpLength[i]);
    assert(numOvl <= dumpLengthMax);

    //  There's no real advantage to saving this file until after we
    //  write it out.  If we crash anywhere during the build, we are
    //  forced to restart from scratch.  I'll argue that removing it
    //  early helps us to not crash from running out of disk space.
    //
    unlink(name);

    fprintf(stderr, "sorting %s (%ld)\n", name, time(NULL) - beginTime);
    qsort_mt(overlapsort, dumpLength[i], sizeof(OVSoverlap), OVSoverlap_sort, nThreads, 16 * 1024 * 1024);

    fprintf(stderr, "writing %s (%ld)\n", name, time(NULL) - beginTime);
    for (uint64 x=0; x<dumpLength[i]; x++)
      AS_OVS_writeOverlapToStore(storeFile, overlapsort + x);
  }

  AS_OVS_closeOverlapStore(storeFile);

  safe_free(overlapsort);

  //  And we have a store.
  //
  exit(0);
}
