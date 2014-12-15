
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

const char *mainid = "$Id$";

#include "AS_global.H"

#include "gkStore.H"
#include "ovStore.H"

//#include "AS_OBT_acceptableOverlap.H"
#warning NOT INCLUDING ACCEPTABLE OVERLAP FROM OBT
bool   AS_OBT_acceptableOverlap(ovsOverlap &ol) {
  return(true);
};

//#include <ctype.h>
//#include <unistd.h>  //  sysconf()

#include <vector>
#include <algorithm>

using namespace std;

static
void
writeToFile(ovsOverlap    *overlap,
            ovFile       **sliceFile,
            uint32         sliceFileMax,
            uint64        *sliceSize,
            uint32         iidPerBucket,
            char          *ovlName,
            uint32         jobIndex,
            bool           useGzip) {

  uint32 df = overlap->a_iid / iidPerBucket + 1;

  if (df > sliceFileMax) {
    char   olapstring[256];
    
    fprintf(stderr, "\n");
    fprintf(stderr, "Too many bucket files when adding overlap:\n");
    fprintf(stderr, "  %s\n", overlap->toString(olapstring));
    fprintf(stderr, "\n");
    fprintf(stderr, "bucket        = "F_U32"\n", df);
    fprintf(stderr, "iidPerBucket  = "F_U32"\n", iidPerBucket);
    fprintf(stderr, "sliceFileMax  = "F_U32"\n", sliceFileMax);
    fprintf(stderr, "\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "This might be a corrupt input file, or maybe you simply need to supply more\n");
    fprintf(stderr, "memory with the runCA option ovlStoreMemory.\n");
    fprintf(stderr, "\n");
    exit(1);
  }

  if (sliceFile[df] == NULL) {
    char name[FILENAME_MAX];

    sprintf(name, "%s/create%04d/slice%03d%s", ovlName, jobIndex, df, (useGzip) ? ".gz" : "");
    sliceFile[df] = new ovFile(name, ovFileFullWrite);
    sliceSize[df] = 0;
  }

  sliceFile[df]->writeOverlap(overlap);
  sliceSize[df]++;
}



//  These are duplicated between ovStoreBucketizer and ovStoreBuild

static
void
markOBT(gkStore *gkp, uint32 maxIID, char *skipFragment) {
  uint64  numMarked = 0;

  if (skipFragment == NULL)
    return;

  fprintf(stderr, "Marking fragments to skip overlap based trimming.\n");

  for (uint64 iid=0; iid<maxIID; iid++) {
    uint32     Lid = gkp->gkStore_getRead(iid)->gkRead_libraryID();
    gkLibrary *L   = gkp->gkStore_getLibrary(Lid);

    if (L == NULL)
      continue;

    if ((L->gkLibrary_removeDuplicateReads()     == false) &&
        (L->gkLibrary_finalTrim()                != FINALTRIM_LARGEST_COVERED) &&
        (L->gkLibrary_finalTrim()                != FINALTRIM_EVIDENCE_BASED) &&
        (L->gkLibrary_removeSpurReads()          == false) &&
        (L->gkLibrary_removeChimericReads()      == false)) {
      numMarked++;
      skipFragment[iid] = true;
    }
  }

  fprintf(stderr, "Marked "F_U64" fragments.\n", numMarked);
}


static
void
markDUP(gkStore *gkp, uint32 maxIID, char *skipFragment) {
  uint64  numMarked = 0;

  if (skipFragment == NULL)
    return;

  fprintf(stderr, "Marking fragments to skip deduplication.\n");

  for (uint64 iid=0; iid<maxIID; iid++) {
    uint32     Lid = gkp->gkStore_getRead(iid)->gkRead_libraryID();
    gkLibrary *L   = gkp->gkStore_getLibrary(Lid);

    if (L == NULL)
      continue;

    if (L->gkLibrary_removeDuplicateReads() == false) {
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

  uint32          doFilterOBT  = 0;

  uint32          jobIndex     = 0;

  double          maxErrorRate = 1.0;
  uint64          maxError     = AS_OVS_encodeQuality(maxErrorRate);

  char           *ovlInput     = NULL;

  bool            useGzip      = true;

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

    } else if (strcmp(argv[arg], "-obt") == 0) {
      doFilterOBT = 1;

    } else if (strcmp(argv[arg], "-dup") == 0) {
      doFilterOBT = 2;

    } else if (strcmp(argv[arg], "-job") == 0) {
      jobIndex = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-i") == 0) {
      ovlInput = argv[++arg];

    } else if (strcmp(argv[arg], "-e") == 0) {
      maxErrorRate = atof(argv[++arg]);
      maxError     = AS_OVS_encodeQuality(maxErrorRate);

    } else if (strcmp(argv[arg], "-raw") == 0) {
      useGzip = false;

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
  if (fileLimit > sysconf(_SC_OPEN_MAX) - 16)
    err++;

  if (err) {
    fprintf(stderr, "usage: %s -c asm.ovlStore -g asm.gkpStore -i file.ovb.gz -job j [opts]\n", argv[0]);
    fprintf(stderr, "  -c asm.ovlStore       path to store to create\n");
    fprintf(stderr, "  -g asm.gkpStore       path to gkpStore for this assembly\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -i file.ovb.gz        input overlaps\n");
    fprintf(stderr, "  -job j                index of this overlap input file\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -F f                  use up to 'f' files for store creation\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -obt                  filter overlaps for OBT\n");
    fprintf(stderr, "  -dup                  filter overlaps for OBT/dedupe\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -e e                  filter overlaps above e fraction error\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -raw                  write uncompressed buckets\n");

    if (ovlName == NULL)
      fprintf(stderr, "ERROR: No overlap store (-o) supplied.\n");
    if (gkpName == NULL)
      fprintf(stderr, "ERROR: No gatekeeper store (-g) supplied.\n");
    if (ovlInput == NULL)
      fprintf(stderr, "ERROR: No input (-i) supplied.\n");
    if (jobIndex == 0)
      fprintf(stderr, "ERROR: No job index (-job) supplied.\n");
    if (fileLimit > sysconf(_SC_OPEN_MAX) - 16)
      fprintf(stderr, "ERROR: Too many jobs (-F); only "F_SIZE_T" supported on this architecture.\n", sysconf(_SC_OPEN_MAX) - 16);

    exit(1);
  }


  {
    if (AS_UTL_fileExists(ovlName, TRUE, FALSE) == false)
      AS_UTL_mkdir(ovlName);
  }


  {
    char name[FILENAME_MAX];

    sprintf(name, "%s/create%04d", ovlName, jobIndex);

    if (AS_UTL_fileExists(name, TRUE, FALSE) == false)
      AS_UTL_mkdir(name);
    else
      fprintf(stderr, "Overwriting previous result; directory '%s' exists.\n", name), exit(0);
  }


  {
    char name[FILENAME_MAX];

    sprintf(name, "%s/bucket%04d/sliceSizes", ovlName, jobIndex);

    if (AS_UTL_fileExists(name, FALSE, FALSE) == true)
      fprintf(stderr, "Job finished; file '%s' exists.\n", name), exit(0);
  }


  gkStore *gkp         = new gkStore(gkpName);

  uint64  maxIID       = gkp->gkStore_getNumReads() + 1;
  uint64  iidPerBucket = (uint64)ceil((double)maxIID / (double)fileLimit);

  ovFile                 **sliceFile     = new ovFile * [fileLimit + 1];
  uint64                  *sliceSize     = new uint64   [fileLimit + 1];

  memset(sliceFile, 0, sizeof(ovFile *) * (fileLimit + 1));
  memset(sliceSize, 0, sizeof(uint64)   * (fileLimit + 1));

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

  uint64   saveTOTAL       = 0;
  uint64   skipERATE       = 0;
  uint64   skipOBT1LQ      = 0;
  uint64   skipOBT2HQ      = 0;
  uint64   skipOBT2LIB     = 0;
  uint64   skipOBT2NODEDUP = 0;

  if (doFilterOBT != 0) {
    skipFragment = new char [maxIID];
    memset(skipFragment, 0, sizeof(char) * maxIID);
  }

  if (doFilterOBT == 1)
    markOBT(gkp, maxIID, skipFragment);

  if (doFilterOBT == 2)
    markDUP(gkp, maxIID, skipFragment);

  ovFile       *inputFile;
  ovsOverlap    foverlap;
  ovsOverlap    roverlap;

  int           df;

  fprintf(stderr, "maxError fraction: %.3f percent: %.3f encoded: "F_U64"\n",
          maxErrorRate, maxErrorRate * 100, maxError);

  fprintf(stderr, "Bucketizing %s\n", ovlInput);

  inputFile = new ovFile(ovlInput);

  //  Do bigger buffers increase performance?  Do small ones hurt?
  //AS_OVS_setBinaryOverlapFileBufferSize(2 * 1024 * 1024);

  while (inputFile->readOverlap(&foverlap)) {

    //  Quick sanity check on IIDs.
    if ((foverlap.a_iid == 0) ||
        (foverlap.b_iid == 0) ||
        (foverlap.a_iid >= maxIID) ||
        (foverlap.b_iid >= maxIID)) {
      char ovlstr[256];

      fprintf(stderr, "Overlap has IDs out of range (maxIID "F_U64"), possibly corrupt input data.\n", maxIID);
      fprintf(stderr, "  %s\n", foverlap.toString(ovlstr));
      exit(1);
    }

    //  Ignore high error overlaps
    if ((foverlap.dat.ovl.erate > maxError)) {
      skipERATE++;
      continue;
    }

    //  If filtering for OBT, skip the crap.
    if ((doFilterOBT == 1) && (AS_OBT_acceptableOverlap(foverlap) == 0)) {
      skipOBT1LQ++;
      continue;
    }

    //  If filtering for OBT, skip overlaps that we're never going to use.
    //  (for now, we allow everything through -- these are used for just about everything)

    //  If filtering for OBTs dedup, skip the good
    if ((doFilterOBT == 2) && (AS_OBT_acceptableOverlap(foverlap) == 1)) {
      skipOBT2HQ++;
      continue;
    }

    //  If filtering for OBTs dedup, skip things we don't dedup, and overlaps between libraries.
    if ((doFilterOBT == 2) &&
        (gkp->gkStore_getRead(foverlap.a_iid)->gkRead_libraryID() !=
         gkp->gkStore_getRead(foverlap.b_iid)->gkRead_libraryID())) {
      skipOBT2LIB++;
      continue;
    }

    if ((doFilterOBT == 2) && (skipFragment[foverlap.a_iid])) {
      skipOBT2NODEDUP++;
      continue;
    }

    writeToFile(&foverlap, sliceFile, fileLimit, sliceSize, iidPerBucket, ovlName, jobIndex, useGzip);
    saveTOTAL++;

    //  flip the overlap -- copy all the dat, then fix whatever
    //  needs to change for the flip.

    roverlap.swapIDs(foverlap);

    writeToFile(&roverlap, sliceFile, fileLimit, sliceSize, iidPerBucket, ovlName, jobIndex, useGzip);
    saveTOTAL++;
  }

  delete inputFile;

  for (uint32 i=0; i<=fileLimit; i++)
    delete sliceFile[i];

  //  Write slice sizes, rename bucket.

  {
    char name[FILENAME_MAX];
    char finl[FILENAME_MAX];

    sprintf(name, "%s/create%04d/sliceSizes", ovlName, jobIndex);

    FILE *F = fopen(name, "w");
    if (errno)
      fprintf(stderr, "ERROR:  Failed to open %s: %s\n", name, strerror(errno)), exit(1);

    AS_UTL_safeWrite(F, sliceSize, "sliceSize", sizeof(uint64), fileLimit + 1);

    fclose(F);

    sprintf(name, "%s/create%04d", ovlName, jobIndex);
    sprintf(finl, "%s/bucket%04d", ovlName, jobIndex);

    errno = 0;
    rename(name, finl);
    if (errno)
      fprintf(stderr, "ERROR:  Failed to rename '%s' to final name '%s': %s\n",
              name, finl, strerror(errno));
  }

  fprintf(stderr, "overlap fate:\n");
  fprintf(stderr, "%16"F_U64P" SAV - overlaps output\n", saveTOTAL);
  fprintf(stderr, "%16"F_U64P" ERR - low quality, more than %.3f fraction error\n", skipERATE, maxErrorRate);
  fprintf(stderr, "%16"F_U64P" OBT - low quality\n", skipOBT1LQ);
  fprintf(stderr, "%16"F_U64P" DUP - non-duplicate overlap\n", skipOBT2HQ);
  fprintf(stderr, "%16"F_U64P" DUP - different library\n", skipOBT2LIB);
  fprintf(stderr, "%16"F_U64P" DUP - dedup not requested\n", skipOBT2NODEDUP);

  delete [] sliceFile;
  delete [] sliceSize;

  delete [] skipFragment;  skipFragment = NULL;
}
