
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

const char *mainid = "$Id: overlapStoreBuild.C,v 1.3 2012-05-09 01:16:54 brianwalenz Exp $";

#include "AS_global.h"

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
//#include <parallel/algorithm>

using namespace std;


#define AS_OVS_CURRENT_VERSION  2





static
uint64
computeIIDperBucket(uint32 fileLimit, uint64 memoryLimit, uint32 maxIID, vector<char *> &fileList) {
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

  for (uint32 i=0; i<fileList.size(); i++) {
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






static
void
writeToDumpFile(OVSoverlap          *overlap,
                BinaryOverlapFile  **dumpFile,
                uint32               dumpFileMax,
                uint64              *dumpLength,
                uint32               iidPerBucket,
                char                *ovlName) {

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
    sprintf(name, "%s/tmp.sort.%03d", ovlName, df);
    dumpFile[df]   = AS_OVS_createBinaryOverlapFile(name, FALSE);
    dumpLength[df] = 0;
  }

  AS_OVS_writeOverlap(dumpFile[df], overlap);
  dumpLength[df]++;
}

















int
main(int argc, char **argv) {
  char           *ovlName      = NULL;
  char           *gkpName      = NULL;
  uint32          fileLimit    = 512;
  uint64          memoryLimit  = 0;

  Ovl_Skip_Type_t ovlSkipOpt   = PLC_ALL;
  uint32          doFilterOBT  = 0;

  double          maxErrorRate = 1.0;
  uint64          maxError     = AS_OVS_encodeQuality(maxErrorRate);

  vector<char *>  fileList;

  uint32          nThreads = 4;

  argc = AS_configure(argc, argv);

  int err=0;
  int arg=1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-o") == 0) {
      ovlName = argv[++arg];

    } else if (strcmp(argv[arg], "-g") == 0) {
      gkpName = argv[++arg];

    } else if (strcmp(argv[arg], "-F") == 0) {
      fileLimit    = atoi(argv[++arg]);
      memoryLimit  = 0;

    } else if (strcmp(argv[arg], "-M") == 0) {
      fileLimit    = 0;
      memoryLimit  = atoi(argv[++arg]);
      memoryLimit *= 1024;
      memoryLimit *= 1024;
      memoryLimit *= 1024;

    } else if (strcmp(argv[arg], "-plc") == 0) {
      //  Former -i option
      //  PLC_NONE, PLC_ALL, PLC_INTERNAL
      ovlSkipOpt = PLC_ALL;

    } else if (strcmp(argv[arg], "-obt") == 0) {
      doFilterOBT = 1;

    } else if (strcmp(argv[arg], "-dup") == 0) {
      doFilterOBT = 2;

    } else if (strcmp(argv[arg], "-e") == 0) {
      maxError = atof(argv[++arg]);
      maxError = AS_OVS_encodeQuality(maxErrorRate);

    } else if (strcmp(argv[arg], "-L") == 0) {
      errno = 0;
      FILE *F = fopen(argv[++arg], "r");
      if (errno)
        fprintf(stderr, "Can't open '%s': %s\n", argv[arg], strerror(errno)), exit(1);

      char *line = new char [FILENAME_MAX];

      fgets(line, FILENAME_MAX, F);

      while (!feof(F)) {
        chomp(line);
        fileList.push_back(line);
        line = new char [FILENAME_MAX];
        fgets(line, FILENAME_MAX, F);
      }

      delete [] line;

      fclose(F);

    } else if ((argv[arg][0] == '-') && (argv[arg][1] != 0)) {
      fprintf(stderr, "%s: unknown option '%s'.\n", argv[0], argv[arg]);
      err++;

    } else {
      //  Assume it's an input file
      fileList.push_back(argv[arg]);
    }

    arg++;
  }
  if (ovlName == NULL)
    err++;
  if (gkpName == NULL)
    err++;
  if (fileList.size() == 0)
    err++;
  if (fileLimit > sysconf(_SC_OPEN_MAX) - 16)
    err++;
  if (err) {
    fprintf(stderr, "usage: %s -o asm.ovlStore -g asm.gkpStore [opts] [-L fileList | *.ovb.gz]\n", argv[0]);
    fprintf(stderr, "  -o asm.ovlStore       path to store to create\n");
    fprintf(stderr, "  -g asm.gkpStore       path to gkpStore for this assembly\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -F f                  use up to 'f' files for store creation\n");
    fprintf(stderr, "  -M m                  use up to 'm' MB memory for store creation\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -plc t                type of filtering for PLC fragments -- NOT SUPPORTED\n");
    fprintf(stderr, "  -obt                  filter overlaps for OBT\n");
    fprintf(stderr, "  -dup                  filter overlaps for OBT/dedupe\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -e e                  filter overlaps above e fraction error\n");
    fprintf(stderr, "  -L fileList           read input filenames from 'flieList'\n");

    if (ovlName == NULL)
      fprintf(stderr, "ERROR: No overlap store (-o) supplied.\n");
    if (gkpName == NULL)
      fprintf(stderr, "ERROR: No gatekeeper store (-g) supplied.\n");
    if (fileList.size() == 0)
      fprintf(stderr, "ERROR: No input overlap files (-L or last on the command line) supplied.\n");
    if (fileLimit > sysconf(_SC_OPEN_MAX) - 16)
      fprintf(stderr, "ERROR: Too many jobs (-F); only "F_SIZE_T" supported on this architecture.\n", sysconf(_SC_OPEN_MAX) - 16);

    exit(1);
  }


  //  We create the store early, allowing it to fail if it already
  //  exists, or just cannot be created.
  //
  OverlapStore    *storeFile = AS_OVS_createOverlapStore(ovlName, TRUE);

  gkStore *gkp         = new gkStore(gkpName, FALSE, FALSE);

  uint64  maxIID       = gkp->gkStore_getNumFragments() + 1;
  //uint64  iidPerBucket = (uint64)ceil((double)maxIID / (double)fileLimit);
  uint64  iidPerBucket = computeIIDperBucket(fileLimit, memoryLimit, maxIID, fileList);

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


  

  for (uint32 i=0; i<fileList.size(); i++) {
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

      writeToDumpFile(&fovrlap, dumpFile, dumpFileMax, dumpLength, iidPerBucket, ovlName);

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

          writeToDumpFile(&rovrlap, dumpFile, dumpFileMax, dumpLength, iidPerBucket, ovlName);
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

          writeToDumpFile(&rovrlap, dumpFile, dumpFileMax, dumpLength, iidPerBucket, ovlName);
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

  OVSoverlap  *overlapsort = new OVSoverlap [dumpLengthMax];

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

    sprintf(name, "%s/tmp.sort.%03d", ovlName, i);
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
    sort(overlapsort, overlapsort + dumpLength[i]);

    fprintf(stderr, "writing %s (%ld)\n", name, time(NULL) - beginTime);
    for (uint64 x=0; x<dumpLength[i]; x++)
      AS_OVS_writeOverlapToStore(storeFile, overlapsort + x);
  }

  AS_OVS_closeOverlapStore(storeFile);

  delete [] overlapsort;

  //  And we have a store.
  //
  exit(0);
}
