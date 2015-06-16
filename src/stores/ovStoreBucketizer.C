
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

#include <vector>
#include <algorithm>

using namespace std;

static
void
writeToFile(ovOverlap    *overlap,
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
    fprintf(stderr, "  Aid "F_U32"  Bid "F_U32"\n",  overlap->a_iid, overlap->b_iid);
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
markOBT(gkStore *gkp, uint32 maxIID, char *skipRead) {
  uint64  numMarked = 0;

  if (skipRead == NULL)
    return;

  fprintf(stderr, "Marking fragments to skip overlap based trimming.\n");

  for (uint64 iid=0; iid<maxIID; iid++) {
    uint32     Lid = gkp->gkStore_getRead(iid)->gkRead_libraryID();
    gkLibrary *L   = gkp->gkStore_getLibrary(Lid);

    if (L == NULL)
      continue;

    if ((L->gkLibrary_removeDuplicateReads()     == false) &&
        (L->gkLibrary_finalTrim()                != FINALTRIM_LARGEST_COVERED) &&
        (L->gkLibrary_removeSpurReads()          == false) &&
        (L->gkLibrary_removeChimericReads()      == false)) {
      numMarked++;
      skipRead[iid] = true;
    }
  }

  fprintf(stderr, "Marked "F_U64" fragments.\n", numMarked);
}


static
void
markDUP(gkStore *gkp, uint32 maxIID, char *skipRead) {
  uint64  numMarked = 0;

  if (skipRead == NULL)
    return;

  fprintf(stderr, "Marking fragments to skip deduplication.\n");

  for (uint64 iid=0; iid<maxIID; iid++) {
    uint32     Lid = gkp->gkStore_getRead(iid)->gkRead_libraryID();
    gkLibrary *L   = gkp->gkStore_getLibrary(Lid);

    if (L == NULL)
      continue;

    if (L->gkLibrary_removeDuplicateReads() == false) {
      numMarked++;
      skipRead[iid] = true;
    }
  }

  fprintf(stderr, "Marked "F_U64" fragments.\n", numMarked);
}





int
main(int argc, char **argv) {
  char           *ovlName      = NULL;
  char           *gkpName      = NULL;
  uint32          fileLimit    = 512;

  uint32          jobIndex     = 0;

  double          maxErrorRate = 1.0;
  uint64          maxError     = AS_OVS_encodeEvalue(maxErrorRate);

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

    } else if (strcmp(argv[arg], "-job") == 0) {
      jobIndex = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-i") == 0) {
      ovlInput = argv[++arg];

    } else if (strcmp(argv[arg], "-e") == 0) {
      maxErrorRate = atof(argv[++arg]);
      maxError     = AS_OVS_encodeEvalue(maxErrorRate);

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

  fprintf(stderr, "maxError fraction: %.3f percent: %.3f encoded: "F_U64"\n",
          maxErrorRate, maxErrorRate * 100, maxError);

  fprintf(stderr, "Bucketizing %s\n", ovlInput);

  ovStoreFilter *filter = new ovStoreFilter(gkp, maxError);
  ovOverlap     foverlap;
  ovOverlap     roverlap;
  ovFile        *inputFile = new ovFile(ovlInput, ovFileFull);

  //  Do bigger buffers increase performance?  Do small ones hurt?
  //AS_OVS_setBinaryOverlapFileBufferSize(2 * 1024 * 1024);

  while (inputFile->readOverlap(&foverlap)) {
    filter->filterOverlap(foverlap, roverlap);  //  The filter copies f into r

    //  If all are skipped, don't bother writing the overlap.

    if ((foverlap.dat.ovl.forUTG == true) ||
        (foverlap.dat.ovl.forOBT == true) ||
        (foverlap.dat.ovl.forDUP == true))
      writeToFile(&foverlap, sliceFile, fileLimit, sliceSize, iidPerBucket, ovlName, jobIndex, useGzip);

    if ((roverlap.dat.ovl.forUTG == true) ||
        (roverlap.dat.ovl.forOBT == true) ||
        (roverlap.dat.ovl.forDUP == true))
      writeToFile(&roverlap, sliceFile, fileLimit, sliceSize, iidPerBucket, ovlName, jobIndex, useGzip);
  }

  delete inputFile;

  filter->reportFate();
  filter->resetCounters();

  delete filter;

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

  delete [] sliceFile;
  delete [] sliceSize;

  return(0);
}
