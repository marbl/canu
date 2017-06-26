
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
 *    src/AS_OVS/overlapStoreSorter.C
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2012-APR-02 to 2013-SEP-09
 *      are Copyright 2012-2013 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2014-DEC-15 to 2015-SEP-21
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

#include "AS_global.H"

#include "gkStore.H"
#include "ovStore.H"

#include <vector>
#include <algorithm>

using namespace std;


//  This is the size of the datastructure that we're using to store overlaps for sorting.
//  At present, with ovOverlap, it is over-allocating a pointer that we don't need, but
//  to make a custom structure, we'd need to duplicate a bunch of code or copy data after
//  loading and before writing.
//
//  Used in both ovStoreSorter.C and ovStoreBuild.C.
//
#define ovOverlapSortSize  (sizeof(ovOverlap))



void
makeSentinel(char *storePath, uint32 fileID, bool forceRun) {
  char name[FILENAME_MAX];

  //  Check if done.

  snprintf(name, FILENAME_MAX, "%s/%04d", storePath, fileID);

  if ((forceRun == false) && (AS_UTL_fileExists(name, FALSE, FALSE)))
    fprintf(stderr, "Job " F_U32 " is finished (remove '%s' or -force to try again).\n", fileID, name), exit(0);

  //  Check if running.

  snprintf(name, FILENAME_MAX, "%s/%04d.ovs", storePath, fileID);

  if ((forceRun == false) && (AS_UTL_fileExists(name, FALSE, FALSE)))
    fprintf(stderr, "Job " F_U32 " is running (remove '%s' or -force to try again).\n", fileID, name), exit(0);

  //  Not done, not running, so create a sentinel to say we're running.

  errno = 0;
  FILE *F = fopen(name, "w");
  if (errno)
    fprintf(stderr, "ERROR: Failed to open '%s' for writing: %s\n", name, strerror(errno)), exit(1);
  fclose(F);
}



void
removeSentinel(char *storePath, uint32 fileID) {
  char name[FILENAME_MAX];
  snprintf(name, FILENAME_MAX, "%s/%04d.ovs", storePath, fileID);
  unlink(name);
}



int
main(int argc, char **argv) {
  char           *storePath      = NULL;
  char           *gkpName        = NULL;

  uint32          fileLimit      = 512;   //  Number of 'slices' from bucketizer
  uint32          fileID         = 0;     //  'slice' that we are going to be sorting
  uint32          jobIdxMax      = 0;     //  Number of 'buckets' from bucketizer

  uint64          maxMemory      = UINT64_MAX;

  bool            deleteIntermediateEarly = false;
  bool            deleteIntermediateLate  = false;

  bool            forceRun = false;

  char            name[FILENAME_MAX];

  argc = AS_configure(argc, argv);

  int err=0;
  int arg=1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-O") == 0) {
      storePath = argv[++arg];

    } else if (strcmp(argv[arg], "-G") == 0) {
      gkpName = argv[++arg];

    } else if (strcmp(argv[arg], "-F") == 0) {
      fileLimit = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-job") == 0) {
      fileID  = atoi(argv[++arg]);
      jobIdxMax = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-M") == 0) {
      maxMemory  = (uint64)ceil(atof(argv[++arg]) * 1024.0 * 1024.0 * 1024.0);

    } else if (strcmp(argv[arg], "-deleteearly") == 0) {
      deleteIntermediateEarly = true;

    } else if (strcmp(argv[arg], "-deletelate") == 0) {
      deleteIntermediateLate  = true;

    } else if (strcmp(argv[arg], "-force") == 0) {
      forceRun = true;

    } else {
      fprintf(stderr, "ERROR: unknown option '%s'\n", argv[arg]);
      err++;
    }

    arg++;
  }
  if (storePath == NULL)
    err++;
  if (fileID == 0)
    err++;
  if (jobIdxMax == 0)
    err++;

  if (err) {
    fprintf(stderr, "usage: %s ...\n", argv[0]);
    fprintf(stderr, "  -G asm.gkpStore  path to gkpStore for this assembly\n");
    fprintf(stderr, "  -O x.ovlStore    path to overlap store to build the final index for\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -F s             number of slices used in bucketizing/sorting\n");
    fprintf(stderr, "  -job j m         index of this overlap input file, and max number of files\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -M m             maximum memory to use, in gigabytes\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -deleteearly     remove intermediates as soon as possible (unsafe)\n");
    fprintf(stderr, "  -deletelate      remove intermediates when outputs exist (safe)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -force           force a recompute, even if the output exists\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    DANGER    DO NOT USE     DO NOT USE     DO NOT USE    DANGER\n");
    fprintf(stderr, "    DANGER                                                DANGER\n");
    fprintf(stderr, "    DANGER   This command is difficult to run by hand.    DANGER\n");
    fprintf(stderr, "    DANGER          Use ovStoreCreate instead.            DANGER\n");
    fprintf(stderr, "    DANGER                                                DANGER\n");
    fprintf(stderr, "    DANGER    DO NOT USE     DO NOT USE     DO NOT USE    DANGER\n");
    fprintf(stderr, "\n");

    if (storePath == NULL)
      fprintf(stderr, "ERROR: No overlap store (-O) supplied.\n");
    if (fileID == 0)
      fprintf(stderr, "ERROR: no slice number (-F) supplied.\n");
    if (jobIdxMax == 0)
      fprintf(stderr, "ERROR: no max job ID (-job) supplied.\n");

    exit(1);
  }

  //  Check if we're running or done (or crashed), then note that we're running.

  makeSentinel(storePath, fileID, forceRun);

  //  Not done.  Let's go!

  gkStore        *gkp    = gkStore::gkStore_open(gkpName);
  ovStoreWriter  *writer = new ovStoreWriter(storePath, gkp, fileLimit, fileID, jobIdxMax);

  //  Get the number of overlaps in each bucket slice.

  fprintf(stderr, "\n");
  fprintf(stderr, "Finding overlaps.\n");

  uint64 *bucketSizes = new uint64 [jobIdxMax + 1];
  uint64  totOvl      = writer->loadBucketSizes(bucketSizes);

  //  Fail if we don't have enough memory to process.

  if (ovOverlapSortSize * totOvl > maxMemory) {
    fprintf(stderr, "ERROR:  Overlaps need %.2f GB memory, but process limited (via -M) to " F_U64 " GB.\n",
            ovOverlapSortSize * totOvl / 1024.0 / 1024.0 / 1024.0, maxMemory >> 30);
    removeSentinel(storePath, fileID);
    exit(1);
  }

  //  Or report that we can process.

  fprintf(stderr, "\n");
  fprintf(stderr, "Loading %10" F_U64P " overlaps using %.2f GB of requested (-M) " F_U64 " GB memory.\n",
          totOvl, ovOverlapSortSize * totOvl / 1024.0 / 1024.0 / 1024.0, maxMemory >> 30);

  //  Load all overlaps - we're guaranteed that either 'name.gz' or 'name' exists (we checked when
  //  we loaded bucket sizes) or funny business is happening with our files.

  ovOverlap *ovls   = ovOverlap::allocateOverlaps(gkp, totOvl);
  uint64    ovlsLen = 0;

  for (uint32 i=0; i<=jobIdxMax; i++)
    writer->loadOverlapsFromSlice(i, bucketSizes[i], ovls, ovlsLen);

  //  Check that we found all the overlaps we were expecting.

  if (ovlsLen != totOvl)
    fprintf(stderr, "ERROR: read " F_U64 " overlaps, expected " F_U64 "\n", ovlsLen, totOvl);
  assert(ovlsLen == totOvl);

  //  Clean up space if told to.

  if (deleteIntermediateEarly)
    writer->removeOverlapSlice();

  //  Sort the overlaps!  Finally!  The parallel STL sort is NOT inplace, and blows up our memory.

  fprintf(stderr, "\n");
  fprintf(stderr, "Sorting.\n");

#ifdef _GLIBCXX_PARALLEL
  __gnu_sequential::sort(ovls, ovls + ovlsLen);
#else
  sort(ovls, ovls + ovlsLen);
#endif

  //  Output to the store.

  fprintf(stderr, "\n");   //  Sorting has no output, so this would generate a distracting extra newline
  fprintf(stderr, "Writing sorted overlaps.\n");

  writer->writeOverlaps(ovls, ovlsLen);

  //  Clean up.  Delete inputs, remove the sentinel, release memory, etc.

  delete [] ovls;
  delete [] bucketSizes;

  removeSentinel(storePath, fileID);

  gkp->gkStore_close();

  if (deleteIntermediateLate) {
    fprintf(stderr, "\n");
    fprintf(stderr, "Removing bucketized overlaps.\n");
    fprintf(stderr, "\n");

    writer->removeOverlapSlice();
  }

  //  Success!

  fprintf(stderr, "Success!\n");

  return(0);
}
