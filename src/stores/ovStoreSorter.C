
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
 *    Brian P. Walenz beginning on 2014-DEC-15
 *      are Copyright 2014-2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

const char *mainid = "$Id$";

#include "AS_global.H"

#include "gkStore.H"
#include "ovStore.H"

#include <vector>
#include <algorithm>

using namespace std;


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

  argc = AS_configure(argc, argv);

  int err=0;
  int arg=1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-o") == 0) {
      storePath = argv[++arg];

    } else if (strcmp(argv[arg], "-g") == 0) {
      gkpName = argv[++arg];

    } else if (strcmp(argv[arg], "-F") == 0) {
      fileLimit = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-job") == 0) {
      fileID  = atoi(argv[++arg]);
      jobIdxMax = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-M") == 0) {
      maxMemory  = atoi(argv[++arg]);
      maxMemory *= 1024;
      maxMemory *= 1024;
      maxMemory *= 1024;

    } else if (strcmp(argv[arg], "-deleteearly") == 0) {
      deleteIntermediateEarly = true;

    } else if (strcmp(argv[arg], "-deletelate") == 0) {
      deleteIntermediateLate  = true;

    } else if (strcmp(argv[arg], "-force") == 0) {
      forceRun = true;

    } else {
      fprintf(stderr, "ERROR: unknown option '%s'\n", argv[arg]);
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
    exit(1);
  }

  //  Check if we're running or done (or crashed), then note that we're running.

  {
    char name[FILENAME_MAX];
    sprintf(name,"%s/%04d.ovs", storePath, fileID);

    if ((forceRun == false) && (AS_UTL_fileExists(name, FALSE, FALSE)))
      fprintf(stderr, "Job "F_U32" is running or finished (remove '%s' or -force to try again).\n", fileID, name), exit(0);

    errno = 0;
    FILE *F = fopen(name, "w");
    if (errno)
      fprintf(stderr, "ERROR: Failed to open '%s' for writing: %s\n", name, strerror(errno)), exit(1);

    fclose(F);
  }

  // Get sizes of each bucket, and the final merge

  uint64   *sliceSizes    = new uint64 [fileLimit + 1];  //  For each overlap job, number of overlaps per bucket
  uint64   *bucketSizes   = new uint64 [jobIdxMax + 1];  //  For each bucket we care about, number of overlaps

  uint64    totOvl        = 0;
  uint64    ovlsLen        = 0;

  for (uint32 i=0; i<=jobIdxMax; i++) {
    bucketSizes[i] = 0;

    char namz[FILENAME_MAX];
    char name[FILENAME_MAX];

    sprintf(namz, "%s/bucket%04d/slice%03d.gz", storePath, i, fileID);
    sprintf(name, "%s/bucket%04d/slice%03d",    storePath, i, fileID);

    if ((AS_UTL_fileExists(namz, FALSE, FALSE) == false) &&
        (AS_UTL_fileExists(name, FALSE, FALSE) == false))
      //  If no file, there are no overlaps.  Skip loading the bucketSizes file.
      //  We expect the gz version to exist (that's the default in bucketizer) more frequently, so
      //  be sure to test for existence of that one first.
      continue;

    sprintf(name, "%s/bucket%04d/sliceSizes", storePath, i);

    FILE *F = fopen(name, "r");
    if (errno)
      fprintf(stderr, "ERROR:  Failed to open %s: %s\n", name, strerror(errno)), exit(1);

    uint64 nr = AS_UTL_safeRead(F, sliceSizes, "sliceSizes", sizeof(uint64), fileLimit + 1);

    fclose(F);

    if (nr != fileLimit + 1) {
      fprintf(stderr, "ERROR: short read on '%s'.\n", name);
      fprintf(stderr, "ERROR: read "F_U64" sizes insteadof "F_U32".\n", nr, fileLimit + 1);
    }
    assert(nr == fileLimit + 1);

    fprintf(stderr, "Found "F_U64" overlaps from '%s'.\n", sliceSizes[fileID], name);

    bucketSizes[i] = sliceSizes[fileID];
    totOvl        += sliceSizes[fileID];
  }

  delete [] sliceSizes;
  sliceSizes = NULL;

  if (sizeof(ovOverlap) * totOvl > maxMemory) {
    fprintf(stderr, "ERROR:  Overlaps need %.2f GB memory, but process limited (via -M) to "F_U64" GB.\n",
            sizeof(ovOverlap) * totOvl / (1024.0 * 1024.0 * 1024.0), maxMemory >> 30);

    char name[FILENAME_MAX];
    sprintf(name,"%s/%04d.ovs", storePath, fileID);

    unlink(name);

    exit(1);
  }

  fprintf(stderr, "Overlaps need %.2f GB memory, allowed to use up to (via -M) "F_U64" GB.\n",
          sizeof(ovOverlap) * totOvl / (1024.0 * 1024.0 * 1024.0), maxMemory >> 30);

  ovOverlap *ovls = ovOverlap::allocateOverlaps(NULL, totOvl);

  //  Load all overlaps - we're guaranteed that either 'name.gz' or 'name' exists (we checked above)
  //  or funny business is happening with our files.

  for (uint32 i=0; i<=jobIdxMax; i++) {
    if (bucketSizes[i] == 0)
      continue;

    char name[FILENAME_MAX];

    sprintf(name, "%s/bucket%04d/slice%03d.gz", storePath, i, fileID);
    if (AS_UTL_fileExists(name, FALSE, FALSE) == false)
      sprintf(name, "%s/bucket%04d/slice%03d", storePath, i, fileID);

    if (AS_UTL_fileExists(name, FALSE, FALSE) == false)
      fprintf(stderr, "ERROR: "F_U64" overlaps claim to exist in bucket '%s', but file not found.\n",
              bucketSizes[i], name);

    fprintf(stderr, "Loading "F_U64" overlaps from '%s'.\n", bucketSizes[i], name);

    ovFile   *bof = new ovFile(name, ovFileFull);
    uint64    num = 0;

    while (bof->readOverlap(ovls + ovlsLen)) {
      ovlsLen++;
      num++;
    }

    if (num != bucketSizes[i])
      fprintf(stderr, "ERROR: expected "F_U64" overlaps, found "F_U64" overlaps.\n", bucketSizes[i], num);
    assert(num == bucketSizes[i]);

    delete bof;
  }

  if (ovlsLen != totOvl)
    fprintf(stderr, "ERROR: read "F_U64" overlaps, expected "F_U64"\n", ovlsLen, totOvl);
  assert(ovlsLen == totOvl);

  if (deleteIntermediateEarly) {
    char name[FILENAME_MAX];

    fprintf(stderr, "Removing inputs.\n");
    for (uint32 i=0; i<=jobIdxMax; i++) {
      if (bucketSizes[i] == 0)
        continue;

      sprintf(name, "%s/bucket%04d/slice%03d.gz", storePath, i, fileID);
      AS_UTL_unlink(name);

      sprintf(name, "%s/bucket%04d/slice%03d", storePath, i, fileID);
      AS_UTL_unlink(name);
    }
  }

  //  Sort the overlaps - at least on FreeBSD 8.2 with gcc46, the parallel STL sort
  //  algorithms are NOT inplace.  Restrict to sequential sorting.
  //
  //  This sort takes at most 2 minutes on 7gb of overlaps.
  //
  fprintf(stderr, "Sorting.\n");

#ifdef _GLIBCXX_PARALLEL
  //  If we have the parallel STL, don't use it!  Sort is not inplace!
  __gnu_sequential::sort(ovls, ovls + ovlsLen);
#else
  sort(ovls, ovls + ovlsLen);
#endif

  //  Output to store format

  fprintf(stderr, "Writing output.\n");
  writeOverlaps(storePath, ovls, ovlsLen, fileID);

  //  Clean up.

  delete [] ovls;

  if (deleteIntermediateLate) {
    char name[FILENAME_MAX];

    fprintf(stderr, "Removing inputs.\n");
    for (uint32 i=0; i<=jobIdxMax; i++) {
      if (bucketSizes[i] == 0)
        continue;

      sprintf(name, "%s/bucket%04d/slice%03d.gz", storePath, i, fileID);
      AS_UTL_unlink(name);

      sprintf(name, "%s/bucket%04d/slice%03d", storePath, i, fileID);
      AS_UTL_unlink(name);
    }
  }

  delete [] bucketSizes;

  //  Remove the sentinel to show we're done.  The output is in "%s/%04d".

  {
    char name[FILENAME_MAX];
    sprintf(name,"%s/%04d.ovs", storePath, fileID);

    unlink(name);
  }

  //  Success!

  return(0);
}
