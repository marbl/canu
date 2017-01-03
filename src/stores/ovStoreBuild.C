
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
 *    src/AS_OVS/overlapStoreBuild.C
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2012-APR-02 to 2013-AUG-01
 *      are Copyright 2012-2013 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2014-JUL-31 to 2015-SEP-21
 *      are Copyright 2014-2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2015-NOV-03
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *    Sergey Koren beginning on 2016-FEB-20
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_global.H"
#include "AS_UTL_decodeRange.H"

#include "gkStore.H"
#include "ovStore.H"

#include <vector>
#include <algorithm>

using namespace std;

#define  MEMORY_OVERHEAD  (256 * 1024 * 1024)

//  This is the size of the datastructure that we're using to store overlaps for sorting.
//  At present, with ovOverlap, it is over-allocating a pointer that we don't need, but
//  to make a custom structure, we'd need to duplicate a bunch of code or copy data after
//  loading and before writing.
//
//  Used in both ovStoreSorter.C and ovStoreBuild.C.
//
#define ovOverlapSortSize  (sizeof(ovOverlap))



static
void
addEvalues(char *ovlName, vector<char *> &fileList) {
  ovStore  *ovs = new ovStore(ovlName, NULL);

  ovs->addEvalues(fileList);

  delete ovs;

  fprintf(stderr, "-  Evalues updated.\n");
}



void
reportConfiguration(char *configOut, uint32 maxIID, uint32 *iidToBucket) {
  char  F[FILENAME_MAX];

  snprintf(F, FILENAME_MAX, "%s.WORKING", configOut);

  errno = 0;
  FILE *C = fopen(configOut, "w");
  if (errno)
    fprintf(stderr, "Failed to open config output file '%s': %s\n", configOut, strerror(errno)), exit(1);

  AS_UTL_safeWrite(C, &maxIID,      "maxIID",      sizeof(uint32), 1);
  AS_UTL_safeWrite(C,  iidToBucket, "iidToBucket", sizeof(uint32), maxIID);

  fclose(C);

  rename(F, configOut);

  delete [] iidToBucket;

  fprintf(stderr, "-  Saved configuration to '%s'.\n", configOut);
}



static
uint32 *
computeIIDperBucket(uint32          fileLimit,
                    uint64          minMemory,
                    uint64          maxMemory,
                    uint32          maxIID,
                    vector<char *> &fileList) {
  uint32  *iidToBucket   = new uint32 [maxIID];
  int64    procMax       = sysconf(_SC_CHILD_MAX);
  int64    openMax       = sysconf(_SC_OPEN_MAX);
  int64    maxFiles      = 0;

  //  As of late August 2016, the ovb files are not gzip compressed, and do not need to use an
  //  external process to decompress.  The support for limiting by number of processes is left in -
  //  but disabled - because it's really just these three lines here, and because the code still
  //  supports gzip inputs.

  if (openMax > 16)
    openMax -= 16;

  if (procMax > 8192)  //  Once saw a case where procMax was 18,446,744,073,709,551,615 (2^64-1)
    procMax = 8192;    //  and openMax was 262,144.  It didn't end well.

  maxFiles    = openMax;  //  MIN(procMax, openMax);    ENABLE THIS TO LIMIT PROCESSES TOO.

  //  If we're reading from stdin, not much we can do but divide the IIDs equally per file.  Note
  //  that the IIDs must be consecutive; the obvious, simple and clean division of 'mod' won't work.

  if (fileList[0][0] == '-') {
    if (maxMemory > 0) {
      minMemory = 0;
      maxMemory = 0;
      fileLimit = maxFiles;

      fprintf(stderr, "WARNING: memory limit (-M) specified, but can't be used with inputs from stdin; using %d files instead.\n", fileLimit);
    } else {
      fprintf(stderr, "Sorting overlaps from stdin using %d files.\n", fileLimit);
    }

    uint32  iidPerBucket   = maxIID / fileLimit;
    uint32  thisBucket     = 1;
    uint32  iidThisBucket  = 0;

    for (uint32 ii=0; ii<maxIID; ii++) {
      iidThisBucket++;
      iidToBucket[ii] = thisBucket;

      if (iidThisBucket > iidPerBucket) {
        iidThisBucket = 0;
        thisBucket++;
      }
    }

    return(iidToBucket);
  }

  //  Otherwise, we have files, and should have counts.  Load them!

  ovStoreHistogram   *hist = new ovStoreHistogram();
  uint32             *oPR = NULL;

  allocateArray(oPR, maxIID);

  for (uint32 i=0; i<fileList.size(); i++)
    hist->loadData(fileList[i]);

  uint64   numOverlaps = hist->getOverlapsPerRead(oPR, maxIID);

  delete hist;  hist = NULL;

  if (numOverlaps == 0)
    fprintf(stderr, "Found no overlaps to sort.\n"), exit(1);

  fprintf(stderr, "Found " F_U64 " (%.2f million) overlaps.\n", numOverlaps, numOverlaps / 1000000.0);

  //  Partition the overlaps into buckets.

  uint64   olapsPerBucketMax = 1;
  double   GBperOlap         = ovOverlapSortSize / 1024.0 / 1024.0 / 1024.0;

  //  If a file limit, distribute the overlaps to equal sized files.
  if (fileLimit > 0) {
    olapsPerBucketMax = (uint64)ceil((double)numOverlaps / (double)fileLimit);
    fprintf(stderr, "Will sort using " F_U32 " files; " F_U64 " (%.2f million) overlaps per bucket; %.2f GB memory per bucket\n",
            fileLimit, olapsPerBucketMax, olapsPerBucketMax / 1000000.0, olapsPerBucketMax * GBperOlap);
  }

  //  If a memory limit, distribute the overlaps to files no larger than the limit.
  //
  //  This will pick the smallest memory size that uses fewer than maxFiles buckets.  Unreasonable
  //  values can break this - either too low memory or too high allowed open files (an OS limit).

  if (maxMemory > 0) {
    fprintf(stderr, "Configuring for %.2f GB to %.2f GB memory and " F_S64 " open files.\n",
            minMemory / 1024.0 / 1024.0 / 1024.0,
            maxMemory / 1024.0 / 1024.0 / 1024.0,
            maxFiles);

    if (minMemory < MEMORY_OVERHEAD + ovOverlapSortSize) {
      fprintf(stderr, "Reset minMemory from " F_U64 " to " F_SIZE_T "\n", minMemory, MEMORY_OVERHEAD + ovOverlapSortSize);
      minMemory  = MEMORY_OVERHEAD + ovOverlapSortSize;
    }

    uint64  incr = (maxMemory - minMemory) / 128;
    if (incr < 1024 * 1024)
      incr = 1024 * 1024;

    uint64  useMemory = minMemory;

    //  Compute the initial number of overlaps per bucket, based on the smallest memory allowed.

    olapsPerBucketMax = (useMemory - MEMORY_OVERHEAD) / ovOverlapSortSize;

    //  Find the smallest memory size that uses fewer files than the OS allows.

    for (;
         ((useMemory <= maxMemory) && (numOverlaps / olapsPerBucketMax + 1 > maxFiles));
         useMemory += incr) {
      olapsPerBucketMax = (useMemory - MEMORY_OVERHEAD) / ovOverlapSortSize;
      fprintf(stderr, "At memory %.3fGB, " F_U64 " olaps per bucket, " F_U64 " buckets (pass 1).\n",
              useMemory / 1024.0 / 1024.0 / 1024.0, olapsPerBucketMax, numOverlaps / olapsPerBucketMax + 1);
    }

    //  If we're at less than half the max, make buckets a little bit bigger to reduce the open file
    //  count.  This helps when multiple bucketizer jobs get scheduled to the same node.

    if (useMemory < minMemory + (maxMemory - minMemory) / 2) {
      for (;
           ((useMemory <= maxMemory) && (numOverlaps / olapsPerBucketMax + 1 > maxFiles / 2));
           useMemory += incr) {
        olapsPerBucketMax = (useMemory - MEMORY_OVERHEAD) / ovOverlapSortSize;
        fprintf(stderr, "At memory %.3fGB, " F_U64 " olaps per bucket, " F_U64 " buckets (pass 2).\n",
                useMemory / 1024.0 / 1024.0 / 1024.0, olapsPerBucketMax, numOverlaps / olapsPerBucketMax + 1);
      }
    }

    //  Give up if we hit our max limit.

    if ((minMemory > maxMemory) ||
        (numOverlaps / olapsPerBucketMax + 1) > maxFiles) {
      fprintf(stderr, "ERROR:  Cannot sort %.2f million overlaps using %.2f GB memory; too few file handles available.\n",
              numOverlaps / 1000000.0,
              maxMemory / 1024.0 / 1024.0 / 1024.0);
      fprintf(stderr, "ERROR:    minMemory      " F_U64 "\n", minMemory);
      fprintf(stderr, "ERROR:    maxMemory      " F_U64 "\n", maxMemory);
      fprintf(stderr, "ERROR:    olapsPerBucket " F_U64 "\n", olapsPerBucketMax);
      fprintf(stderr, "ERROR:    buckets        " F_U64 "\n", numOverlaps / olapsPerBucketMax + 1);
      fprintf(stderr, "ERROR:    SC_CHILD_MAX   " F_S64 "\n", (int64)sysconf(_SC_CHILD_MAX));
      fprintf(stderr, "ERROR:    SC_OPEN_MAX    " F_S64 "\n", (int64)sysconf(_SC_OPEN_MAX));
      fprintf(stderr, "ERROR:  Increase memory size (in canu, ovsMemory; in ovStoreBuild, -M)\n");
      exit(1);
    }

    fprintf(stderr, "Will sort using " F_U64 " files; " F_U64 " (%.2f million) overlaps per bucket; %.2f GB memory per bucket\n",
            numOverlaps / olapsPerBucketMax + 1,
            olapsPerBucketMax,
            olapsPerBucketMax / 1000000.0,
            olapsPerBucketMax * GBperOlap + MEMORY_OVERHEAD / 1024.0 / 1024.0 / 1024.0);
  }

  //  Given the limit on each bucket, count the number of buckets needed, then reset the limit on
  //  each bucket to have the same number of overlaps for every bucket.

  {
    uint64  olaps  = 0;
    uint32  bucket = 1;

    for (uint32 ii=0; ii<maxIID; ii++) {
      olaps            += oPR[ii];
      iidToBucket[ii]   = bucket;

      if (olaps >= olapsPerBucketMax) {
        olaps = 0;
        bucket++;
      }
    }

    olapsPerBucketMax = (uint64)ceil((double)numOverlaps / (double)bucket);
  }

  //  And, finally, assign IIDs to buckets.

  {
    uint64  olaps  = 0;
    uint32  bucket = 1;

    for (uint32 ii=0; ii<maxIID; ii++) {
      olaps            += oPR[ii];
      iidToBucket[ii]   = bucket;

      if (olaps >= olapsPerBucketMax) {
        fprintf(stderr, "  bucket %3d has " F_U64 " olaps.\n", bucket, olaps);
        olaps = 0;
        bucket++;
      }
    }

    fprintf(stderr, "  bucket %3d has " F_U64 " olaps.\n", bucket, olaps);
  }

  fprintf(stderr, "Will sort %.3f million overlaps per bucket, using %u buckets %.2f GB per bucket.\n",
          olapsPerBucketMax / 1000000.0,
          iidToBucket[maxIID-1],
          olapsPerBucketMax * GBperOlap + MEMORY_OVERHEAD / 1024.0 / 1024.0 / 1024.0);
  fprintf(stderr, "\n");

  delete hist;

  return(iidToBucket);
}



static
void
writeToDumpFile(gkStore          *gkp,
                ovOverlap       *overlap,
                ovFile          **dumpFile,
                uint64           *dumpLength,
                uint32           *iidToBucket,
                char             *ovlName) {
  uint32 df = iidToBucket[overlap->a_iid];

  //  If the dump file isn't opened, open it.

  if (dumpFile[df] == NULL) {
    char name[FILENAME_MAX];
    snprintf(name, FILENAME_MAX, "%s/tmp.sort.%03d", ovlName, df);
    fprintf(stderr, "-- Create bucket '%s'\n", name);
    dumpFile[df]   = new ovFile(gkp, name, ovFileFullWriteNoCounts);
    dumpLength[df] = 0;
  }

  //  And write the overlap.

  dumpFile[df]->writeOverlap(overlap);
  dumpLength[df]++;
}



int
main(int argc, char **argv) {
  char           *ovlName        = NULL;
  char           *gkpName        = NULL;
  uint32          fileLimit      = 0;
  uint64          minMemory      = (uint64)1 * 1024 * 1024 * 1024;
  uint64          maxMemory      = (uint64)4 * 1024 * 1024 * 1024;

  double          maxError     = 1.0;
  uint32          minOverlap   = 0;

  vector<char *>  fileList;

  uint32          nThreads     = 4;

  bool            eValues      = false;
  char           *configOut    = NULL;

  argc = AS_configure(argc, argv);

  int err=0;
  int arg=1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-O") == 0) {
      ovlName = argv[++arg];

    } else if (strcmp(argv[arg], "-G") == 0) {
      gkpName = argv[++arg];

    } else if (strcmp(argv[arg], "-F") == 0) {
      fileLimit    = atoi(argv[++arg]);
      minMemory    = 0;
      maxMemory    = 0;

    } else if (strcmp(argv[arg], "-M") == 0) {
      double lo=0.0, hi=0.0;

      AS_UTL_decodeRange(argv[++arg], lo, hi);

      minMemory = (uint64)ceil(lo * 1024.0 * 1024.0 * 1024.0);
      maxMemory = (uint64)ceil(hi * 1024.0 * 1024.0 * 1024.0);
      fileLimit = 0;

    } else if (strcmp(argv[arg], "-e") == 0) {
      maxError = atof(argv[++arg]);

    } else if (strcmp(argv[arg], "-e") == 0) {
      minOverlap = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-L") == 0) {
      AS_UTL_loadFileList(argv[++arg], fileList);

    } else if (strcmp(argv[arg], "-evalues") == 0) {
      eValues = true;

    } else if (strcmp(argv[arg], "-config") == 0) {
      configOut = argv[++arg];

    } else if (((argv[arg][0] == '-') && (argv[arg][1] == 0)) ||
               (AS_UTL_fileExists(argv[arg]))) {
      //  Assume it's an input file
      fileList.push_back(argv[arg]);

    } else {
      fprintf(stderr, "%s: unknown option '%s'.\n", argv[0], argv[arg]);
      err++;
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
  if (maxMemory < MEMORY_OVERHEAD)
    err++;
  if (err) {
    fprintf(stderr, "usage: %s -O asm.ovlStore -G asm.gkpStore [opts] [-L fileList | *.ovb.gz]\n", argv[0]);
    fprintf(stderr, "  -O asm.ovlStore       path to store to create\n");
    fprintf(stderr, "  -G asm.gkpStore       path to gkpStore for this assembly\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -L fileList           read input filenames from 'flieList'\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -F f                  use up to 'f' files for store creation\n");
    fprintf(stderr, "  -M g                  use up to 'g' gigabytes memory for sorting overlaps\n");
    fprintf(stderr, "                          default 4; g-0.25 gb is available for sorting overlaps\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -e e                  filter overlaps above e fraction error\n");
    fprintf(stderr, "  -l l                  filter overlaps below l bases overlap length (needs gkpStore to get read lengths!)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Non-building options:\n");
    fprintf(stderr, "  -evalues              input files are evalue updates from overlap error adjustment\n");
    fprintf(stderr, "  -config out.dat       don't build a store, just dump a binary partitioning file for ovStoreBucketizer\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Sizes and Limits:\n");
    fprintf(stderr, "  ovOverlap             " F_S32 " words of " F_S32 " bits each.\n", (int32)ovOverlapNWORDS, (int32)ovOverlapWORDSZ);
    fprintf(stderr, "  ovOverlapSortSize     " F_S32 " bits\n",      (int32)ovOverlapSortSize * 8);
    fprintf(stderr, "  SC_CHILD_MAX          " F_S32 " processes\n", (int32)sysconf(_SC_CHILD_MAX));
    fprintf(stderr, "  SC_OPEN_MAX           " F_S32 " files\n",     (int32)sysconf(_SC_OPEN_MAX));
    fprintf(stderr, "\n");

    if (ovlName == NULL)
      fprintf(stderr, "ERROR: No overlap store (-O) supplied.\n");
    if (gkpName == NULL)
      fprintf(stderr, "ERROR: No gatekeeper store (-G) supplied.\n");
    if (fileList.size() == 0)
      fprintf(stderr, "ERROR: No input overlap files (-L or last on the command line) supplied.\n");
    if (fileLimit > sysconf(_SC_OPEN_MAX) - 16)
      fprintf(stderr, "ERROR: Too many jobs (-F); only " F_SIZE_T " supported on this architecture.\n", sysconf(_SC_OPEN_MAX) - 16);
    if (maxMemory < MEMORY_OVERHEAD)
      fprintf(stderr, "ERROR: Memory (-M) must be at least %.3f GB to account for overhead.\n", MEMORY_OVERHEAD / 1024.0 / 1024.0 / 1024.0);

    exit(1);
  }

  //  If only updating evalues, do it and quit.

  if (eValues)
    addEvalues(ovlName, fileList), exit(0);

  //  Open reads, figure out a partitioning scheme.

  gkStore  *gkp         = gkStore::gkStore_open(gkpName);
  uint32    maxIID      = gkp->gkStore_getNumReads() + 1;
  uint32   *iidToBucket = computeIIDperBucket(fileLimit, minMemory, maxMemory, maxIID, fileList);

  uint32    maxFiles    = sysconf(_SC_OPEN_MAX);

  if (iidToBucket[maxIID-1] > maxFiles - 8) {
    fprintf(stderr, "ERROR:\n");
    fprintf(stderr, "ERROR:  Operating system limit of " F_U32 " open files.  The current -F/-M settings\n", maxFiles);
    fprintf(stderr, "ERROR:  will need to create " F_U32 " files to construct the store.\n", iidToBucket[maxIID-1]);
    fprintf(stderr, "ERROR:\n");
    exit(1);
  }

  //  But if only asked to report the configuration, do it and quit.

  if (configOut)
    reportConfiguration(configOut, maxIID, iidToBucket), gkp->gkStore_close(), exit(0);

  //  Read the gkStore to determine which fragments we care about.

  ovStoreFilter *filter = new ovStoreFilter(gkp, maxError);

  //

  fprintf(stderr, "\n");
  fprintf(stderr, "-- BUCKETIZING --\n");
  fprintf(stderr, "\n");

  //  And load reads into the store!  We used to create the store before filtering, so it could fail
  //  quicker, but the filter should be much faster with the mmap()'d gkpStore in canu.

  ovStoreWriter  *store   = new ovStoreWriter(ovlName, gkp);

  uint32          dumpFileMax  = iidToBucket[maxIID-1] + 1;
  ovFile        **dumpFile     = new ovFile * [dumpFileMax];
  uint64         *dumpLength   = new uint64   [dumpFileMax];

  memset(dumpFile,   0, sizeof(ovFile *) * dumpFileMax);
  memset(dumpLength, 0, sizeof(uint64)   * dumpFileMax);

  for (uint32 i=0; i<fileList.size(); i++) {
    ovOverlap    foverlap(gkp);
    ovOverlap    roverlap(gkp);

    fprintf(stderr, "-  Bucketizing '%s'\n", fileList[i]);

    ovFile *inputFile = new ovFile(gkp, fileList[i], ovFileFull);

    while (inputFile->readOverlap(&foverlap)) {
      filter->filterOverlap(foverlap, roverlap);  //  The filter copies f into r

      //  Check that overlap IDs are valid.
#warning not checking overlap IDs for validity

      //  If all are skipped, don't bother writing the overlap.

      if ((foverlap.dat.ovl.forUTG == true) ||
          (foverlap.dat.ovl.forOBT == true) ||
          (foverlap.dat.ovl.forDUP == true))
        writeToDumpFile(gkp, &foverlap, dumpFile, dumpLength, iidToBucket, ovlName);

      if ((roverlap.dat.ovl.forUTG == true) ||
          (roverlap.dat.ovl.forOBT == true) ||
          (roverlap.dat.ovl.forDUP == true))
        writeToDumpFile(gkp, &roverlap, dumpFile, dumpLength, iidToBucket, ovlName);
    }

    delete inputFile;
  }

  for (uint32 i=0; i<dumpFileMax; i++)
    delete dumpFile[i];

  //  Report the fate of filtering

  fprintf(stderr, "-  Bucketizing finished:\n");

  if (filter->savedDedupe() > 0) {
    fprintf(stderr, "-- Saved      " F_U64 " dedupe overlaps\n", filter->savedDedupe());
    fprintf(stderr, "-- Discarded  " F_U64 " don't care " F_U64 " different library " F_U64 " obviously not duplicates\n", filter->filteredNoDedupe(), filter->filteredNotDupe(), filter->filteredDiffLib());
  }

  if (filter->savedTrimming() > 0) {
    fprintf(stderr, "-- Saved      " F_U64 " trimming overlaps\n", filter->savedTrimming());
    fprintf(stderr, "-- Discarded  " F_U64 " don't care " F_U64 " too similar " F_U64 " too short\n", filter->filteredNoTrim(), filter->filteredBadTrim(), filter->filteredShortTrim());
  }

  if (filter->savedUnitigging() > 0) {
    fprintf(stderr, "-- Saved      " F_U64 " unitigging overlaps\n", filter->savedUnitigging());
  }

  if (filter->filteredErate() > 0)
    fprintf(stderr, "-- Discarded  " F_U64 " low quality, more than %.4f fraction error\n", filter->filteredErate(), maxError);

  delete filter;

  //
  //  Read each bucket, sort it, and dump it to the store
  //

  fprintf(stderr, "\n");
  fprintf(stderr, "-- SORTING --\n");
  fprintf(stderr, "\n");

  uint64 dumpLengthMax = 0;
  for (uint32 i=0; i<dumpFileMax; i++)
    if (dumpLengthMax < dumpLength[i])
      dumpLengthMax = dumpLength[i];

  ovStoreHistogram   *histogram = new ovStoreHistogram;

  ovOverlap  *overlapsort = ovOverlap::allocateOverlaps(gkp, dumpLengthMax);

  for (uint32 i=0; i<dumpFileMax; i++) {
    char      name[FILENAME_MAX];
    ovFile   *bof = NULL;

    if (dumpLength[i] == 0)
      continue;

    //  We're vastly more efficient if we skip the AS_OVS interface and just suck in the whole file
    //  directly....BUT....we can't do that because the AS_OVS interface is rearranging the data to
    //  make sure the store is cross-platform compatible.

    snprintf(name, FILENAME_MAX, "%s/tmp.sort.%03d", ovlName, i);
    fprintf(stderr, "-  Loading '%s'\n", name);

    bof = new ovFile(gkp, name, ovFileFull);

    uint64 numOvl = 0;
    while (bof->readOverlap(overlapsort + numOvl)) {

      //  Quick sanity check on IIDs.

      if ((overlapsort[numOvl].a_iid == 0) ||
          (overlapsort[numOvl].b_iid == 0) ||
          (overlapsort[numOvl].a_iid >= maxIID) ||
          (overlapsort[numOvl].b_iid >= maxIID)) {
        char ovlstr[256];

        fprintf(stderr, "Overlap has IDs out of range (maxIID " F_U32 "), possibly corrupt input data.\n", maxIID);
        fprintf(stderr, "  Aid " F_U32 "  Bid " F_U32 "\n",  overlapsort[numOvl].a_iid, overlapsort[numOvl].b_iid);
        exit(1);
      }

      numOvl++;
    }

    delete bof;

    assert(numOvl == dumpLength[i]);
    assert(numOvl <= dumpLengthMax);

    //  There's no real advantage to saving this file until after we write it out.  If we crash
    //  anywhere during the build, we are forced to restart from scratch.  I'll argue that removing
    //  it early helps us to not crash from running out of disk space.

    unlink(name);

    fprintf(stderr, "-  Sorting\n");

#ifdef _GLIBCXX_PARALLEL
    //  If we have the parallel STL, don't use it!  Sort is not inplace!
    __gnu_sequential::sort(overlapsort, overlapsort + dumpLength[i]);
#else
    sort(overlapsort, overlapsort + dumpLength[i]);
#endif

    fprintf(stderr, "-  Writing\n");

    for (uint64 x=0; x<dumpLength[i]; x++)
      store->writeOverlap(overlapsort + x);
  }

  fprintf(stderr, "\n");
  fprintf(stderr, "-- FINISHING --\n");
  fprintf(stderr, "\n");

  delete    store;
  delete [] overlapsort;

  gkp->gkStore_close();

  //  And we have a store.

  exit(0);
}
