
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
 *  Modifications by:
 *
 *    Brian P. Walenz beginning on 2018-MAR-15
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
#include "ovStoreConfig.H"

#include <vector>
#include <algorithm>

using namespace std;





void
ovStoreConfig::assignReadsToSlices(gkStore        *gkp,
                                   uint64          minMemory,
                                   uint64          maxMemory) {

  int64    procMax       = sysconf(_SC_CHILD_MAX);
  int64    openMax       = sysconf(_SC_OPEN_MAX) - 16;

  //
  //  Load the number of overlaps per read.
  //

  fprintf(stderr, "\n");
  fprintf(stderr, "Finding number of overlaps per read and per file.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "      Molaps\n");
  fprintf(stderr, "------------ ----------------------------------------\n");

  uint64             *oPF         = new uint64 [_numInputs];
  uint32             *oPR         = new uint32 [_maxID + 1];
  uint64              numOverlaps = 0;

  memset(oPR, 0, sizeof(uint32) * (_maxID + 1));

  for (uint32 ii=0; ii<_numInputs; ii++) {
    ovFile            *inputFile = new ovFile(gkp, _inputNames[ii], ovFileFull);

    oPF[ii]      = inputFile->getCounts()->numOverlaps();       //  Used for load balancing.
    numOverlaps += inputFile->getCounts()->numOverlaps() * 2;   //  Because we symmetrize!

    for (uint32 rr=0; rr<_maxID + 1; rr++)
      oPR[rr] += inputFile->getCounts()->numOverlaps(rr);

    delete inputFile;

    fprintf(stderr, "%12.3f %40s\n", oPF[ii] / 1000000.0, _inputNames[ii]);
  }

  fprintf(stderr, "------------ ----------------------------------------\n");
  fprintf(stderr, "%12.3f overlaps in inputs\n", numOverlaps / 2 / 1000000.0);
  fprintf(stderr, "%12.3f overlaps to sort\n",   numOverlaps     / 1000000.0);
  fprintf(stderr, "\n");

  if (numOverlaps == 0)
    fprintf(stderr, "Found no overlaps to sort.\n"), exit(1);


  //
  //  Partition the overlaps into buckets.
  //
  //  This will pick the smallest memory size that uses fewer than maxFiles buckets.  Unreasonable
  //  values can break this - either too low memory or too high allowed open files (an OS limit).
  //

  uint64  olapsPerSliceMin = (minMemory - OVSTORE_MEMORY_OVERHEAD) / ovOverlapSortSize;
  uint64  olapsPerSliceMax = (maxMemory - OVSTORE_MEMORY_OVERHEAD) / ovOverlapSortSize;

  fprintf(stderr, "Configuring for %.2f GB to %.2f GB memory.  Not to exceed " F_S64 " processes or " F_S64 " open files.\n",
          minMemory / 1024.0 / 1024.0 / 1024.0,
          maxMemory / 1024.0 / 1024.0 / 1024.0,
          procMax, openMax);
  fprintf(stderr, "  At minimum memory, %4" F_U64P " sort processes with %.2f Molaps each.\n",
          numOverlaps / olapsPerSliceMin + 1,
          olapsPerSliceMin / 1000000.0);
  fprintf(stderr, "  At maximum memory, %4" F_U64P " sort processes with %.2f Molaps each.\n",
          numOverlaps / olapsPerSliceMax + 1,
          olapsPerSliceMax / 1000000.0);
  fprintf(stderr, "\n");

  //  More processes -> higher bandwidth utilization, but can also thrash the FS
  //  More memory    -> harder to get machines to run on
  //
  //  So arbitrarily pick a memory size 75% of the maximum.

  uint64  sortMemory       = minMemory + 3 * (maxMemory - minMemory) / 4;

  uint64  olapsPerSlice    = (sortMemory - OVSTORE_MEMORY_OVERHEAD) / ovOverlapSortSize;
  uint64  maxOlapsPerSlice = 0;

  //  With that upper limit on the number of overlaps per slice, count how many slices
  //  we need to make.

  _numSlices = 1;

  for (uint64 olaps=0, ii=0; ii<_maxID+1; ii++) {
    olaps += oPR[ii];

    if (olaps >= olapsPerSlice) {
      olaps = 0;
      _numSlices++;
    }
  }

  //  Now, divide those overlaps evenly among the slices and find the maximum number of overlaps in any slice.

  olapsPerSlice = (uint64)ceil((double)numOverlaps / (double)_numSlices) + 1;

  for (uint64 olaps=0, ii=0; ii<_maxID+1; ii++) {
    olaps += oPR[ii];

    if (maxOlapsPerSlice < olaps)
      maxOlapsPerSlice = olaps;

    if (olaps >= olapsPerSlice)
      olaps = 0;
  }

  //  Assign inputs to each bucketizer.  Greedy load balancing.

  //  Essentially a free parameter - lower makes bigger buckets and fewer files.
  _numBuckets = min(_numInputs, _numSlices);

  uint64  *olapsPerBucket = new uint64 [_numBuckets];

  for (uint32 ii=0; ii<_numBuckets; ii++)
    olapsPerBucket[ii] = 0;

  for (uint32 ss=0, ii=0; ii<_numInputs; ii++) {
    uint32  mb = 0;

    for (uint32 bb=0; bb<_numBuckets; bb++)
      if (olapsPerBucket[bb] < olapsPerBucket[mb])
        mb = bb;

    _inputToBucket[ii]  = mb;
    olapsPerBucket[mb] += oPF[ii] * 2;
  }

  delete [] oPF;

  //  Report results that nobody will read.

  _sortMemory = (maxOlapsPerSlice * ovOverlapSortSize + OVSTORE_MEMORY_OVERHEAD) / 1024.0 / 1024.0 / 1024.0;

  fprintf(stderr, "Will bucketize using " F_U32 " processes.\n", _numBuckets);
  fprintf(stderr, "\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Will sort using " F_U32 " processes.\n",  _numSlices);
  fprintf(stderr, "  Up to %7.2f M overlaps\n", maxOlapsPerSlice / 1000000.0);
  fprintf(stderr, "        %7.2f GB memory\n", _sortMemory);
  fprintf(stderr, "\n");


  //  Report ovb to bucket mapping.

  fprintf(stderr, "         number    number of\n");
  fprintf(stderr, " bucket  inputs     overlaps\n");
  fprintf(stderr, "------- ------- ------------\n");

  uint64  totOlaps = 0;

  for (uint32 ii=0; ii<_numBuckets; ii++) {
    uint32  ni = 0;

    for (uint32 xx=0; xx<_numInputs; xx++) {
      if (_inputToBucket[xx] == ii)
        ni++;
    }

    fprintf(stderr, "%7" F_U32P " %7" F_U32P " %12" F_U64P "\n",
            ii, ni, olapsPerBucket[ii]);

    totOlaps += olapsPerBucket[ii];
  }

  delete [] olapsPerBucket;

  fprintf(stderr, "------- ------- ------------\n");
  fprintf(stderr, "                %12" F_U64P "\n", totOlaps);
  fprintf(stderr, "\n");

  //  Assign reads to slices.

  fprintf(stderr, "          number of\n");
  fprintf(stderr, " slice     overlaps\n");
  fprintf(stderr, "------ ------------\n");

  totOlaps = 0;

  {
    uint64  olaps = 0;
    uint32  slice = 0;

    for (uint32 ii=0; ii<_maxID+1; ii++) {
      olaps            += oPR[ii];
      _readToSlice[ii]   = slice;

      if (olaps >= olapsPerSlice) {
        fprintf(stderr, "%6" F_U32P " %12" F_U64P "\n", slice, olaps);
        totOlaps += olaps;
        olaps = 0;
        slice++;
      }
    }

    fprintf(stderr, "%6" F_U32P " %12" F_U64P "\n", slice, olaps);
    totOlaps += olaps;
    olaps = 0;
    slice++;
  }

  fprintf(stderr, "------ ------------\n");
  fprintf(stderr, "       %12" F_U64P "\n", totOlaps);

  fprintf(stderr, "\n");

  delete [] oPR;
}




int
main(int argc, char **argv) {
  char           *gkpName    = NULL;
  uint64          minMemory  = (uint64)1 * 1024 * 1024 * 1024;
  uint64          maxMemory  = (uint64)4 * 1024 * 1024 * 1024;

  vector<char *>  fileList;

  char           *configOut  = NULL;
  char           *configIn   = NULL;

  argc = AS_configure(argc, argv);

  vector<char *>  err;
  int             arg=1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-G") == 0) {
      gkpName = argv[++arg];

    } else if (strcmp(argv[arg], "-M") == 0) {
      double lo=0.0, hi=0.0;

      AS_UTL_decodeRange(argv[++arg], lo, hi);

      minMemory = (uint64)ceil(lo * 1024.0 * 1024.0 * 1024.0);
      maxMemory = (uint64)ceil(hi * 1024.0 * 1024.0 * 1024.0);

    } else if (strcmp(argv[arg], "-L") == 0) {
      AS_UTL_loadFileList(argv[++arg], fileList);

    } else if (strcmp(argv[arg], "-create") == 0) {
      configOut = argv[++arg];

    } else if (strcmp(argv[arg], "-describe") == 0) {
      configIn = argv[++arg];

    } else if (((argv[arg][0] == '-') && (argv[arg][1] == 0)) ||
               (AS_UTL_fileExists(argv[arg]))) {
      //  Assume it's an input file
      fileList.push_back(argv[arg]);

    } else {
      char *s = new char [1024];
      snprintf(s, 1024, "%s: unknown option '%s'.\n", argv[0], argv[arg]);
      err.push_back(s);
    }

    arg++;
  }

  if ((gkpName == NULL) && (configIn == NULL))
    err.push_back("ERROR: No gatekeeper store (-G) supplied.\n");

  if ((fileList.size() == 0) && (configIn == NULL))
    err.push_back("ERROR: No input overlap files (-L or last on the command line) supplied.\n");

  if ((configOut != NULL) && (configIn != NULL))
    err.push_back("ERROR: Can't both -create -describe a config.\n");

  if ((configOut == NULL) && (configIn == NULL))
    err.push_back("ERROR: Must supply one of -create or -describe.\n");

  if (maxMemory < OVSTORE_MEMORY_OVERHEAD + ovOverlapSortSize)
    fprintf(stderr, "ERROR: Memory (-M) must be at least 0.25 GB to account for overhead.\n");  //  , OVSTORE_MEMORY_OVERHEAD / 1024.0 / 1024.0 / 1024.0

  if (err.size() > 0) {
    fprintf(stderr, "usage: %s -G asm.gkpStore -create out.config [opts] [-L fileList | *.ovb]\n", argv[0]);
    fprintf(stderr, "  -G asm.gkpStore       path to gkpStore for this assembly\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -L fileList           a list of ovb files in 'fileList'\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -M g                  use up to 'g' gigabytes memory for sorting overlaps\n");
    fprintf(stderr, "                          default 4; g-0.25 gb is available for sorting overlaps\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -create config        write overlap store configuration to file 'config'\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -describe config      write a description of the config in 'config' to the screen\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Sizes and Limits:\n");
    fprintf(stderr, "  ovOverlap             " F_S32 " words of " F_S32 " bits each.\n", (int32)ovOverlapNWORDS, (int32)ovOverlapWORDSZ);
    fprintf(stderr, "  ovOverlapSortSize     " F_S32 " bits\n",      (int32)ovOverlapSortSize * 8);
    fprintf(stderr, "  SC_CHILD_MAX          " F_S32 " processes\n", (int32)sysconf(_SC_CHILD_MAX));
    fprintf(stderr, "  SC_OPEN_MAX           " F_S32 " files\n",     (int32)sysconf(_SC_OPEN_MAX));
    fprintf(stderr, "\n");

    for (uint32 ii=0; ii<err.size(); ii++)
      if (err[ii])
        fputs(err[ii], stderr);

    exit(1);
  }

  ovStoreConfig  *config = NULL;

  //  If describing, load and describe.

  if (configIn) {
    config = new ovStoreConfig(configIn);
  }

  //  Check parameters, reset some of them.

  else {
    if (minMemory < OVSTORE_MEMORY_OVERHEAD + ovOverlapSortSize) {
      fprintf(stderr, "Reset minMemory from " F_U64 " to " F_SIZE_T "\n", minMemory, OVSTORE_MEMORY_OVERHEAD + ovOverlapSortSize);
      minMemory  = OVSTORE_MEMORY_OVERHEAD + ovOverlapSortSize;
    }

    gkStore        *gkp    = gkStore::gkStore_open(gkpName);
    uint32          maxID  = gkp->gkStore_getNumReads();

    config = new ovStoreConfig(fileList, maxID);

    config->assignReadsToSlices(gkp, minMemory, maxMemory);
    config->writeConfig(configOut);
 
    gkp->gkStore_close();
  }

  //  If we have a config, report parameters.

  if (config) {
    fprintf(stdout, "\n");
    fprintf(stdout, "Configured for:\n");
    fprintf(stdout, "  numBuckets %8" F_U32P "\n", config->numBuckets());
    fprintf(stdout, "  numSlices  %8" F_U32P "\n", config->numSlices());
    fprintf(stdout, "  sortMemory %8" F_U32P " GB (%5.3f GB)\n",
            (uint32)(config->sortMemory()) + 1,  //  Adds an extra 0.5 to 1.4 gb.
            config->sortMemory());
  }

  delete config;

  exit(0);
}
