
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

#include "runtime.H"

#include "strings.H"

#include "sqStore.H"
#include "ovStore.H"
#include "ovStoreConfig.H"

#include <vector>
#include <algorithm>

using namespace std;





void
ovStoreConfig::assignReadsToSlices(sqStore        *seq,
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
  fprintf(stderr, "   Moverlaps\n");
  fprintf(stderr, "------------ ----------------------------------------\n");

  uint64             *oPF                = new uint64 [_numInputs];
  uint32             *oPR                = new uint32 [_maxID + 1];
  uint64              numOverlaps        = 0;

  memset(oPF, 0, sizeof(uint64) * (_numInputs));
  memset(oPR, 0, sizeof(uint32) * (_maxID + 1));

  for (uint32 ii=0; ii<_numInputs; ii++) {
    ovFile            *inputFile = new ovFile(seq, _inputNames[ii], ovFileFullCounts);

    for (uint32 rr=0; rr<_maxID + 1; rr++) {
      oPF[ii] += inputFile->getCounts()->numOverlaps(rr) / 2;   //  Reports counts as if they were
      oPR[rr] += inputFile->getCounts()->numOverlaps(rr);       //  already symmetrized.
    }

    numOverlaps += oPF[ii] * 2;

    delete inputFile;

    fprintf(stderr, "%12.3f %40s\n", oPF[ii] / 1000000.0, _inputNames[ii]);
  }

  fprintf(stderr, "------------ ----------------------------------------\n");
  fprintf(stderr, "%12.3f Moverlaps in inputs\n", numOverlaps / 2 / 1000000.0);
  fprintf(stderr, "%12.3f Moverlaps to sort\n",   numOverlaps     / 1000000.0);
  fprintf(stderr, "\n");

  if (numOverlaps == 0)
    fprintf(stderr, "Found no overlaps to sort.\n");


  //
  //  Partition the overlaps into buckets.
  //
  //  This will pick the smallest memory size that uses fewer than maxFiles buckets.  Unreasonable
  //  values can break this - either too low memory or too high allowed open files (an OS limit).
  //

  uint64  olapsPerSliceMin = (minMemory - OVSTORE_MEMORY_OVERHEAD) / ovOverlapSortSize;
  uint64  olapsPerSliceMax = (maxMemory - OVSTORE_MEMORY_OVERHEAD) / ovOverlapSortSize;

  //  Reset the limits so that the maximum number of overlaps per read can be held in one slice.

  uint64  maxOverlapsPerRead = 0;

  for (uint64 olaps=0, ii=0; ii<_maxID+1; ii++)
    if (maxOverlapsPerRead < oPR[ii])
      maxOverlapsPerRead = oPR[ii];

  if ((olapsPerSliceMin < maxOverlapsPerRead) ||
      (olapsPerSliceMax < maxOverlapsPerRead)) {
    fprintf(stderr, "WARNING:\n");

    if (olapsPerSliceMin < maxOverlapsPerRead) {
      fprintf(stderr, "WARNING:  Increasing minimum memory to handle " F_U64 " overlaps per read.\n", maxOverlapsPerRead);
      olapsPerSliceMin = maxOverlapsPerRead;
      minMemory        = maxOverlapsPerRead * ovOverlapSortSize + OVSTORE_MEMORY_OVERHEAD;
    }

    if (olapsPerSliceMax < maxOverlapsPerRead) {
      fprintf(stderr, "WARNING:  Increasing maximum memory to handle " F_U64 " overlaps per read.\n", maxOverlapsPerRead);
      olapsPerSliceMax = maxOverlapsPerRead;
      maxMemory        = maxOverlapsPerRead * ovOverlapSortSize + OVSTORE_MEMORY_OVERHEAD;
    }

    fprintf(stderr, "WARNING:\n");
    fprintf(stderr, "\n");
  }

  fprintf(stderr, "Configuring for:\n");
  fprintf(stderr, "  Up to " F_S64 " processes.\n", procMax);
  fprintf(stderr, "  Up to " F_S64 " open files.\n", openMax);

  if (minMemory < maxMemory)
    fprintf(stderr, "  At least %5.2f GB memory -> %4" F_U64P " sort processes with %.2f Moverlaps each.\n",
            minMemory / 1024.0 / 1024.0 / 1024.0,
            numOverlaps / olapsPerSliceMin + 1,
            olapsPerSliceMin / 1000000.0);

  fprintf(stderr, "  At most  %5.2f GB memory -> %4" F_U64P " sort processes with %.2f Moverlaps each.\n",
          maxMemory / 1024.0 / 1024.0 / 1024.0,
          numOverlaps / olapsPerSliceMax + 1,
          olapsPerSliceMax / 1000000.0);
  fprintf(stderr, "\n");



  //  More processes -> higher bandwidth utilization, but can also thrash the FS
  //  More memory    -> harder to get machines to run on
  //
  //  So arbitrarily pick a memory size 75% of the maximum.

  uint64  sortMemory       = minMemory + 3 * (maxMemory - minMemory) / 4;

  uint64  olapsPerSlice    = (sortMemory - OVSTORE_MEMORY_OVERHEAD) / ovOverlapSortSize;

  //  With that upper limit on the number of overlaps per slice, count how many slices
  //  we need to make.

  _numSlices = 1;

  uint64 total = 0;

  for (uint64 olaps=0, ii=0; ii<_maxID+1; ii++) {
    if (olaps + oPR[ii] > olapsPerSlice) {
      olaps = 0;
      _numSlices++;
    }

    olaps += oPR[ii];
    total += oPR[ii];
  }

  //  Divide those overlaps evenly among the slices, but no smaller than
  //  our minimum (oddly called 'maxOverlapsPerRead').

  olapsPerSlice = (uint64)ceil((double)numOverlaps / (double)_numSlices) + 1;

  if (olapsPerSlice < maxOverlapsPerRead)
    olapsPerSlice = maxOverlapsPerRead;

  _sortMemory = (olapsPerSlice * ovOverlapSortSize + OVSTORE_MEMORY_OVERHEAD) / 1024.0 / 1024.0 / 1024.0;

  //  One more time, just to count the number of slices we're making.

  _numSlices = 1;

  for (uint64 olaps=0, ii=0; ii<_maxID+1; ii++) {
    if (olaps + oPR[ii] > olapsPerSlice) {
      olaps = 0;
      _numSlices++;
    }

    olaps += oPR[ii];
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


  //  Report ovb to bucket mapping.

  fprintf(stderr, "------------------------------------------------------------\n");
  fprintf(stderr, "Will bucketize using " F_U32 " processes.\n", _numBuckets);
  fprintf(stderr, "\n");
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

  fprintf(stderr, "\n");
  fprintf(stderr, "------------------------------------------------------------\n");
  fprintf(stderr, "Will sort using " F_U32 " processes.\n",  _numSlices);
  fprintf(stderr, "  Up to %7.2f M overlaps\n", olapsPerSlice / 1000000.0);
  fprintf(stderr, "        %7.2f GB memory\n", _sortMemory);
  fprintf(stderr, "\n");
  fprintf(stderr, "          number of\n");
  fprintf(stderr, " slice     overlaps       read range\n");
  fprintf(stderr, "------ ------------ ---------------------\n");

  totOlaps = 0;

  {
    uint32  first = 1;
    uint64  olaps = 0;
    uint32  slice = 0;

    for (uint32 ii=0; ii<_maxID+1; ii++) {
      if (olaps + oPR[ii] > olapsPerSlice) {
        fprintf(stderr, "%6" F_U32P " %12" F_U64P " %10" F_U32P "-%-10" F_U32P "\n", slice, olaps, first, ii-1);
        totOlaps += olaps;
        olaps = 0;
        slice++;
        first = ii;
      }

      olaps            += oPR[ii];
      _readToSlice[ii]  = slice;
    }

    fprintf(stderr, "%6" F_U32P " %12" F_U64P " %10" F_U32P "-%-10" F_U32P "\n", slice, olaps, first, _maxID);
    totOlaps += olaps;

    if (slice + 1 != _numSlices)
      fprintf(stderr, "ERROR: expected %u slices, generated %u instead.  %lu overlaps per slice allowed.\n",
              _numSlices, slice, olapsPerSlice);
    assert(slice + 1 == _numSlices);
  }

  fprintf(stderr, "------ ------------\n");
  fprintf(stderr, "       %12" F_U64P "\n", totOlaps);

  fprintf(stderr, "\n");

  delete [] oPR;
}




int
main(int argc, char **argv) {
  char                 *seqName         = NULL;
  uint64                minMemory       = (uint64)1 * 1024 * 1024 * 1024;
  uint64                maxMemory       = (uint64)4 * 1024 * 1024 * 1024;

  vector<char const *>  fileList;

  char const           *configOut       = NULL;
  char const           *configIn        = NULL;

  bool                  writeNumBuckets = false;
  bool                  writeNumSlices  = false;
  bool                  writeMemory     = false;
  uint32                writeInputs     = 0;
  uint32                writeSlices     = 0;

  argc = AS_configure(argc, argv);

  vector<char const *>  err;
  int                   arg=1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-S") == 0) {
      seqName = argv[++arg];

    } else if (strcmp(argv[arg], "-M") == 0) {
      double lo=0.0, hi=0.0;

      decodeRange(argv[++arg], lo, hi);

      minMemory = (uint64)ceil(lo * 1024.0 * 1024.0 * 1024.0);
      maxMemory = (uint64)ceil(hi * 1024.0 * 1024.0 * 1024.0);

    } else if (strcmp(argv[arg], "-L") == 0) {
      AS_UTL_loadFileList(argv[++arg], fileList);

    } else if (strcmp(argv[arg], "-create") == 0) {
      configOut = argv[++arg];

    } else if (strcmp(argv[arg], "-describe") == 0) {
      configIn = argv[++arg];

    } else if (strcmp(argv[arg], "-numbuckets") == 0) {
      writeNumBuckets = true;
    } else if (strcmp(argv[arg], "-numslices") == 0) {
      writeNumSlices = true;
    } else if (strcmp(argv[arg], "-sortmemory") == 0) {
      writeMemory = true;
    } else if (strcmp(argv[arg], "-listinputs") == 0) {
      writeInputs = strtouint32(argv[++arg]);
    } else if (strcmp(argv[arg], "-listslices") == 0) {
      writeSlices = strtouint32(argv[++arg]);

    } else if (((argv[arg][0] == '-') && (argv[arg][1] == 0)) ||
               (fileExists(argv[arg]))) {
      //  Assume it's an input file
      fileList.push_back(argv[arg]);

    } else {
      char *s = new char [1024];
      snprintf(s, 1024, "%s: unknown option '%s'.\n", argv[0], argv[arg]);
      err.push_back(s);
    }

    arg++;
  }

  if ((seqName == NULL) && (configIn == NULL))
    err.push_back("ERROR: No sequence store (-S) supplied.\n");

  if ((fileList.size() == 0) && (configIn == NULL))
    err.push_back("ERROR: No input overlap files (-L or last on the command line) supplied.\n");

  if ((configOut != NULL) && (configIn != NULL))
    err.push_back("ERROR: Can't both -create -describe a config.\n");

  if ((configOut == NULL) && (configIn == NULL))
    err.push_back("ERROR: Must supply one of -create or -describe.\n");

  if ((minMemory <= OVSTORE_MEMORY_OVERHEAD) ||
      (maxMemory <= OVSTORE_MEMORY_OVERHEAD + ovOverlapSortSize))
    err.push_back("ERROR: Memory (-M) must be at least 0.25 GB to account for overhead.\n");  //  , OVSTORE_MEMORY_OVERHEAD / 1024.0 / 1024.0 / 1024.0

  if (err.size() > 0) {
    fprintf(stderr, "usage: %s -S asm.seqStore -create out.config [opts] [-L fileList | *.ovb]\n", argv[0]);
    fprintf(stderr, "  -S asm.seqStore       path to seqStore for this assembly\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -L fileList           a list of ovb files in 'fileList'\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -M g                  use up to 'g' gigabytes memory for sorting overlaps\n");
    fprintf(stderr, "                          default 4; g-0.25 gb is available for sorting overlaps\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -create config        write overlap store configuration to file 'config'\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -describe config      write a readable description of the config in 'config' to the screen\n");
    fprintf(stderr, "  -numbuckets           write the number of buckets to the screen\n");
    fprintf(stderr, "  -numslices            write the number of slices to the screen\n");
    fprintf(stderr, "  -sortmemory           write the memory needed (in GB) for a sort job to the screen\n");
    fprintf(stderr, "  -listinputs n         write a list of the input ovb files needed for bucketizer job 'n'");
    fprintf(stderr, "  -listslices n         write a list of the input slice files needed for sorter job 'n'\n");
    fprintf(stderr, "\n");
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

    sqStore        *seq    = new sqStore(seqName);
    uint32          maxID  = seq->sqStore_lastReadID();

    config = new ovStoreConfig(fileList, maxID);

    config->assignReadsToSlices(seq, minMemory, maxMemory);
    config->writeConfig(configOut);

    delete seq;
  }

  //  If we have a config, report parameters.

  if (config) {
    uint32  memGB = (uint32)ceil(config->sortMemory() + 0.5);  //  Adds an extra 0.5 to 1.49 gb.

    if (writeNumBuckets) {
      fprintf(stdout, F_U32 "\n", config->numBuckets());
    }

    else if (writeNumSlices) {
      fprintf(stdout, F_U32 "\n", config->numSlices());
    }

    else if (writeMemory) {
      fprintf(stdout, F_U32 "\n", memGB);
    }

    else if (writeInputs) {
      for (uint32 ff=0; ff<config->numInputs(writeInputs); ff++)
        fprintf(stdout, "%s\n", config->getInput(writeInputs, ff));
    }

    else if (writeSlices) {
      for (uint32 bb=1; bb<=config->numBuckets(); bb++) {
        fprintf(stdout, "bucket%04" F_U32P "/slice%04" F_U32P "\n", bb, writeSlices);
        fprintf(stdout, "bucket%04" F_U32P "/sliceSizes\n", bb);
      }
    }

    else {
      fprintf(stdout, "\n");
      fprintf(stdout, "Configured for:\n");
      fprintf(stdout, "  numBuckets %8" F_U32P "\n", config->numBuckets());
      fprintf(stdout, "  numSlices  %8" F_U32P "\n", config->numSlices());
      fprintf(stdout, "  sortMemory %8" F_U32P " GB (%5.3f GB)\n", memGB, config->sortMemory());
    }
  }

  //  Cleanup and quit!

  delete config;

  exit(0);
}
