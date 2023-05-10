
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

#include "system.H"
#include "strings.H"

#include "sqStore.H"
#include "ovStore.H"
#include "ovStoreConfig.H"

#include <vector>
#include <algorithm>



void
adjustMemoryLimits(uint64 maxPerRead,
                   uint64 olapsPerSliceMin, uint64 &minMemory,
                   uint64 olapsPerSliceMax, uint64 &maxMemory) {

  if ((maxPerRead <= olapsPerSliceMin) &&
      (maxPerRead <= olapsPerSliceMax))
    return;

  fprintf(stderr, "\n");
  fprintf(stderr, "WARNING:\n");

  if (olapsPerSliceMin < maxPerRead) {
    fprintf(stderr, "WARNING:  Increasing minimum memory to handle " F_U64 " overlaps per read.\n", maxPerRead);
    olapsPerSliceMin = maxPerRead;
    minMemory        = maxPerRead * ovOverlapSortSize + OVSTORE_MEMORY_OVERHEAD;
  }

  if (olapsPerSliceMax < maxPerRead) {
    fprintf(stderr, "WARNING:  Increasing maximum memory to handle " F_U64 " overlaps per read.\n", maxPerRead);
    olapsPerSliceMax = maxPerRead;
    maxMemory        = maxPerRead * ovOverlapSortSize + OVSTORE_MEMORY_OVERHEAD;
  }

  fprintf(stderr, "WARNING:\n");
  fprintf(stderr, "\n");
}


double
GBforOlaps(uint64 no) {
  return((OVSTORE_MEMORY_OVERHEAD + no * ovOverlapSortSize) / 1024.0 / 1024.0 / 1024.0);
}


void
ovStoreConfig::assignReadsToSlices(sqStore        *seq,
                                   uint64          minMemory,
                                   uint64          maxMemory) {

  //
  //  Load the number of overlaps per read.
  //

  uint64    *oPF         = new uint64 [_numInputs];   //  Overlaps per file.
  uint64    *oPR         = new uint64 [_maxID + 1];   //  Overlaps per read.
  uint64     numOverlaps = 0;                         //  Total overlaps.
  uint64     maxPerRead  = 0;

  memset(oPF, 0, sizeof(uint64) * (_numInputs));
  memset(oPR, 0, sizeof(uint64) * (_maxID + 1));

  for (uint32 ii=0; ii<_numInputs; ii++) {
    ovFile  *inputFile = new ovFile(seq, _inputNames[ii], ovFileFullCounts);

    for (uint32 rr=0; rr<_maxID + 1; rr++) {
      uint32 no = inputFile->getCounts()->numOverlaps(rr);
      oPF[ii]     += no;
      oPR[rr]     += no;   //  So much negativity here, sorry.
      numOverlaps += no;
    }

    delete inputFile;
  }

  for (uint64 ii=0; ii<_maxID+1; ii++)
    maxPerRead = std::max(maxPerRead, oPR[ii]);

  fprintf(stderr, "\n");
  fprintf(stderr, "   Moverlaps Input file\n");
  fprintf(stderr, "------------ ----------\n");

  for (uint32 ii=0; ii<_numInputs; ii++)
    fprintf(stderr, "%12.3f %s\n", oPF[ii] / 1000000.0 / 2, _inputNames[ii]);   //  Report number of canonical overlaps.

  fprintf(stderr, "------------\n");
  fprintf(stderr, "%12.3f Moverlaps in inputs\n", numOverlaps / 2 / 1000000.0);
  fprintf(stderr, "%12.3f Moverlaps to sort\n",   numOverlaps     / 1000000.0);
  fprintf(stderr, "\n");

  //
  //  Decide how many overlaps per slice we want to allow.  This directly
  //  affects how much memory we need for sorting.  The user has specified a
  //  minimum and maximum sort size, which could be the same.  We'll adjust
  //  this as needed to fit overlaps for extremely repetitive in one slice,
  //  but then pick a value 3/4 between min and max memory allowed.
  //

  uint64  olapsPerSliceMin   = (minMemory - OVSTORE_MEMORY_OVERHEAD) / ovOverlapSortSize;
  uint64  olapsPerSlice      = 0;
  uint64  olapsPerSliceExtra = 0;
  uint64  olapsPerSliceMax   = (maxMemory - OVSTORE_MEMORY_OVERHEAD) / ovOverlapSortSize;

  //  Reset the memory limits so that the maximum number of overlaps per read
  //  can be held in one slice.

  adjustMemoryLimits(maxPerRead,
                     olapsPerSliceMin, minMemory,
                     olapsPerSliceMax, maxMemory);

  //  Compute how many overlaps we can fit in the largest memory size allowed,
  //  then count how many slices we'd need to do that.
  //
  //  Now knowing how many slices we need, recompute olapsPerSlice to make
  //  jobs more equal (but never smaller than the max number of overlaps per
  //  read).
  //
  //  But if thise size is close to maxMemory, add another (and another and
  //  another) slice until it is below the max.
  //
  //  With that, compute the final (estimated) memory size for sorting jobs,
  //  and a limit on how many additional overlaps we can allow in a slice
  //  before exceeding the maximum limit.

  olapsPerSlice = (maxMemory - OVSTORE_MEMORY_OVERHEAD) / ovOverlapSortSize;

  _numSlices = 1;

  for (uint64 olaps=0, ii=0; ii<_maxID+1; ii++) {
    if (olaps + oPR[ii] > olapsPerSlice) {
      olaps = 0;
      _numSlices++;
    }

    olaps += oPR[ii];
  }

  olapsPerSlice = (uint64)ceil((double)numOverlaps / (double)_numSlices) + 1;

  while (GBforOlaps(olapsPerSlice) + 1.0 > maxMemory)
    _numSlices++;

  if (olapsPerSlice < maxPerRead)   //  Slices cannot be smaller than the maximum
    olapsPerSlice = maxPerRead;     //  number of overlaps per read.

  _sortMemory = (olapsPerSlice * ovOverlapSortSize + OVSTORE_MEMORY_OVERHEAD) / 1024.0 / 1024.0 / 1024.0;

  olapsPerSliceExtra = ((maxMemory - _sortMemory) * 1024.0 * 1024.0 * 1024.0) / ovOverlapSortSize;

  //  Greedily assign inputs to each bucketizer task.
  //
  //  The number of buckets to make is essentially a free parameter - lower
  //  makes bigger buckets and fewer files.

  _numBuckets = std::min(_numInputs, _numSlices);

  uint64  *olapsPerBucket = new uint64 [_numBuckets];
  uint32  *inputPerBucket = new uint32 [_numBuckets];

  for (uint32 ii=0; ii<_numBuckets; ii++) {
    olapsPerBucket[ii] = 0;
    inputPerBucket[ii] = 0;
  }

  for (uint32 ss=0, ii=0; ii<_numInputs; ii++) {
    uint32  mb = 0;

    for (uint32 bb=0; bb<_numBuckets; bb++)          //  Find the bucket with the
      if (olapsPerBucket[bb] < olapsPerBucket[mb])   //  fewest number of overlaps
        mb = bb;                                     //  assigned to it.

    _inputToBucket[ii]  = mb;                        //  And add this input to
    olapsPerBucket[mb] += oPF[ii];                   //  that bucket.
    inputPerBucket[mb] += 1;
  }

  delete [] oPF;

  //  Report ovb to bucket mapping.

  {
    fprintf(stderr, "\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "         number    number of\n");
    fprintf(stderr, " bucket  inputs     overlaps\n");
    fprintf(stderr, "------- ------- ------------\n");

    uint64 totOlaps = 0;
    for (uint32 ii=0; ii<_numBuckets; ii++) {
      fprintf(stderr, "%7" F_U32P " %7" F_U32P " %12" F_U64P "\n",
              ii, inputPerBucket[ii], olapsPerBucket[ii]);
      totOlaps += olapsPerBucket[ii];
    }

    fprintf(stderr, "------- ------- ------------\n");
    fprintf(stderr, "                %12" F_U64P "\n", totOlaps);
    fprintf(stderr, "\n");

    assert(totOlaps == numOverlaps);
  }

  delete [] olapsPerBucket;
  delete [] inputPerBucket;

  //
  //  Assign reads to slices.
  //

  fprintf(stderr, "\n");
  fprintf(stderr, "          number of    expected\n");
  fprintf(stderr, " slice     overlaps memory (GB)       read range\n");
  fprintf(stderr, "------ ------------ ----------- ---------------------\n");

  {
    uint32  first    = 1;
    uint64  olaps    = 0;
    uint32  slice    = 0;
    uint64  totOlaps = 0;

    for (uint32 ii=0; ii<_maxID+1; ii++) {
      //  All overlaps fit easily into this slice.
      if (olaps + oPR[ii] < olapsPerSlice) {
        _readToSlice[ii]  =   slice;
        olaps            += oPR[ii];
      }

      //  All overlaps fit mostly in this slice; add them to it and move to the next slice.
      else if (olaps + oPR[ii] < olapsPerSlice + olapsPerSliceExtra) {
        fprintf(stderr, "%6u %12lu %11.3f %10u-%-10u\n", slice, olaps + oPR[ii], GBforOlaps(olaps + oPR[ii]), first, ii);
        _readToSlice[ii]  = slice++;
        first             = ii+1;
        olaps             = 0;
      }

      //  All overlaps do not fit at all in this slice; add them to the next slice.
      else {
        fprintf(stderr, "%6u %12lu %11.3f %10u-%-10u\n", slice, olaps, GBforOlaps(olaps), first, ii-1);
        _readToSlice[ii]  = ++slice;
        first             = ii;
        olaps             = oPR[ii];
      }

      totOlaps += oPR[ii];
    }

    if (olaps > 0)
      fprintf(stderr, "%6u %12lu %11.3f %10u-%-10u\n", slice, olaps, GBforOlaps(olaps), first, _maxID);

    fprintf(stderr, "------ ------------ ----------- ---------------------\n");
    fprintf(stderr, "       %12" F_U64P "\n", totOlaps);
    fprintf(stderr, "\n");

    assert(totOlaps == numOverlaps);
  }

  delete [] oPR;
}




int
main(int argc, char **argv) {
  char                      *seqName         = NULL;
  uint64                     minMemory       = (uint64)1 * 1024 * 1024 * 1024;
  uint64                     maxMemory       = (uint64)4 * 1024 * 1024 * 1024;

  stringList                 fileList;

  char const                *configOut       = NULL;
  char const                *configIn        = NULL;

  bool                       writeNumBuckets = false;
  bool                       writeNumSlices  = false;
  bool                       writeMemory     = false;
  uint32                     writeInputs     = 0;
  uint32                     writeSlices     = 0;

  argc = AS_configure(argc, argv, 1);

  std::vector<char const *>  err;
  for (int32 arg=1; arg < argc; arg++) {
    if        (strcmp(argv[arg], "-S") == 0) {
      seqName = argv[++arg];

    } else if (strcmp(argv[arg], "-M") == 0) {
      double lo=0.0, hi=0.0;

      decodeRange<double>(argv[++arg], lo, hi);

      minMemory = (uint64)ceil(lo * 1024.0 * 1024.0 * 1024.0);
      maxMemory = (uint64)ceil(hi * 1024.0 * 1024.0 * 1024.0);

    } else if (strcmp(argv[arg], "-L") == 0) {
      fileList.load(argv[++arg]);

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
      fileList.add(argv[arg]);

    } else {
      char *s = new char [1024];
      snprintf(s, 1024, "%s: unknown option '%s'.\n", argv[0], argv[arg]);
      err.push_back(s);
    }
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
    sqStore        *seq    = new sqStore(seqName);
    uint32          maxID  = seq->sqStore_lastReadID();

    config = new ovStoreConfig(fileList.getVector(), maxID);

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
