
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
 *    Brian P. Walenz beginning on 2015-NOV-03
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



class ovStoreConfig {
public:
  ovStoreConfig() {
    _maxID         = 0;
    _numInputs     = 0;
    _inputNames    = NULL;
    _inputToBucket = NULL;
    _readToSlice   = NULL;
  };

  ovStoreConfig(vector<char *> &names, uint32 maxID) {
    _maxID         = maxID;

    _numInputs  = names.size();
    _inputNames = new char * [_numInputs];

    for (uint32 ii=0; ii<names.size(); ii++)
      _inputNames[ii] = duplicateString(names[ii]);

    _inputToBucket = new uint32 [_numInputs];
    _readToSlice   = new uint16 [_maxID+1];
  };

  ~ovStoreConfig() {
    for (uint32 ii=0; ii<_numInputs; ii++)
      delete [] _inputNames[ii];

    delete [] _inputNames;

    delete [] _inputToBucket;
    delete [] _readToSlice;
  };

  void    loadConfig(const char *configName) {
    FILE *C = AS_UTL_openInputFile(configName);

    AS_UTL_safeRead(C, &_maxID,     "maxID",     sizeof(uint32), 1);
    AS_UTL_safeRead(C, &_numInputs, "numInputs", sizeof(uint32), 1);

    _inputNames = new char * [_numInputs];

    for (uint32 ii=0; ii<_numInputs; ii++) {
      uint32  nl = 0;

      AS_UTL_safeRead(C, &nl, "nameLen", sizeof(uint32), 1);

      _inputNames[ii] = new char [nl+1];

      AS_UTL_safeRead(C, _inputNames[ii], "name", sizeof(uint32), nl+1);
    }

    _inputToBucket = new uint32 [_numInputs];
    _readToSlice   = new uint16 [_maxID+1];

    AS_UTL_safeRead(C, _inputToBucket, "inputToBucket", sizeof(uint32), _numInputs);
    AS_UTL_safeRead(C, _readToSlice,   "readToSlice",   sizeof(uint16), _maxID+1);

    AS_UTL_closeFile(C, configName);
  };

  void    writeConfig(const char *configName) {
    FILE *C = AS_UTL_openOutputFile(configName);

    AS_UTL_safeWrite(C, &_maxID,     "maxID",     sizeof(uint32), 1);
    AS_UTL_safeWrite(C, &_numInputs, "numInputs", sizeof(uint32), 1);

    for (uint32 ii=0; ii<_numInputs; ii++) {
      uint32 nl = strlen(_inputNames[ii]);

      AS_UTL_safeWrite(C, &nl,             "nameLen", sizeof(uint32), 1);
      AS_UTL_safeWrite(C, _inputNames[ii], "name",    sizeof(char),   nl+1);
    }

    AS_UTL_safeWrite(C, _inputToBucket, "inputToBucket", sizeof(uint32), _numInputs);
    AS_UTL_safeWrite(C, _readToSlice,   "readToSlice",   sizeof(uint16), _maxID+1);

    AS_UTL_closeFile(C, configName);

    fprintf(stderr, "Saved configuration to '%s'.\n", configName);
  };

  void    getNumberOfInputs(void);
  uint32  getInputBucket(uint32 bucketNumber);
  char   *getInput(uint32 bucketNumber, uint32 fileNumber);

  void    assignReadsToSlices(uint64  minMemory,
                              uint64  maxMemory);

private:
  uint32     _maxID;

  uint32     _numInputs;       //  Number of input ovb files.
  char     **_inputNames;      //  Input ovb files.

  uint32    *_inputToBucket;   //  Maps an input name to a bucket.
  uint16    *_readToSlice;      //  Map each read ID to a slice.
};





void
ovStoreConfig::assignReadsToSlices(uint64          minMemory,
                                   uint64          maxMemory) {

  int64    procMax       = sysconf(_SC_CHILD_MAX);
  int64    openMax       = sysconf(_SC_OPEN_MAX) - 16;

  //
  //  Load the number of overlaps per read.
  //
  //  WE NEED TO KNOW OVERLAPS PER INPUT FILE, SO WE CAN BALANCE BUCKETIZER JOBS

  fprintf(stderr, "Finding number of overlaps per read and per file.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "      Molaps\n");
  fprintf(stderr, "------------ ----------------------------------------\n");

  uint64             *oPF  = new uint64 [_numInputs];
  uint32             *oPR  = new uint32 [_maxID+1];
  ovStoreHistogram   *hist = new ovStoreHistogram();

  for (uint32 ii=0; ii<_numInputs; ii++) {
    oPF[ii] = hist->loadData(_inputNames[ii], _maxID+1);
    fprintf(stderr, "%12.3f %40s\n", oPF[ii] / 1024.0 / 1024.0, _inputNames[ii]);
  }

  uint64  numOverlaps = hist->getOverlapsPerRead(oPR, _maxID+1);

  delete hist;

  fprintf(stderr, "------------ ----------------------------------------\n");
  fprintf(stderr, "%12.3f %40" F_U32P "\n", numOverlaps / 1024.0 / 1024.0, _numInputs);
  fprintf(stderr, "\n");

  if (numOverlaps == 0)
    fprintf(stderr, "Found no overlaps to sort.\n"), exit(1);


  //
  //  Partition the overlaps into buckets.
  //
  //  This will pick the smallest memory size that uses fewer than maxFiles buckets.  Unreasonable
  //  values can break this - either too low memory or too high allowed open files (an OS limit).
  //

  uint64  olapsPerSliceMin = (minMemory - MEMORY_OVERHEAD) / ovOverlapSortSize;
  uint64  olapsPerSliceMax = (maxMemory - MEMORY_OVERHEAD) / ovOverlapSortSize;

  fprintf(stderr, "Configuring for %.2f GB to %.2f GB memory.  Not to exceed " F_S64 " processes or " F_S64 " open files.\n",
          minMemory / 1024.0 / 1024.0 / 1024.0,
          maxMemory / 1024.0 / 1024.0 / 1024.0,
          procMax, openMax);
  fprintf(stderr, "  At minimum memory, %4" F_U64P " sort processes with %.2f Molaps each.\n",
          numOverlaps / olapsPerSliceMin + 1,
          olapsPerSliceMin / 1024.0 / 1024.0);
  fprintf(stderr, "  At maximum memory, %4" F_U64P " sort processes with %.2f Molaps each.\n",
          numOverlaps / olapsPerSliceMax + 1,
          olapsPerSliceMax / 1024.0 / 1024.0);
  fprintf(stderr, "\n");

  //  More processes -> higher bandwidth utilization, but can also thrash the FS
  //  More memory    -> harder to get machines to run on

  uint64  useMemory     = minMemory + 3 * (maxMemory - minMemory) / 4;
  uint64  olapsPerSlice = (useMemory - MEMORY_OVERHEAD) / ovOverlapSortSize;

  //  With that upper limit on the number of overlaps per slice, count how many slices
  //  we need to make.

  uint32  numSlices        = 1;
  //uint32  numBuckers       = 1;
  uint64  maxOlapsPerSlice = 0;

  for (uint64 olaps=0, ii=0; ii<_maxID+1; ii++) {
    olaps += oPR[ii];

    if (olaps >= olapsPerSlice) {
      olaps = 0;
      numSlices++;
    }
  }

  //  Now, divide those overlaps evenly among the slices and find the maximum number of overlaps in any slice.

  olapsPerSlice = (uint64)ceil((double)numOverlaps / (double)numSlices) + 1;

  for (uint64 olaps=0, ii=0; ii<_maxID+1; ii++) {
    olaps += oPR[ii];

    if (maxOlapsPerSlice < olaps)
      maxOlapsPerSlice = olaps;

    if (olaps >= olapsPerSlice)
      olaps = 0;
  }

  //  Assign inputs to each bucketizer.  Greedy load balancing.

  uint64  *olapsPerBucket = new uint64 [numSlices];

  for (uint32 ii=0; ii<numSlices; ii++)
    olapsPerBucket[ii] = 0;

  for (uint32 ss=0, ii=0; ii<_numInputs; ii++) {
    uint32  mb = 0;

    for (uint32 bb=0; bb<numSlices; bb++)
      if (olapsPerBucket[bb] < olapsPerBucket[mb])
        mb = bb;

    _inputToBucket[ii]  = mb;
    olapsPerBucket[mb] += oPF[ii];
  }

  //  Report results to nobody.

  fprintf(stderr, "Will bucketize using " F_U32 " processes.\n", numSlices);
  fprintf(stderr, "\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Will sort using " F_U32 " slices.\n",  numSlices);
  fprintf(stderr, "  Up to %7.2f M overlaps per slice\n", maxOlapsPerSlice / 1000000.0);
  fprintf(stderr, "        %7.2f GB memory  per slice\n", (maxOlapsPerSlice * ovOverlapSortSize + MEMORY_OVERHEAD) / 1024.0 / 1024.0 / 1024.0);
  fprintf(stderr, "\n");


  //  Report ovb to bucket mapping.

  fprintf(stderr, "         number    number of\n");
  fprintf(stderr, " bucket  inputs     overlaps\n");
  fprintf(stderr, "------- ------- ------------\n");

  for (uint32 ii=0; ii<numSlices; ii++) {
    uint32  ni = 0;

    for (uint32 xx=0; xx<_numInputs; xx++) {
      if (_inputToBucket[xx] == ii)
        ni++;
    }

    fprintf(stderr, "%7" F_U32P " %7" F_U32P " %12" F_U64P "\n",
            ii, ni, olapsPerBucket[ii]);
  }

  fprintf(stderr, "\n");

  //  Assign reads to slices.

  fprintf(stderr, "          number of\n");
  fprintf(stderr, " slice     overlaps\n");
  fprintf(stderr, "------ ------------\n");

  {
    uint64  olaps = 0;
    uint32  slice = 0;

    for (uint32 ii=0; ii<_maxID+1; ii++) {
      olaps            += oPR[ii];
      _readToSlice[ii]   = slice;

      if (olaps >= olapsPerSlice) {
        fprintf(stderr, "%6" F_U32P " %12" F_U64P "\n", slice, olaps);
        olaps = 0;
        slice++;
      }
    }

    fprintf(stderr, "%6" F_U32P " %12" F_U64P "\n", slice, olaps);
  }

  fprintf(stderr, "\n");

  delete [] oPR;
}




int
main(int argc, char **argv) {
  char           *gkpName        = NULL;
  uint64          minMemory      = (uint64)1 * 1024 * 1024 * 1024;
  uint64          maxMemory      = (uint64)4 * 1024 * 1024 * 1024;

  vector<char *>  fileList;

  char           *configOut    = NULL;

  argc = AS_configure(argc, argv);

  int err=0;
  int arg=1;
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
  if (gkpName == NULL)
    err++;
  if (fileList.size() == 0)
    err++;
  if (maxMemory < MEMORY_OVERHEAD + ovOverlapSortSize)
    err++;
  if (err) {
    fprintf(stderr, "usage: %s -G asm.gkpStore -config out.config [opts] [-L fileList | *.ovb]\n", argv[0]);
    fprintf(stderr, "  -G asm.gkpStore       path to gkpStore for this assembly\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -L fileList           a list of ovb files in 'fileList'\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -M g                  use up to 'g' gigabytes memory for sorting overlaps\n");
    fprintf(stderr, "                          default 4; g-0.25 gb is available for sorting overlaps\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -config out.config    write overlap store configuration to file 'out.config'\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Sizes and Limits:\n");
    fprintf(stderr, "  ovOverlap             " F_S32 " words of " F_S32 " bits each.\n", (int32)ovOverlapNWORDS, (int32)ovOverlapWORDSZ);
    fprintf(stderr, "  ovOverlapSortSize     " F_S32 " bits\n",      (int32)ovOverlapSortSize * 8);
    fprintf(stderr, "  SC_CHILD_MAX          " F_S32 " processes\n", (int32)sysconf(_SC_CHILD_MAX));
    fprintf(stderr, "  SC_OPEN_MAX           " F_S32 " files\n",     (int32)sysconf(_SC_OPEN_MAX));
    fprintf(stderr, "\n");

    if (gkpName == NULL)
      fprintf(stderr, "ERROR: No gatekeeper store (-G) supplied.\n");
    if (fileList.size() == 0)
      fprintf(stderr, "ERROR: No input overlap files (-L or last on the command line) supplied.\n");

    if (maxMemory < MEMORY_OVERHEAD + ovOverlapSortSize)
      fprintf(stderr, "ERROR: Memory (-M) must be at least %.3f GB to account for overhead.\n", MEMORY_OVERHEAD / 1024.0 / 1024.0 / 1024.0);

    exit(1);
  }

  //  Check parameters, reset some of them.

  if (minMemory < MEMORY_OVERHEAD + ovOverlapSortSize) {
    fprintf(stderr, "Reset minMemory from " F_U64 " to " F_SIZE_T "\n", minMemory, MEMORY_OVERHEAD + ovOverlapSortSize);
    minMemory  = MEMORY_OVERHEAD + ovOverlapSortSize;
  }

  //

  gkStore        *gkp    = gkStore::gkStore_open(gkpName);
  uint32          maxID = gkp->gkStore_getNumReads();

  gkp->gkStore_close();

  ovStoreConfig  *config = new ovStoreConfig(fileList, maxID);

  config->assignReadsToSlices(minMemory, maxMemory);

  config->writeConfig(configOut);

  delete config;

  fprintf(stderr, "Bye.\n");

  exit(0);
}
