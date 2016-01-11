
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
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_global.H"

#include "gkStore.H"
#include "ovStore.H"

#include <vector>
#include <algorithm>

using namespace std;


uint32  lastLibFirstIID = 0;
uint32  lastLibLastIID  = 0;


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

    numOverlaps += 2 * no / sizeof(ovOverlap);
  }

  fprintf(stderr, "Found %.3f million overlaps.\n", numOverlaps / 1000000.0);
  assert(numOverlaps > 0);

  //  Why the +1 below?  Consider the case when the number of overlaps is less than the number of
  //  fragments.  This value is used to figure out how many IIDs we can fit into a single bucket,
  //  and making it too large means we'll get maybe one more bucket and the buckets will be smaller.
  //  Yeah, we probably could have just used ceil.
  //
  double  overlapsPerBucket   = (double)memoryLimit / (double)sizeof(ovOverlap);
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
writeToDumpFile(ovOverlap       *overlap,
                ovFile          **dumpFile,
                uint32            dumpFileMax,
                uint64           *dumpLength,
                uint32            iidPerBucket,
                char             *ovlName) {

  uint32 df = overlap->a_iid / iidPerBucket;

  if (lastLibFirstIID > 0) {
    uint32  firstHighDensity = lastLibFirstIID;      //  IID of first pacBio read
    uint32  lastHighDensity  = lastLibLastIID + 1;  //  IID of last pacBio read, plus 1
    uint32  numHighDensity   = lastHighDensity - firstHighDensity;

    uint32  lowDensity       = firstHighDensity /  64;  //  64 buckets for illumina overlaps
    uint32  highDensity      = numHighDensity   / 128;  //  128 buckets for dense overlaps

    if (overlap->a_iid < firstHighDensity)
      df = overlap->a_iid / lowDensity;
    else
      df = (overlap->a_iid - firstHighDensity) / highDensity + 64;  //  plus 64 buckets from above
  }

  //fprintf(stderr, "IID %u DF %u\n", overlap->a_iid, df);

  if (df >= dumpFileMax) {
    char   olapstring[256];

    fprintf(stderr, "\n");
    fprintf(stderr, "Too many bucket files when adding overlap:\n");
    fprintf(stderr, "  Aid "F_U32"  Bid "F_U32"\n",  overlap->a_iid, overlap->b_iid);
    fprintf(stderr, "\n");
    fprintf(stderr, "bucket       = "F_U32"\n", df);
    fprintf(stderr, "iidPerBucket = "F_U32"\n", iidPerBucket);
    fprintf(stderr, "dumpFileMax  = "F_U32"\n", dumpFileMax);
    fprintf(stderr, "\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "This might be a corrupt input file, or maybe you simply need to supply more\n");
    fprintf(stderr, "memory with the canu option ovlStoreMemory.\n");
    fprintf(stderr, "\n");
    exit(1);
  }

  if (dumpFile[df] == NULL) {
    char name[FILENAME_MAX];
    sprintf(name, "%s/tmp.sort.%03d", ovlName, df);
    fprintf(stderr, "CREATE bucket '%s'\n", name);
    dumpFile[df]   = new ovFile(name, ovFileFullWrite);
    dumpLength[df] = 0;
  }

  dumpFile[df]->writeOverlap(overlap);
  dumpLength[df]++;
}



int
main(int argc, char **argv) {
  char           *ovlName      = NULL;
  char           *gkpName      = NULL;
  uint32          fileLimit    = sysconf(_SC_OPEN_MAX) - 16;
  uint64          memoryLimit  = 0;

  double          maxError     = 1.0;
  uint32          minOverlap   = 0;

  vector<char *>  fileList;

  uint32          nThreads = 4;

  bool            eValues = false;


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
      memoryLimit  = 0;

    } else if (strcmp(argv[arg], "-M") == 0) {
      fileLimit    = 0;
      memoryLimit  = atoi(argv[++arg]);
      memoryLimit *= 1024;
      memoryLimit *= 1024;
      memoryLimit *= 1024;

    } else if (strcmp(argv[arg], "-e") == 0) {
      maxError = atof(argv[++arg]);

    } else if (strcmp(argv[arg], "-e") == 0) {
      minOverlap = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-L") == 0) {
      AS_UTL_loadFileList(argv[++arg], fileList);

    } else if (strcmp(argv[arg], "-big") == 0) {
      lastLibFirstIID = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-evalues") == 0) {
      eValues = true;

    } else if ((argv[arg][0] == '-') && (argv[arg][1] != 0)) {
      fprintf(stderr, "%s: unknown option '%s'.\n", argv[0], argv[arg]);
      err++;

    } else if (AS_UTL_fileExists(argv[arg])) {
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
  if (err) {
    fprintf(stderr, "usage: %s -O asm.ovlStore -G asm.gkpStore [opts] [-L fileList | *.ovb.gz]\n", argv[0]);
    fprintf(stderr, "  -O asm.ovlStore       path to store to create\n");
    fprintf(stderr, "  -G asm.gkpStore       path to gkpStore for this assembly\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -F f                  use up to 'f' files for store creation\n");
    fprintf(stderr, "  -M m                  use up to 'm' gigabytes memory for store creation\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -e e                  filter overlaps above e fraction error\n");
    fprintf(stderr, "  -l l                  filter overlaps below l bases overlap length (needs gkpStore to get read lengths!)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -L fileList           read input filenames from 'flieList'\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -big iid              handle a large number of overlaps in the last library\n");
    fprintf(stderr, "                        iid is the first read iid in the last library, from\n");
    fprintf(stderr, "                        'gatekeeper -dumpinfo *gkpStore'\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -evalues              Input files are evalue updates from overlap error adjustment\n");
    fprintf(stderr, "\n");

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





  if (eValues) {
    ovStore  *ovs = new ovStore(ovlName, NULL);

    for (uint32 i=0; i<fileList.size(); i++) {
      errno = 0;
      FILE  *fp = fopen(fileList[i], "r");
      if (errno)
        fprintf(stderr, "Failed to open evalues file '%s': %s\n", fileList[i], strerror(errno));

      uint32        bgnID = 0;
      uint32        endID = 0;
      uint64        len   = 0;

      fprintf(stderr, "loading evalues from '%s'\n", fileList[i]);

      AS_UTL_safeRead(fp, &bgnID, "loid",   sizeof(uint32), 1);
      AS_UTL_safeRead(fp, &endID, "hiid",   sizeof(uint32), 1);
      AS_UTL_safeRead(fp, &len,   "len",    sizeof(uint64), 1);

      uint16 *evalues = new uint16 [len];

      AS_UTL_safeRead(fp, evalues, "evalues", sizeof(uint16), len);

      fclose(fp);

      fprintf(stderr, "loading evalues from '%s' -- ID range "F_U32"-"F_U32" with "F_U64" overlaps\n",
              fileList[i], bgnID, endID, len);

      ovs->addEvalues(bgnID, endID, evalues, len);

      delete [] evalues;
    }

    delete ovs;

    exit(0);
  }






  //  We create the store early, allowing it to fail if it already
  //  exists, or just cannot be created.

  gkStore  *gkp         = gkStore::gkStore_open(gkpName);
  ovStore  *storeFile   = new ovStore(ovlName, gkp, ovStoreWrite);

  uint64    maxIID       = gkp->gkStore_getNumReads() + 1;
  uint64    iidPerBucket = computeIIDperBucket(fileLimit, memoryLimit, maxIID, fileList);

  lastLibLastIID       = gkp->gkStore_getNumReads();

  uint32    dumpFileMax  = sysconf(_SC_OPEN_MAX) + 1;
  ovFile  **dumpFile     = new ovFile * [dumpFileMax];
  uint64   *dumpLength   = new uint64   [dumpFileMax];

  memset(dumpFile,   0, sizeof(ovFile *) * dumpFileMax);
  memset(dumpLength, 0, sizeof(uint64)   * dumpFileMax);

  if (maxIID / iidPerBucket + 1 > dumpFileMax - 16) {
    fprintf(stderr, "ERROR:\n");
    fprintf(stderr, "ERROR:  Operating system limit of %d open files.  The current -F setting\n", dumpFileMax);
    fprintf(stderr, "ERROR:  will need to create "F_U64" files to construct the store.\n", maxIID / iidPerBucket + 1);
    exit(1);
  }



  //  Read the gkStore to determine which fragments we care about.

  ovStoreFilter *filter = new ovStoreFilter(gkp, maxError);

  for (uint32 i=0; i<fileList.size(); i++) {
    ovOverlap    foverlap(gkp);
    ovOverlap    roverlap(gkp);

    fprintf(stderr, "bucketizing %s\n", fileList[i]);

    ovFile *inputFile = new ovFile(fileList[i], ovFileFull);

    while (inputFile->readOverlap(&foverlap)) {
      filter->filterOverlap(foverlap, roverlap);  //  The filter copies f into r

      //  If all are skipped, don't bother writing the overlap.

      if ((foverlap.dat.ovl.forUTG == true) ||
          (foverlap.dat.ovl.forOBT == true) ||
          (foverlap.dat.ovl.forDUP == true))
        writeToDumpFile(&foverlap, dumpFile, dumpFileMax, dumpLength, iidPerBucket, ovlName);

      if ((roverlap.dat.ovl.forUTG == true) ||
          (roverlap.dat.ovl.forOBT == true) ||
          (roverlap.dat.ovl.forDUP == true))
        writeToDumpFile(&roverlap, dumpFile, dumpFileMax, dumpLength, iidPerBucket, ovlName);
    }

    delete inputFile;

    //  AFTER EVERY FILE

    filter->reportFate();
    filter->resetCounters();
  }

  delete filter;

  for (uint32 i=0; i<dumpFileMax; i++)
    delete dumpFile[i];

  fprintf(stderr, "bucketizing DONE!\n");



  //
  //  Read each bucket, sort it, and dump it to the store
  //

  uint64 dumpLengthMax = 0;
  for (uint32 i=0; i<dumpFileMax; i++)
    if (dumpLengthMax < dumpLength[i])
      dumpLengthMax = dumpLength[i];

  ovOverlap  *overlapsort = ovOverlap::allocateOverlaps(gkp, dumpLengthMax);

  time_t  beginTime = time(NULL);

  for (uint32 i=0; i<dumpFileMax; i++) {
    char      name[FILENAME_MAX];
    ovFile   *bof = NULL;

    if (dumpLength[i] == 0)
      continue;

    //  We're vastly more efficient if we skip the AS_OVS interface
    //  and just suck in the whole file directly....BUT....we can't do
    //  that because the AS_OVS interface is rearranging the data to
    //  make sure the store is cross-platform compatible.

    sprintf(name, "%s/tmp.sort.%03d", ovlName, i);
    fprintf(stderr, "reading %s (%ld)\n", name, time(NULL) - beginTime);

    bof = new ovFile(name, ovFileFull);

    uint64 numOvl = 0;
    while (bof->readOverlap(overlapsort + numOvl)) {

      //  Quick sanity check on IIDs.

      if ((overlapsort[numOvl].a_iid == 0) ||
          (overlapsort[numOvl].b_iid == 0) ||
          (overlapsort[numOvl].a_iid >= maxIID) ||
          (overlapsort[numOvl].b_iid >= maxIID)) {
        char ovlstr[256];

        fprintf(stderr, "Overlap has IDs out of range (maxIID "F_U64"), possibly corrupt input data.\n", maxIID);
        fprintf(stderr, "  Aid "F_U32"  Bid "F_U32"\n",  overlapsort[numOvl].a_iid, overlapsort[numOvl].b_iid);
        exit(1);
      }

      numOvl++;
    }

    delete bof;


    assert(numOvl == dumpLength[i]);
    assert(numOvl <= dumpLengthMax);

    //  There's no real advantage to saving this file until after we
    //  write it out.  If we crash anywhere during the build, we are
    //  forced to restart from scratch.  I'll argue that removing it
    //  early helps us to not crash from running out of disk space.
    //
    unlink(name);

    fprintf(stderr, "sorting %s (%ld)\n", name, time(NULL) - beginTime);

#ifdef _GLIBCXX_PARALLEL
    //  If we have the parallel STL, don't use it!  Sort is not inplace!
    __gnu_sequential::sort(overlapsort, overlapsort + dumpLength[i]);
#else
    sort(overlapsort, overlapsort + dumpLength[i]);
#endif

    fprintf(stderr, "writing %s (%ld)\n", name, time(NULL) - beginTime);
    for (uint64 x=0; x<dumpLength[i]; x++)
      storeFile->writeOverlap(overlapsort + x);
  }

  delete    storeFile;
  delete [] overlapsort;

  //  And we have a store.

  exit(0);
}
