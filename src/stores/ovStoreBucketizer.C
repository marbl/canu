
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
 *    src/AS_OVS/overlapStoreBucketizer.C
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2012-APR-02 to 2013-AUG-01
 *      are Copyright 2012-2013 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2014-AUG-22 to 2015-SEP-21
 *      are Copyright 2014-2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2015-NOV-08
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *    Sergey Koren beginning on 2016-FEB-29
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_global.H"

#include "gkStore.H"
#include "ovStore.H"


static
void
writeToFile(gkStore       *gkp,
            ovOverlap    *overlap,
            ovFile       **sliceFile,
            uint32         sliceFileMax,
            uint64        *sliceSize,
            uint32        *iidToBucket,
            char          *ovlName,
            uint32         jobIndex,
            bool           useGzip) {

  uint32 df = iidToBucket[overlap->a_iid];

  if (sliceFile[df] == NULL) {
    char name[FILENAME_MAX];

    snprintf(name, FILENAME_MAX, "%s/create%04d/slice%04d%s", ovlName, jobIndex, df, (useGzip) ? ".gz" : "");
    sliceFile[df] = new ovFile(gkp, name, ovFileFullWriteNoCounts);
    sliceSize[df] = 0;
  }

  sliceFile[df]->writeOverlap(overlap);
  sliceSize[df]++;
}



int
main(int argc, char **argv) {
  char           *ovlName      = NULL;
  char           *gkpName      = NULL;
  char           *cfgName      = NULL;
  uint32          maxFiles     = sysconf(_SC_OPEN_MAX) - 16;
  uint32          fileLimit    = maxFiles;

  uint32          jobIndex     = 0;

  double          maxErrorRate = 1.0;
  uint64          maxError     = AS_OVS_encodeEvalue(maxErrorRate);

  char           *ovlInput     = NULL;

  bool            useGzip      = false;

  argc = AS_configure(argc, argv);

  int err=0;
  int arg=1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-O") == 0) {
      ovlName = argv[++arg];

    } else if (strcmp(argv[arg], "-G") == 0) {
      gkpName = argv[++arg];

    } else if (strcmp(argv[arg], "-C") == 0) {
      cfgName = argv[++arg];

    } else if (strcmp(argv[arg], "-F") == 0) {
      fileLimit = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-job") == 0) {
      jobIndex = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-i") == 0) {
      ovlInput = argv[++arg];

    } else if (strcmp(argv[arg], "-e") == 0) {
      maxErrorRate = atof(argv[++arg]);
      maxError     = AS_OVS_encodeEvalue(maxErrorRate);

    } else if (strcmp(argv[arg], "-gzip") == 0) {
      useGzip = true;

    } else {
      fprintf(stderr, "ERROR: unknown option '%s'\n", argv[arg]);
      err++;
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
  if (fileLimit > maxFiles)
    err++;

  if (err) {
    fprintf(stderr, "usage: %s -O asm.ovlStore -G asm.gkpStore -i file.ovb -job j [opts]\n", argv[0]);
    fprintf(stderr, "  -O asm.ovlStore       path to store to create\n");
    fprintf(stderr, "  -G asm.gkpStore       path to gkpStore for this assembly\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -C config             path to previously created ovStoreBuild config data file\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -i file.ovb[.gz]      input overlaps\n");
    fprintf(stderr, "  -job j                index of this overlap input file\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -F f                  use up to 'f' files for store creation\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -obt                  filter overlaps for OBT\n");
    fprintf(stderr, "  -dup                  filter overlaps for OBT/dedupe\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -e e                  filter overlaps above e fraction error\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -gzip                 compress buckets even more\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    DANGER    DO NOT USE     DO NOT USE     DO NOT USE    DANGER\n");
    fprintf(stderr, "    DANGER                                                DANGER\n");
    fprintf(stderr, "    DANGER   This command is difficult to run by hand.    DANGER\n");
    fprintf(stderr, "    DANGER          Use ovStoreCreate instead.            DANGER\n");
    fprintf(stderr, "    DANGER                                                DANGER\n");
    fprintf(stderr, "    DANGER    DO NOT USE     DO NOT USE     DO NOT USE    DANGER\n");
    fprintf(stderr, "\n");

    if (ovlName == NULL)
      fprintf(stderr, "ERROR: No overlap store (-O) supplied.\n");
    if (gkpName == NULL)
      fprintf(stderr, "ERROR: No gatekeeper store (-G) supplied.\n");
    if (ovlInput == NULL)
      fprintf(stderr, "ERROR: No input (-i) supplied.\n");
    if (jobIndex == 0)
      fprintf(stderr, "ERROR: No job index (-job) supplied.\n");
    if (fileLimit > maxFiles)
      fprintf(stderr, "ERROR: Too many jobs (-F); only " F_U32 " supported on this architecture.\n", maxFiles);

    exit(1);
  }


  {
    if (AS_UTL_fileExists(ovlName, TRUE, FALSE) == false)
      AS_UTL_mkdir(ovlName);
  }


  {
    char name[FILENAME_MAX];

    snprintf(name, FILENAME_MAX, "%s/create%04d", ovlName, jobIndex);

    if (AS_UTL_fileExists(name, TRUE, FALSE) == false)
      AS_UTL_mkdir(name);
    else
      fprintf(stderr, "Overwriting previous result; directory '%s' exists.\n", name), exit(0);
  }


  {
    char name[FILENAME_MAX];

    snprintf(name, FILENAME_MAX, "%s/bucket%04d/sliceSizes", ovlName, jobIndex);

    if (AS_UTL_fileExists(name, FALSE, FALSE) == true)
      fprintf(stderr, "Job finished; file '%s' exists.\n", name), exit(0);
  }


  gkStore *gkp         = gkStore::gkStore_open(gkpName);

  uint32  maxIID       = gkp->gkStore_getNumReads() + 1;
  uint32 *iidToBucket  = new uint32 [maxIID];

  {
    errno = 0;
    FILE *C = fopen(cfgName, "r");
    if (errno)
      fprintf(stderr, "ERROR: failed to open config file '%s' for reading: %s\n", cfgName, strerror(errno)), exit(1);

    uint32  maxIIDtest  = 0;

    AS_UTL_safeRead(C, &maxIIDtest,  "maxIID",      sizeof(uint32), 1);
    AS_UTL_safeRead(C,  iidToBucket, "iidToBucket", sizeof(uint32), maxIID);

    if (maxIIDtest != maxIID)
      fprintf(stderr, "ERROR: maxIID in store (" F_U32 ") differs from maxIID in config file (" F_U32 ").\n",
              maxIID, maxIIDtest), exit(1);
  }


  ovFile       **sliceFile = new ovFile * [fileLimit + 1];
  uint64        *sliceSize = new uint64   [fileLimit + 1];

  memset(sliceFile, 0, sizeof(ovFile *) * (fileLimit + 1));
  memset(sliceSize, 0, sizeof(uint64)   * (fileLimit + 1));

  fprintf(stderr, "maxError fraction: %.3f percent: %.3f encoded: " F_U64 "\n",
          maxErrorRate, maxErrorRate * 100, maxError);

  fprintf(stderr, "Bucketizing %s\n", ovlInput);

  ovStoreFilter *filter = new ovStoreFilter(gkp, maxError);
  ovOverlap      foverlap(gkp);
  ovOverlap      roverlap(gkp);
  ovFile         *inputFile = new ovFile(gkp, ovlInput, ovFileFull);

  //  Do bigger buffers increase performance?  Do small ones hurt?
  //AS_OVS_setBinaryOverlapFileBufferSize(2 * 1024 * 1024);

  while (inputFile->readOverlap(&foverlap)) {
    filter->filterOverlap(foverlap, roverlap);  //  The filter copies f into r, and checks IDs

    //  If all are skipped, don't bother writing the overlap.

    if ((foverlap.dat.ovl.forUTG == true) ||
        (foverlap.dat.ovl.forOBT == true) ||
        (foverlap.dat.ovl.forDUP == true))
      writeToFile(gkp, &foverlap, sliceFile, fileLimit, sliceSize, iidToBucket, ovlName, jobIndex, useGzip);

    if ((roverlap.dat.ovl.forUTG == true) ||
        (roverlap.dat.ovl.forOBT == true) ||
        (roverlap.dat.ovl.forDUP == true))
      writeToFile(gkp, &roverlap, sliceFile, fileLimit, sliceSize, iidToBucket, ovlName, jobIndex, useGzip);
  }

  delete inputFile;

#warning not reporting fate
  //filter->reportFate();
  //filter->resetCounters();

  delete filter;

  for (uint32 i=0; i<=fileLimit; i++)
    delete sliceFile[i];

  //  Write slice sizes, rename bucket.

  {
    char name[FILENAME_MAX];
    char finl[FILENAME_MAX];

    snprintf(name, FILENAME_MAX, "%s/create%04d/sliceSizes", ovlName, jobIndex);

    FILE *F = fopen(name, "w");
    if (errno)
      fprintf(stderr, "ERROR:  Failed to open %s: %s\n", name, strerror(errno)), exit(1);

    AS_UTL_safeWrite(F, sliceSize, "sliceSize", sizeof(uint64), fileLimit + 1);

    fclose(F);

    snprintf(name, FILENAME_MAX, "%s/create%04d", ovlName, jobIndex);
    snprintf(finl, FILENAME_MAX, "%s/bucket%04d", ovlName, jobIndex);

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
