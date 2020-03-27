
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

#include "sqStore.H"
#include "ovStore.H"
#include "ovStoreConfig.H"


static
void
writeToFile(sqStore          *seq,
            ovOverlap        *overlap,
            ovFile          **sliceFile,
            uint64           *sliceSize,
            ovStoreConfig    *config,
            char const       *ovlName,
            uint32            bucketNum) {

  uint32 df = config->getAssignedSlice(overlap->a_iid);

  if (sliceFile[df] == NULL) {
    char name[FILENAME_MAX];

    snprintf(name, FILENAME_MAX, "%s/create%04d/slice%04d", ovlName, bucketNum, df);
    sliceFile[df] = new ovFile(seq, name, ovFileFullWriteNoCounts);
    sliceSize[df] = 0;
  }

  if ((df < 1) ||
      (df > config->numSlices() + 1)) {
    char ovlstr[256];

    fprintf(stderr, "Invalid slice file %u in overlap %s\n",
            df, overlap->toString(ovlstr, ovOverlapAsUnaligned, false));
  }

  sliceFile[df]->writeOverlap(overlap);
  sliceSize[df]++;
}



int
main(int argc, char **argv) {
  char const     *ovlName        = NULL;
  char const     *seqName        = NULL;
  char const     *cfgName        = NULL;
  uint32          bucketNum      = UINT32_MAX;

  double          maxErrorRate   = 1.0;

  bool            deleteInputs   = false;
  bool            forceOverwrite = false;

  char            createName[FILENAME_MAX+1];
  char            sliceSName[FILENAME_MAX+1];
  char            bucketName[FILENAME_MAX+1];

  argc = AS_configure(argc, argv);

  vector<char const *>  err;
  int                   arg=1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-O") == 0) {
      ovlName = argv[++arg];

    } else if (strcmp(argv[arg], "-S") == 0) {
      seqName = argv[++arg];

    } else if (strcmp(argv[arg], "-C") == 0) {
      cfgName = argv[++arg];

    } else if (strcmp(argv[arg], "-b") == 0) {
      bucketNum = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-e") == 0) {
      maxErrorRate = atof(argv[++arg]);

    } else if (strcmp(argv[arg], "-delete") == 0) {
      deleteInputs = true;

    } else if (strcmp(argv[arg], "-f") == 0) {
      forceOverwrite = true;

    } else {
      char *s = new char [1024];
      snprintf(s, 1024, "%s: unknown option '%s'.\n", argv[0], argv[arg]);
      err.push_back(s);
    }

    arg++;
  }

  if (ovlName == NULL)
    err.push_back("ERROR: No overlap store (-O) supplied.\n");

  if (seqName == NULL)
    err.push_back("ERROR: No sequence store (-S) supplied.\n");

  if (cfgName == NULL)
    err.push_back("ERROR: No store configuration (-C) supplied.\n");

  if (bucketNum == UINT32_MAX)
    err.push_back("ERROR: Invalid or no bucket (-b) supplied.\n");

  if (err.size() > 0) {
    fprintf(stderr, "usage: %s -O asm.ovlStore -S asm.seqStore -C ovStoreConfig -b bucket [opts]\n", argv[0]);
    fprintf(stderr, "  -O asm.ovlStore       path to overlap store to create\n");
    fprintf(stderr, "  -S asm.seqStore       path to a sequence store\n");
    fprintf(stderr, "  -C config             path to ovStoreConfig configuration file\n");
    fprintf(stderr, "  -b bucket             bucket to create (1 ... N)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -e e                  filter overlaps above e fraction error\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -f                    force overwriting existing data\n");
    fprintf(stderr, "\n");

    for (uint32 ii=0; ii<err.size(); ii++)
      if (err[ii])
        fputs(err[ii], stderr);

    exit(1);
  }

  //  Load the config.

  ovStoreConfig  *config = new ovStoreConfig(cfgName);

  //  Create the output directory names and check if we're running or done.  Or if the user is a moron.

  snprintf(createName, FILENAME_MAX, "%s/create%04u",            ovlName, bucketNum);
  snprintf(sliceSName, FILENAME_MAX, "%s/create%04u/sliceSizes", ovlName, bucketNum);  //  sliceSizes name before rename.
  snprintf(bucketName, FILENAME_MAX, "%s/bucket%04u",            ovlName, bucketNum);

  if (directoryExists(createName) == true) {
    if (forceOverwrite) {
      fprintf(stderr, "Overwriting incomplete result from presumed crashed job in directory '%s'.\n", createName);
    } else {
      fprintf(stderr, "Job (appears to be) in progress (or crashed); directory '%s' exists.\n", createName);
      exit(1);
    }
  }

  if (directoryExists(bucketName) == true) {
    fprintf(stderr, "Job finished; directory '%s' exists.\n", bucketName);
    exit(0);
  }

  if ((bucketNum == 0) ||
      (bucketNum > config->numBuckets())) {
    fprintf(stderr, "No bucket " F_U32 " exists; only buckets 1-" F_U32 " exist.\n", bucketNum, config->numBuckets());
    exit(1);
  }

  //  Open inputs.

  sqStore        *seq    = new sqStore(seqName);

  fprintf(stderr, "\n");
  fprintf(stderr, "Opened '%s' with %u reads.\n", seqName, seq->sqStore_lastReadID());
  fprintf(stderr, "\n");

  //  Report options.

  fprintf(stderr, "Constructing slice " F_U32 " for store '%s'.\n", bucketNum, ovlName);
  fprintf(stderr, " - Filtering overlaps over %.4f fraction error.\n", maxErrorRate);
  fprintf(stderr, "\n");

  //  Make directories.

  AS_UTL_mkdir(ovlName);
  AS_UTL_mkdir(createName);

  //  Allocate stuff.

  ovFile        **sliceFile = new ovFile * [config->numSlices() + 1];
  uint64         *sliceSize = new uint64   [config->numSlices() + 1];

  memset(sliceFile, 0, sizeof(ovFile *) * (config->numSlices() + 1));
  memset(sliceSize, 0, sizeof(uint64)   * (config->numSlices() + 1));

  ovStoreFilter *filter = new ovStoreFilter(seq, maxErrorRate);
  ovOverlap      foverlap;
  ovOverlap      roverlap;

  //  And process each input!

  for (uint32 ff=0; ff<config->numInputs(bucketNum); ff++) {
    fprintf(stderr, "Bucketizing input %4" F_U32P " out of %4" F_U32P " - '%s'\n",
            ff+1, config->numInputs(bucketNum), config->getInput(bucketNum, ff));

    ovFile  *inputFile = new ovFile(seq, config->getInput(bucketNum, ff), ovFileFull);

    //  Do bigger buffers increase performance?  Do small ones hurt?
    //AS_OVS_setBinaryOverlapFileBufferSize(2 * 1024 * 1024);

    while (inputFile->readOverlap(&foverlap)) {
      filter->filterOverlap(foverlap, roverlap);  //  The filter copies f into r, and checks IDs

      //  Write the overlap if anything requests it.  These can be non-symmetric; e.g., if
      //  we only want to trim reads 1-1000, we'll not output any overlaps for a_iid > 1000.

      if ((foverlap.dat.ovl.forUTG == true) ||
          (foverlap.dat.ovl.forOBT == true) ||
          (foverlap.dat.ovl.forDUP == true))
        writeToFile(seq, &foverlap, sliceFile, sliceSize, config, ovlName, bucketNum);

      if ((roverlap.dat.ovl.forUTG == true) ||
          (roverlap.dat.ovl.forOBT == true) ||
          (roverlap.dat.ovl.forDUP == true))
        writeToFile(seq, &roverlap, sliceFile, sliceSize, config, ovlName, bucketNum);
    }

    delete inputFile;
  }

  //  Report what we've filtered.



  //  Write the outputs.

  for (uint32 i=0; i<config->numSlices() + 1; i++)
    delete sliceFile[i];

  AS_UTL_saveFile(sliceSName, sliceSize, config->numSlices() + 1);

  //  Rename the bucket to show we're done.

  AS_UTL_rename(createName, bucketName);

  //  Delete the inputs, if requested.

  if (deleteInputs) {
    fprintf(stderr, "\n");
    fprintf(stderr, "Deleting inputs.\n");
    fprintf(stderr, "\n");

    for (uint32 ff=0; ff<config->numInputs(bucketNum); ff++) {
      fprintf(stderr, "Deleting input %4" F_U32P " out of %4" F_U32P " - '%s'\n",
              ff+1, config->numInputs(bucketNum), config->getInput(bucketNum, ff));

      ovFile::deleteDiskFiles(config->getInput(bucketNum, ff));
    }

    fprintf(stderr, "\n");
  }

  //  Cleanup and be done.

  delete [] sliceFile;
  delete [] sliceSize;

  delete seq;

  delete    filter;
  delete    config;

  fprintf(stderr, "Success!\n");

  return(0);
}
