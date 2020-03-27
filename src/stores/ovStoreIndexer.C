
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



int
main(int argc, char **argv) {
  char const     *ovlName     = NULL;
  char const     *seqName     = NULL;
  char const     *cfgName     = NULL;
  bool            deleteInter = false;

  argc = AS_configure(argc, argv);

  vector<char const *>  err;
  int                   arg=1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-O") == 0) {
      ovlName = argv[++arg];

    } else if (strcmp(argv[arg], "-S") == 0) {    //  Yup, not used, but left in
      seqName = argv[++arg];                      //  so it's the same as the others.

    } else if (strcmp(argv[arg], "-C") == 0) {
      cfgName = argv[++arg];

    } else if (strcmp(argv[arg], "-delete") == 0) {
      deleteInter = true;

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
    err.push_back("ERROR: No config (-C) supplied.\n");

  if (err.size() > 0) {
    fprintf(stderr, "usage: %s -O asm.ovlStore -S asm.seqStore -C ovStoreConfig [options]\n", argv[0]);
    fprintf(stderr, "  -O asm.ovlStore    path to overlap store to create\n");
    fprintf(stderr, "  -S asm.seqStore    path to sequence store\n");
    fprintf(stderr, "  -C config          path to ovStoreConfig configuration file\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -delete          remove intermediate files when the index is\n");
    fprintf(stderr, "                   successfully created\n");
    fprintf(stderr, "\n");

    for (uint32 ii=0; ii<err.size(); ii++)
      if (err[ii])
        fputs(err[ii], stderr);

    exit(1);
  }

  sqStore             *seq    = new sqStore(seqName);
  ovStoreConfig       *config = new ovStoreConfig(cfgName);
  ovStoreSliceWriter  *writer = new ovStoreSliceWriter(ovlName, seq, 0, config->numSlices(), config->numBuckets());

  writer->checkSortingIsComplete();
  writer->mergeInfoFiles();
  writer->mergeHistogram();

  if (deleteInter == true)
    writer->removeAllIntermediateFiles();

  delete writer;
  delete config;

  //  Test.  Open the store and get the number of overlaps per read.

  ovStore *tester = new ovStore(ovlName, seq);
  tester->testStore();
  delete    tester;

  //  And we have a store.  Cleanup and success!

  delete seq;

  fprintf(stderr, "\n");
  fprintf(stderr, "Success!\n");

  exit(0);
}

