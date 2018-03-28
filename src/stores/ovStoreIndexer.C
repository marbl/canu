
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
 *    src/AS_OVS/overlapStoreIndexer.C
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2012-APR-02 to 2013-AUG-01
 *      are Copyright 2012-2013 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2014-DEC-15 to 2015-SEP-21
 *      are Copyright 2014-2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2016-JAN-11
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_global.H"

#include "gkStore.H"
#include "ovStore.H"
#include "ovStoreConfig.H"



int
main(int argc, char **argv) {
  char           *ovlName     = NULL;
  char           *gkpName     = NULL;
  char           *cfgName     = NULL;
  bool            deleteInter = false;
  bool            testFailed  = false;

  argc = AS_configure(argc, argv);

  vector<char *>  err;
  int             arg=1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-O") == 0) {
      ovlName = argv[++arg];

    } else if (strcmp(argv[arg], "-G") == 0) {    //  Yup, not used, but left in
      gkpName = argv[++arg];                      //  so it's the same as the others.

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

  if (gkpName == NULL)
    err.push_back("ERROR: No gatekeeper store (-G) supplied.\n");

  if (cfgName == NULL)
    err.push_back("ERROR: No config (-C) supplied.\n");

  if (err.size()) {
    fprintf(stderr, "usage: %s -O asm.ovlStore -G asm.gkpStore -C ovStoreConfig [options]\n", argv[0]);
    fprintf(stderr, "  -O asm.ovlStore    path to overlap store to create\n");
    fprintf(stderr, "  -G asm.gkpStore       path to gatekeeper store\n");
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

  //  Load the config, make a new writer.

  gkStore             *gkp    = gkStore::gkStore_open(gkpName);
  ovStoreConfig       *config = new ovStoreConfig(cfgName);
  ovStoreSliceWriter  *writer = new ovStoreSliceWriter(ovlName, gkp, 0, config->numSlices(), config->numBuckets());

  //  Finalize the store if not testing.

  fprintf(stderr, "\n");
  fprintf(stderr, "Checking that sorting is complete (every slice should have an info file).\n");

  writer->checkSortingIsComplete();

  fprintf(stderr, "\n");
  fprintf(stderr, "Merging indices.\n");

  writer->mergeInfoFiles();

  fprintf(stderr, "\n");
  fprintf(stderr, "Merging histograms.\n");

  writer->mergeHistogram();

  //  Do some tests.

#if 0
  testFailed = (writer->testIndex(doFixes) == true) ? false : true;

  if (testFailed == false)
    fprintf(stderr, "Index looks correct.\n");
  else
    deleteInter = false;
#endif

  //  Remove the intermediate files.

  if (deleteInter == true) {
    fprintf(stderr, "\n");
    fprintf(stderr, "Removing intermediate files.\n");

    writer->removeAllIntermediateFiles();
  }

  delete writer;

  gkp->gkStore_close();

  if (testFailed)
    fprintf(stderr, "Index not valid.  Intermediate files not removed.\n");
  else
    fprintf(stderr, "Success!\n");

  exit(testFailed);
}

