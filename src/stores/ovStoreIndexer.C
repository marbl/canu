
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



int
main(int argc, char **argv) {
  char           *storePath    = NULL;
  uint32          fileLimit    = 0;         //  Number of 'slices' from bucketizer

  bool            deleteIntermediates = true;

  bool            doExplicitTest = false;
  bool            doFixes        = false;

  char            name[FILENAME_MAX];

  argc = AS_configure(argc, argv);

  int err=0;
  int arg=1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-O") == 0) {
      storePath = argv[++arg];

    } else if (strcmp(argv[arg], "-F") == 0) {
      fileLimit = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-f") == 0) {
      doFixes = true;

    } else if (strcmp(argv[arg], "-t") == 0) {
      doExplicitTest = true;
      storePath = argv[++arg];

    } else if (strcmp(argv[arg], "-nodelete") == 0) {
      deleteIntermediates = false;

    } else {
      fprintf(stderr, "ERROR: unknown option '%s'\n", argv[arg]);
      err++;
    }

    arg++;
  }
  if (storePath == NULL)
    err++;
  if ((fileLimit == 0) && (doExplicitTest == false))
    err++;

  if (err) {
    fprintf(stderr, "usage: %s ...\n", argv[0]);
    fprintf(stderr, "  -O x.ovlStore    path to overlap store to build the final index for\n");
    fprintf(stderr, "  -F s             number of slices used in bucketizing/sorting\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -t x.ovlStore    explicitly test a previously constructed index\n");
    fprintf(stderr, "  -f               when testing, also create a new 'idx.fixed' which might\n");
    fprintf(stderr, "                   resolve rare problems\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -nodelete        do not remove intermediate files when the index is\n");
    fprintf(stderr, "                   successfully created\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    DANGER    DO NOT USE     DO NOT USE     DO NOT USE    DANGER\n");
    fprintf(stderr, "    DANGER                                                DANGER\n");
    fprintf(stderr, "    DANGER   This command is difficult to run by hand.    DANGER\n");
    fprintf(stderr, "    DANGER          Use ovStoreCreate instead.            DANGER\n");
    fprintf(stderr, "    DANGER                                                DANGER\n");
    fprintf(stderr, "    DANGER    DO NOT USE     DO NOT USE     DO NOT USE    DANGER\n");
    fprintf(stderr, "\n");

    if (storePath == NULL)
      fprintf(stderr, "ERROR: No overlap store (-O) supplied.\n");
    if ((fileLimit == 0) && (doExplicitTest == false))
      fprintf(stderr, "ERROR: One of -F (number of slices) or -t (test a store) must be supplied.\n");

    exit(1);
  }

  //  Do the test, and maybe fix things up.

  //gkStore        *gkp    = gkStore::gkStore_open(gkpName);
  ovStoreWriter  *writer = new ovStoreWriter(storePath, NULL, fileLimit, 0, 0);

  if (doExplicitTest == true) {
    bool  passed = writer->testIndex(doFixes);
    if (passed == true)
      fprintf(stderr, "Index looks correct.\n");
    delete writer;
    exit(passed == false);
  }

  //  Check that all segments are present.  Every segment should have an info file.

  writer->checkSortingIsComplete();

  //  Merge the indices and histogram data.

  writer->mergeInfoFiles();
  writer->mergeHistogram();

  //  Diagnostics.

  if (writer->testIndex(false) == false) {
    fprintf(stderr, "ERROR: index failed tests.\n");
    delete writer;
    exit(1);
  }

  //  Remove intermediates.  For the buckets, we keep going until there are 10 in a row not present.
  //  During testing, on a microbe using 2850 buckets, some buckets were empty.

  if (deleteIntermediates == false) {
    fprintf(stderr, "\n");
    fprintf(stderr, "Not removing intermediate files.  Finished.\n");
    exit(0);
  }

  fprintf(stderr, "\n");
  fprintf(stderr, "Removing intermediate files.\n");

  writer->removeAllIntermediateFiles();

  fprintf(stderr, "Success!\n");

  exit(0);
}

