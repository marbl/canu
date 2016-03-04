
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
  char           *ovlName      = NULL;
  uint32          maxJob       = 0;

  bool            deleteIntermediates = true;

  bool            doExplicitTest = false;
  bool            doFixes        = false;

  char            name[FILENAME_MAX];

  argc = AS_configure(argc, argv);

  int err=0;
  int arg=1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-O") == 0) {
      ovlName = argv[++arg];

    } else if (strcmp(argv[arg], "-F") == 0) {
      maxJob = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-f") == 0) {
      doFixes = true;

    } else if (strcmp(argv[arg], "-t") == 0) {
      doExplicitTest = true;
      ovlName = argv[++arg];

    } else if (strcmp(argv[arg], "-nodelete") == 0) {
      deleteIntermediates = false;

    } else {
      fprintf(stderr, "ERROR: unknown option '%s'\n", argv[arg]);
    }

    arg++;
  }
  if (ovlName == NULL)
    err++;
  if ((maxJob == 0) && (doExplicitTest == false))
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

    if (ovlName == NULL)
      fprintf(stderr, "ERROR: No overlap store (-O) supplied.\n");
    if ((maxJob == 0) && (doExplicitTest == false))
      fprintf(stderr, "ERROR: One of -F (number of slices) or -t (test a store) must be supplied.\n");

    exit(1);
  }

  //  Do the test, and maybe fix things up.

  if (doExplicitTest == true) {
    bool passed = testIndex(ovlName, doFixes);

    exit((passed == true) ? 0 : 1);
  }

  //  Check that all segments are present.  Every segment should have an info file.

  uint32  cntJob = 0;

  for (uint32 i=1; i<=maxJob; i++) {
    uint32  complete = 0;

    sprintf(name, "%s/%04d", ovlName, i);
    if (AS_UTL_fileExists(name, FALSE, FALSE) == true)
      complete++;
    else
      fprintf(stderr, "ERROR: Segment "F_U32" data not present  (%s)\n", i, name);

    sprintf(name, "%s/%04d.info", ovlName, i);
    if (AS_UTL_fileExists(name, FALSE, FALSE) == true)
      complete++;
    else
      fprintf(stderr, "ERROR: Segment "F_U32" info not present (%s)\n", i, name);

    sprintf(name, "%s/%04d.index", ovlName, i);
    if (AS_UTL_fileExists(name, FALSE, FALSE) == true)
      complete++;
    else
      fprintf(stderr, "ERROR: Segment "F_U32" index not present (%s)\n", i, name);

    if (complete == 3)
      cntJob++;
  }

  if (cntJob != maxJob) {
    fprintf(stderr, "ERROR: Expected "F_U32" segments, only found "F_U32".\n", maxJob, cntJob);
    exit(1);
  }

  //  Merge the stuff.

  mergeInfoFiles(ovlName, maxJob);

  //  Diagnostics.

  if (testIndex(ovlName, false) == false) {
    fprintf(stderr, "ERROR: index failed tests.\n");
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

  //  Removing indices is easy, beacuse we know how many there are.

  for (uint32 i=1; i<=maxJob; i++) {
    sprintf(name, "%s/%04u.index", ovlName, i);   AS_UTL_unlink(name);
    sprintf(name, "%s/%04u.info",  ovlName, i);   AS_UTL_unlink(name);
  }

  //  We don't know how many buckets there are, so we remove until we fail to find ten
  //  buckets in a row.

  for (uint32 missing=0, i=1; missing<10; i++) {
    sprintf(name, "%s/bucket%04d", ovlName, i);

    if (AS_UTL_fileExists(name, TRUE, FALSE) == FALSE) {
      missing++;
      continue;
    }

    missing = 0;

    sprintf(name, "%s/bucket%04d/sliceSizes", ovlName, i);
    AS_UTL_unlink(name);

    sprintf(name, "%s/bucket%04d", ovlName, i);
    rmdir(name);
  }

  fprintf(stderr, "Finished.\n");

  exit(0);
}

