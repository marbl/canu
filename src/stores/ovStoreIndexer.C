
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 2007, J. Craig Venter Institute. All rights reserved.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received (LICENSE.txt) a copy of the GNU General Public
 * License along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *************************************************************************/

const char *mainid = "$Id$";

#include "AS_global.H"

#include "gkStore.H"
#include "ovStore.H"



int
main(int argc, char **argv) {
  char           *storePath    = NULL;
  uint32          maxJob       = 0;

  bool            deleteIntermediates = true;

  bool            doExplicitTest = false;
  bool            doFixes        = false;

  char            name[FILENAME_MAX];

  argc = AS_configure(argc, argv);

  int err=0;
  int arg=1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-o") == 0) {
      storePath = argv[++arg];

    } else if (strcmp(argv[arg], "-F") == 0) {
      maxJob = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-f") == 0) {
      doFixes = true;

    } else if (strcmp(argv[arg], "-t") == 0) {
      doExplicitTest = true;
      storePath = argv[++arg];

    } else if (strcmp(argv[arg], "-nodelete") == 0) {
      deleteIntermediates = false;

    } else {
      fprintf(stderr, "ERROR: unknown option '%s'\n", argv[arg]);
    }

    arg++;
  }
  if (storePath == NULL)
    err++;
  if ((maxJob == 0) && (doExplicitTest == false))
    err++;
  if (err) {
    fprintf(stderr, "usage: %s ...\n", argv[0]);
    fprintf(stderr, "  -o x.ovlStore    path to overlap store to build the final index for\n");
    fprintf(stderr, "  -F s             number of slices used in bucketizing/sorting\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -t x.ovlStore    explicitly test a previously constructed index\n");
    fprintf(stderr, "  -f               when testing, also create a new 'idx.fixed' which might\n");
    fprintf(stderr, "                   resolve rare problems\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -nodelete        do not remove intermediate files when the index is\n");
    fprintf(stderr, "                   successfully created\n");
    exit(1);
  }

  //  Do the test, and maybe fix things up.

  if (doExplicitTest == true) {
    bool passed = testIndex(storePath, doFixes);

    exit((passed == true) ? 0 : 1);
  }

  //  Check that all segments are present.  Every segment should have an info file.

  uint32  cntJob = 0;

  for (uint32 i=1; i<=maxJob; i++) {
    uint32  complete = 0;

    sprintf(name, "%s/%04d", storePath, i);
    if (AS_UTL_fileExists(name, FALSE, FALSE) == true)
      complete++;
    else
      fprintf(stderr, "ERROR: Segment "F_U32" data not present  (%s)\n", i, name);

    sprintf(name, "%s/%04d.info", storePath, i);
    if (AS_UTL_fileExists(name, FALSE, FALSE) == true)
      complete++;
    else
      fprintf(stderr, "ERROR: Segment "F_U32" info not present (%s)\n", i, name);

    sprintf(name, "%s/%04d.index", storePath, i);
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

  mergeInfoFiles(storePath, maxJob);

  //  Diagnostics.

  if (testIndex(storePath, false) == false) {
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
    sprintf(name, "%s/%04u.index", storePath, i);   AS_UTL_unlink(name);
    sprintf(name, "%s/%04u.info",  storePath, i);   AS_UTL_unlink(name);
  }

  //  We don't know how many buckets there are, so we remove until we fail to find ten
  //  buckets in a row.

  for (uint32 missing=0, i=1; missing<10; i++) {
    sprintf(name, "%s/bucket%04d", storePath, i);

    if (AS_UTL_fileExists(name, TRUE, FALSE) == FALSE) {
      missing++;
      continue;
    }

    missing = 0;

    sprintf(name, "%s/bucket%04d/sliceSizes", storePath, i);
    AS_UTL_unlink(name);

    sprintf(name, "%s/bucket%04d", storePath, i);
    rmdir(name);
  }

  fprintf(stderr, "Finished.\n");

  exit(0);
}

