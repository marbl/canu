
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 1999-2004, Applera Corporation. All rights reserved.
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

const char *mainid = "$Id: AS_GKP_bench.c,v 1.13 2011-12-29 06:11:02 brianwalenz Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fcntl.h>
#include <string.h>
#include <unistd.h>

#include <sys/types.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <sys/utsname.h>
#include <sys/param.h>
#include <sys/mount.h>

#include "AS_global.h"
#include "AS_PER_genericStore.h"
#include "AS_PER_gkpStore.h"
#include "AS_PER_encodeSequenceQuality.h"


double  startTime = 0.0;


static
double
getTime(void) {
  struct timeval  tp;
  gettimeofday(&tp, NULL);
  return(tp.tv_sec + (double)tp.tv_usec / 1000000.0);
}


static
void
printrusage(char *gkpName, double startTime) {
#if defined(__FreeBSD__) || defined(__APPLE__)
  struct   statfs  sf;
#endif
  struct   utsname un;
  struct   rusage  ru;

#if defined(__FreeBSD__) || defined(__APPLE__)
  errno = 0;
  if (statfs(gkpName, &sf) == -1)
    fprintf(stdout, "statfs() call failed: %s\n", strerror(errno));
#endif

  errno = 0;
  if (uname(&un) == -1)
    fprintf(stdout, "uname() call failed: %s\n", strerror(errno));

  errno = 0;
  if (getrusage(RUSAGE_SELF, &ru) == -1)
    fprintf(stdout, "getrusage() call failed: %s\n", strerror(errno));

#if defined(__FreeBSD__) || defined(__APPLE__)
  fprintf(stdout, "%s (%s|%s) ",
          gkpName,
          sf.f_fstypename,
          sf.f_mntfromname);
#else
  fprintf(stdout, "%s ",
          gkpName);
#endif
  fprintf(stdout, "%s (%s/%s %s) ",
          un.nodename,
          un.sysname,
          un.machine,
          un.release);
  fprintf(stdout, "%.3fu %.3fs %.3fw\n",
          ru.ru_utime.tv_sec + (double)ru.ru_utime.tv_usec / 1000000.0,
          ru.ru_stime.tv_sec + (double)ru.ru_stime.tv_usec / 1000000.0,
          getTime() - startTime);
}


static
void
addRandomFrags(char *gkpName, uint32 numFrags) {
  gkStore                  *gkp = new gkStore(gkpName, TRUE, TRUE);
  gkFragment                frag;
  char                      acgt[10] = {'a', 'c', 'g', 't', 'n', 'A', 'C', 'G', 'T', 'N'};
  int                       len;
  int                       i, j;
  AS_UID                    uid;

  //  This is a special case for gatekeeper; we never call
  //  gkStore_getFragment() and so we never set up the gkFragment.
  //
  frag.gkFragment_enableGatekeeperMode(gkp);

  char *seq = frag.gkFragment_getSequence();
  char *qlt = frag.gkFragment_getQuality();

  for (i=0; i<numFrags; i++) {
    uid.isString = 0;
    uid.UID      = i + 1 + 2000000000;

    frag.gkFragment_setReadUID(uid);
    frag.gkFragment_setType(GKFRAGMENT_NORMAL);

    len = 600 + lrand48() % 600;

    //  Generate a bogus sequence and quality
    for (j=0; j<len; j++) {
      seq[j] = acgt[lrand48() % 10];
      qlt[j] = '0' + lrand48() % 60;
    }
    seq[len] = 0;
    qlt[len] = 0;

    frag.gkFragment_setLength(len);

    if ((i > 50000) && (i % 10000) == 0)
      fprintf(stderr, "%.3f ops/sec %.3f%% complete\n", i / (getTime() - startTime), 100.0 * i / numFrags);

    frag.clrBgn = frag.vecBgn = 0;
    frag.clrEnd = frag.vecEnd = len;

    frag.maxBgn = frag.tntBgn = 1;
    frag.maxEnd = frag.tntEnd = 0;

    gkp->gkStore_addFragment(&frag);
  }

  fprintf(stderr, "%.3f ops/sec\n", i / (getTime() - startTime));

  delete gkp;
}



static
void
updateRandomMates(char *gkpName, uint32 numMates) {
  gkStore           *gkp   = new gkStore(gkpName, FALSE, TRUE);
  gkFragment         frag1;
  gkFragment         frag2;
  AS_IID             frag1IID;
  AS_IID             frag2IID;
  int                i, lo, hi;
  int                totalFrags = gkp->gkStore_getNumFragments();

  frag1.gkFragment_enableGatekeeperMode(gkp);
  frag2.gkFragment_enableGatekeeperMode(gkp);

  for (i=0; i<numMates; i++) {
    frag1IID = (lrand48() % totalFrags) + 1;

    lo = 1;
    hi = totalFrags;

    if (frag1IID > 50000)
      lo = frag1IID - 50000;
    if (frag1IID <= totalFrags - 50000)
      hi = frag1IID + 50000;

    frag2IID = lo + (lrand48() % (hi - lo));

    gkp->gkStore_getFragment(frag1IID, &frag1, GKFRAGMENT_INF);
    gkp->gkStore_getFragment(frag2IID, &frag2, GKFRAGMENT_INF);

    frag1.gkFragment_setMateIID(frag2IID);
    frag2.gkFragment_setMateIID(frag1IID);

    gkp->gkStore_setFragment(&frag1);
    gkp->gkStore_setFragment(&frag2);

    if ((i > 50000) && (i % 10000) == 0)
      fprintf(stderr, "%.3f ops/sec %.3f%% complete\n", i / (getTime() - startTime), 100.0 * i / numMates);
  }

  fprintf(stderr, "%.3f ops/sec\n", i / (getTime() - startTime));

  delete gkp;
}



static
void
readRandomFragments(char *gkpName, uint32 numReads) {
  gkStore                   *gkp   = new gkStore(gkpName, FALSE, FALSE);
  gkFragment                 frg;
  AS_IID                     fragIID;
  int                        i;
  int                        totalFrags = gkp->gkStore_getNumFragments();

  for (i=0; i<numReads; i++) {
    fragIID = (lrand48() % totalFrags) + 1;

    gkp->gkStore_getFragment(fragIID, &frg, GKFRAGMENT_QLT);

    if ((i > 50000) && (i % 10000) == 0)
      fprintf(stderr, "%.3f ops/sec %.3f%% complete\n", i / (getTime() - startTime), 100.0 * i / numReads);
  }

  fprintf(stderr, "%.3f ops/sec\n", i / (getTime() - startTime));

  delete gkp;
}




int
main(int argc, char **argv) {
  char      gkpName[FILENAME_MAX] = {0};
  uint32    numFrags   = 0;  //  Create a store with numFrags bogus frags in it
  uint32    numMates   = 0;  //  Add mates to random frags
  uint32    numReads   = 0;  //  Read random frags

  srand48(time(NULL));

  argc = AS_configure(argc, argv);

  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-g") == 0) {
      strcpy(gkpName, argv[++arg]);
    } else if (strcmp(argv[arg], "-create") == 0) {
      numFrags = atoi(argv[++arg]);
    } else if (strcmp(argv[arg], "-mates") == 0) {
      numMates = atoi(argv[++arg]);
    } else if (strcmp(argv[arg], "-reads") == 0) {
      numReads = atoi(argv[++arg]);
    } else if (strcmp(argv[arg], "-seed") == 0) {
      srand48(atoi(argv[++arg]));
    } else {
      err++;
    }
    arg++;
  }
  if (((numFrags > 0) + (numMates > 0) + (numReads > 0)) != 1) {
    fprintf(stderr, "Exactly one of -n, -m and -r must be supplied.\n\n");
  }
  if ((err) || (gkpName[0] == 0)) {
    fprintf(stderr, "usage: %s -g gkpStoreName [opts]\n", argv[0]);
    fprintf(stderr, "  -g      gkpStoreName    create/read/write the store called 'gkpStoreName'\n");
    fprintf(stderr, "  -seed   s               use random seed s\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -create numFrags        add numFrags random fragments\n");
    fprintf(stderr, "  -mates  numMates        update numMates random mated fragments\n");
    fprintf(stderr, "  -reads  numReads        read numReads random fragments\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "-n is not a very useful benchmark.  It is somewhat CPU bound, and simply writes\n");
    fprintf(stderr, "sequentially to a handful of files.  This isn't the primary task of this benchmark,\n");
    fprintf(stderr, "we just need to create the files somehow.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "-m is possibly the most brutal test.  It reads and writes randomly to a moderately\n");
    fprintf(stderr, "large file.  Record size is 104 bytes.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "-r is (presumed to be) the majority of accesses made by the assembler.  It reads a\n");
    fprintf(stderr, "random fragment from the store.  It reads the 104 byte record from one file, and\n");
    fprintf(stderr, "a variable length (800 to 1200 bytes) sequence from a larger file.\n");
    fprintf(stderr, "\n");
    exit(1);
  }

  startTime = getTime();

  if (numFrags > 0) {
    addRandomFrags(gkpName, numFrags);
  }

  if (numMates > 0) {
    updateRandomMates(gkpName, numMates);
  }

  if (numReads > 0) {
    readRandomFragments(gkpName, numReads);
  }

  printrusage(gkpName, startTime);

  exit(0);
}
