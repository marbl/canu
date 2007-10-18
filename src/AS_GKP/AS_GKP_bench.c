
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

static char const *rcsid = "$Id: AS_GKP_bench.c,v 1.1 2007-10-18 07:44:34 brianwalenz Exp $";

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
  struct   statfs  sf = {0};
  struct   utsname un = {0};
  struct   rusage  ru = {0};

  errno = 0;
  if (statfs(gkpName, &sf) == -1)
    fprintf(stdout, "statfs() call failed: %s\n", strerror(errno));

  errno = 0;
  if (uname(&un) == -1)
    fprintf(stdout, "uname() call failed: %s\n", strerror(errno));

  errno = 0;
  if (getrusage(RUSAGE_SELF, &ru) == -1)
    fprintf(stdout, "getrusage() call failed: %s\n", strerror(errno));

  fprintf(stdout, "%s (%s|%s) ",
          gkpName,
          sf.f_fstypename,
          sf.f_mntfromname);
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
  GateKeeperStore          *gkp   = createGateKeeperStore(gkpName);
  GateKeeperFragmentRecord  gkf   = {0};
  char                      seq[AS_FRAG_MAX_LEN];
  char                      qlt[AS_FRAG_MAX_LEN];
  char                      enc[AS_FRAG_MAX_LEN];
  char                      acgt[10] = {'a', 'c', 'g', 't', 'n', 'A', 'C', 'G', 'T', 'N'};
  int                       i, j;

  for (i=0; i<numFrags; i++) {
    clearGateKeeperFragmentRecord(&gkf);

    gkf.readUID = getLastElemStore(gkp->frg) + 1 + 2000000000;
    gkf.readIID = getLastElemStore(gkp->frg) + 1;

    gkf.seqLen = 600 + lrand48() % 600;
    gkf.hpsLen = 0;
    gkf.srcLen = 0;

    //  Generate a bogus sequence and quality
    for (j=0; j<gkf.seqLen; j++) {
      seq[j] = acgt[lrand48() % 10];
      qlt[j] = '0' + lrand48() % 60;
    }

    {
      StoreStat   stats;

      statsStore(gkp->seq, &stats);
      gkf.seqOffset = stats.lastElem;

      statsStore(gkp->qlt, &stats);
      gkf.qltOffset = stats.lastElem;

      statsStore(gkp->hps, &stats);
      gkf.hpsOffset = stats.lastElem;

      statsStore(gkp->src, &stats);
      gkf.srcOffset = stats.lastElem;
    }

    setGatekeeperUIDtoIID(gkp, gkf.readUID, gkf.readIID, AS_IID_FRG);
    appendIndexStore(gkp->frg, &gkf);

    appendVLRecordStore(gkp->seq, seq, gkf.seqLen);

    encodeSequenceQuality(enc, seq, qlt);
    appendVLRecordStore(gkp->qlt, enc, gkf.seqLen);

    appendVLRecordStore(gkp->hps, NULL, gkf.hpsLen);
    appendVLRecordStore(gkp->src, NULL, gkf.srcLen);

    if ((i > 50000) && (i % 10000) == 0)
      fprintf(stderr, "%.3f ops/sec %.3f%% complete\n", i / (getTime() - startTime), 100.0 * i / numFrags);
  }

  fprintf(stderr, "%.3f ops/sec\n", i / (getTime() - startTime));

  closeGateKeeperStore(gkp);
}



static
void
updateRandomMates(char *gkpName, uint32 numMates) {
  GateKeeperStore           *gkp   = openGateKeeperStore(gkpName, TRUE);
  GateKeeperFragmentRecord   gkFrag1;
  GateKeeperFragmentRecord   gkFrag2;
  CDS_IID_t                  frag1IID;
  CDS_IID_t                  frag2IID;
  int                        i, lo, hi;
  int                        totalFrags = getLastElemStore(gkp->frg);

  for (i=0; i<numMates; i++) {
    frag1IID = (lrand48() % totalFrags) + 1;

    lo = 1;
    hi = totalFrags;

    if (frag1IID > 50000)
      lo = frag1IID - 50000;
    if (frag1IID <= totalFrags - 50000)
      hi = frag1IID + 50000;

    frag2IID = lo + (lrand48() % (hi - lo));

    getGateKeeperFragment(gkp, frag1IID, &gkFrag1);
    getGateKeeperFragment(gkp, frag2IID, &gkFrag2);

    gkFrag1.mateIID    = frag2IID;
    gkFrag2.mateIID    = frag1IID;

    setIndexStore(gkp->frg, frag1IID, &gkFrag1);
    setIndexStore(gkp->frg, frag2IID, &gkFrag2);

    if ((i > 50000) && (i % 10000) == 0)
      fprintf(stderr, "%.3f ops/sec %.3f%% complete\n", i / (getTime() - startTime), 100.0 * i / numMates);
  }

  fprintf(stderr, "%.3f ops/sec\n", i / (getTime() - startTime));

  closeGateKeeperStore(gkp);
}



static
void
readRandomFragments(char *gkpName, uint32 numReads) {
  GateKeeperStore           *gkp   = openGateKeeperStore(gkpName, FALSE);
  fragRecord                 frg;
  CDS_IID_t                  fragIID;
  int                        i;
  int                        totalFrags = getLastElemStore(gkp->frg);

  for (i=0; i<numReads; i++) {
    fragIID = (lrand48() % totalFrags) + 1;

    getFrag(gkp, fragIID, &frg, FRAG_S_ALL);

    if ((i > 50000) && (i % 10000) == 0)
      fprintf(stderr, "%.3f ops/sec %.3f%% complete\n", i / (getTime() - startTime), 100.0 * i / numReads);
  }

  fprintf(stderr, "%.3f ops/sec\n", i / (getTime() - startTime));

  closeGateKeeperStore(gkp);
}




int
main(int argc, char **argv) {
  char      gkpName[FILENAME_MAX] = {0};
  uint32    numFrags   = 0;  //  Create a store with numFrags bogus frags in it
  uint32    numMates   = 0;  //  Add mates to random frags
  uint32    numReads   = 0;  //  Read random frags

  srand48(time(NULL));

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
