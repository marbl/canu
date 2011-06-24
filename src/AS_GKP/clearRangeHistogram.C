
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 2011, J. Craig Venter Institute.
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

const char *mainid = "$Id: clearRangeHistogram.C,v 1.1 2011-06-24 12:44:35 brianwalenz Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fcntl.h>
#include <string.h>
#include <unistd.h>
#include <ctype.h>

#include "AS_global.h"
#include "AS_PER_gkpStore.h"




int
main(int argc, char **argv) {
  char     *gkpName   = NULL;
  gkStore  *gkpStore  = NULL;

  uint32    clearRegion = AS_READ_CLEAR_LATEST;

  char     *outPrefix = NULL;
  char      outName[FILENAME_MAX];

  uint32    lib[1024] = { 1, 0 };
  uint32   *bgn[1024] = { NULL };
  uint32   *end[1024] = { NULL };
  uint32    max       = 0;

  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-g") == 0) {
      gkpName = argv[++arg];

    } else if (strcmp(argv[arg], "-c") == 0) {
      clearRegion = gkStore_decodeClearRegionLabel(argv[++arg]);

    } else if (strcmp(argv[arg], "-l") == 0) {
      lib[0] = 0;
      lib[atoi(argv[++arg])]++;

    } else if (strcmp(argv[arg], "-o") == 0) {
      outPrefix = argv[++arg];
    } else {
      err++;
    }

    arg++;
  }
  if (gkpName == NULL)
    err++;
  if (outPrefix == NULL)
    err++;
  if (err) {
    fprintf(stderr, "usage: %s -g gkpStore -l lib [-l lib] -o output-prefix\n", argv[0]);
    if (gkpName == NULL)
      fprintf(stderr, "ERROR:  No gkpStore supplied (-g).\n");
    if (outPrefix == NULL)
      fprintf(stderr, "ERROR:  No output-prefix supplied (-o).\n");
    exit(1);
  }

  gkpStore = new gkStore(gkpName, FALSE, FALSE, true);

  uint32     numLibs  = gkpStore->gkStore_getNumLibraries();
  uint32     numFrags = gkpStore->gkStore_getNumFragments();

  if (lib[0] == 1)
    for (uint32 l=1; l<=numLibs; l++)
      lib[l] = 1;

  for (uint32 l=1; l<=numLibs; l++) {
    if (lib[l] == 0)
      continue;

    bgn[l] = new uint32 [AS_READ_MAX_NORMAL_LEN];
    end[l] = new uint32 [AS_READ_MAX_NORMAL_LEN];

    memset(bgn[l], 0, sizeof(uint32) * AS_READ_MAX_NORMAL_LEN);
    memset(end[l], 0, sizeof(uint32) * AS_READ_MAX_NORMAL_LEN);
  }

  gkFragment  fr;
  gkStream   *fs = new gkStream(gkpStore, 0, 0, GKFRAGMENT_INF);

  while (fs->next(&fr) == true) {
    uint32  l = fr.gkFragment_getLibraryIID();

    if (lib[l] == 0)
      continue;

    uint32  b = fr.gkFragment_getClearRegionBegin(clearRegion);
    uint32  e = fr.gkFragment_getClearRegionEnd(clearRegion);

    if (max < b)
      max = b+1;
    if (max < e)
      max = e+1;

    bgn[l][b]++;
    end[l][e]++;
  }


  sprintf(outName, "%s.bgn.dat", outPrefix);
  FILE *datbgn = fopen(outName, "w");

  sprintf(outName, "%s.end.dat", outPrefix);
  FILE *datend = fopen(outName, "w");

  uint32  numLibsActual = 0;
  for (uint32 l=1; l<=numLibs; l++)
    numLibsActual += lib[l];


  for (uint32 p=0; p<=max; p++) {
    fprintf(datbgn, "%u", p);
    fprintf(datend, "%u", p);

    for (uint32 l=1; l<=numLibs; l++) {
      if (lib[l] == 0)
        continue;

      fprintf(datbgn, "\t%u", bgn[l][p]);
      fprintf(datend, "\t%u", end[l][p]);
    }

    fprintf(datbgn, "\n");
    fprintf(datend, "\n");
  }

  fclose(datbgn);
  fclose(datend);



  sprintf(outName, "%s.bgn.gp", outPrefix);
  FILE *gpbgn = fopen(outName, "w");

  fprintf(gpbgn, "set terminal png\n");
  fprintf(gpbgn, "set output \"%s.bgn.png\"\n", outPrefix);
  fprintf(gpbgn, "\n");
  fprintf(gpbgn, "set title \"%s begin trim point\"\n", AS_READ_CLEAR_NAMES[clearRegion]);
  fprintf(gpbgn, "\n");
  fprintf(gpbgn, "set logscale y\n");
  fprintf(gpbgn, "\n");
  fprintf(gpbgn, "plot [-10:] \\\n");

  for (uint32 l=1, a=1; l<=numLibs; l++) {
    if (lib[l] == 0)
      continue;

    if (++a <= numLibsActual)
      fprintf(gpbgn, "     \"%s.bgn.dat\" using 1:%u with lines title \"%u\", \\\n", outPrefix, a, l);
    else
      fprintf(gpbgn, "     \"%s.bgn.dat\" using 1:%u with lines title \"%u\"\n", outPrefix, a, l);
  }

  fclose(gpbgn);



  sprintf(outName, "%s.end.gp", outPrefix);
  FILE *gpend = fopen(outName, "w");

  fprintf(gpend, "set terminal png\n");
  fprintf(gpend, "set output \"%s.end.png\"\n", outPrefix);
  fprintf(gpend, "\n");
  fprintf(gpend, "set title \"%s end trim point\"\n", AS_READ_CLEAR_NAMES[clearRegion]);
  fprintf(gpend, "\n");
  fprintf(gpbgn, "set logscale y\n");
  fprintf(gpend, "\n");
  fprintf(gpend, "plot [-10:] \\\n");

  for (uint32 l=1, a=1; l<=numLibs; l++) {
    if (lib[l] == 0)
      continue;

    if (++a <= numLibsActual)
      fprintf(gpend, "     \"%s.end.dat\" using 1:%u with lines title \"%u\", \\\n", outPrefix, a, l);
    else
      fprintf(gpend, "     \"%s.end.dat\" using 1:%u with lines title \"%u\"\n", outPrefix, a, l);
  }

  fclose(gpend);

  sprintf(outName, "gnuplot < %s.bgn.gp", outPrefix);
  system(outName);

  sprintf(outName, "gnuplot < %s.end.gp", outPrefix);
  system(outName);

  exit(0);
}
