
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//  A quick hack to compute a histogram of coverage depth using
//  the runCA-OBT posmap files.

#define DEPTHSIZE   (128 * 1024 * 1024)

int
main(int argc, char **argv) {
  int              i = 0;

  unsigned short  *depth = NULL;

  unsigned long    uid = 0;
  int              beg = 0;
  int              end = 0;

  unsigned long    lastuid = 0;
  int              lastend = 0;

  int              histogram[256] = { 0 };
  int              histmax = 0;

  int              minSize = 0;
  int              maxSize = DEPTHSIZE;

  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-min") == 0) {
      minSize = atoi(argv[++arg]);
    } else if (strcmp(argv[arg], "-max") == 0) {
      maxSize = atoi(argv[++arg]);
    } else {
      fprintf(stderr, "unknown option '%s'\n", argv[arg]);
      err++;
    }
    arg++;
  }
  if (err) {
    fprintf(stderr, "usage: %s [-min N] [-max N] < x.posmap.frgscf", argv[0]);
    fprintf(stderr, "  -min N     count scaffolds at least N bases long.\n");
    fprintf(stderr, "  -max N     count scaffolds at most N bases long.\n");
  }

  depth = (unsigned short *)malloc(sizeof(unsigned short) * DEPTHSIZE);

  while (3 == fscanf(stdin, " %*ld %ld %d %d %*d ", &uid, &beg, &end)) {
    if (uid != lastuid) {
      if ((minSize <= lastend) && (lastend <= maxSize)) {
        for (i=0; i<lastend; i++) {
          if (histmax < depth[i])
            histmax = depth[i];
          histogram[depth[i]]++;
          depth[i] = 0;
        }
      }

      lastuid = uid;
      lastend = 0;
    }

    for (i=beg; i<end; i++)
      depth[i]++;

    if (lastend < end)
      lastend = end;
  }

  for (i=0; i<=histmax; i++)
    fprintf(stdout, "%d\t%d\n", i, histogram[i]);

  exit(0);
}
