
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 2010, The Venter Institute. All rights reserved.
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

const char *mainid = "$Id: analyzeBest.C,v 1.1 2010-10-09 14:05:54 brianwalenz Exp $";

#include "AS_global.h"
#include "AS_PER_gkpStore.h"
#include "AS_UTL_splitToWords.H"

//  Read the 'best.edges' and 'best.contains' outputs from BOG, compute
//  how many singletons, spurs and contains are present per library.
//
//  Assumes it is run from 4-unitigger.

int
main(int argc, char **argv) {
  char  *gkpName = 0L;
  char  *bEdges  = "best.edges";
  char  *bCont   = "best.contains";

  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-g") == 0) {
      gkpName = argv[++arg];

    } else if (strcmp(argv[arg], "-e") == 0) {
      bEdges = argv[++arg];

    } else if (strcmp(argv[arg], "-c") == 0) {
      bCont = argv[++arg];

    } else {
      err++;
    }

    arg++;
  }
  if ((err) || (gkpName == 0L)) {
    fprintf(stderr, "usage: %s -g gkpName [-e best.edges] [-c best.contains]\n", argv[0]);
    exit(1);
  }

  fprintf(stderr, "Opening best edges and contains.\n");

  errno = 0;
  FILE *be = fopen(bEdges, "r");
  if (errno)
    fprintf(stderr, "Failed to open '%s' for reading: %s\n", bEdges, strerror(errno)), exit(1);

  FILE *bc = fopen(bCont,  "r");
  if (errno)
    fprintf(stderr, "Failed to open '%s' for reading: %s\n", bCont, strerror(errno)), exit(1);

  fprintf(stderr, "Loading fragment to library mapping.\n");

  gkStore    *gkp = new gkStore(gkpName, false, false);
  gkStream   *str = new gkStream(gkp, 0, 0, GKFRAGMENT_INF);
  gkFragment  fr;

  uint32   numFrg = gkp->gkStore_getNumFragments();
  uint32   numLib = gkp->gkStore_getNumLibraries();

  uint32  *frgToLib  = new uint32 [numFrg + 1];  memset(frgToLib, 0, sizeof(uint32) * (numFrg + 1));

  uint64  *fragPerLib = new uint64 [numLib + 1];  memset(fragPerLib, 0, sizeof(uint64) * (numLib + 1));
  uint64  *deldPerLib = new uint64 [numLib + 1];  memset(deldPerLib, 0, sizeof(uint64) * (numLib + 1));

  while (str->next(&fr)) {
    frgToLib[fr.gkFragment_getReadIID()] = fr.gkFragment_getLibraryIID();

    if (fr.gkFragment_getIsDeleted())
      deldPerLib[fr.gkFragment_getLibraryIID()]++;
    else
      fragPerLib[fr.gkFragment_getLibraryIID()]++;
  }

  delete str;

  uint64  *cntdPerLib = new uint64 [numLib + 1];  memset(cntdPerLib, 0, sizeof(uint64) * (numLib + 1));
  uint64  *cntrPerLib = new uint64 [numLib + 1];  memset(cntrPerLib, 0, sizeof(uint64) * (numLib + 1));

  uint64  *singPerLib = new uint64 [numLib + 1];  memset(singPerLib, 0, sizeof(uint64) * (numLib + 1));
  uint64  *spu5PerLib = new uint64 [numLib + 1];  memset(spu5PerLib, 0, sizeof(uint64) * (numLib + 1));
  uint64  *spu3PerLib = new uint64 [numLib + 1];  memset(spu3PerLib, 0, sizeof(uint64) * (numLib + 1));
  uint64  *dovePerLib = new uint64 [numLib + 1];  memset(dovePerLib, 0, sizeof(uint64) * (numLib + 1));

  int32          lineMax = 1048576;
  char          *line    = new char [lineMax];
  splitToWords   W;

  fprintf(stderr, "Processing best.contains.\n");

  while (!feof(bc)) {
    fgets(line, lineMax, bc);  chomp(line);
    W.split(line);

    if (line[0] == '#')
      continue;

    cntdPerLib[frgToLib[W(0)]]++;
    cntrPerLib[frgToLib[W(3)]]++;
  }


  fprintf(stderr, "Processing best.edges.\n");

  while (!feof(be)) {
    fgets(line, lineMax, be);  chomp(line);
    W.split(line);

    if (line[0] == '#')
      continue;

    if      ((W(2) == 0) && (W(4) == 0))
      singPerLib[frgToLib[W(0)]]++;
    else if (W(2) == 0)
      spu5PerLib[frgToLib[W(0)]]++;
    else if (W(4) == 0)
      spu3PerLib[frgToLib[W(0)]]++;
    else
      dovePerLib[frgToLib[W(0)]]++;
  }

  fprintf(stderr, "libIID   libUID                             #frg     #del             cnt'd             cnt'r              sing             spur5             spur3              dove\n");

  for (uint32 i=0; i<numLib+1; i++) {
    double tot = fragPerLib[i] + deldPerLib[i];

    fprintf(stderr, "%-8u %-30s %8lu %8lu (%4.1f%%)  %8lu (%4.1f%%)  %8lu (%4.1f%%)  %8lu (%4.1f%%)  %8lu (%4.1f%%)  %8lu (%4.1f%%)  %8lu (%4.1f%%)\n",
            i,
            (i == 0) ? "(none)" : AS_UID_toString(gkp->gkStore_getLibrary(i)->libraryUID),
            fragPerLib[i],
            deldPerLib[i], (deldPerLib[i] == 0) ? 0.0 : 100.0 * deldPerLib[i] / tot,
            cntdPerLib[i], (cntdPerLib[i] == 0) ? 0.0 : 100.0 * cntdPerLib[i] / tot,
            cntrPerLib[i], (cntrPerLib[i] == 0) ? 0.0 : 100.0 * cntrPerLib[i] / tot,
            singPerLib[i], (singPerLib[i] == 0) ? 0.0 : 100.0 * singPerLib[i] / tot,
            spu5PerLib[i], (spu5PerLib[i] == 0) ? 0.0 : 100.0 * spu5PerLib[i] / tot,
            spu3PerLib[i], (spu3PerLib[i] == 0) ? 0.0 : 100.0 * spu3PerLib[i] / tot,
            dovePerLib[i], (dovePerLib[i] == 0) ? 0.0 : 100.0 * dovePerLib[i] / tot);
  }

  delete gkp;
  fclose(be);
  fclose(bc);

  delete [] frgToLib;
  delete [] fragPerLib;
  delete [] deldPerLib;

  delete [] cntdPerLib;
  delete [] cntrPerLib;

  delete [] singPerLib;
  delete [] spu5PerLib;
  delete [] spu3PerLib;
  delete [] dovePerLib;

  delete [] line;

  return(0);
}

