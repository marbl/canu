
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
 *    src/AS_BOG/analyzeBest.C
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2010-OCT-09 to 2013-AUG-01
 *      are Copyright 2010-2011,2013 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2014-OCT-09 to 2015-APR-10
 *      are Copyright 2014-2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2015-DEC-07
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_global.H"
#include "AS_PER_gkpStore.H"
#include "splitToWords.H"


//  Read the 'best.edges' and 'best.contains' outputs from BOG, compute
//  how many singletons, spurs and contains are present per library.
//
//  Assumes it is run from 4-unitigger.

int
main(int argc, char **argv) {
  char  *gkpName = 0L;
  char  *bEdge   = "best.edges";
  char  *bCont   = "best.contains";
  char  *bSing   = "best.singletons";

  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-g") == 0) {
      gkpName = argv[++arg];

    } else if (strcmp(argv[arg], "-e") == 0) {
      bEdge = argv[++arg];

    } else if (strcmp(argv[arg], "-c") == 0) {
      bCont = argv[++arg];

    } else if (strcmp(argv[arg], "-s") == 0) {
      bSing = argv[++arg];

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
  FILE *be = fopen(bEdge, "r");
  if (errno)
    fprintf(stderr, "Failed to open '%s' for reading: %s\n", bEdge, strerror(errno)), exit(1);

  FILE *bc = fopen(bCont,  "r");
  if (errno)
    fprintf(stderr, "Failed to open '%s' for reading: %s\n", bCont, strerror(errno)), exit(1);

  FILE *bs = fopen(bSing,  "r");
  if (errno)
    fprintf(stderr, "Failed to open '%s' for reading: %s\n", bSing, strerror(errno)), exit(1);

  fprintf(stderr, "Loading fragment to library mapping.\n");

  gkStore    *gkp = gkStore::gkStore_open(gkpName, false, false);
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

  uint64      *cntdPerLib = new uint64 [numLib + 1];  memset(cntdPerLib, 0, sizeof(uint64) * (numLib + 1));
  uint64      *cntrPerLib = new uint64 [numLib + 1];  memset(cntrPerLib, 0, sizeof(uint64) * (numLib + 1));
  uint32      *cntr       = new uint32 [numFrg + 1];  memset(cntr,       0, sizeof(uint32) * (numFrg + 1));

  uint64      *singPerLib = new uint64 [numLib + 1];  memset(singPerLib, 0, sizeof(uint64) * (numLib + 1));
  uint64      *spu5PerLib = new uint64 [numLib + 1];  memset(spu5PerLib, 0, sizeof(uint64) * (numLib + 1));
  uint64      *spu3PerLib = new uint64 [numLib + 1];  memset(spu3PerLib, 0, sizeof(uint64) * (numLib + 1));
  uint64      *dovePerLib = new uint64 [numLib + 1];  memset(dovePerLib, 0, sizeof(uint64) * (numLib + 1));

  int32          lineMax = 1048576;
  char          *line    = new char [lineMax];
  splitToWords   W;


  fprintf(stderr, "Processing best.singletons.\n");

  while (!feof(bs)) {
    fgets(line, lineMax, bs);  chomp(line);
    W.split(line);

    if (line[0] == '#')
      continue;

    singPerLib[frgToLib[W(0)]]++;
  }


  fprintf(stderr, "Processing best.contains.\n");

  while (!feof(bc)) {
    fgets(line, lineMax, bc);  chomp(line);
    W.split(line);

    if (line[0] == '#')
      continue;

    cntdPerLib[frgToLib[W(0)]]++;

    cntr[W(3)]++;
  }

  for (uint32 i=1; i<=numFrg; i++) {
    if (cntr[i] > 0)
      cntrPerLib[frgToLib[i]]++;
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

    fprintf(stderr, "%-8u %-30s %8"F_U64P" %8"F_U64P" (%4.1f%%)  %8"F_U64P" (%4.1f%%)  %8"F_U64P" (%4.1f%%)  %8"F_U64P" (%4.1f%%)  %8"F_U64P" (%4.1f%%)  %8"F_U64P" (%4.1f%%)  %8"F_U64P" (%4.1f%%)\n",
            i,
            gkp->gkStore_getLibrary(i)->libraryName,
            fragPerLib[i],
            deldPerLib[i], (deldPerLib[i] == 0) ? 0.0 : 100.0 * deldPerLib[i] / tot,
            cntdPerLib[i], (cntdPerLib[i] == 0) ? 0.0 : 100.0 * cntdPerLib[i] / tot,
            cntrPerLib[i], (cntrPerLib[i] == 0) ? 0.0 : 100.0 * cntrPerLib[i] / tot,
            singPerLib[i], (singPerLib[i] == 0) ? 0.0 : 100.0 * singPerLib[i] / tot,
            spu5PerLib[i], (spu5PerLib[i] == 0) ? 0.0 : 100.0 * spu5PerLib[i] / tot,
            spu3PerLib[i], (spu3PerLib[i] == 0) ? 0.0 : 100.0 * spu3PerLib[i] / tot,
            dovePerLib[i], (dovePerLib[i] == 0) ? 0.0 : 100.0 * dovePerLib[i] / tot);
  }

  gkp->gkStore_close();

  fclose(be);
  fclose(bc);

  delete [] frgToLib;

  delete [] fragPerLib;
  delete [] deldPerLib;

  delete [] cntdPerLib;
  delete [] cntrPerLib;
  delete [] cntr;

  delete [] singPerLib;
  delete [] spu5PerLib;
  delete [] spu3PerLib;
  delete [] dovePerLib;

  delete [] line;

  return(0);
}

