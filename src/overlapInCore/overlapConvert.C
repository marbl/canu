
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
 *  Modifications by:
 *
 *    Brian P. Walenz from 2015-FEB-11 to 2015-JUN-25
 *      are Copyright 2015 Battelle National Biodefense Institute, and
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
#include "gkStore.H"
#include "ovStore.H"

#include <vector>

using namespace std;


int
main(int argc, char **argv) {
  char                  *gkpStoreName = NULL;
  gkStore               *gkpStore = NULL;

  ovOverlapDisplayType   dt = ovOverlapAsCoords;

  vector<char *>         files;


  int32     arg = 1;
  int32     err = 0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-G") == 0) {
      gkpStoreName = argv[++arg];

    } else if (strcmp(argv[arg], "-coords") == 0) {
      dt = ovOverlapAsCoords;

    } else if (strcmp(argv[arg], "-hangs") == 0) {
      dt = ovOverlapAsHangs;

    } else if (strcmp(argv[arg], "-raw") == 0) {
      dt = ovOverlapAsRaw;

    } else if (AS_UTL_fileExists(argv[arg])) {
      files.push_back(argv[arg]);

    } else {
      fprintf(stderr, "ERROR:  invalid arg '%s'\n", argv[arg]);
      err++;
    }

    arg++;
  }

  if ((gkpStoreName == NULL) && (dt == ovOverlapAsCoords))
    err++;

  if ((err) || (files.size() == 0)) {
    fprintf(stderr, "usage: %s [options] file.ovb[.gz]\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "  -G             gkpStore (needed for -coords, the default)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -coords        output coordiantes on reads\n");
    fprintf(stderr, "  -hangs         output hangs on reads\n");
    fprintf(stderr, "  -raw           output raw hangs on reads\n");

    if ((gkpStoreName == NULL) && (dt == ovOverlapAsCoords))
      fprintf(stderr, "ERROR:  -coords mode requires a gkpStore (-G)\n");

    if (files.size() == 0)
      fprintf(stderr, "ERROR:  no overlap files supplied\n");

    exit(1);
  }

  if (gkpStoreName)
    gkpStore = gkStore::gkStore_open(gkpStoreName);

  char  *ovStr = new char [1024];

  for (uint32 ff=0; ff<files.size(); ff++) {
    ovFile      *of = new ovFile(files[ff], ovFileFull);
    ovOverlap   ov(gkpStore);

    while (of->readOverlap(&ov))
      fputs(ov.toString(ovStr, dt, true), stdout);

    delete of;
  }

  delete [] ovStr;

  gkpStore->gkStore_close();

  exit(0);
}
