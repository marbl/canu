
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
 *    Brian P. Walenz from 2015-JAN-26 to 2015-JUL-01
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
#include "tgStore.H"


int
main(int argc, char **argv) {
  tgTig  tig;
  char  *gkpName     = NULL;
  char  *tigFileName = NULL;

  argc = AS_configure(argc, argv);

  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-G") == 0) {
      gkpName = argv[++arg];

    } else if (strcmp(argv[arg], "-t") == 0) {
      tigFileName = argv[++arg];

    } else {
      err++;
    }

    arg++;
  }
  if (gkpName == NULL)
    err++;
  if (tigFileName == NULL)
    err++;
  if (err) {
    fprintf(stderr, "usage: %s -G gkpStore -t tigFile\n", argv[0]);
    exit(1);
  }

  gkStore  *gkpStore = gkStore::gkStore_open(gkpName);

  FILE *F = fopen(tigFileName, "r");

  tig.loadFromStream(F);

  fclose(F);

  uint32  displayWidth    = 250;
  uint32  displaySpacing  = 10;
  bool    withQV          = false;
  bool    withDots        = true;

  tig.display(stdout, gkpStore, displayWidth, displaySpacing, withQV, withDots);

  exit(0);
}
