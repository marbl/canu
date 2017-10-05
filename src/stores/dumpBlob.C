
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
 *    src/AS_GKP/gkpStoreDumpFASTQ.C
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2012-FEB-06 to 2013-AUG-01
 *      are Copyright 2012-2013 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2014-NOV-21 to 2015-JUN-03
 *      are Copyright 2014-2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2015-NOV-23
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_global.H"
#include "gkStore.H"

#include "AS_UTL_decodeRange.H"
#include "AS_UTL_fileIO.H"

#include "memoryMappedFile.H"



int
main(int argc, char **argv) {
  char            *blobName  = NULL;
  uint64           offset    = 0;

  argc = AS_configure(argc, argv);

  int arg = 1;
  int err = 0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-b") == 0) {
      blobName = argv[++arg];

    } else if (strcmp(argv[arg], "-o") == 0) {
      offset = strtouint64(argv[++arg]);

    } else {
      err++;
      fprintf(stderr, "ERROR: unknown option '%s'\n", argv[arg]);
    }
    arg++;
  }

  if (blobName == NULL)
    err++;
  if (err) {
    fprintf(stderr, "usage: %s -b blobFile\n", argv[0]);
    fprintf(stderr, "\n");

    if (blobName == NULL)
      fprintf(stderr, "ERROR: no blob file (-b) supplied.\n");

    exit(1);
  }

  memoryMappedFile  *blobMap = new memoryMappedFile(blobName, memoryMappedFile_readOnly);
  uint8             *blob    = (uint8 *)blobMap->get(0) + offset;
  uint8             *blobMax = (uint8 *)blobMap->get(0) + blobMap->length();
  uint64             blobPos = offset;

  char               chunk[5] = {0};
  uint32             chunkLen = 0;

  while (blob < blobMax) {
    chunk[0] = blob[0];
    chunk[1] = blob[1];
    chunk[2] = blob[2];
    chunk[3] = blob[3];
    chunk[4] = 0;

    chunkLen = *((uint32 *)blob + 1); 
    blob    += 8;
    blobPos += 8;

    if ((chunk[0] == 'B') &&
        (chunk[1] == 'L') &&
        (chunk[2] == 'O') &&
        (chunk[3] == 'B')) { 
      fprintf(stdout, "START %s pos %8lu max %8lu length %8u\n", chunk, blobPos, blobMap->length(), chunkLen);
    }

    else {
      blob    += chunkLen;
      blobPos += chunkLen;

      fprintf(stdout, "      %s pos %8lu max %8lu length %8u\n", chunk, blobPos, blobMap->length(), chunkLen);
    }
  }

  delete blobMap;

  exit(0);
}
