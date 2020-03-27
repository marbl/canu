
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' r4587 (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' r1994 (http://kmer.sourceforge.net)
 *
 *  Except as indicated otherwise, this is a 'United States Government Work',
 *  and is released in the public domain.
 *
 *  File 'README.licenses' in the root directory of this distribution
 *  contains full conditions and disclaimers.
 */

#include "runtime.H"
#include "sqStore.H"

#include "files.H"
#include "strings.H"




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
      fprintf(stdout, "START %s pos %8" F_U64P " max %8" F_SIZE_TP " length %8" F_U32P "\n", chunk, blobPos, blobMap->length(), chunkLen);
    }

    else {
      blob    += chunkLen;
      blobPos += chunkLen;

      fprintf(stdout, "      %s pos %8" F_U64P " max %8" F_SIZE_TP " length %8" F_U32P "\n", chunk, blobPos, blobMap->length(), chunkLen);
    }
  }

  delete blobMap;

  exit(0);
}
