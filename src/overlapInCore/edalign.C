
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
 *    Brian P. Walenz beginning on 2018-JAN-14
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_global.H"
#include "files.H"

#include "strings.H"

#include "edlib.H"


bool
readLine(FILE *file, char *line, int32 lineMax, int32 &len, splitToWords &s) {

  fgets(line, lineMax, file);

  if (feof(file))
    return(false);

  chomp(line);

  s.split(line);

  len = strlen(s[1]);

  return(true);
}



int
main(int argc, char **argv) {
  char    *nameA           = NULL;
  char    *nameB           = NULL;

  argc = AS_configure(argc, argv);

  int err=0;
  int arg=1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-a") == 0) {
      nameA = argv[++arg];

    } else if (strcmp(argv[arg], "-b") == 0) {
      nameB = argv[++arg];

    } else {
      err++;
    }

    arg++;
  }

  if (nameA == NULL)
    err++;
  if (nameB == NULL)
    err++;

  if (err) {
    fprintf(stderr, "usage: %s -a fileA -b fileB ...\n", argv[0]);
    fprintf(stderr, "  -a fileA     Mandatory, path to first input file\n");
    fprintf(stderr, "  -b fileB     Mandatory, path to second input file\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  Aligns corresponding lines from fileA and B, reporting cigar string.\n");
    fprintf(stderr, "  Lines are currently limited to 1 Mbp.\n");
    exit(1);
  }

  FILE *fileA = AS_UTL_openInputFile(nameA);
  FILE *fileB = AS_UTL_openInputFile(nameB);

  int32 lineMax = 1024 * 1024;
  int32 lenA    = 0;
  int32 lenB    = 0;

  char *lineA = new char [lineMax + 1];
  char *lineB = new char [lineMax + 1];

  splitToWords  sA;
  splitToWords  sB;

  readLine(fileA, lineA, lineMax, lenA, sA);
  readLine(fileB, lineB, lineMax, lenB, sB);

  while (1) {

    EdlibAlignResult result = edlibAlign(sA[0], lenA,
                                         sB[0], lenB,
                                         edlibNewAlignConfig(lenA + lenB, EDLIB_MODE_NW, EDLIB_TASK_PATH));

    assert(result.numLocations > 0);

    if (strcmp(sA[0], sB[0]) != 0)
      fprintf(stdout, "lost sync A %s B %s\n", sA[0], sB[0]), exit(1);

    char *cigar = edlibAlignmentToCigar(result.alignment,
                                        result.alignmentLength, (1) ? EDLIB_CIGAR_STANDARD : EDLIB_CIGAR_EXTENDED);

    edlibFreeAlignResult(result);

    if (strlen(cigar) > 50) {
      cigar[47] = '.';
      cigar[48] = '.';
      cigar[49] = '.';
      cigar[50] = 0;
    }

    fprintf(stdout, "%s ident %6.2f len(a-b) %6d cigar %s\n",
            sA[0],
            100.0 - 100.0 * result.editDistance / result.alignmentLength,
            lenA - lenB,
            cigar);

    delete [] cigar;

    //  The B file is allowed to have duplicate sequences.

    if (readLine(fileB, lineB, lineMax, lenB, sB) == false)
      break;

    if (strcmp(sA[0], sB[0]) != 0)
      if (readLine(fileA, lineA, lineMax, lenA, sA) == false)
        break;

    while (strcmp(sA[0], sB[0]) < 0) {
      fprintf(stdout, "B (at %s) lost sequence %s\n", sB[0], sA[0]);
      if (readLine(fileA, lineA, lineMax, lenA, sA) == false)
        break;
    }

    while (strcmp(sA[0], sB[0]) > 0) {
      fprintf(stdout, "A (at %s) lost sequence %s\n", sA[0], sB[0]);
      if (readLine(fileB, lineB, lineMax, lenB, sB) == false)
        break;
    }
  }

  delete [] lineA;
  delete [] lineB;

  AS_UTL_closeFile(fileA, nameA);
  AS_UTL_closeFile(fileB, nameB);

  //fprintf(stderr, "\n");
  //fprintf(stderr, "Bye.\n");

  return(0);
}
