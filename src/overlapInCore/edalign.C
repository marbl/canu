
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

#include "sequence.H"
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



void
pairAlign(char *nameA, char *nameB) {
  dnaSeqFile  *fileA = new dnaSeqFile(nameA);
  dnaSeqFile  *fileB = new dnaSeqFile(nameB);
  dnaSeq       seqA;
  dnaSeq       seqB;

  while (fileA->loadSequence(seqA) &&
         fileB->loadSequence(seqB)) {

    EdlibAlignResult result = edlibAlign(seqA.bases(), seqA.length(),
                                         seqB.bases(), seqB.length(),
                                         edlibNewAlignConfig(0.001 * (seqA.length() + seqB.length()), EDLIB_MODE_HW, EDLIB_TASK_PATH));

    assert(result.numLocations > 0);

    //if (strcmp(sA[0], sB[0]) != 0)
    //  fprintf(stdout, "lost sync A %s B %s\n", sA[0], sB[0]), exit(1);

    char *cigar = edlibAlignmentToCigar(result.alignment,
                                        result.alignmentLength, (1) ? EDLIB_CIGAR_STANDARD : EDLIB_CIGAR_EXTENDED);

    edlibFreeAlignResult(result);

    //if (strlen(cigar) > 50) {
    //  cigar[47] = '.';
    //  cigar[48] = '.';
    //  cigar[49] = '.';
    //  cigar[50] = 0;
    //}

    fprintf(stdout, "%s ident %6.2f len(a-b) %6d cigar %s\n",
            seqA.name(),
            100.0 - 100.0 * result.editDistance / result.alignmentLength,
            (int32)seqA.length() - (int32)seqB.length(),
            cigar);

    delete [] cigar;

    //  The B file is allowed to have duplicate sequences.

#if 0
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
#endif
  }

  delete fileA;
  delete fileB;
}



void
refAlign(char *nameA, char *nameB) {
  dnaSeq        seqA;
  dnaSeqFile   *fileA = new dnaSeqFile(nameA);
  dnaSeq        seqB;
  dnaSeqFile   *fileB = new dnaSeqFile(nameB);

  fileB->loadSequence(seqB);

  while (fileA->loadSequence(seqA) == true) {
    EdlibAlignResult result = edlibAlign(seqA.bases(), seqA.length(),   //  Free end gaps!
                                         seqB.bases(), seqB.length(),
                                         edlibNewAlignConfig(0.25 * seqA.length(), EDLIB_MODE_HW, EDLIB_TASK_PATH));

    char *cigar = edlibAlignmentToCigar(result.alignment,
                                        result.alignmentLength, (1) ? EDLIB_CIGAR_STANDARD : EDLIB_CIGAR_EXTENDED);

    if (strlen(cigar) > 50) {
      cigar[47] = '.';
      cigar[48] = '.';
      cigar[49] = '.';
      cigar[50] = 0;
    }

    uint32   nMatch, nMismatch, nInsertOpen, nInsert, nDeleteOpen, nDelete;

    edlibAlignmentAnalyze(result.alignment, result.alignmentLength, nMatch, nMismatch, nInsertOpen, nInsert, nDeleteOpen, nDelete);

    if (result.numLocations > 0) {
      fprintf(stdout, "%s %8d-%-8d alignLen %6d %6.2f%% gap %6.2f %6.2f match %7u mismatch %7u ins %7u %7u del %7u %7u cigar %s\n",
              nameA,
              result.startLocations[0],
              result.endLocations[0] + 1,
              result.alignmentLength,
              100.0 - 100.0 * result.editDistance / result.alignmentLength,
              100.0 * (nInsertOpen + nDeleteOpen) / result.alignmentLength,
              100.0 * (nInsert     + nDelete)     / result.alignmentLength,
              nMatch, nMismatch, nInsertOpen, nInsert, nDeleteOpen, nDelete,
              cigar);
    }

    else {
      fprintf(stdout, "%s %8d-%-8d %6.2f%% %7u %7u ins %7u %7u del %7u %7u cigar %s\n",
              nameA, 0, 0, 0.0, 0, 0, 0, 0, 0, 0, "0M");
    }

    delete [] cigar;

    edlibFreeAlignResult(result);
  }

  delete fileA;
  delete fileB;
}



int
main(int argc, char **argv) {
  char    *nameA           = NULL;
  char    *nameB           = NULL;
  bool     reference       = false;

  argc = AS_configure(argc, argv);

  int err=0;
  int arg=1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-a") == 0) {
      nameA = argv[++arg];

    } else if (strcmp(argv[arg], "-b") == 0) {
      nameB = argv[++arg];

    } else if (strcmp(argv[arg], "-ref") == 0) {
      nameB = argv[++arg];
      reference = true;
      arg++;
      break;

    } else {
      err++;
    }

    arg++;
  }

  if ((reference == false) && ((nameA == NULL) || (nameB == NULL)))  err++;
  if ((reference == true)  && (arg == argc))                         err++;

  if (err) {
    fprintf(stderr, "usage: %s [-a file -b file] [-ref file file ...]\n", argv[0]);
    fprintf(stderr, "  PAIRWISE MODE:  Align two sequences globally.\n");
    //fprintf(stderr, "  Aligns corresponding lines from fileA and B, reporting cigar string.\n");
    //fprintf(stderr, "  Lines are currently limited to 1 Mbp.\n");
    fprintf(stderr, "    -a fileA     Mandatory, path to first input file\n");
    fprintf(stderr, "    -b fileB     Mandatory, path to second input file\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  REFERENCE MODE:  Align multiple sequences to reference.\n");
    fprintf(stderr, "  -ref R.fasta\n");
    fprintf(stderr, "\n");
    exit(1);
  }

  if (reference == false)
    pairAlign(nameA, nameB);
  else
    while (arg < argc)
      refAlign(argv[arg++], nameB);

  //fprintf(stderr, "\n");
  //fprintf(stderr, "Bye.\n");

  return(0);
}
