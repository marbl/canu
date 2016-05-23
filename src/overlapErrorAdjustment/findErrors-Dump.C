
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
 *    Brian P. Walenz beginning on 2016-MAY-20
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "findErrors.H"





int
main(int argc, char **argv) {
  char       *redName = NULL;

  argc = AS_configure(argc, argv);

  int arg = 1;
  int err = 0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-r") == 0) {
      redName = argv[++arg];

    } else {
      fprintf(stderr, "Unknown option '%s'\n", argv[arg]);
      err++;
    }

    arg++;
  }

  if (redName == NULL)
    err++;

  if (err > 0) {
    fprintf(stderr, "usage: %s -r file.red\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "Dumps, as ASCII, the results from findErrors *.red files.\n");

    if (redName == NULL)
      fprintf(stderr, "ERROR: no *.red file (-r) supplied.\n");

    exit(1);
  }

  char  *typeName[13] = { "IDENT",
                          "DELETE",
                          "A_SUBST",
                          "C_SUBST",
                          "G_SUBST",
                          "T_SUBST",
                          "A_INSERT",
                          "C_INSERT",
                          "G_INSERT",
                          "T_INSERT",
                          "NO_VOTE",
                          "EXTENSION",
                          NULL };

  memoryMappedFile     *Cfile    = new memoryMappedFile(redName);
  Correction_Output_t  *C        = (Correction_Output_t *)Cfile->get();
  uint64                Cpos     = 0;
  uint64                Clen     = Cfile->length() / sizeof(Correction_Output_t);

  for (uint32 ii=0; ii<Clen; ii++) {
    fprintf(stdout, "%8u %12s %8u %c %c\n",
            C[ii].readID,
            typeName[C[ii].type],
            C[ii].pos,
            C[ii].keep_left  ? 't' : 'f',
            C[ii].keep_right ? 't' : 'f');
  }

  exit(0);
}


