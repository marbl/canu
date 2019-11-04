
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
 *    Brian P. Walenz beginning on 2018-NOV-07
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

//  g++6 -o filesTest -I.. -I. filesTest.C files.C

#include "files.H"

typedef  uint8   TYPE;


int32
main(int32 argc, char **argv) {
  uint64   nObj  = (uint64)16 * 1024 * 1024;
  TYPE    *array = new TYPE [nObj];
  TYPE     value = 0;


  if (1) {
    fprintf(stderr, "Initializing.\n");

    for (uint64 ii=0; ii<nObj; ii++)
      array[ii] = ii;

    fprintf(stderr, "Writing.\n");

    FILE *OUT = AS_UTL_openOutputFile("./filesTest.dat");

    writeToFile(array, "array", nObj, OUT);

    AS_UTL_closeFile(OUT);
  }


  if (1) {
    fprintf(stderr, "Reading - as one block.\n");

    FILE *IN = AS_UTL_openInputFile("./filesTest.dat");
    loadFromFile(array, "array", nObj, IN);
    AS_UTL_closeFile(IN);

    for (uint64 ii=0; ii<nObj; ii++)
      assert(array[ii] == (TYPE)ii);
  }


  if (1) {
    fprintf(stderr, "Reading.\n");

    FILE *IN = AS_UTL_openInputFile("./filesTest.dat");

    for (uint64 ii=0; ii<nObj; ii++) {
      loadFromFile(value, "value", IN);

      assert(value == (TYPE)ii);
    }

    fprintf(stderr, "Reading - one after eof.\n");
    loadFromFile(value, "value", IN, false);
    loadFromFile(value, "value", IN, true);

    AS_UTL_closeFile(IN);
  }


  exit(0);
}

