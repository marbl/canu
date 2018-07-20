
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
 *    src/AS_BAT/memoryMappedFileTest.C
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2014-NOV-26 to 2014-NOV-27
 *      are Copyright 2014 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "memoryMappedFile.H"

int
main(int argc, char **argv) {

  memoryMappedFileRW<uint32>  *write = new memoryMappedFileRW<uint32>("testWrite");

  for (uint32 ii=0; ii<1000000000; ii++)
    (*write)[ii] = ii;

  delete write;

  exit(0);
}
