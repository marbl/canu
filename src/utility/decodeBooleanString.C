
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
 *    Brian P. Walenz on 2014-NOV-26
 *      are Copyright 2014 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

static
bool
decodeBoolean(char *feature, char *value) {
  bool ret = false;

  //  Decodes a string with 0/1, false/true, no/yes into an integer flag.

  switch (value[0]) {
    case '0':
    case 'f':
    case 'F':
    case 'n':
    case 'N':
      ret = false;
      break;
    case '1':
    case 't':
    case 'T':
    case 'y':
    case 'Y':
      ret = true;
      break;
    default:
      fprintf(stderr, "decodeBoolean()-- feature '%s' has unknown boolean value '%s'\n",
              feature, value);
      break;
  }

  return(ret);
}
