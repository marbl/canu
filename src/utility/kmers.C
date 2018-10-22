
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
 *    Brian P. Walenz beginning on 2018-JUL-21
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "kmers.H"
#include "bits.H"

#include "files.H"


uint32 kmerTiny::_merSize   = 0;
uint64 kmerTiny::_fullMask  = 0;
uint64 kmerTiny::_leftMask  = 0;
uint32 kmerTiny::_leftShift = 0;


char *
constructBlockName(char   *prefix,
                   uint64  outIndex,
                   uint32  numFiles,
                   uint32  iteration,
                   bool    isIndex) {
  char *name = new char [FILENAME_MAX+1];
  char  bits[67] = { 0 };

  bits[0] = '0';
  bits[1] = 'x';

  uint32 bp = 2;

  for (uint32 mask=1; mask < numFiles; mask <<= 1)   //  Count the number of digits we need.
    bp++;

  for (uint32 mask=1; mask < numFiles; mask <<= 1)   //  Then make the name from right to left.
    bits[--bp] = (outIndex & mask) ? '1' : '0';

  if (iteration == 0)
    snprintf(name, FILENAME_MAX, "%s/%s.%s", prefix, bits, (isIndex == false) ? "merylData" : "merylIndex");
  else
    snprintf(name, FILENAME_MAX, "%s/%s[%03u].%s", prefix, bits, iteration, (isIndex == false) ? "merylData" : "merylIndex");

  return(name);
}



FILE *
openOutputBlock(char   *prefix,
                uint64  fileIndex,
                uint32  numFiles,
                uint32  iteration) {
  char    *name = constructBlockName(prefix, fileIndex, numFiles, iteration, false);

  FILE *F = AS_UTL_openOutputFile(name);

  delete [] name;

  return(F);
}



FILE *
openInputBlock(char   *prefix,
               uint64  fileIndex,
               uint32  numFiles,
               uint32  iteration) {
  char    *name = constructBlockName(prefix, fileIndex, numFiles, iteration, false);

  FILE *F = AS_UTL_openInputFile(name);

  delete [] name;

  return(F);
}
