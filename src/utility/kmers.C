
/******************************************************************************
 *
 *  This file is part of 'sequence' and/or 'meryl', software programs for
 *  working with DNA sequence files and k-mers contained in them.
 *
 *  Modifications by:
 *
 *    Brian P. Walenz beginning on 2018-FEB-26
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.license' in the root directory of this distribution contains
 *  full conditions and disclaimers.
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
  char  bits[67];

  bits[0] = '0';
  bits[1] = 'x';

  uint32 mask = 1;
  uint32 bp   = 2;

  for (; mask < numFiles; mask <<= 1)            //  Technically, writes the reverse of the
    bits[bp++] = (outIndex & mask) ? '1' : '0';  //  prefix, but who cares?

  bits[bp] = 0;

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



void
removeBlock(char   *prefix,
            uint64  fileIndex,
            uint32  numFiles,
            uint32  iteration) {
  char    *name = constructBlockName(prefix, fileIndex, numFiles, iteration, false);

  AS_UTL_unlink(name);

  delete [] name;
}

