
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
 *    src/stores/gkStore.C
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2014-NOV-26 to 2015-AUG-10
 *      are Copyright 2014-2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2015-OCT-09
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *    Sergey Koren beginning on 2015-DEC-09
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "sqStore.H"
#include "sequence.H"
#include "files.H"



//  Lowest level function to load data into a read.
//
void
sqRead::sqRead_decodeBlob(void) {

  //  Resize strings.

  uint32  rawLength = 0;
  uint32  corLength = 0;

  if (_rawU)   rawLength = _rawU->sqReadSeq_length();
  if (_corU)   corLength = _corU->sqReadSeq_length();

  resizeArray(_rawBases, 0, _rawBasesAlloc, rawLength+1, resizeArray_doNothing);
  resizeArray(_corBases, 0, _corBasesAlloc, corLength+1, resizeArray_doNothing);

  //  Forget what sequence we previously returned to the user.

  _retFlags = 0;

  //  Decode the blob data until there is no more data.

  for (uint32 blobPos=0; blobPos < _blobLen; ) {
    char   *chunkName =  (char *)  (_blob + blobPos + 0);
    uint32  chunkLen  = *(uint32 *)(_blob + blobPos + 4);
    uint8  *chunk     =            (_blob + blobPos + 8);

    //  Decode the NAME?

    if      (strncmp(chunkName, "NAME", 4) == 0) {
      resizeArray(_name, 0, _nameAlloc, chunkLen + 1, resizeArray_doNothing);
      memcpy(_name, chunk, chunkLen);
      _name[chunkLen] = 0;
    }

    //  Decode the raw bases?

    else if ((rawLength > 0) && (strncmp(chunkName, "2SQR", 4) == 0))
      decode2bitSequence(chunk, chunkLen, _rawBases, rawLength);

    else if ((rawLength > 0) && (strncmp(chunkName, "3SQR", 4) == 0))
      decode3bitSequence(chunk, chunkLen, _rawBases, rawLength);

    else if ((rawLength > 0) && (strncmp(chunkName, "USQR", 4) == 0))
      decode8bitSequence(chunk, chunkLen, _rawBases, rawLength);


    //  Decode the corrected bases

    else if ((corLength > 0) && (strncmp(chunkName, "2SQC", 4) == 0))
      decode2bitSequence(chunk, chunkLen, _corBases, corLength);

    else if ((corLength > 0) && (strncmp(chunkName, "3SQC", 4) == 0))
      decode3bitSequence(chunk, chunkLen, _corBases, corLength);

    else if ((corLength > 0) && (strncmp(chunkName, "USQC", 4) == 0))
      decode8bitSequence(chunk, chunkLen, _corBases, corLength);


    //  No idea what this is then.

    else
      fprintf(stderr, "sqRead::sqRead_loadDataFromBlob()--  unknown chunk type %02x %02x %02x %02x '%c%c%c%c' skipped\n",
              chunkName[0], chunkName[1], chunkName[2], chunkName[3],
              chunkName[0], chunkName[1], chunkName[2], chunkName[3]);

    //  All done.  Move to the next block.

    blobPos += 4 + 4 + chunkLen;
  }
}
