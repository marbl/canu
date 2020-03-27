
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

  assert((_rawU != NULL) || (_corU != NULL));

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

    else {
      fprintf(stderr, "sqRead::sqRead_loadDataFromBlob()--  unknown chunk type 0x%02x%02x%02x%02x '%c%c%c%c' skipped (lengths %u %u)\n",
              chunkName[0], chunkName[1], chunkName[2], chunkName[3],
              chunkName[0], chunkName[1], chunkName[2], chunkName[3],
              rawLength,
              corLength);
      assert(0);
    }

    //  All done.  Move to the next block.

    blobPos += 4 + 4 + chunkLen;
  }
}
