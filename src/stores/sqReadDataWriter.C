
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



void
sqReadDataWriter::sqReadDataWriter_importData(sqRead *read) {

  _meta = read->_meta;
  _rawU = read->_rawU;
  _rawC = read->_rawC;
  _corU = read->_corU;
  _corC = read->_corC;

  sqReadDataWriter_setName          (read->_name);
  sqReadDataWriter_setRawBases      (read->_rawBases, read->sqRead_length(sqRead_raw));
  sqReadDataWriter_setCorrectedBases(read->_corBases, read->sqRead_length(sqRead_corrected));
}



void
sqReadDataWriter::sqReadDataWriter_writeBlob(writeBuffer *buffer) {

  //  Encode the data, so we know (approximately) how much data we're storing.
  //  If that is over the limit for each file, move to the next file.
  //
  //  Note that during encoding, the read metadata is updated,

  uint8   *null     = NULL;

  uint8   *rseq     = NULL;
  uint32   rseq2Len = 0;
  uint32   rseq3Len = 0;
  uint32   rseqULen = 0;

  uint8   *cseq     = NULL;
  uint32   cseq2Len = 0;
  uint32   cseq3Len = 0;
  uint32   cseqULen = 0;

  //  The sqReadSeq pointers are NULL when we're writing to a non-store file.
  //  But if we're writing to the store, they all need to be present.

  if ((_rawU != NULL) ||
      (_rawC != NULL) ||
      (_corU != NULL) ||
      (_corC != NULL))
    assert((_rawU != NULL) && (_rawC != NULL) && (_corU != NULL) && (_corC != NULL));

  if ((_rawU == NULL) ||
      (_rawC == NULL) ||
      (_corU == NULL) ||
      (_corC == NULL))
    assert((_rawU == NULL) && (_rawC == NULL) && (_corU == NULL) && (_corC == NULL));

  if ((_rawBases != NULL) && (_rawBases[0] != 0)) {
    assert(_rawBasesLen > 0);

    if ((_rawU) && (_rawU->sqReadSeq_valid() == false))   _rawU->sqReadSeq_setLength(_rawBases, _rawBasesLen-1, false);
    if ((_rawC) && (_rawC->sqReadSeq_valid() == false))   _rawC->sqReadSeq_setLength(_rawBases, _rawBasesLen-1, true);

    rseq     = NULL;
    rseq2Len =                                      encode2bitSequence(rseq, _rawBases, _rawU->sqReadSeq_length());
    rseq3Len = (rseq2Len == 0)                    ? encode3bitSequence(rseq, _rawBases, _rawU->sqReadSeq_length()) : 0;
    rseqULen = (rseq2Len == 0) && (rseq3Len == 0) ? encode8bitSequence(rseq, _rawBases, _rawU->sqReadSeq_length()) : 0;
  }

  if ((_corBases != NULL) && (_corBases[0] != 0)) {
    assert(_corBasesLen > 0);

    if ((_corU) && (_corU->sqReadSeq_valid() == false))   _corU->sqReadSeq_setLength(_corBases, _corBasesLen-1, false);
    if ((_corC) && (_corC->sqReadSeq_valid() == false))   _corC->sqReadSeq_setLength(_corBases, _corBasesLen-1, true);

    cseq     = NULL;
    cseq2Len =                                      encode2bitSequence(cseq, _corBases, _corU->sqReadSeq_length());
    cseq3Len = (cseq2Len == 0)                    ? encode3bitSequence(cseq, _corBases, _corU->sqReadSeq_length()) : 0;
    cseqULen = (cseq2Len == 0) && (cseq3Len == 0) ? encode8bitSequence(cseq, _corBases, _corU->sqReadSeq_length()) : 0;
  }

  //  Write the header and name.

  buffer->writeIFFchunk("BLOB");
  buffer->writeIFFchunk("NAME", _name, _nameLen);

  //  Write raw bases.

  if (rseq2Len > 0)
    buffer->writeIFFchunk("2SQR", rseq, rseq2Len);    //  Two-bit encoded sequence (ACGT only)
  if (rseq3Len > 0)
    buffer->writeIFFchunk("3SQR", rseq, rseq3Len);    //  Three-bit encoded sequence (ACGTN)
  if (rseqULen > 0)
    buffer->writeIFFchunk("USQR", rseq, rseqULen);    //  Unencoded sequence

  //  Write corrected bases.

  if (cseq2Len > 0)
    buffer->writeIFFchunk("2SQC", cseq, cseq2Len);    //  Two-bit encoded sequence (ACGT only)
  if (cseq3Len > 0)
    buffer->writeIFFchunk("3SQC", cseq, cseq3Len);    //  Three-bit encoded sequence (ACGTN)
  if (cseqULen > 0)
    buffer->writeIFFchunk("USQC", cseq, cseqULen);    //  Unencoded sequence

  //  And terminate the blob.

  buffer->closeIFFchunk("BLOB");

  delete [] rseq;
  delete [] cseq;
}


