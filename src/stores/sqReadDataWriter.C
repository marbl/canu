
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



void
sqReadDataWriter::sqReadDataWriter_importData(sqRead *read) {

  _meta = read->_meta;   //  These point back to the sqStore.
  _rawU = read->_rawU;
  _rawC = read->_rawC;
  _corU = read->_corU;
  _corC = read->_corC;

  //  Copy the name, raw bases and corrected bases.  We cannot use
  //  setRawBases/setCorrectedBases because those want to change the metadata
  //  in sqStore -- even though it should be "changing" the value to the same
  //  thing.
  //
  //  Note that the lengths of these arrays include the NUL terminating byte,
  //  where the *Len variables below do not.

  uint32   namLen = strlen(read->_name);
  uint32   rawLen = read->sqRead_length(sqRead_raw);
  uint32   corLen = read->sqRead_length(sqRead_corrected);

  assert(0 == read->_name[namLen]);
  assert(0 == read->_rawBases[rawLen]);
  assert(0 == read->_corBases[corLen]);

  duplicateArray(_name,     _nameLen,     _nameAlloc,     read->_name,     namLen + 1);
  duplicateArray(_rawBases, _rawBasesLen, _rawBasesAlloc, read->_rawBases, rawLen + 1);
  duplicateArray(_corBases, _corBasesLen, _corBasesAlloc, read->_corBases, corLen + 1);

  assert(0 == _name[namLen]);
  assert(0 == _rawBases[rawLen]);
  assert(0 == _corBases[corLen]);
}



void
sqReadDataWriter::sqReadDataWriter_setName(const char *N) {
  duplicateArray(_name, _nameLen, _nameAlloc, N, (uint32)strlen(N) + 1);
}



void
sqReadDataWriter::sqReadDataWriter_setRawBases(const char *S, uint32 Slen) {

  setArraySize(_rawBases, _rawBasesLen, _rawBasesAlloc, Slen+1, _raAct::doNothing);

  for (uint32 ii=0; ii<Slen; ii++)
    _rawBases[ii] = _charMap[S[ii]];
  _rawBases[Slen] = 0;

  _rawBasesLen = Slen + 1;   //  Length INCLUDING NUL, remember?

  assert(_rawU->sqReadSeq_valid() == false);
  assert(_rawC->sqReadSeq_valid() == false);

  _rawU->sqReadSeq_setLength(_rawBases, _rawBasesLen-1, false);
  _rawC->sqReadSeq_setLength(_rawBases, _rawBasesLen-1, true);
}



void
sqReadDataWriter::sqReadDataWriter_setCorrectedBases(const char *S, uint32 Slen) {

  setArraySize(_corBases, _corBasesLen, _corBasesAlloc, Slen+1, _raAct::doNothing);

  for (uint32 ii=0; ii<Slen; ii++)
    _corBases[ii] = _charMap[S[ii]];
  _corBases[Slen] = 0;

  _corBasesLen = Slen + 1;   //  Length INCLUDING NUL, remember?

  assert(_corU->sqReadSeq_valid() == false);
  assert(_corC->sqReadSeq_valid() == false);

  _corU->sqReadSeq_setLength(_corBases, _corBasesLen-1, false);
  _corC->sqReadSeq_setLength(_corBases, _corBasesLen-1, true);
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

    assert(_rawU->sqReadSeq_valid() == true);   //  setRawBases() should be setting this true.
    assert(_rawC->sqReadSeq_valid() == true);

    rseq     = NULL;
    rseq2Len =                                      encode2bitSequence(rseq, _rawBases, _rawU->sqReadSeq_length());
    rseq3Len = (rseq2Len == 0)                    ? encode3bitSequence(rseq, _rawBases, _rawU->sqReadSeq_length()) : 0;
    rseqULen = (rseq2Len == 0) && (rseq3Len == 0) ? encode8bitSequence(rseq, _rawBases, _rawU->sqReadSeq_length()) : 0;
  }

  if ((_corBases != NULL) && (_corBases[0] != 0)) {
    assert(_corBasesLen > 0);

    assert(_corU->sqReadSeq_valid() == true);   //  setCorrectedBases should be setting this true.
    assert(_corC->sqReadSeq_valid() == true);

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


