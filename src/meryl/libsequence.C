
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
 *    Brian P. Walenz beginning on 2018-JAN-29
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "libsequence.H"


dnaSeqFile::dnaSeqFile(const char *filename, bool indexed) {

  _file     = new compressedFileReader(filename);
  _buffer   = new readBuffer(_file->file());

  _index    = NULL;
  _indexLen = 0;
  _indexMax = 0;

  if (indexed == false)
    return;

  if (_file->isCompressed() == true)
    fprintf(stderr, "ERROR: cannot index compressed input '%s'.\n", filename), exit(1);

  if (_file->isNormal() == false)
    fprintf(stderr, "ERROR: cannot index pipe input.\n"), exit(1);

  fprintf(stderr, "ERROR: indexing not supported (yet).\n"), exit(1);
}



dnaSeqFile::~dnaSeqFile() {
  delete    _file;
  delete    _buffer;
  delete [] _index;
}



uint64
dnaSeqFile::loadFASTA(char   *&name,     uint32   nameMax,
                      char   *&seq,
                      uint8  *&qlt,      uint64   seqMax) {
  uint64  nameLen = 0;
  uint64  seqLen  = 0;
  char    ch      = _buffer->read();

  assert(ch == '>');

  //  Read the header line into the name string.

  for (ch=_buffer->read(); (ch != '\n') && (ch != 0); ch=_buffer->read()) {
    if (nameLen+1 >= nameMax)
      resizeArray(name, nameLen, nameMax, 3 * nameMax / 2);
    name[nameLen++] = ch;
  }

  //  Read sequence, skipping whitespace, until we hit a new sequence (or eof).

  for (ch=_buffer->readuntil('>'); (ch != '>') && (ch != 0); ch=_buffer->readuntil('>')) {
    if (ch == '\n')
      continue;

    assert(_buffer->eof() == false);

    if (seqLen+1 >= seqMax)
      resizeArrayPair(seq, qlt, seqLen, seqMax, 3 * seqMax / 2);

    seq[seqLen] = ch;
    qlt[seqLen] = 0;

    seqLen++;
  }

  name[nameLen] = 0;
  seq[seqLen] = 0;
  qlt[seqLen] = 0;

  assert(nameLen < nameMax);
  assert(seqLen  < seqMax);

  return(seqLen);
}



uint64
dnaSeqFile::loadFASTQ(char   *&name,     uint32   nameMax,
                      char   *&seq,
                      uint8  *&qlt,      uint64   seqMax) {
  uint32  nameLen = 0;
  uint64  seqLen  = 0;
  uint64  qltLen  = 0;
  char    ch      = _buffer->read();  

  assert(ch == '@');

  //  Read the header line into the name string.

  for (ch=_buffer->read(); (ch != '\n') && (ch != 0); ch=_buffer->read()) {
    if (nameLen+1 >= nameMax)
      resizeArray(name, nameLen, nameMax, 3 * nameMax / 2);
    name[nameLen++] = ch;
  }

  //  Read sequence.

  for (ch=_buffer->read(); (ch != '\n') && (ch != 0); ch=_buffer->read()) {
    if (seqLen+1 >= seqMax)
      resizeArrayPair(seq, qlt, seqLen, seqMax, 3 * seqMax / 2);
    seq[seqLen++] = ch;
  }

  //  Skip header line

  for (ch=_buffer->read(); (ch != '\n') && (ch != 0); ch=_buffer->read()) {
    ;
  }

  //  Read qualities.

  for (ch=_buffer->read(); (ch != '\n') && (ch != 0); ch=_buffer->read()) {
    if (qltLen+1 >= seqMax)
      resizeArrayPair(seq, qlt, qltLen, seqMax, 3 * seqMax / 2);
    qlt[qltLen++] = ch;
  }

  //fprintf(stderr, "READ FASTQ name %u seq %lu qlt %lu\n", nameLen, seqLen, qltLen);

  name[nameLen] = 0;
  seq[seqLen] = 0;
  qlt[qltLen] = 0;

  assert(nameLen < nameMax);
  assert(seqLen  < seqMax);
  assert(qltLen  < seqMax);
  assert(seqLen == qltLen);

  return(seqLen);
}



bool
dnaSeqFile::loadSequence(char   *&name,     uint32   nameMax,
                         char   *&seq,
                         uint8  *&qlt,      uint64   seqMax,
                         uint64  &seqLen) {

  if (nameMax == 0)
    resizeArray(name, 0, nameMax, (uint32)1024);

  if (seqMax == 0)
    resizeArrayPair(seq, qlt, 0, seqMax, (uint64)65536);

  while (_buffer->peek() == '\n')
    _buffer->read();

  if      (_buffer->peek() == '>')
    seqLen = loadFASTA(name, nameMax,
                       seq,
                       qlt, seqMax);

  else if (_buffer->peek() == '@')
    seqLen = loadFASTQ(name, nameMax,
                       seq,
                       qlt, seqMax);

  else
    return(false);

  return(true);
}



bool
dnaSeqFile::loadBases(char    *seq,
                      uint64   maxLength,
                      uint64  &seqLength) {

  seqLength = 0;

  if (_buffer->eof() == true)
    return(false);

  //  If this is a new file, skip the first name line.

  if (_buffer->tell() == 0) {
    while (_buffer->peek() == '\n')    //  Skip whitespace before the first name line.
      _buffer->read();

    for (char ch = _buffer->read(); (ch != '\n') && (ch != 0); ch = _buffer->read())
      ;
  }

  //  Skip whitespace.

  while (_buffer->peek() == '\n')
    _buffer->read();

  //  We're now at sequence, so load until we're not in sequence or out of space.
  //  Read the first letter.

  char    ch = _buffer->read();

  while (_buffer->eof() == false) {

    //  End of a FASTA sequence.  Return the bases loaded so far, even if there are none.
    if (ch == '>') {
      for (ch = _buffer->read(); (ch != '\n') && (ch != 0); ch = _buffer->read())   //  Skip the name of the next sequence
        ;
      //fprintf(stderr, "loadBases()-- hit new FASTA, return seq length %lu pos %lu\n", seqLength, _buffer->tell());
      return(true);
    }

    //  End of a FASTQ sequence.  Skip the qualities.  Then return the bases loaded so far, even if there are none.
    if (ch == '+') {
      for (ch = _buffer->read(); (ch != '\n') && (ch != 0); ch = _buffer->read())   //  Skip the quality name line
        ;
      for (ch = _buffer->read(); (ch != '\n') && (ch != 0); ch = _buffer->read())   //  Skip the qualities
        ;
      for (ch = _buffer->read(); (ch != '\n') && (ch != 0); ch = _buffer->read())   //  Skip the name of the next sequence
        ;
      //fprintf(stderr, "loadBases()-- hit FASTQ qlt, return seq length %lu pos %lu\n", seqLength, _buffer->tell());
      return(true);
    }

    //  Otherwise, add the base and move ahead.

    if (ch != '\n') {
      //fprintf(stderr, "loadBases()-- LOAD %d %c into position %d pos %lu\n", ch, ch, seqLength, _buffer->tell());
      seq[seqLength++] = ch;
    }

    if (seqLength == maxLength) {
      //fprintf(stderr, "loadBases()-- buffer full, return seq length %lu pos %lu.  next: %d %c\n", seqLength, _buffer->tell(), _buffer->peek(), _buffer->peek());
      return(true);
    }

    //fprintf(stderr, "loadBases()-- READ %d %c at buffer pos %lu\n", _buffer->peek(), _buffer->peek(), _buffer->tell());
    ch = _buffer->read();
  }

  //  We've hit EOF.  If bases were loaded, return them, otherwise, just say we're done.

  //fprintf(stderr, "loadBases()-- file empty, return seq length %lu pos %lu.\n", seqLength, _buffer->tell());
  return(seqLength > 0);
}
