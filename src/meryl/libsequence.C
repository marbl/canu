
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

  while ((ch = _buffer->readuntil('\n')) != 0) {
    if (nameLen+1 >= nameMax)
      resizeArray(name, nameLen, nameMax, 3 * nameMax / 2);
    name[nameLen++] = ch;
  }

  //  Read sequence, skipping whitespace, into seq until we hit a new sequence.

  while ((ch = _buffer->readuntil('>')) != 0) {
    if (ch == '\n')
      continue;

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
  uint64  nameLen = 0;
  uint64  seqLen  = 0;
  uint64  qltLen  = 0;
  char    ch      = _buffer->read();

  assert(ch == '@');

  //  Read the header line into the name string.

  while ((ch = _buffer->readuntil('\n')) != 0) {
    if (nameLen+1 >= nameMax)
      resizeArray(name, nameLen, nameMax, 3 * nameMax / 2);
    name[nameLen++] = ch;
  }

  //  Read sequence.

  while ((ch = _buffer->readuntil('\n')) != 0) {
    if (seqLen+1 >= seqMax)
      resizeArrayPair(seq, qlt, seqLen, seqMax, 3 * seqMax / 2);
    seq[seqLen++] = ch;
  }

  //  Skip header line

  while ((ch = _buffer->readuntil('\n')) != 0)
    ;

  //  Read qualities.

  while ((ch = _buffer->readuntil('\n')) != 0) {
    if (seqLen+1 >= seqMax)
      resizeArrayPair(seq, qlt, seqLen, seqMax, 3 * seqMax / 2);
    qlt[qltLen++] = ch;
  }

  name[nameLen] = 0;
  seq[seqLen] = 0;
  qlt[qltLen] = 0;

  assert(nameLen < nameMax);
  assert(seqLen  < seqMax);
  assert(qltLen  < seqMax);
  assert(seqLen == qltLen);

  return(seqLen);
}



uint64
dnaSeqFile::loadSequence(char   *&name,     uint32   nameMax,
                         char   *&seq,
                         uint8  *&qlt,      uint64   seqMax) {
  char  ch = _buffer->peek();

  if (nameMax == 0)
    resizeArray(name, 0, nameMax, (uint32)1024);

  if (seqMax == 0)
    resizeArrayPair(seq, qlt, 0, seqMax, (uint64)65536);

  if        (ch == '>') {
    return(loadFASTA(name, nameMax,
                     seq,
                     qlt, seqMax));

  } else if (ch == '@') {
    return(loadFASTQ(name, nameMax,
                     seq,
                     qlt, seqMax));

  } else {
    return(0);
  }
}



uint64
dnaSeqFile::loadBases(char   *&seq,
                      uint32   length) {
  uint64  seqLen  = 0;
  char    ch      = _buffer->peek();

  //  If at the start of a sequence, skip the header

  if ((ch == '>') ||
      (ch == '@')) {
    while ((ch = _buffer->readuntil('\n')) != 0)
      ;
  }

  //  We're at sequence, so load until we're not in sequence or out of space.

  while ((_buffer->eof() == false) &&
         (seqLen < length)) {

    //  End of a FASTA sequence, just return the length.
    if (ch == '>') {
      return(seqLen);
    }

    //  End of a FASTQ sequence, skip ahead to the next sequence, then return the length.
    if (ch == '+') {
      while ((ch = _buffer->readuntil('\n')) != 0)
        ;
      while ((ch = _buffer->readuntil('\n')) != 0)
        ;
      return(seqLen);
    }

    //  Otherwise, add the base and move ahead.

    if (ch != '\n')
      seq[seqLen++] = ch;

    ch = _buffer->read();
  }

  //  Out of space in the user buffer, return the length.

  return(seqLen);
}









