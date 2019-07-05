
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
 *    Sergey Koren beginning on 2018-AUG-30
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "sequence.H"




static
const
char
inv[256] = {
   0,  0,  0,  0,  0,  0,  0,  0,  //  0x00 -
   0,  0,  0,  0,  0,  0,  0,  0,  //  0x08 -
   0,  0,  0,  0,  0,  0,  0,  0,  //  0x10 -
   0,  0,  0,  0,  0,  0,  0,  0,  //  0x18 -
   0,  0,  0,  0,  0,  0,  0,  0,  //  0x20 -  !"#$%&'
   0,  0,  0,  0,  0, '-', 0,  0,  //  0x28 - ()*+,-./
   0,  0,  0,  0,  0,  0,  0,  0,  //  0x30 - 01234567
   0,  0,  0,  0,  0,  0,  0,  0,  //  0x38 - 89:;<=>?
   0, 'T', 0, 'G', 0,  0,  0, 'C', //  0x40 - @ABCDEFG
   0,  0,  0,  0,  0,  0, 'N', 0,  //  0x48 - HIJKLMNO
   0,  0,  0,  0, 'A', 0,  0,  0,  //  0x50 - PQRSTUVW
   0,  0,  0,  0,  0,  0,  0,  0,  //  0x58 - XYZ[\]^_
   0, 't', 0, 'g', 0,  0,  0, 'c', //  0x60 - `abcdefg
   0,  0,  0,  0,  0,  0, 'n', 0,  //  0x68 - hijklmno
   0,  0,  0,  0, 'a', 0,  0,  0,  //  0x70 - pqrstuvw
   0,  0,  0,  0,  0,  0,  0,  0,  //  0x78 - xyz{|}~
   0,  0,  0,  0,  0,  0,  0,  0,  //  0x80 -
   0,  0,  0,  0,  0,  0,  0,  0,  //  0x88 -
   0,  0,  0,  0,  0,  0,  0,  0,  //  0x90 -
   0,  0,  0,  0,  0,  0,  0,  0,  //  0x98 -
   0,  0,  0,  0,  0,  0,  0,  0,  //  0xa0 -
   0,  0,  0,  0,  0,  0,  0,  0,  //  0xa8 -
   0,  0,  0,  0,  0,  0,  0,  0,  //  0xb0 -
   0,  0,  0,  0,  0,  0,  0,  0,  //  0xb8 -
   0,  0,  0,  0,  0,  0,  0,  0,  //  0xc0 -
   0,  0,  0,  0,  0,  0,  0,  0,  //  0xc8 -
   0,  0,  0,  0,  0,  0,  0,  0,  //  0xd0 -
   0,  0,  0,  0,  0,  0,  0,  0,  //  0xd8 -
   0,  0,  0,  0,  0,  0,  0,  0,  //  0xe0 -
   0,  0,  0,  0,  0,  0,  0,  0,  //  0xe8 -
   0,  0,  0,  0,  0,  0,  0,  0,  //  0xf0 -
   0,  0,  0,  0,  0,  0,  0,  0   //  0xf8 -
};

static
const
char
Dacgtn[5] = { 'A',
              'C',
              'G',
              'T',
              'N' };

static
const
uint8
Eacgtn[256] = {
   0,    0,    0,    0,    0,    0,    0,    0,    //  0x00 -
   0,    0,    0,    0,    0,    0,    0,    0,    //  0x08 -
   0,    0,    0,    0,    0,    0,    0,    0,    //  0x10 -
   0,    0,    0,    0,    0,    0,    0,    0,    //  0x18 -
   0,    0,    0,    0,    0,    0,    0,    0,    //  0x20 -  !"#$%&'
   0,    0,    0,    0,    0,    0,    0,    0,    //  0x28 - ()*+,-./
   0,    0,    0,    0,    0,    0,    0,    0,    //  0x30 - 01234567
   0,    0,    0,    0,    0,    0,    0,    0,    //  0x38 - 89:;<=>?
   0, 0x00,    0, 0x01,    0,    0,    0, 0x02,    //  0x40 - @ABCDEFG
   0,    0,    0,    0,    0,    0, 0x04,    0,    //  0x48 - HIJKLMNO
   0,    0,    0,    0, 0x03,    0,    0,    0,    //  0x50 - PQRSTUVW
   0,    0,    0,    0,    0,    0,    0,    0,    //  0x58 - XYZ[\]^_
   0, 0x00,    0, 0x01,    0,    0,    0, 0x02,    //  0x60 - `abcdefg
   0,    0,    0,    0,    0,    0, 0x04,    0,    //  0x68 - hijklmno
   0,    0,    0,    0, 0x03,    0,    0,    0,    //  0x70 - pqrstuvw
   0,    0,    0,    0,    0,    0,    0,    0,    //  0x78 - xyz{|}~
   0,    0,    0,    0,    0,    0,    0,    0,    //  0x80 -
   0,    0,    0,    0,    0,    0,    0,    0,    //  0x88 -
   0,    0,    0,    0,    0,    0,    0,    0,    //  0x90 -
   0,    0,    0,    0,    0,    0,    0,    0,    //  0x98 -
   0,    0,    0,    0,    0,    0,    0,    0,    //  0xa0 -
   0,    0,    0,    0,    0,    0,    0,    0,    //  0xa8 -
   0,    0,    0,    0,    0,    0,    0,    0,    //  0xb0 -
   0,    0,    0,    0,    0,    0,    0,    0,    //  0xb8 -
   0,    0,    0,    0,    0,    0,    0,    0,    //  0xc0 -
   0,    0,    0,    0,    0,    0,    0,    0,    //  0xc8 -
   0,    0,    0,    0,    0,    0,    0,    0,    //  0xd0 -
   0,    0,    0,    0,    0,    0,    0,    0,    //  0xd8 -
   0,    0,    0,    0,    0,    0,    0,    0,    //  0xe0 -
   0,    0,    0,    0,    0,    0,    0,    0,    //  0xe8 -
   0,    0,    0,    0,    0,    0,    0,    0,    //  0xf0 -
   0,    0,    0,    0,    0,    0,    0,    0     //  0xf8 -
};




void
reverseComplementSequence(char *seq, int len) {
  char   c=0;
  char  *s=seq,  *S=seq+len-1;

  if (len == 0) {
    len = strlen(seq);
    S = seq + len - 1;
  }

  while (s < S) {
    c    = *s;
    *s++ =  inv[*S];
    *S-- =  inv[c];
  }

  if (s == S)
    *s = inv[*s];
}



char *
reverseComplementCopy(char *seq, int len) {
  char  *rev = new char [len+1];

  assert(len > 0);

  for (int32 p=len, q=0; p>0; )
    rev[q++] = inv[seq[--p]];

  rev[len] = 0;

  return(rev);
}



template<typename qvType>
void
reverseComplement(char *seq, qvType *qlt, int len) {
  char    c=0;
  char   *s=seq,  *S=seq+len-1;
  qvType *q=qlt,  *Q=qlt+len-1;

  if (qlt == NULL) {
    reverseComplementSequence(seq, len);
    return;
  }

  if (len == 0) {
    len = strlen(seq);
    S = seq + len - 1;
    Q = qlt + len - 1;
  }

  while (s < S) {
    c    = *s;
    *s++ =  inv[*S];
    *S-- =  inv[c];

    c    = *q;
    *q++ = *Q;
    *Q-- =  c;
  }

  if (s == S)
    *s = inv[*s];
}

template void reverseComplement<char> (char *seq, char  *qlt, int len);   //  Give the linker
template void reverseComplement<uint8>(char *seq, uint8 *qlt, int len);   //  something to link



uint32
homopolyCompress(char *bases, char *compr) {
  uint32  cc = 0;  //  position of the start of the run
  uint32  rr = 1;  //  position of the scan head
  uint32  sl = 0;  //  length of the compressed sequence

  //  Assume bases is never empty.

  if (compr)
    compr[sl] = bases[cc];

  if (bases[0] == 0)
    return(0);

  sl += 1;

  while (bases[rr] != 0) {
    //  In a run, move the scan head one position.
    if (bases[cc] == bases[rr]) {
      rr += 1;
    }

    //  Just ended a run.  Set the start of the (next) run
    //  to the current position, move the current position
    //  to the next base, and increase the length of the
    //  compressed sequence.
    else {
      if (compr)
        compr[sl] = bases[rr];
      cc  = rr;
      rr += 1;
      sl += 1;
    }
  }

  if (compr)
    compr[sl] = 0;

  return(sl);
}



void
decode2bitSequence(uint8 *chunk, uint32 chunkLen, char *seq, uint32 seqLen) {
  uint32       chunkPos = 0;

  if (seq == NULL)
    seq = new char [seqLen + 1];

  for (uint32 ii=0; ii<seqLen; ) {
    if (chunkPos == chunkLen) {
      fprintf(stderr, "decode2bit()-- ran out of chunk (length %u) before end of sequence (at %u out of %u)\n",
              chunkLen, ii, seqLen);
    }
    assert(chunkPos < chunkLen);

    uint8  byte = chunk[chunkPos++];

    if (ii + 4 < seqLen) {
      seq[ii++] = Dacgtn[((byte >> 6) & 0x03)];
      seq[ii++] = Dacgtn[((byte >> 4) & 0x03)];
      seq[ii++] = Dacgtn[((byte >> 2) & 0x03)];
      seq[ii++] = Dacgtn[((byte >> 0) & 0x03)];
    }

    else {
      if (ii < seqLen)  seq[ii++] = Dacgtn[((byte >> 6) & 0x03)];  //  This if is redundant,
      if (ii < seqLen)  seq[ii++] = Dacgtn[((byte >> 4) & 0x03)];  //    but pretty.
      if (ii < seqLen)  seq[ii++] = Dacgtn[((byte >> 2) & 0x03)];
      if (ii < seqLen)  seq[ii++] = Dacgtn[((byte >> 0) & 0x03)];
    }
  }

  seq[seqLen] = 0;
}



uint32
encode2bitSequence(uint8 *&chunk, char *seq, uint32 seqLen) {

  for (uint32 ii=0; ii<seqLen; ii++) {       //  If non-ACGT present, return
    char  base = seq[ii];                    //  0 to indicate we can't encode.

    if ((base != 'a') && (base != 'A') &&
        (base != 'c') && (base != 'C') &&
        (base != 'g') && (base != 'G') &&
        (base != 't') && (base != 'T')) {
      fprintf(stderr, "Invalid base %c detected at position %u\n", base, ii);
      return(0);
    }
  }

  uint32 chunkLen = 0;

  if (chunk == NULL)
    chunk = new uint8 [ seqLen / 4 + 1];

  for (uint32 ii=0; ii<seqLen; ) {
    uint8  byte = 0;

    if (ii + 4 < seqLen) {
      byte  = Eacgtn[seq[ii++]];  byte <<= 2;
      byte |= Eacgtn[seq[ii++]];  byte <<= 2;
      byte |= Eacgtn[seq[ii++]];  byte <<= 2;
      byte |= Eacgtn[seq[ii++]];
    }

    else {
      if (ii < seqLen)  { byte |= Eacgtn[seq[ii++]]; }   byte <<= 2;  //  This if is redundant, but pretty.
      if (ii < seqLen)  { byte |= Eacgtn[seq[ii++]]; }   byte <<= 2;  //  Yes, all three always shift,
      if (ii < seqLen)  { byte |= Eacgtn[seq[ii++]]; }   byte <<= 2;  //  not conditionally.
      if (ii < seqLen)  { byte |= Eacgtn[seq[ii++]]; }
    }

    chunk[chunkLen++] = byte;
  }

  return(chunkLen);
}



void
decode3bitSequence(uint8 *chunk, uint32 chunkLen, char *seq, uint32 seqLen) {
  uint32       chunkPos = 0;

  if (seq == NULL)
    seq = new char [seqLen + 1];

  for (uint32 ii=0; ii<seqLen; ) {
    if (chunkPos == chunkLen) {
      fprintf(stderr, "decode3bit()-- ran out of chunk (length %u) before end of sequence (at %u out of %u)\n",
              chunkLen, ii, seqLen);
    }
    assert(chunkPos < chunkLen);

    uint8  byte = chunk[chunkPos++];
    uint8  c1   = 0;
    uint8  c2   = 0;
    uint8  c3   = 0;

    if (ii + 3 < seqLen) {
      uint8 c1 = byte / 5 / 5;    byte -= c1 * 5 * 5;
      uint8 c2 = byte / 5;        byte -= c2 * 5;
      uint8 c3 = byte;

      seq[ii++] = Dacgtn[c1];
      seq[ii++] = Dacgtn[c2];
      seq[ii++] = Dacgtn[c3];
    }

    else {
      uint8 c1 = byte / 5 / 5;    byte -= c1 * 5 * 5;
      uint8 c2 = byte / 5;        byte -= c2 * 5;
      uint8 c3 = byte;

      if (ii < seqLen)  seq[ii++] = Dacgtn[c1];  //  This if is redundant,
      if (ii < seqLen)  seq[ii++] = Dacgtn[c2];  //    but pretty.
      if (ii < seqLen)  seq[ii++] = Dacgtn[c3];
    }
  }

  seq[seqLen] = 0;
}



uint32
encode3bitSequence(uint8 *&chunk, char *seq, uint32 seqLen) {

  for (uint32 ii=0; ii<seqLen; ii++) {       //  If non-ACGTN present, return
    char  base = seq[ii];                    //  0 to indicate we can't encode.

    if ((base != 'a') && (base != 'A') &&
        (base != 'c') && (base != 'C') &&
        (base != 'g') && (base != 'G') &&
        (base != 't') && (base != 'T') &&
        (base != 'n') && (base != 'N')) {
      fprintf(stderr, "Invalid base %c detected at position %u\n", base, ii);
      return(0);
    }
  }

  uint32 chunkLen = 0;

  if (chunk == NULL)
    chunk = new uint8 [ seqLen / 3 + 1];

  for (uint32 ii=0; ii<seqLen; ) {
    uint8  byte = 0;

    if (ii + 3 < seqLen) {
      byte += Eacgtn[seq[ii++]] * 5 * 5;
      byte += Eacgtn[seq[ii++]] * 5;
      byte += Eacgtn[seq[ii++]];
    }

    else {
      if (ii < seqLen)   byte += Eacgtn[seq[ii++]] * 5 * 5;
      if (ii < seqLen)   byte += Eacgtn[seq[ii++]] * 5;
      if (ii < seqLen)   byte += Eacgtn[seq[ii++]];
    }

    chunk[chunkLen++] = byte;
  }

  return(chunkLen);
}



void
decode8bitSequence(uint8 *chunk, uint32 chunkLen, char *seq, uint32 seqLen) {

  if (seq == NULL)
    seq = new char [seqLen + 1];

  for (uint32 ii=0; ii<seqLen; ii++)
    seq[ii] = chunk[ii];

  seq[seqLen] = 0;
}



uint32
encode8bitSequence(uint8 *&chunk, char *seq, uint32 seqLen) {

  if (chunk == NULL)
    chunk = new uint8 [ seqLen ];

  for (uint32 ii=0; ii<seqLen; ii++)
    chunk[ii] = seq[ii];

  return(seqLen);
}




//  Saves the file offset of the first byte in the record:
//    for FASTA, the '>'
//    for FASTQ, the '@'.

class dnaSeqIndexEntry {
public:
  dnaSeqIndexEntry() {
    _fileOffset     = UINT64_MAX;
    _sequenceLength = 0;
  };
  ~dnaSeqIndexEntry() {
  };

  uint64   _fileOffset;
  uint64   _sequenceLength;
};



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

  generateIndex();
}



dnaSeqFile::~dnaSeqFile() {
  delete    _file;
  delete    _buffer;
  delete [] _index;
}



bool
dnaSeqFile::findSequence(uint64 i) {

  if (_indexLen == 0)   return(false);
  if (_indexLen <= i)   return(false);

  _buffer->seek(_index[i]._fileOffset);

  return(true);
}



uint64
dnaSeqFile::sequenceLength(uint64 i) {

  if (_indexLen == 0)   return(UINT64_MAX);
  if (_indexLen <= i)   return(UINT64_MAX);

  return(_index[i]._sequenceLength);
}




bool
dnaSeqFile::findSequence(const char *name) {
  fprintf(stderr, "dnaSeqFile::findSequence(const char *) not supported.\n");
  exit(1);
  return(false);
}



bool
dnaSeqFile::loadIndex(void) {
  char   indexName[FILENAME_MAX+1];

  snprintf(indexName, FILENAME_MAX, "%s.index", _file->filename());

  if (fileExists(indexName) == false)
    return(false);

  FILE   *indexFile = AS_UTL_openInputFile(indexName);

  loadFromFile(_indexLen, "dnaSeqFile::indexLen", indexFile);

  _index = new dnaSeqIndexEntry [_indexLen];

  loadFromFile(_index, "dnaSeqFile::index", _indexLen, indexFile);

  AS_UTL_closeFile(indexFile, indexName);

  return(true);
}



void
dnaSeqFile::saveIndex(void) {
  char   indexName[FILENAME_MAX+1];

  snprintf(indexName, FILENAME_MAX, "%s.index", _file->filename());

  FILE   *indexFile = AS_UTL_openOutputFile(indexName);

  writeToFile(_indexLen, "dnaSeqFile::indexLen",            indexFile);
  writeToFile(_index,    "dnaSeqFile::index",    _indexLen, indexFile);

  AS_UTL_closeFile(indexFile, indexName);
}



void
dnaSeqFile::generateIndex(void) {
  uint32          nameMax = 0;
  char           *name    = NULL;
  uint64          seqMax  = 0;
  char           *seq     = NULL;
  uint8          *qlt     = NULL;
  uint64          seqLen  = 0;

  if (loadIndex() == true)
    return;

  _indexLen = 0;
  _indexMax = 1048576;
  _index    = new dnaSeqIndexEntry [_indexMax];

  _index[_indexLen]._fileOffset     = _buffer->tell();
  _index[_indexLen]._sequenceLength = 0;

  //  While we read sequences:
  //    update the length of the sequence (we've already save the position)
  //    make space for more sequences
  //    save the position of the next sequence
  //
  while (loadSequence(name, nameMax, seq, qlt, seqMax, seqLen) == true) {
    _index[_indexLen]._sequenceLength = seqLen;

    increaseArray(_index, _indexLen, _indexMax, 1048576);

    _indexLen++;

    _index[_indexLen]._fileOffset     = _buffer->tell();
    _index[_indexLen]._sequenceLength = 0;
  }

  //for (uint32 ii=0; ii<_indexLen; ii++)
  //  fprintf(stderr, "%u offset %lu length %lu\n", ii, _index[ii]._fileOffset, _index[ii]._sequenceLength);

  if (_indexLen > 0)
    saveIndex();
}



uint64
dnaSeqFile::loadFASTA(char   *&name,     uint32  &nameMax,
                      char   *&seq,
                      uint8  *&qlt,      uint64  &seqMax) {
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
dnaSeqFile::loadFASTQ(char   *&name,     uint32  &nameMax,
                      char   *&seq,
                      uint8  *&qlt,      uint64  &seqMax) {
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
dnaSeqFile::loadSequence(char   *&name,     uint32  &nameMax,
                         char   *&seq,
                         uint8  *&qlt,      uint64  &seqMax,
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
                      uint64  &seqLength,
                      bool    &endOfSequence) {

  seqLength     = 0;
  endOfSequence = false;

  if (_buffer->eof() == true)
    return(false);

  //  If this is a new file, skip the first name line.

  if (_buffer->tell() == 0) {
    while (_buffer->peek() == '\n')    //  Skip whitespace before the first name line.
      _buffer->read();

    _buffer->skipAhead('\n', true);
  }

  //  Skip whitespace.

  while (_buffer->peek() == '\n')
    _buffer->read();

  //  If now at EOF, that's it.

  if  (_buffer->eof() == true)
    return(false);

  //  Otherwise, we must be in the middle of sequence, so load
  //  until we're not in sequence or out of space.

  while (_buffer->eof() == false) {

    //  If we're at the start of a new sequence, skip over any QV's and
    //  the next name line, set endOfSequence and return.

    if (_buffer->peek() == '>') {
      _buffer->skipAhead('\n', true);      //  Skip the name line.
      endOfSequence = true;
      return(true);
    }

    if (_buffer->peek() == '+') {
      _buffer->skipAhead('\n', true);      //  Skip the + line.
      _buffer->skipAhead('\n', true);      //  Skip the QV line.
      _buffer->skipAhead('\n', true);      //  Skip the @ line for the next sequence.
      endOfSequence = true;
      return(true);
    }

    //  Read some bases.

    seqLength += _buffer->copyUntil('\n', seq + seqLength, maxLength - seqLength);

    if (seqLength == maxLength)
      return(true);

    //  We're at a newline (or end of file), either way, suck in the next letter
    //  (or nothing) and keep going.

    _buffer->read();
  }

  //  We hit EOF.  If there are bases loaded, then we're at the end of
  //  a sequence, and should return that we loaded bases.

  endOfSequence = (seqLength > 0);

  return(endOfSequence);
}
