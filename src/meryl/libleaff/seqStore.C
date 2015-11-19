
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
 *    Brian P. Walenz on 2014-DEC-08
 *      are Copyright 2014 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2015-OCT-29
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "seqStore.H"
#include "seqCache.H"
#include "dnaAlphabets.H"
#include "speedCounter.H"

//  Says 'kmerSeqStoreFile'
#define SEQSTORE_MAGICNUMBER1  0x5371655372656d6bULL
#define SEQSTORE_MAGICNUMBER2  0x656c694665726f74ULL


seqStore::seqStore(const char *filename) {
  clear();

  strcpy(_filename, filename);

  errno = 0;
  FILE *F = fopen(_filename, "r");
  if (errno)
    fprintf(stderr, "seqStore::seqStore()--  Failed to open '%s': %s\n",
            _filename, strerror(errno)), exit(1);
  fread(&_header,   sizeof(seqStoreHeader), 1, F);
  fclose(F);

  //_indexBPF = new bitPackedFile(_filename, _header._indexStart);
  //_blockBPF = new bitPackedFile(_filename, _header._blockStart);
  //_namesBPF = new bitPackedFile(_filename, _header._namesStart);

  _bpf      = new bitPackedFile(_filename, sizeof(seqStoreHeader));

  _numberOfSequences = _header._numberOfSequences;
}



seqStore::seqStore() {
  clear();
}



seqStore::~seqStore() {
  //if ((_filename) && (_filename[0] != 0))
  //  fprintf(stderr, "Closing seqStore '%s'\n", _filename);
  delete    _bpf;
  delete [] _index;
  delete [] _block;
  delete [] _names;
  delete    _indexBPF;
  delete    _blockBPF;
  delete    _namesBPF;
}



seqFile *
seqStore::openFile(const char *filename) {
  uint64         magic1, magic2;
  struct stat    st;

  errno = 0;
  stat(filename, &st);
  if (errno)
    return(0L);
  if ((st.st_mode & S_IFREG) == 0)
    return(0L);

  //  Check the magic.  Fail if not correct.

  errno = 0;
  FILE *F = fopen(filename, "r");
  if (errno)
    return(0L);
  fread(&magic1, sizeof(uint64), 1, F);
  fread(&magic2, sizeof(uint64), 1, F);
  fclose(F);
  if ((magic1 != SEQSTORE_MAGICNUMBER1) || (magic2 != SEQSTORE_MAGICNUMBER2))
    return(0L);

  return(new seqStore(filename));
}



//  If this proves far too slow, rewrite the _names string to separate IDs with 0xff, then use
//  strstr on the whole thing.  To find the ID, scan down the string counting the number of 0xff's.
//
//  Similar code is used for fastaFile::find()
//
uint32
seqStore::find(const char *sequencename) {

  if (_names == NULL)
    loadIndex();

  char   *ptr = _names;

  for (uint32 iid=0; iid < _header._numberOfSequences; iid++) {
    if (strcmp(sequencename, ptr) == 0)
      return(iid);

    while (*ptr)
      ptr++;
    ptr++;
  }

  return(~uint32ZERO);
}



uint32
seqStore::getSequenceLength(uint32 iid) {
  if (_index == NULL)
    loadIndex();
  return((iid < _header._numberOfSequences) ? _index[iid]._seqLength : 0);
}



bool
seqStore::getSequence(uint32 iid,
                      char *&h, uint32 &hLen, uint32 &hMax,
                      char *&s, uint32 &sLen, uint32 &sMax) {

  if (_index == NULL)
    loadIndex();

  if (iid >= _header._numberOfSequences) {
    fprintf(stderr, "seqStore::getSequence(full)--  iid "F_U32" more than number of sequences "F_U32"\n",
            iid, _header._numberOfSequences);
    return(false);
  }

  if (sMax == 0)  s = 0L;  //  So the delete below doesn't bomb
  if (hMax == 0)  h = 0L;

  if (sMax < _index[iid]._seqLength + 1) {
    sMax = _index[iid]._seqLength + 1024;
    delete [] s;
    s = new char [sMax];
  }

  if (hMax < _index[iid]._hdrLength + 1) {
    hMax = _index[iid]._hdrLength + 1024;
    delete [] h;
    h = new char [hMax];
  }

  hLen = 0;
  sLen = 0;

  //  Copy the defline into h

  memcpy(h, _names + _index[iid]._hdrPosition, _index[iid]._hdrLength);

  h[_index[iid]._hdrLength] = 0;

  //  Decode and copy the sequence into s

  uint32 seqLen  = _index[iid]._seqLength;
  uint32 block   = _index[iid]._block;
  uint64 seekpos = _index[iid]._seqPosition * 2;

  _bpf->seek(seekpos);

  while (sLen < seqLen) {
    assert(_bpf->tell() == _block[block]._bpf * 2);
    assert(sLen == _block[block]._pos);

    if (_block[block]._isACGT == 0) {
      memset(s + sLen, 'N', _block[block]._len);
      sLen += _block[block]._len;
    } else {
      for (uint32 xx=0; xx<_block[block]._len; xx++) {
        s[sLen++] = alphabet.bitsToLetter(_bpf->getBits(2));
      }
    }

    block++;
  }

  s[sLen] = 0;

  return(true);
}



bool
seqStore::getSequence(uint32 iid,
                      uint32 bgn, uint32 end, char *s) {

  if (_index == NULL)
    loadIndex();

  if (iid >= _header._numberOfSequences) {
    fprintf(stderr, "seqStore::getSequence(part)--  iid "F_U32" more than number of sequences "F_U32"\n",
            iid, _header._numberOfSequences);
    return(false);
  }

  if (bgn >= end) {
    fprintf(stderr, "seqStore::getSequence(part)--  for iid "F_U32"; invalid bgn="F_U32" end="F_U32"; seqLen="F_U32"\n",
            iid, bgn, end, _index[iid]._seqLength);
    return(false);
  }

  //  Decode and copy the sequence into s

  uint32 block  = _index[iid]._block;
  uint32 sLen   = 0;  //  length of sequence we've copied
  uint32 sPos   = 0;  //  position in the sequence

  //  Skip blocks before we care.
  //
  while (sPos + _block[block]._len < bgn) {
    sPos += _block[block]._len;
    block++;
  }

  assert(sPos == _block[block]._pos);

  //  Move into the block (we could just set sPos = bgn...).
  sPos += bgn - _block[block]._pos;

  //  Handle the partial block.  Copy what is left in the block, or
  //  the requested size, whichever is smaller.

  uint32 partLen = MIN((_block[block]._pos + _block[block]._len - bgn),
                       (end - bgn));

  if (_block[block]._isACGT == 0) {
    memset(s, 'N', partLen);
    sLen += partLen;
    _bpf->seek(_block[block+1]._bpf * 2);
  } else {
    _bpf->seek((_block[block]._bpf + bgn - _block[block]._pos) * 2);

    for (uint32 xx=0; xx<partLen; xx++)
      s[sLen++] = alphabet.bitsToLetter(_bpf->getBits(2));
  }

  sPos += partLen;

  block++;

  while (sPos < end) {
    assert(_bpf->tell() == _block[block]._bpf * 2);
    assert(sPos == _block[block]._pos);

    //  Like the partial block above, pick how much to copy as the
    //  smaller of the block size and what is left to fill.

    partLen = MIN((_block[block]._len), (end - sPos));

    if (_block[block]._isACGT == 0) {
      memset(s + sLen, 'N', partLen);
      sLen += partLen;
    } else {
      for (uint32 xx=0; xx<partLen; xx++)
        s[sLen++] = alphabet.bitsToLetter(_bpf->getBits(2));
    }

    sPos += partLen;

    block++;
  }

  s[sLen] = 0;

  return(true);
}



void
seqStore::clear(void) {
  memset(_filename, 0, FILENAME_MAX);
  memset(_typename, 0, FILENAME_MAX);

  strcpy(_typename, "seqStore");

  _numberOfSequences = 0;

  _bpf = 0L;

  memset(&_header, 0, sizeof(seqStoreHeader));

  _index = 0L;
  _block = 0L;
  _names = 0L;

  _indexBPF = 0L;
  _blockBPF = 0L;
  _namesBPF = 0L;

  _lastIIDloaded = ~uint32ZERO;
}



void
seqStore::loadIndex(void) {

  if (_index)
    return;

  delete _indexBPF;  _indexBPF = 0L;
  delete _blockBPF;  _blockBPF = 0L;
  delete _namesBPF;  _namesBPF = 0L;

  errno = 0;
  FILE *F = fopen(_filename, "r");
  if (errno)
    fprintf(stderr, "seqStore::seqStore()--  Failed to open '%s': %s\n",
            _filename, strerror(errno)), exit(1);

  fread(&_header,   sizeof(seqStoreHeader), 1, F);

  //fprintf(stderr, "seqStore::seqStore()--  Allocating space for "F_U32" sequences ("F_U64"MB)\n", _header._numberOfSequences, _header._numberOfSequences * sizeof(seqStoreIndex) / 1024 / 1024);
  //fprintf(stderr, "seqStore::seqStore()--  Allocating space for "F_U32" blocks    ("F_U64"MB)\n", _header._numberOfBlocks,    _header._numberOfBlocks    * sizeof(seqStoreBlock) / 1024 / 1024);
  //fprintf(stderr, "seqStore::seqStore()--  Allocating space for "F_U32" labels    ("F_U64"MB)\n", _header._namesLength,       _header._namesLength       * sizeof(char)          / 1024 / 1024);

  _index = new seqStoreIndex [_header._numberOfSequences];
  _block = new seqStoreBlock [_header._numberOfBlocks];
  _names = new char          [_header._namesLength];

  fseeko(F, _header._indexStart, SEEK_SET);
  fread( _index,   sizeof(seqStoreIndex), _header._numberOfSequences, F);

#if 0
  for (uint32 i=0; i<_header._numberOfSequences; i++)
    fprintf(stderr, "IDX[%4u] hdrPos=%u hdrLen=%u seqPos=%llu seqLen=%u block=%u\n",
            i,
            _index[i]._hdrPosition,
            _index[i]._hdrLength,
            _index[i]._seqPosition,
            _index[i]._seqLength,
            _index[i]._block);
#endif

  fseeko(F, _header._blockStart, SEEK_SET);
  fread( _block,   sizeof(seqStoreBlock), _header._numberOfBlocks, F);

  fseeko(F, _header._namesStart, SEEK_SET);

  fread( _names,   sizeof(char),          _header._namesLength, F);
  if (errno)
    fprintf(stderr, "seqStore::seqStore()--  Failed to read index from '%s': %s\n",
            _filename, strerror(errno)), exit(1);

  fclose(F);
}













static
void
addSeqStoreBlock(uint32          &BLOKmax,
                 uint32          &BLOKlen,
                 seqStoreBlock*  &BLOK,
                 seqStoreBlock   &b,
                 uint32          &nBlockACGT,
                 uint32          &nBlockGAP,
                 uint64          &nACGT) {

  //fprintf(stderr, "addSeqStoreBlock()-- BLOK max=%u len=%u ACGT=%u GAP=%u nACGT=%lu\n",
  //        BLOKmax, BLOKlen, nBlockACGT, nBlockGAP, nACGT);

  if (b._len == 0)
    return;

  if (b._isACGT == 1) {
    nBlockACGT++;
    nACGT += b._len;
  } else {
    nBlockGAP++;
  }

  BLOK[BLOKlen++] = b;

  if (BLOKlen >= BLOKmax) {
    BLOKmax *= 2;
    seqStoreBlock *nb = new seqStoreBlock [BLOKmax];
    memcpy(nb, BLOK, BLOKlen * sizeof(seqStoreBlock));
    delete [] BLOK;
    BLOK = nb;
  }
}



void
constructSeqStore(char *filename, seqCache *inputseq) {

  fprintf(stderr, "constructSeqStore()-- constructing seqStore '%s' from seqCache '%s' of type '%s'.\n",
          filename, inputseq->getSourceName(), inputseq->getFileTypeName());

  seqStoreHeader    HEAD;
  memset(&HEAD, sizeof(seqStoreHeader), 0);

  bitPackedFile    *DATA    = new bitPackedFile(filename, sizeof(seqStoreHeader), true);

  uint32            INDXmax = 1048576;
  seqStoreIndex    *INDX    = new seqStoreIndex [INDXmax];

  uint32            BLOKmax = 1048576;
  uint32            BLOKlen = 0;
  seqStoreBlock    *BLOK    = new seqStoreBlock [BLOKmax];

  uint32            NAMEmax = 32 * 1024 * 1024;
  uint32            NAMElen = 0;
  char             *NAME    = new char [NAMEmax];

  seqInCore        *sic        = inputseq->getSequenceInCore();

  uint64            nACGT      = 0;
  uint32            nBlockACGT = 0;
  uint32            nBlockGAP  = 0;
  uint32            nSequences = 0;

  speedCounter      C(" reading sequences %7.0f sequences -- %5.0f sequences/second\r", 1.0, 0x1ffff, true);

  while (sic != NULL) {
    if (sic->sequence()) {
      char          *seq = sic->sequence();
      seqStoreBlock  b;

      if (nSequences >= INDXmax) {
        seqStoreIndex *I = new seqStoreIndex[INDXmax * 2];
        memcpy(I, INDX, sizeof(seqStoreIndex) * nSequences);
        delete [] INDX;
        INDXmax *= 2;
        INDX     = I;
      }

      INDX[nSequences]._hdrPosition = NAMElen;
      INDX[nSequences]._hdrLength   = sic->headerLength();
      INDX[nSequences]._seqPosition = DATA->tell() / 2;
      INDX[nSequences]._seqLength   = sic->sequenceLength();
      INDX[nSequences]._block       = BLOKlen;

#if 0
      fprintf(stderr, "ADD SEQUENCE hdr pos=%u len=%u seq pos=%u len=%u blok=%u\n",
              INDX[nSequences]._hdrPosition,
              INDX[nSequences]._hdrLength,
              INDX[nSequences]._seqPosition,
              INDX[nSequences]._seqLength,
              INDX[nSequences]._block);
#endif

#if SEQSTOREBLOCK_MAXPOS < uint64MASK(32)
      if (sic->sequenceLength() > SEQSTOREBLOCK_MAXPOS)
        fprintf(stderr, "constructSeqStore()-- sequence %s too long, must be shorter than "F_U64" Gbp.\n",
                sic->header(), SEQSTOREBLOCK_MAXPOS / 1024 / 1024 / 1024), exit(1);
#endif

#if SEQSTOREBLOCK_MAXIID < uint64MASK(32)
      if (sic->getIID() > SEQSTOREBLOCK_MAXPOS)
        fprintf(stderr, "constructSeqStore()-- too many sequences, must be fewer than "F_U64".\n",
                SEQSTOREBLOCK_MAXIID), exit(1);
#endif

      if (NAMElen + sic->headerLength() + 1 > NAMEmax) {
        NAMEmax += 32 * 1024 * 1024;
        char *nm = new char [NAMEmax];
        memcpy(nm, NAME, sizeof(char) * NAMElen);
        delete [] NAME;
        NAME = nm;
      }
      strcpy(NAME + NAMElen, sic->header());
      NAMElen += sic->headerLength() + 1;

      b._isACGT = 0;
      b._iid    = sic->getIID();
      b._pos    = 0;
      b._len    = 0;
      b._bpf    = DATA->tell() / 2;

      for (uint32 p=0; p<sic->sequenceLength(); p++) {
        uint64   bits = alphabet.letterToBits(seq[p]);

        //  If the length of the current block is too big (which would
        //  soon overflow the bit field storing length) write out a
        //  block and reset the length.
        //
        if (b._len == SEQSTOREBLOCK_MAXLEN) {
          addSeqStoreBlock(BLOKmax, BLOKlen, BLOK, b, nBlockACGT, nBlockGAP, nACGT);

          b._pos    = p;
          b._len    = 0;
          b._bpf    = DATA->tell() / 2;
        }


        if (bits == 0xff) {
          //  This letter is NOT ACGT.  If the current block is an ACGT block, write it
          //  and reset.
          //
          if (b._isACGT == 1) {
            addSeqStoreBlock(BLOKmax, BLOKlen, BLOK, b, nBlockACGT, nBlockGAP, nACGT);

            b._isACGT = 0;
            b._iid    = sic->getIID();
            b._pos    = p;
            b._len    = 0;
            b._bpf    = DATA->tell() / 2;
          }

        } else {

          //  This letter is ACGT.  If the current block is NOT an ACGT block, write it
          //  and reset.
          //
          if (b._isACGT == 0) {
            addSeqStoreBlock(BLOKmax, BLOKlen, BLOK, b, nBlockACGT, nBlockGAP, nACGT);

            b._isACGT = 1;
            b._iid    = sic->getIID();
            b._pos    = p;
            b._len    = 0;
            b._bpf    = DATA->tell() / 2;
          }
        }

        //  Always add one to the length of the current block, and
        //  write out the base if the letter is ACGT.
        //
        b._len++;

        if (bits != 0xff)
          DATA->putBits(bits, 2);
      }

      //  Emit the last block
      //
      addSeqStoreBlock(BLOKmax, BLOKlen, BLOK, b, nBlockACGT, nBlockGAP, nACGT);
    }

    //  If there is no sequence, the index record for this sequence is left blank.
    //
    nSequences++;

    C.tick();

    delete sic;
    sic = inputseq->getSequenceInCore();
  }

  //  And a sentinel EOF block -- gets the last position in the file,
  //  useful for the binary search.  We always have a space block at
  //  the end of the list, but we don't care if we just used the last
  //  block (and so we don't bother to reallocate the array if it is
  //  full).

  BLOK[BLOKlen]._isACGT = 0;
  BLOK[BLOKlen]._iid    = uint32MASK(32);
  BLOK[BLOKlen]._pos    = uint32MASK(31);
  BLOK[BLOKlen]._len    = 0;
  BLOK[BLOKlen]._bpf    = DATA->tell() / 2;

  BLOKlen++;

  //  Update the header, assemble the final file.

  delete DATA;

  HEAD._magic[0]           = SEQSTORE_MAGICNUMBER1;
  HEAD._magic[1]           = SEQSTORE_MAGICNUMBER2;
  HEAD._pad                = uint32ZERO;
  HEAD._numberOfSequences  = nSequences;
  HEAD._numberOfACGT       = nACGT;
  HEAD._numberOfBlocksACGT = nBlockACGT;
  HEAD._numberOfBlocksGAP  = nBlockGAP;
  HEAD._numberOfBlocks     = BLOKlen;
  HEAD._namesLength        = NAMElen;
  HEAD._indexStart         = uint64ZERO;
  HEAD._blockStart         = uint64ZERO;
  HEAD._namesStart         = uint64ZERO;

  errno = 0;
  FILE *F = fopen(filename, "r+");
  if (errno)
    fprintf(stderr, "constructSeqStore()--  Failed to reopen '%s' to write data: %s\n",
            filename, strerror(errno)), exit(1);

  fseeko(F, 0, SEEK_END);
  HEAD._indexStart = ftello(F);
  fwrite(INDX, sizeof(seqStoreIndex), HEAD._numberOfSequences, F);

  fseeko(F, 0, SEEK_END);
  HEAD._blockStart = ftello(F);
  fwrite(BLOK, sizeof(seqStoreBlock), HEAD._numberOfBlocks, F);

  fseeko(F, 0, SEEK_END);
  HEAD._namesStart = ftello(F);
  fwrite(NAME, sizeof(char), HEAD._namesLength, F);

  fseeko(F, 0, SEEK_SET);
  fwrite(&HEAD, sizeof(seqStoreHeader), 1, F);

  fclose(F);

  if (errno)
    fprintf(stderr, "constructSeqStore()--  Failed to write data to '%s': %s\n",
            filename, strerror(errno)), exit(1);

  delete [] INDX;
  delete [] BLOK;
  delete [] NAME;

  //  ESTmapper depends on this output.

  fprintf(stderr, "constructSeqStore()-- seqStore '%s' constructed ("F_U32" sequences, "F_U64" ACGT letters, "F_U32" ACGT blocks, "F_U32" GAP blocks).\n",
          filename, HEAD._numberOfSequences, HEAD._numberOfACGT, HEAD._numberOfBlocksACGT, HEAD._numberOfBlocksGAP);
}
