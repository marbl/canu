
#include "seqStore.H"
#include "seqCache.H"
#include "alphabet.h"

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

  _index = new seqStoreIndex [_header._numberOfSequences];
  _block = new seqStoreBlock [_header._numberOfBlocks];
  _names = new char          [_header._namesLength];

  fseeko(F, _header._indexStart, SEEK_SET);
  fread( _index,   sizeof(seqStoreIndex), _header._numberOfSequences, F);

  fseeko(F, _header._blockStart, SEEK_SET);
  fread( _block,   sizeof(seqStoreBlock), _header._numberOfBlocks, F);

  fseeko(F, _header._namesStart, SEEK_SET);
  fread( _names,   sizeof(char),          _header._namesLength, F);

  if (errno)
    fprintf(stderr, "seqStore::seqStore()--  Failed to read index from '%s': %s\n",
            _filename, strerror(errno)), exit(1);

  fclose(F);

  _bpf = new bitPackedFile(_filename, sizeof(seqStoreHeader));

  _numberOfSequences = _header._numberOfSequences;
}



seqStore::seqStore() {
  clear();
}



seqStore::~seqStore() {
  delete    _bpf;
  delete [] _index;
  delete [] _names;
}



seqFile *
seqStore::openFile(const char *filename) {
  u64bit         magic1, magic2;
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
  fread(&magic1, sizeof(u64bit), 1, F);
  fread(&magic2, sizeof(u64bit), 1, F);
  fclose(F);
  if ((magic1 != SEQSTORE_MAGICNUMBER1) || (magic2 != SEQSTORE_MAGICNUMBER2))
    return(0L);

  return(new seqStore(filename));
}



u32bit
seqStore::find(const char *sequencename) {
  char   *ptr = _names;

  //  If this proves far too slow, rewrite the _names string to
  //  separate IDs with 0xff, then use strstr on the whole thing.  To
  //  find the ID, scan down the string counting the number of 0xff's.
  //
  //  Similar code is used for fastaFile::find()
  //
  for (u32bit iid=0; iid < _header._numberOfSequences; iid++) {
    //fprintf(stderr, "seqStore::find()-- '%s' vs '%s'\n", sequencename, ptr);
    if (strcmp(sequencename, ptr) == 0)
      return(iid);

    while (*ptr)
      ptr++;
    ptr++;
  }

  return(~u32bitZERO);
}



u32bit
seqStore::getSequenceLength(u32bit iid) {
  return((iid < _header._numberOfSequences) ? _index[iid]._seqLength : 0);
}



bool
seqStore::getSequence(u32bit iid,
                      char *&h, u32bit &hLen, u32bit &hMax,
                      char *&s, u32bit &sLen, u32bit &sMax) {

  if (iid >= _header._numberOfSequences) {
    fprintf(stderr, "seqStore::getSequence(full)--  iid "u32bitFMT" more than number of sequences "u32bitFMT"\n",
            iid, _header._numberOfSequences);
    return(false);
  }

  if (sMax == 0)
    s = 0L;

  if (hMax == 0)
    h = 0L;

  if (sMax < _index[iid]._seqLength) {
    sMax = _index[iid]._seqLength + 1024;
    delete [] s;
    s = new char [sMax];
  }

  if (hMax < _index[iid]._hdrLength) {
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

  _bpf->seek(_index[iid]._seqPosition * 2);

  u32bit seqLen = _index[iid]._seqLength;
  u32bit block  = _index[iid]._block;

  while (sLen < seqLen) {
    assert(_bpf->tell() == _block[block]._bpf * 2);
    assert(sLen == _block[block]._pos);

    if (_block[block]._isACGT == 0) {
      memset(s + sLen, 'N', _block[block]._len);
      sLen += _block[block]._len;
    } else {
      for (u32bit xx=0; xx<_block[block]._len; xx++) {
        s[sLen++] = bitsToLetter[_bpf->getBits(2)];
      }
    }

    block++;
  }

  s[sLen] = 0;
  
  return(true);
}



bool
seqStore::getSequence(u32bit iid,
                      u32bit bgn, u32bit end, char *s) {

  if (iid >= _header._numberOfSequences) {
    fprintf(stderr, "seqStore::getSequence(part)--  iid "u32bitFMT" more than number of sequences "u32bitFMT"\n",
            iid, _header._numberOfSequences);
    return(false);
  }

  if (bgn >= end) {
    fprintf(stderr, "seqStore::getSequence(part)--  for iid "u32bitFMT"; invalid bgn="u32bitFMT" end="u32bitFMT"; seqLen="u32bitFMT"\n",
            iid, bgn, end, _index[iid]._seqLength);
    return(false);
  }

  //  Copy the defline into h

  //memcpy(h, _names + _index[iid]._hdrPosition, _index[iid]._hdrLength);
  //h[_index[iid]._hdrLength] = 0;

  //  Decode and copy the sequence into s

  u32bit block  = _index[iid]._block;
  u32bit sLen   = 0;  //  length of sequence we've handled so far
  u32bit sPos   = 0;  //  length of sequence we've copied (== sLen - bgn)

  //  Skip blocks before we care.
  //
  while (sLen + _block[block]._len < bgn) {
    sLen += _block[block]._len;
    block++;
  }

  assert(sLen == _block[block]._pos);

  //  Handle the partial block.  Make sure the stupid user didn't ask
  //  for something less than a block size.

  u32bit partLen = _block[block]._pos + _block[block]._len - bgn;

  if (partLen < end - bgn)
    partLen = end - bgn;

  if (_block[block]._isACGT == 0) {
    memset(s, 'N', partLen);
    sPos  = partLen;
    sLen += partLen;

    _bpf->seek(_block[block+1]._bpf * 2);
  } else {
    _bpf->seek((_block[block]._bpf + bgn - _block[block]._pos) * 2);

    for (u32bit xx=0; xx<partLen; xx++)
      s[sLen++] = bitsToLetter[_bpf->getBits(2)];
  }

  block++;

  while (sLen < end) {
    assert(_bpf->tell() == _block[block]._bpf * 2);
    assert(sLen == _block[block]._pos);

    if (_block[block]._isACGT == 0) {
      memset(s + sLen, 'N', _block[block]._len);
      sLen += _block[block]._len;
    } else {
      for (u32bit xx=0; xx<_block[block]._len; xx++)
        s[sLen++] = bitsToLetter[_bpf->getBits(2)];
    }

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
  _names = 0L;
}



void
constructSeqStore(char *filename, seqCache *inputseq) {

  fprintf(stderr, "constructSeqStore()-- constructing seqStore '%s' from seqCache '%s' of type '%s'.\n",
          filename, inputseq->getSourceName(), inputseq->getFileTypeName());

  seqStoreHeader    HEAD;
  memset(&HEAD, sizeof(seqStoreHeader), 0);

  bitPackedFile    *DATA    = new bitPackedFile(filename, sizeof(seqStoreHeader), true);

  seqStoreIndex    *INDX    = new seqStoreIndex [inputseq->getNumberOfSequences()];

  u32bit            BLOKmax = 4 * inputseq->getNumberOfSequences();
  u32bit            BLOKlen = 0;
  seqStoreBlock    *BLOK    = new seqStoreBlock [BLOKmax];

  u32bit            NAMEmax = 32 * 1024 * 1024;
  u32bit            NAMElen = 0;
  char             *NAME    = new char [NAMEmax];

  u64bit            nACGT      = 0;
  u32bit            nBlockACGT = 0;
  u32bit            nBlockGAP  = 0;

  for (u32bit iid=0; iid<inputseq->getNumberOfSequences(); iid++) {
    seqInCore     *sic = inputseq->getSequenceInCore(iid);
    char          *seq = sic->sequence();

    seqStoreBlock  b;

    if (seq) {
      INDX[iid]._hdrPosition = NAMElen;
      INDX[iid]._hdrLength   = sic->headerLength();
      INDX[iid]._seqPosition = DATA->tell() / 2;
      INDX[iid]._seqLength   = sic->sequenceLength();
      INDX[iid]._block       = BLOKlen;

      if (NAMElen + sic->headerLength() + 1 > NAMEmax) {
        NAMEmax += 32 * 1024 * 1024;
        char *nm = new char [NAMEmax];
        memcpy(nm, NAME, sizeof(char) * NAMElen);
        delete [] NAME;
        NAME = nm;
      }
      strcpy(NAME + NAMElen, sic->header());
      NAMElen += sic->headerLength() + 1;

      //fprintf(stderr, "name: '%s'\n", sic->header());

      b._isACGT = 0;
      b._iid    = sic->getIID();
      b._pos    = 0;
      b._len    = 0;
      b._bpf    = DATA->tell() / 2;

      for (u32bit p=0; p<sic->sequenceLength(); p++) {
        u64bit   bits = letterToBits[seq[p]];

        if (bits == 0xff) {

          //  Letter is NOT ACGI, write out the last block if it was
          //  ACGT, then increment our length.

          if (b._isACGT == 1) {
            if (b._len > 0) {
              nBlockACGT++;
              nACGT += b._len;
              BLOK[BLOKlen++] = b;
            }

            b._isACGT = 0;
            b._iid    = sic->getIID();
            b._pos    = p;
            b._len    = 0;
            b._bpf    = DATA->tell() / 2;
          }

          b._len++;

        } else {

          //  Letter is ACGT.  Write out last block if it was a gap,
          //  then emit this letter.

          if (b._isACGT == 0) {
            if (b._len > 0) {
              nBlockGAP++;
              BLOK[BLOKlen++] = b;
            }

            b._isACGT = 1;
            b._iid    = sic->getIID();
            b._pos    = p;
            b._len    = 0;
            b._bpf    = DATA->tell() / 2;
          }

          b._len++;

          DATA->putBits(bits, 2);
        }
      }

      //  Emit the last block

      if (b._isACGT == 1) {
        nBlockACGT++;
        nACGT += b._len;
      } else {
        nBlockGAP++;
      }

      BLOK[BLOKlen++] = b;
    }

    delete seq;
  }

  //  And a sentinel EOF block -- gets the last position in the file,
  //  useful for the binary search.

  BLOK[BLOKlen]._isACGT = 0;
  BLOK[BLOKlen]._iid    = u32bitMASK(32);
  BLOK[BLOKlen]._pos    = u32bitMASK(31);
  BLOK[BLOKlen]._len    = 0;
  BLOK[BLOKlen]._bpf    = DATA->tell() / 2;

  BLOKlen++;

  //  Update the header, assemble the final file.

  delete DATA;

  HEAD._magic[0]           = SEQSTORE_MAGICNUMBER1;
  HEAD._magic[1]           = SEQSTORE_MAGICNUMBER2;
  HEAD._pad                = u32bitZERO;
  HEAD._numberOfSequences  = inputseq->getNumberOfSequences();
  HEAD._numberOfACGT       = nACGT;
  HEAD._numberOfBlocksACGT = nBlockACGT;
  HEAD._numberOfBlocksGAP  = nBlockGAP;
  HEAD._numberOfBlocks     = BLOKlen;
  HEAD._namesLength        = NAMElen;
  HEAD._indexStart         = u64bitZERO;
  HEAD._blockStart         = u64bitZERO;
  HEAD._namesStart         = u64bitZERO;

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

  fprintf(stderr, "constructSeqStore()-- seqStore '%s' constructed ("u32bitFMT" sequences, "u64bitFMT" ACGT letters, "u32bitFMT" ACGT blocks, "u32bitFMT" GAP blocks).\n",
          filename, HEAD._numberOfSequences, HEAD._numberOfACGT, HEAD._numberOfBlocksACGT, HEAD._numberOfBlocksGAP);
}
