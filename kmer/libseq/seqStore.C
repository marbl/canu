
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

  fread(&_header,   sizeof(seqStoreIndex), 1, F);

  _index = new seqStoreIndex [_header._numberOfSequences];
  _names = new char          [_header._namesLength];

  fseeko(F, _header._indexStart, SEEK_SET);
  fread( _index,   sizeof(seqStoreIndex), _header._numberOfSequences, F);

  fseeko(F, _header._blockStart, SEEK_SET);
  fread( _block,   sizeof(seqStoreBlock), _header._numberOfBlocks, F);

  fseeko(F, _header._namesStart, SEEK_SET);
  fread( _names,   sizeof(char),          _header._namesLength,       F);

  if (errno)
    fprintf(stderr, "seqStore::seqStore()--  Failed to read index from '%s': %s\n",
            _filename, strerror(errno)), exit(1);

  fclose(F);

  _rb = new readBuffer(_filename);

  _numberOfSequences = _header._numberOfSequences;
}



seqStore::seqStore() {
  clear();
}



seqStore::~seqStore() {
  delete    _rb;
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
#warning unimplemented
  fprintf(stderr, "seqStore::find()--  Not implemented.\n");
  exit(1);
}



u32bit
seqStore::getSequenceLength(u32bit iid) {
  return((iid < _header._numberOfSequences) ? _index[iid]._seqLength : 0);
}



bool
seqStore::getSequence(u32bit id,
                      char *&h, u32bit &hLen, u32bit &hMax,
                      char *&s, u32bit &sLen, u32bit &sMax) {
#warning unimplemented
  fprintf(stderr, "seqStore::getSequence()-- Not implemented.\n");
  exit(1);
  return(false);
}



bool
seqStore::getSequence(u32bit iid,
                      u32bit bgn, u32bit end, char *s) {
#warning unimplemented
  fprintf(stderr, "seqStore::getSequence()-- Not implemented.\n");
  exit(1);
  return(false);
}



void
seqStore::clear(void) {
  memset(_filename, 0, FILENAME_MAX);
  memset(_typename, 0, FILENAME_MAX);

  strcpy(_typename, "seqStore");

  _numberOfSequences = 0;

  _rb = 0L;

  memset(&_header, 0, sizeof(seqStoreIndex));

  _index = 0L;
  _names = 0L;
}



void
seqStore::construct(char *filename, seqCache *inputseq) {

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

  for (u32bit iid=0; iid<inputseq->getNumberOfSequences(); iid++) {
    seqInCore     *sic = inputseq->getSequenceInCore(iid);
    char          *seq = sic->sequence();

    seqStoreBlock *b   = BLOK + BLOKlen;

    if (seq) {
      INDX[iid]._hdrPosition = NAMElen;
      INDX[iid]._hdrLength   = sic->headerLength();
      INDX[iid]._seqPosition = DATA->tell();
      INDX[iid]._seqLength   = sic->sequenceLength();

      if (NAMElen + sic->headerLength() + 1 > NAMEmax) {
        NAMEmax += 32 * 1024 * 1024;
        char *nm = new char [NAMEmax];
        memcpy(nm, NAME, sizeof(char) * NAMElen);
        delete [] NAME;
        NAME = nm;
      }
      strcpy(NAME + NAMElen, sic->header());
      NAMElen += sic->headerLength() + 1;

      b->_isACGT = 0;
      b->_iid    = sic->getIID();
      b->_pos    = 0;
      b->_len    = 0;
      b->_bpf    = DATA->tell();

      for (u32bit p=0; p<sic->sequenceLength(); p++) {
        u64bit   bits = letterToBits[seq[p]];

        if (bits == 0xff) {

          //  Letter is NOT ACGI, write out the last block if it was
          //  ACGT, then increment our length.
          if (b->_isACGT == 1) {
            if (b->_len > 0)
              BLOK[BLOKlen++] = *b;

            b->_isACGT = 0;
            b->_iid    = sic->getIID();
            b->_pos    = p;
            b->_len    = 0;
            b->_bpf    = DATA->tell() / 2;
          }

          b->_len++;

        } else {

          //  Letter is ACGT.  Write out last block if it was a gap,
          //  then emit this letter.
          if (b->_isACGT == 0) {
            if (b->_len > 0)
              BLOK[BLOKlen++] = *b;

            b->_isACGT = 1;
            b->_iid    = sic->getIID();
            b->_pos    = p;
            b->_len    = 0;
            b->_bpf    = DATA->tell() / 2;
          }

          b->_len++;

          DATA->putBits(bits, 2);
        }
      }

      //  Emit the last block
      BLOK[BLOKlen++] = *b;
    }

    delete seq;
  }

  //  And a sentinel EOF block -- gets the last position in the file,
  //  useful for the binary search.

  BLOK[BLOKlen]._isACGT = 0;
  BLOK[BLOKlen]._iid    = ~u32bitZERO;
  BLOK[BLOKlen]._pos    = ~u32bitZERO;
  BLOK[BLOKlen]._len    = 0;
  BLOK[BLOKlen]._bpf    = DATA->tell();

  BLOKlen++;

  //  Update the header, assembly the final file.

  delete DATA;

  HEAD._magic[0]          = SEQSTORE_MAGICNUMBER1;
  HEAD._magic[1]          = SEQSTORE_MAGICNUMBER2;
  HEAD._numberOfSequences = inputseq->getNumberOfSequences();
  HEAD._numberOfBlocks    = BLOKlen;
  HEAD._namesLength       = NAMElen;

  errno = 0;
  FILE *F = fopen(filename, "r+");
  if (errno)
    fprintf(stderr, "seqStore::construct()--  Failed to reopen '%s' to write data: %s\n",
            _filename, strerror(errno)), exit(1);

  HEAD._indexStart = fseeko(F, 0, SEEK_END);
  fwrite(INDX, sizeof(seqStoreIndex), HEAD._numberOfSequences, F);

  HEAD._blockStart = fseeko(F, 0, SEEK_END);
  fwrite(BLOK, sizeof(seqStoreBlock), HEAD._numberOfBlocks, F);

  HEAD._namesStart = fseeko(F, 0, SEEK_END);
  fwrite(NAME, sizeof(char), HEAD._namesLength, F);

  fseeko(F, 0, SEEK_SET);
  fwrite(&HEAD, sizeof(seqStoreIndex), 1, F);

  fclose(F);

  if (errno)
    fprintf(stderr, "seqStore::construct()--  Failed to write data to '%s': %s\n",
            _filename, strerror(errno)), exit(1);

  delete [] INDX;
  delete [] BLOK;
  delete [] NAME;
}
