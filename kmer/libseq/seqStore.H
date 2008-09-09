#ifndef SEQSTORE_H
#define SEQSTORE_H

#include "util++.H"
#include "seqCache.H"

//  A binary fasta file.
//
//  HEADER
//    magic number
//    number of sequences
//    optional - alphabet size
//    optional - alphabet map (0x00 -> 'a', etc)
//    position of index start
//    position of data start
//  DATA
//  INDEX
//    position of sequence start in DATA
//    header length
//    sequence length
//  MAP
//    name to IID mapping

struct seqStoreHeader {
  u64bit  _magic[2];
  u32bit  _numberOfSequences;
  u32bit  _numberOfBlocks;
  u32bit  _namesLength;

  u64bit  _indexStart;
  u64bit  _blockStart;
  u64bit  _namesStart;
};


//  This index allows us to return a complete sequence
//
struct seqStoreIndex {
  u32bit  _hdrPosition;  //  Offset into _names for the defline
  u32bit  _hdrLength;    //  Length of the defline
  u64bit  _seqPosition;  //  Offset into _rb for the sequence data
  u32bit  _seqLength;    //  Length, in bases, of the sequence
};


//  This index allows us to seek to a specific base in the
//  file of sequences.  Each block is either:
//    ACGT - and has data
//    N    - no data
//  It will map a specific ACGT location to the sequence, and the ID
//  of that sequence (seq ID and location in that sequence).
//
struct seqStoreBlock {
  u64bit      _isACGT:1;    // block is acgt
  u64bit      _pos:31;      // position in sequence
  u64bit      _iid:32;      // iid of the sequence we are in
  u64bit      _len:24;      // length of block
  u64bit      _bpf:40;      // position in the bit file of sequence
};


class seqStore : public seqFile {
protected:
  seqStore(const char *filename);
public:
  seqStore();
  ~seqStore();

  seqFile            *openFile(const char *filename);

  u32bit              find(const char *sequencename);

  u32bit              getSequenceLength(u32bit iid);
  bool                getSequence(u32bit id,
                                  char *&h, u32bit &hLen, u32bit &hMax,
                                  char *&s, u32bit &sLen, u32bit &sMax);
  bool                getSequence(u32bit iid,
                                  u32bit bgn, u32bit end, char *s);

private:
  void                clear(void);
  void                construct(char *filename, seqCache *inputseq);

  readBuffer        *_rb;

  seqStoreHeader     _header;
  seqStoreIndex     *_index;
  seqStoreBlock     *_block;
  char              *_names;
};



#endif  //  SEQSTORE_H
