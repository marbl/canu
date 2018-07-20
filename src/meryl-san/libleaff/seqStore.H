
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
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#ifndef SEQSTORE_H
#define SEQSTORE_H

#include "seqCache.H"
#include "bitPackedFile.H"

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
  uint64  _magic[2];
  uint32  _pad;
  uint32  _numberOfSequences;
  uint64  _numberOfACGT;
  uint32  _numberOfBlocksACGT;
  uint32  _numberOfBlocksGAP;
  uint32  _numberOfBlocks;
  uint32  _namesLength;

  uint64  _indexStart;
  uint64  _blockStart;
  uint64  _namesStart;
};


//  This index allows us to return a complete sequence
//
struct seqStoreIndex {
  uint32  _hdrPosition;  //  Offset into _names for the defline
  uint32  _hdrLength;    //  Length of the defline
  uint64  _seqPosition;  //  Offset into _bpf for the sequence data
  uint32  _seqLength;    //  Length, in bases, of the sequence
  uint32  _block;        //  The seqStoreBlock that starts this sequence
};


//  This index allows us to seek to a specific base in the
//  file of sequences.  Each block is either:
//    ACGT - and has data
//    N    - no data
//  It will map a specific ACGT location to the sequence, and the ID
//  of that sequence (seq ID and location in that sequence).
//
struct seqStoreBlock {
  uint64      _isACGT:1;    // block is acgt
  uint64      _pos:32;      // position in sequence
  uint64      _iid:32;      // iid of the sequence we are in
  uint64      _len:23;      // length of block
  uint64      _bpf:40;      // position in the bit file of sequence
};

#define SEQSTOREBLOCK_MAXPOS uint64MASK(32)
#define SEQSTOREBLOCK_MAXIID uint64MASK(32)
#define SEQSTOREBLOCK_MAXLEN uint64MASK(23)

class seqStore : public seqFile {
protected:
  seqStore(const char *filename);
  seqStore();

public:
  ~seqStore();

protected:
  seqFile            *openFile(const char *filename);

public:
  uint32              find(const char *sequencename);

  uint32              getSequenceLength(uint32 iid);
  bool                getSequence(uint32 iid,
                                  char *&h, uint32 &hLen, uint32 &hMax,
                                  char *&s, uint32 &sLen, uint32 &sMax);
  bool                getSequence(uint32 iid,
                                  uint32 bgn, uint32 end, char *s);

private:
  void                clear(void);
  void                loadIndex(void);

  bitPackedFile     *_bpf;

  seqStoreHeader     _header;

  seqStoreIndex     *_index;
  seqStoreBlock     *_block;
  char              *_names;

  bitPackedFile     *_indexBPF;
  bitPackedFile     *_blockBPF;
  bitPackedFile     *_namesBPF;

  uint32             _lastIIDloaded;

  friend class seqFactory;
};


//  Construct a new seqStore 'filename' from input file 'inputseq'.
//
void
constructSeqStore(char *filename,
                  seqCache *inputseq);


#endif  //  SEQSTORE_H
