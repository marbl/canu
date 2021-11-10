
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

#include "runtime.H"

#include "sqStore.H"
#include "ovStore.H"
#include "tgStore.H"

#include <map>
#include <set>
#include <string>

//
//  Stores read sequence, compressed, in memory.
//


//  An entry in the cache.
//
//  _basesLength is the length of the sequence stored in the blob.  It
//  is NOT the length of the read we will eventually return.
//
//  _bgn and _end tell what bases we will be returning, _end - _bgn is the
//  length of the read we are storing.
//    If not trimmed, these must be 0 and the actual length of the read.
//    If homopoly compressed, this will be the compressed length of the read.
//
struct sqCacheEntry {
  uint32  _nLen           = 0;           //  Length of unencoded data, not
  uint32  _sLen           = 0;           //  including any terminating nul.

  uint32  _bgn            = 0;           //  Portion of the read we'll return.
  uint32  _end            = 0;

  uint8  *_nData          = nullptr;
  uint8  *_sData          = nullptr;
};



//  Loads part/all of a sqStore, or all of a fasta/fastq file, into memory.
//
class sqCache {
private:
  void         loadMetadata(void);

public:
  sqCache(sqStore       *seqStore     = nullptr,
          sqRead_which   which        = sqRead_defaultVersion);
  ~sqCache();

private:
  void         loadRead(uint32 id);
  void         loadRead(dnaSeq &seq);

private:

public:
  //  Read accessors.
  uint32       sqCache_getLength(uint32      id) {
    return(_reads[id]._end - _reads[id]._bgn);
  };

  char        *sqCache_getSequence(uint32    id);

  char        *sqCache_getSequence(uint32    id,
                                   char    *&seq,
                                   uint32   &seqLen,
                                   uint32   &seqMax);

  uint32       sqCache_mapNameToID(char const *readName);

public:
  //  Data loaders that load from a seqStore supplied in the constructor.

  void         sqCache_loadReads(bool verbose=false);
  void         sqCache_loadReads(uint32 bgnID, uint32 endID, bool verbose=false);
  void         sqCache_loadReads(std::set<uint32> reads, bool verbose=false);
  void         sqCache_loadReads(std::map<uint32, uint32> reads, bool verbose=false);
  void         sqCache_loadReads(ovOverlap *ovl, uint32 nOvl, bool verbose=false);
  void         sqCache_loadReads(tgTig *tig, bool verbose=false);

  //  Data loader from a single fasta/fastq file.  The 'id' of the read is
  //  assigned incrementally starting from 1.  ALL reads are loaded by this call.

  void         sqCache_loadReads(char const *filename);

public:
  void         sqCache_saveReadToBuffer(writeBuffer *B, uint32 id, sqRead *rd, sqReadDataWriter *wr);


private:
  sqStore         *_seqStore        = nullptr;

  sqRead_which     _which           = sqRead_unset;
  bool             _compressed      = false;
  bool             _trimmed         = false;

  sqCacheEntry    *_reads           = nullptr;
  uint64           _readsLen        = 0;
  uint64           _readsMax        = 0;

  typedef std::map<std::string, uint32>    nameMap;
  nameMap          _nameToID;

  void            allocateNewBlock(void) {
    increaseArray(_dataBlocks, _dataBlocksLen, _dataBlocksMax, 16);

    assert(_dataBlocksLen < _dataBlocksMax);

    _dataLen = 0;
    _data    = _dataBlocks[_dataBlocksLen++] = new uint8 [_dataMax];
  };

  uint32           _dataBlocksLen = 0;   //  Pointers to allocated blocks.
  uint32           _dataBlocksMax = 0;
  uint8          **_dataBlocks    = nullptr;

  uint64           _dataLen = 0;         //  The active block, current
  uint64           _dataMax = 33554432;  //  and maximum length (32MB).
  uint8           *_data    = nullptr;

  sqRead           _read;                //  Used mostly as a buffer for blob data.
};

