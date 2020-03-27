
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

#include <set>
using namespace std;


//  Stores read sequence, compressed, in memory.
//
//   - load all reads.
//   - load all reads in a list.
//   - load all reads in a list of overlaps.
//   - load all reads in a tig.
//


class sqCacheEntry {
public:
  sqCacheEntry() {
    _basesLength    = 0;
    _bgn            = 0;
    _end            = 0;
    //_dataAge        = 0;
    _dataExpiration = UINT32_MAX;
    _data           = NULL;
  };

  ~sqCacheEntry() {
    //if (_data)
    //  fprintf(stderr, "DELETE READ\n");

    delete [] _data;
  };

  //  _sequenceLength is the length of the sequence stored in the blob.  It
  //  is NOT the length of the read we will eventually return.
  //
  //  _bgn and _end tell what bases we will be returning, _end - _bgn is the
  //  length of the read we are storing.
  //    If not trimmed, these must be 0 and the actual length of the read.
  //    If homopoly compressed, this will be the compressed length of the read.

  uint32  _basesLength;

  uint32  _bgn;
  uint32  _end;

  //  For expiring data from the cache, two possibilities:
  //   - We know ahead of time how many times we're going to request
  //     each read, and can remove the read from the cache when
  //     _dataAge == _dataExpiration
  //
  //   - We want to keep only the most recently used reads in the
  //     cache; if we run out of memory, throw out the least recently
  //     used reads, those with the largest _dataAge.

  //uint32  _dataAge;
  uint32  _dataExpiration;

  uint8  *_data;
};



class sqCache {
public:
  sqCache(sqStore *seqStore, sqRead_which which=sqRead_defaultVersion, uint64 memoryLimit=0);
  ~sqCache();

private:
  void         loadRead(uint32 id, uint32 expiration=1);
  void         removeRead(uint32 id);
  void         increaseAge(void);

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

public:
  //  Data loaders.
  void         sqCache_loadReads(bool verbose=false);
  void         sqCache_loadReads(uint32 bgnID, uint32 endID, bool verbose=false);
  void         sqCache_loadReads(set<uint32> reads, bool verbose=false);
  void         sqCache_loadReads(map<uint32, uint32> reads, bool verbose=false);
  void         sqCache_loadReads(ovOverlap *ovl, uint32 nOvl, bool verbose=false);
  void         sqCache_loadReads(tgTig *tig, bool verbose=false);

  void         sqCache_purgeReads(void);


private:
  sqStore         *_seqStore;
  uint32           _nReads;

  bool             _trackAge;
  bool             _trackExpiration;
  bool             _noMoreLoads;

  sqRead_which     _which;
  bool             _compressed;
  bool             _trimmed;

  uint64           _memoryLimit;

  sqCacheEntry    *_reads;

  void            allocateNewBlock(void) {
    increaseArray(_dataBlocks, _dataBlocksLen, _dataBlocksMax, 16);

    assert(_dataBlocksLen < _dataBlocksMax);

    _dataBlocks[_dataBlocksLen++] = new uint8 [_dataMax];

    _dataLen = 0;
    _data    = _dataBlocks[_dataBlocksLen - 1];
  };

  uint32           _dataBlocksLen;   //  Pointers to allocated blocks.
  uint32           _dataBlocksMax;
  uint8          **_dataBlocks;

  uint64           _dataLen;         //  The active block, current
  uint64           _dataMax;         //  and maximum length.
  uint8           *_data;

  sqRead           _read;            //  Used mostly as a buffer for blob data.
};

