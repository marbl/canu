
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
 *  This file is derived from:
 *
 *    src/overlapInCore/overlapPair-readCache.H
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2015-JUN-16 to 2015-JUN-23
 *      are Copyright 2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2017-SEP-26
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_global.H"

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
    _readLength     = 0;
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

  uint32  _readLength;   //  Length of the untrimmed read.

  uint32  _bgn;          //  Trim points.  If not trimmed, set
  uint32  _end;          //  to 0 and _readLength.

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
  sqCache(sqStore *seqStore, sqRead_version version=sqRead_latest, uint64 memoryLimit=0);
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
  void         sqCache_loadReads(void);
  void         sqCache_loadReads(set<uint32> reads);
  void         sqCache_loadReads(map<uint32, uint32> reads);
  void         sqCache_loadReads(ovOverlap *ovl, uint32 nOvl);
  void         sqCache_loadReads(tgTig *tig);

  void         sqCache_purgeReads(void);


private:
  sqStore         *_seqStore;
  uint32           _nReads;

  bool             _trackAge;
  bool             _trackExpiration;

  sqRead_version   _version;

  uint64           _memoryLimit;

  sqCacheEntry    *_reads;

  void            allocateNewBlock(void) {
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

  sqReadData       _readData;
};


