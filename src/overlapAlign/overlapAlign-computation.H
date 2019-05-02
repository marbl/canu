
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
 *    Brian P. Walenz beginning on 2019-APR-22
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */


//  This object is created by the loader thread, overlapReader().
//  During construction of this object, overlaps are loaded,
//  and a pointer to the cache is saved so we can retrieve the
//  reads on demand during the compute.

#include <tuple>

using namespace std;


class maComputation {
public:
  maComputation(uint32 id, trReadData *readData, sqCache *seqCache, ovStore *ovlStore) {

    _readData    = readData;

    _seqCache    = seqCache;

    _overlaps    = NULL;
    _overlapsLen = 0;
    _overlapsMax = 0;

    _aID   = id;
    _bID   = UINT32_MAX;

    _aRead = _bRead = NULL;
    _aMax  = _bMax  = 0;

    //  Load overlaps.

    _overlapsLen = ovlStore->loadOverlapsForRead(id, _overlaps, _overlapsMax);

#if 0
    if (_overlapsLen > 0)
      fprintf(stderr, "READ %8u found %8u overlaps; first IDs %8u %8u\n",
              _aID, _overlapsLen, _overlaps[0].a_iid, _overlaps[0].b_iid);
    else
      fprintf(stderr, "READ %8u found   NONE   overlaps; first IDs -------- --------\n",
              _aID);
#endif

    //  Allocate space for the a read.

    _aMax  = _seqCache->sqCache_getLength(id) + 1;
    _aRead = new char [_aMax];

    //  Find the longest b read and allocate space for b read sequences.

    for (uint32 ii=0; ii<_overlapsLen; ii++)
      _bMax = max(_bMax, _seqCache->sqCache_getLength(_overlaps[ii].b_iid) + 1);

    _bRead = new char [_bMax];
  };

  ~maComputation() {
    delete [] _overlaps;
    delete [] _aRead;
    delete [] _bRead;
  };


public:
  void                 trimRead(uint32  minOverlapLength,
                                double  maxErate);

private:
  bool                 isWellContained(int32 maxEdge);
  tuple<int32, int32>  findMinThickestEdge(double coverage);

  void                 fetchUntrimmedRead(uint32 id_, bool isA_=true, bool revComp_=false);
  void                 fetchTrimmedRead  (uint32 id_, bool isA_=true, bool revComp_=false);
  bool                 trimOverlap_Normal(ovOverlap *ovl);
  bool                 trimOverlap_Flipped(ovOverlap *ovl);
  bool                 trimOverlap(ovOverlap *ovl);

public:
  void                 computeAlignments(uint32 minOverlapLength,
                                         double maxErate);

private:
  sqCache    *_seqCache;

public:

  //  A pointer to the output object.

  trReadData *_readData;

  //  Overlaps are loaded by the loader thread during creation
  //  of this object.

  ovOverlap  *_overlaps;
  uint32      _overlapsLen;
  uint32      _overlapsMax;

  //  Sequence storage during alignment.  We assume the read
  //  sequences are cached in core, and that we can get them
  //  on demand, for free, during the computation.

  uint32      _aID,      _bID;
  uint32      _aMax,     _bMax;    //  Allocated length of aRead.
  char       *_aRead,   *_bRead;   //
};


