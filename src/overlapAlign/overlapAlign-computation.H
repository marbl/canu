
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

class maComputation {
public:
  maComputation(uint32 id, sqCache *seqCache, ovStore *ovlStore) {

    _seqCache    = seqCache;

    _overlaps    = NULL;
    _overlapsLen = 0;
    _overlapsMax = 0;

    _aID   = id;
    _aRead = _bRead = NULL;
    _aLen  = _bLen  = 0;
    _aMax  = _bMax  = 0;

    _bgn = 0;
    _end = 0;

    //  Load overlaps.

    _overlapsLen = ovlStore->loadOverlapsForRead(id, _overlaps, _overlapsMax);

    //  Allocate space for the a read.

    _aLen  = 0;
    _aMax  = _seqCache->sqCache_getLength(id) + 1;
    _aRead = new char [_aMax];

    //  Find the longest b read and allocate space for b read sequences.

    _bLen  = 0;
    _bMax  = 0;

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
  void         trimRead(uint32  minOverlapLength, double  maxErate, int32 &aClrBgn, int32 &aClrEnd);
  void         computeAlignments(uint32 minOverlapLength, double maxErate);

private:
  sqCache    *_seqCache;

public:

  //  Overlaps are loaded by the loader thread during creation
  //  of this object.

  ovOverlap  *_overlaps;
  uint32      _overlapsLen;
  uint32      _overlapsMax;

  //  Sequence storage during alignment.  We assume the read
  //  sequences are cached in core, and that we can get them
  //  on demand, for free, during the computation.

  uint32      _aID;
  char       *_aRead,   *_bRead;
  uint32      _aLen,     _bLen;
  uint32      _aMax,     _bMax;

  //  Errors, evidence, alignments.

  

  //  Final trim points.

  uint32      _bgn;
  uint32      _end;
};


