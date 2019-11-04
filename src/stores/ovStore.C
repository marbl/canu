
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
 *    src/AS_OVS/AS_OVS_overlapStore.C
 *    src/AS_OVS/AS_OVS_overlapStore.c
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2007-MAR-08 to 2013-AUG-01
 *      are Copyright 2007-2013 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Sergey Koren on 2007-MAY-08
 *      are Copyright 2007 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Sergey Koren from 2011-JUN-02 to 2011-JUN-03
 *      are Copyright 2011 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Gregory Sims from 2012-FEB-01 to 2012-FEB-14
 *      are Copyright 2012 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2014-DEC-09 to 2015-AUG-14
 *      are Copyright 2014-2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2015-OCT-12
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *    Sergey Koren beginning on 2015-DEC-15
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *    Sergey Nurk beginning on 2019-SEP-18
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "ovStore.H"



ovStore::ovStore(const char *path, sqStore *seq) {
  char  name[FILENAME_MAX];

  //  Save the path name.

  if (path == NULL)
    fprintf(stderr, "ovStore::ovStore()-- ERROR: no name supplied.\n"), exit(1);

  if ((path[0] == '-') &&
      (path[1] == 0))
    fprintf(stderr, "ovStore::ovStore()-- ERROR: name cannot be '-' (stdin).\n"), exit(1);

  memset(_storePath, 0, FILENAME_MAX+1);
  strncpy(_storePath, path, FILENAME_MAX);

  //  Load the info file, or die trying.

  _info.load(_storePath);

  //  Initialize.

  _seq              = seq;

  _curID            = 1;
  _bgnID            = 1;
  _endID            = _info.maxID();

  _curOlap          = 0;

  _index            = NULL;

  _evaluesMap       = NULL;
  _evalues          = NULL;

  _bof              = NULL;
  _bofSlice         = 0;
  _bofPiece         = 0;

  //  Open the index

  _index = new ovStoreOfft [_info.maxID()+1];

  AS_UTL_loadFile(_storePath, '/', "index", _index, _info.maxID()+1);

  //  Open and load erates

  snprintf(name, FILENAME_MAX, "%s/evalues", _storePath);

  if (fileExists(name)) {
    _evaluesMap  = new memoryMappedFile(name, memoryMappedFile_readOnly);
    _evalues     = (uint16 *)_evaluesMap->get(0);
  }
}




ovStore::~ovStore() {
  delete [] _index;
  delete    _evaluesMap;
  delete    _bof;
}



uint32
ovStore::readOverlap(ovOverlap *overlap) {

  //  If we've finished reading overlaps for the current read, find the next read.

  if (_curOlap == _index[_curID]._numOlaps) {
    _curOlap  = 0;
    _curID   += 1;

    while ((_curID <= _endID) &&
           (_index[_curID]._numOlaps == 0))
      _curID++;

    if (_curID > _endID)   //  Out of reads to return overlaps for.
      return(0);

    assert(_index[_curID]._slice > 0);
    assert(_index[_curID]._piece > 0);

    if ((_bofSlice != _index[_curID]._slice) ||     //  Make sure we're in the correct file.
        (_bofPiece != _index[_curID]._piece)) {
      delete _bof;

      assert(_index[_curID]._slice > 0);
      assert(_index[_curID]._piece > 0);

      _bofSlice = _index[_curID]._slice;
      _bofPiece = _index[_curID]._piece;

      _bof = new ovFile(_seq, _storePath, _bofSlice, _bofPiece, ovFileNormal);
      _bof->seekOverlap(_index[_curID]._offset);
    }
  }

  //  If we can read the next overlap, return it.

  if (_bof->readOverlap(overlap) == true) {
    overlap->a_iid = _curID;
    overlap->g     = _seq;

    if (_evalues)
      overlap->evalue(_evalues[_index[_curID]._overlapID++]);

    _curOlap++;

    return(1);
  }

  //  Otherwise....what?

  assert(0);
  return(0);
}



//  Bulk loads overlaps into ovl.
//  Doesn't split overlaps for a single read across multiple blocks.
//
uint32
ovStore::loadBlockOfOverlaps(ovOverlap *ovl,
                             uint32     ovlMax) {
  uint32  ovlLen = 0;

  while ((ovlLen + _index[_curID]._numOlaps < ovlMax) &&
         (_curID <= _endID)) {

    //  Open a new file if the file changed (but only if this read actually HAS overlaps, otherwise,
    //  the slice/piece it claims to be in is invalid).

    if ((_index[_curID]._numOlaps > 0) &&
        ((_bofSlice != _index[_curID]._slice) ||
         (_bofPiece != _index[_curID]._piece))) {
      delete _bof;

      assert(_index[_curID]._slice > 0);
      assert(_index[_curID]._piece > 0);

      _bofSlice = _index[_curID]._slice;
      _bofPiece = _index[_curID]._piece;

      _bof = new ovFile(_seq, _storePath, _bofSlice, _bofPiece, ovFileNormal);
      _bof->seekOverlap(_index[_curID]._offset);
    }

    //  Load all overlaps for this read.  No need to check anything; we're guaranteed
    //  all these overlaps exist in this file.

    for (uint32 oo=0; oo<_index[_curID]._numOlaps; oo++) {
      if (_bof->readOverlap(ovl + ovlLen) == false) {
        fprintf(stderr, "ovStore::loadBLockOfOverlaps()-- Failed to load overlap %u out of %u for read %u.\n", oo, _index[_curID]._numOlaps, _curID);
        exit(1);
      }

      ovl[ovlLen].a_iid = _curID;
      ovl[ovlLen].g     = _seq;

      if (_evalues)
        ovl[ovlLen].evalue(_evalues[_index[_curID]._overlapID++]);

      ovlLen++;
    }

    _curID   += 1;     //  Advance to the next read.
    _curOlap  = 0;     //  We've read no overlaps for this read.
  }

  return(ovlLen);
}



uint32
ovStore::loadOverlapsForRead(uint32       id,
                             ovOverlap  *&ovl,
                             uint32      &ovlMax) {

  _curID   = id;
  _curOlap = 0;

  //  Not a requested overlap?  Do nothing.

  if (_curID < _bgnID)
    return(0);
  if (_endID < _curID)
    return(0);

  //  Nothing there?  Do nothing.

  if (_index[_curID]._numOlaps == 0) {
    _curID++;
    return(0);
  }

  //  Make more space if needed.

  if (ovlMax < _index[_curID]._numOlaps) {
    if (ovlMax > 0)
      delete [] ovl;

    ovlMax = _index[_curID]._numOlaps * 1.2;
    ovl    = new ovOverlap [ovlMax];
  }

  //  If we're not in the correct file, open the correct file.

  if ((_index[_curID]._numOlaps > 0) &&
      ((_bofSlice != _index[_curID]._slice) ||
       (_bofPiece != _index[_curID]._piece))) {

    assert(_index[_curID]._slice > 0);
    assert(_index[_curID]._piece > 0);

    _bofSlice = _index[_curID]._slice;
    _bofPiece = _index[_curID]._piece;

    delete _bof;

    _bof = new ovFile(_seq, _storePath, _index[_curID]._slice, _index[_curID]._piece, ovFileNormal);
  }

  //  Always reposition (unless there are no overlaps).
  //  I assume this will do nothing if not needed.

  if (_index[_curID]._numOlaps > 0)
    _bof->seekOverlap(_index[_curID]._offset);

  //  Load the overlaps.  By the construction of the store, we're guaranteed
  //  all overlaps will be in this ovFile, so can just load load load.

  for (uint32 oo=0; oo<_index[_curID]._numOlaps; oo++) {
    if (_bof->readOverlap(ovl + oo) == false) {
      fprintf(stderr, "ovStore::loadOverlapsForRead()-- Failed to load overlap %u out of %u for read %u.\n", oo, _index[_curID]._numOlaps, _curID);
      exit(1);
    }

    ovl[oo].a_iid = _curID;
    ovl[oo].g     = _seq;

    if (_evalues)
      ovl[oo].evalue(_evalues[_index[_curID]._overlapID++]);
  }

  _curID   += 1;     //  Advance to the next read.
  _curOlap  = 0;     //  We've read no overlaps for this read.

  //  Done!

  return(_index[id]._numOlaps);
}




void
ovStore::setRange(uint32 bgnID, uint32 endID) {

  //  Remove the old file.

  delete _bof;
  _bof = NULL;

  //  Set ranges, limiting them to the last read (possibly last read with overlaps).

  _bgnID = min(bgnID, _info.maxID());
  _curID = min(bgnID, _info.maxID());
  _endID = min(endID, _info.maxID());

  //  Skip reads with no overlaps.

  while ((_curID <= _endID) &&
         (_index[_curID]._numOlaps == 0))
    _curID++;

  //  If no overlaps, the range is already exhausted and we can just return.

  if (_curID > _endID)
    return;

  //  If no slice or piece, that's kind of bad and we blow ourself up.

  if ((_index[_curID]._slice == 0) ||
      (_index[_curID]._piece == 0))
    fprintf(stderr, "ovStore::setRange(%u, %u)-- invalid slice/piece %u/%u for curID %u\n",
            _bgnID, _endID, _index[_curID]._slice, _index[_curID]._piece, _curID);

  assert(_index[_curID]._slice != 0);
  assert(_index[_curID]._piece != 0);

  //  Open new file, and position at the correct spot.

  _bof = new ovFile(_seq, _storePath, _index[_curID]._slice, _index[_curID]._piece, ovFileNormal);
  _bof->seekOverlap(_index[_curID]._offset);
}



void
ovStore::restartIteration(void) {
  _curID   = _bgnID;
  _curOlap = 0;
}



void
ovStore::endIteration(void) {
  _curID   = _endID;
  _curOlap = _index[_curID]._numOlaps;
}



uint64
ovStore::numOverlapsInRange(void) {
  uint64    numOlaps = 0;

  for (uint32 ii=_bgnID; ii<=_endID; ii++)
    numOlaps += _index[ii]._numOlaps;

  return(numOlaps);
}



//  Return an array with the number of overlaps per read.
//  If numReads is more than zero, only those reads will be loaded.
//
uint32 *
ovStore::numOverlapsPerRead(void) {
  uint32  *olapsPerRead = new uint32 [_info.maxID() + 1];

  for (uint32 ii=0; ii <= _info.maxID(); ii++)
    olapsPerRead[ii] = _index[ii]._numOlaps;

  return(olapsPerRead);
}



void
ovStore::dumpMetaData(uint32 bgnID, uint32 endID) {

  fprintf(stdout, "   readID slice piece    offset  numOlaps overlapID\n");
  fprintf(stdout, "--------- ----- ----- --------- --------- ---------\n");

  for (uint32 ii=bgnID; ii<=endID; ii++)
    fprintf(stdout, "%9u %5u %5u %9u %9u %9lu\n", ii,
            _index[ii]._slice,
            _index[ii]._piece,
            _index[ii]._offset,
            _index[ii]._numOlaps,
            _index[ii]._overlapID);

  fprintf(stdout, "--------- ----- ----- --------- --------- ---------\n");
}

