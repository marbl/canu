
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



//  Test that the store can be accessed.  This is not testing the implementation
//  of ovStore, just that the data on disk can be accessed successfully.
void
ovStore::testStore(bool verbose) {
  uint32    lid   = _seq->sqStore_lastReadID();

  //  Step 1: load the number of overlaps per read.

  uint32   *nopr  = numOverlapsPerRead();
  uint32    nrds  = 0;
  uint64    novl  = 0;

  for (uint32 xx=0; xx <= lid; xx++)
    if (nopr[xx] > 0) {
      nrds += 1;
      novl += nopr[xx];
    }

  //  Step 2: load overlaps for some reads.

  uint32      t[5];
  uint32      ovlLen = 0;
  uint32      ovlMax = 0;
  ovOverlap  *ovl    = NULL;

  for (uint32 ii=0; ii<5; ii++) {
    uint32 tid = ii * lid / 4;

    while ((nopr[tid] == 0) && (tid < lid))   //  Search forward for the next
      tid++;                                  //  read with overlaps.

    while ((nopr[tid] == 0) && (tid > 1))     //  Search backward.
      tid--;

    if (nopr[tid] > 0) {
      ovlLen = loadOverlapsForRead(tid, ovl, ovlMax);

      if (ovlLen != nopr[tid]) {
        fprintf(stderr, "ERRROR:  Failed loading overlaps for read %u: expected %u overlaps, loaded %u.\n",
                tid, nopr[tid], ovlLen);
        exit(1);
      }

      if (verbose)
        fprintf(stderr, "Loaded %6u overlaps for read %6u\n", ovlLen, tid);
    }
  }

  //  Passed!  Cleanup and celebrate.

  delete [] nopr;
  delete [] ovl;

  if (verbose)
    fprintf(stderr, "Success!  Found %u reads with %lu overlaps out of %u reads total.\n",
            nrds, novl, lid);
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

    if (_seq)                         //  Is this needed anymore?  (29 Jan 2020)
      overlap->sqStoreAttach(_seq);   //  (there are three in this file)

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
ovStore::loadBlockOfOverlaps(ovOverlap *&ovl,
                             uint32     &ovlMax) {
  uint32  ovlLen = 0;

  //  If we don't have space for the overlaps from the next read,
  //  reallocate space for just those.

  if (ovlMax < _index[_curID]._numOlaps) {
    delete [] ovl;

    ovlMax = _index[_curID]._numOlaps;
    ovl    = new ovOverlap [ovlMax];
  }

  //  Now load overlaps for reads until we run out of space.

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

      if (_seq)                            //  Is this needed anymore?  (29 Jan 2020)
        ovl[ovlLen].sqStoreAttach(_seq);   //  (there are three in this file)

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

    if (_seq)                        //  Is this needed anymore?  (29 Jan 2020)
      ovl[oo].sqStoreAttach(_seq);   //  (there are three in this file)

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

