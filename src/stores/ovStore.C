
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
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "ovStore.H"

const uint64 ovStoreVersion         = 2;
const uint64 ovStoreMagic           = 0x53564f3a756e6163;   //  == "canu:OVS - store complete
const uint64 ovStoreMagicIncomplete = 0x50564f3a756e6163;   //  == "canu:OVP - store under construction


void
ovStore::ovStore_write(void) {
  AS_UTL_mkdir(_storePath);

  char name[FILENAME_MAX];

  sprintf(name, "%s/info", _storePath);

  //  If the ovs file exists, AND has a valid magic number, then the store is complete and we should
  //  abort before the valid store is destroyed.

  if (AS_UTL_fileExists(name, false, false)) {
    errno = 0;
    FILE *ovsinfo = fopen(name, "r");
    if (errno) {
      fprintf(stderr, "ERROR: failed to read store metadata from '%s': %s\n", name, strerror(errno));
      exit(1);
    }

    AS_UTL_safeRead(ovsinfo, &_info, "ovStore::ovStore::testinfo", sizeof(ovStoreInfo), 1);

    fclose(ovsinfo);

    if (_info._ovsMagic == ovStoreMagic)
      fprintf(stderr, "ERROR:  overlapStore '%s' is a valid overlap store, will not overwrite.\n",
              _storePath), exit(1);
  }

  //  Create a new incomplete info file.

  errno = 0;
  FILE *ovsinfo = fopen(name, "w");

  if (errno)
    fprintf(stderr, "failed to create overlap store '%s': %s\n", _storePath, strerror(errno)), exit(1);

  AS_UTL_safeWrite(ovsinfo, &_info, "ovStore::ovStore::saveinfo", sizeof(ovStoreInfo), 1);

  fclose(ovsinfo);

  sprintf(name, "%s/index", _storePath);

  errno = 0;
  _offtFile = fopen(name, "w");
  if (errno)
    fprintf(stderr, "AS_OVS_createOverlapStore()-- failed to open offset file '%s': %s\n", name, strerror(errno)), exit(1);

  _overlapsThisFile = 0;
  _currentFileIndex = 0;
  _bof              = NULL;
}



void
ovStore::ovStore_read(void) {
  char  name[FILENAME_MAX];

  sprintf(name, "%s/info", _storePath);
  errno = 0;
  FILE *ovsinfo = fopen(name, "r");
  if (errno)
    fprintf(stderr, "ERROR: directory '%s' is not an ovelrapStore; failed to open info file '%s': %s\n",
            _storePath, name, strerror(errno)), exit(1);

  AS_UTL_safeRead(ovsinfo, &_info, "ovStore::ovStore::info", sizeof(ovStoreInfo), 1);

  fclose(ovsinfo);

  if ((_info._ovsMagic != ovStoreMagic) && (_info._ovsMagic != ovStoreMagicIncomplete))
    fprintf(stderr, "ERROR:  directory '%s' is not an overlapStore; magic number 0x%016"F_X64P" incorrect.\n",
            _storePath, _info._ovsMagic), exit(1);

  if ((_info._ovsMagic != ovStoreMagic) && (_info._ovsMagic != ovStoreMagicIncomplete))
    fprintf(stderr, "ERROR:  overlapStore '%s' is incomplate; creation crashed?\n",
            _storePath), exit(1);

  if (_info._ovsVersion != ovStoreVersion)
    fprintf(stderr, "ERROR:  overlapStore '%s' is version "F_U64"; this code supports only version "F_U64".\n",
            _storePath, _info._ovsVersion, ovStoreVersion), exit(1);

  if (_info._maxReadLenInBits != AS_MAX_READLEN_BITS)
    fprintf(stderr, "ERROR:  overlapStore '%s' is for AS_MAX_READLEN_BITS="F_U64"; this code supports only %d bits.\n",
            _storePath, _info._maxReadLenInBits, AS_MAX_READLEN_BITS), exit(1);

  //  Load stats

#if 0
  sprintf(name, "%s/statistics", _storePath);
  errno = 0;
  FILE *ost = fopen(name, "r");
  if (errno)
    fprintf(stderr, "failed to open the stats file '%s': %s\n", name, strerror(errno)), exit(1);
  AS_UTL_safeRead(ost, &_stats, "ovStore::ovStore::stats", sizeof(OverlapStoreStats), 1);
  fclose(ost);
#endif

  //  Open the index

  sprintf(name, "%s/index", _storePath);

  errno = 0;
  _offtFile = fopen(name, "r");
  if (errno)
    fprintf(stderr, "ERROR:  failed to open offset file '%s': %s\n", name, strerror(errno)), exit(1);

  //  Open erates

  sprintf(name, "%s/evalues", _storePath);

  if (AS_UTL_fileExists(name)) {
    _evaluesMap  = new memoryMappedFile(name, memoryMappedFile_readOnly);
    _evalues     = (uint16 *)_evaluesMap->get(0);
  }
}




ovStore::ovStore(const char *path, gkStore *gkp, ovStoreType cType) {

  if (path == NULL)
    fprintf(stderr, "ovStore::ovStore()-- ERROR: no name supplied.\n"), exit(1);

  if ((path[0] == '-') &&
      (path[1] == 0))
    fprintf(stderr, "ovStore::ovStore()-- ERROR: name cannot be '-' (stdin).\n"), exit(1);

  memset(_storePath, 0, FILENAME_MAX);
  strncpy(_storePath, path, FILENAME_MAX-1);

  _isOutput  = (cType & ovStoreWrite)   ? true : false;

  _info._ovsMagic         = ovStoreMagicIncomplete;  //  Appropriate for a new store.
  _info._ovsVersion       = ovStoreVersion;
  _info._smallestIID      = UINT64_MAX;
  _info._largestIID       = 0;
  _info._numOverlapsTotal = 0;
  _info._highestFileIndex = 0;
  _info._maxReadLenInBits = AS_MAX_READLEN_BITS;

  _offtFile        = NULL;
  _offt.clear();
  _offm.clear();

  _evaluesMap         = NULL;
  _evalues            = NULL;

  _overlapsThisFile  = 0;
  _currentFileIndex  = 0;
  _bof               = NULL;

  //  Now open an existing store, or a create a new store.

  if (_isOutput == false)
    ovStore_read();
  else
    ovStore_write();

  //  AFTER the info is loaded, set the ranges.

  _firstIIDrequested      = _info._smallestIID;
  _lastIIDrequested       = _info._largestIID;

  _gkp = gkp;
}




ovStore::~ovStore() {

  //  If output, write the last index element (don't forget to fill in gaps);
  //             update the info, using the final magic number

  if (_isOutput) {
    if (_offt._numOlaps > 0) {
      for (; _offm._a_iid < _offt._a_iid; _offm._a_iid++) {
        _offm._fileno   = _offt._fileno;
        _offm._offset   = _offt._offset;
        _offm._numOlaps = 0;

        AS_UTL_safeWrite(_offtFile, &_offm, "ovStore::~ovStore::offm", sizeof(ovStoreOfft), 1);
      }

      AS_UTL_safeWrite(_offtFile, &_offt, "ovStore::~ovStore::offt", sizeof(ovStoreOfft), 1);
    }

    _info._ovsMagic         = ovStoreMagic;
    _info._ovsVersion       = ovStoreVersion;
    _info._highestFileIndex = _currentFileIndex;

    char name[FILENAME_MAX];

    sprintf(name, "%s/info", _storePath);
    errno = 0;
    FILE *ovsinfo = fopen(name, "w");
    if (errno)
      fprintf(stderr, "failed to create overlap store '%s': %s\n", _storePath, strerror(errno)), exit(1);

    AS_UTL_safeWrite(ovsinfo, &_info, "ovStore::~ovStore::info", sizeof(ovStoreInfo), 1);

    fclose(ovsinfo);

    fprintf(stderr, "Closing the new store:\n");
    fprintf(stderr, "  info._ovsMagic           = 0x%016"F_X64P"\n", _info._ovsMagic);
    fprintf(stderr, "  info._ovsVersion         = "F_U64"\n", _info._ovsVersion);
    fprintf(stderr, "  info._smallestIID        = "F_U64"\n", _info._smallestIID);
    fprintf(stderr, "  info._largestIID         = "F_U64"\n", _info._largestIID);
    fprintf(stderr, "  info._numOverlapsTotal   = "F_U64"\n", _info._numOverlapsTotal);
    fprintf(stderr, "  info._highestFileIndex   = "F_U64"\n", _info._highestFileIndex);
    fprintf(stderr, "  info._maxReadLenInBits   = "F_U64"\n", _info._maxReadLenInBits);
  }

  if (_evaluesMap) {
    delete _evaluesMap;

    _evaluesMap = NULL;
    _evalues    = NULL;
  }

#if 0
  if (_statsUpdated) {
    fprintf(stderr, "Writing new stats.\n");

    char name [FILENAME_MAX];

    sprintf(name, "%s/ost", _storePath);
    errno = 0;
    FILE *ost = fopen(name, "w");
    if (errno)
      fprintf(stderr, "failed to write overlap stats '%s': %s\n", name, strerror(errno)), exit(1);

    AS_UTL_safeWrite(ost, &_stats, "AS_OVS_closeOverlapStore", sizeof(OverlapStoreStats), 1);

    fclose(ost);
  }
#endif

  delete _bof;

  fclose(_offtFile);
}



uint32
ovStore::readOverlap(ovOverlap *overlap) {

  assert(_isOutput == FALSE);

  //  If we've finished reading overlaps for the current a_iid, get
  //  another a_iid.  If we hit EOF here, we're all done, no more
  //  overlaps.

  while (_offt._numOlaps == 0)
    if (0 == AS_UTL_safeRead(_offtFile, &_offt, "ovStore::readOverlap::offset",
                             sizeof(ovStoreOfft), 1))
      return(0);

  //  And if we've exited the range of overlaps requested, return.

  if (_offt._a_iid > _lastIIDrequested)
    return(0);

  while ((_bof == NULL) ||
         (_bof->readOverlap(overlap) == FALSE)) {
    char name[FILENAME_MAX];

    //  We read no overlap, open the next file and try again.

    if (_bof)
      delete _bof;

    _currentFileIndex++;

    sprintf(name, "%s/%04d", _storePath, _currentFileIndex);
    _bof = new ovFile(name, ovFileNormal);
  }

  overlap->a_iid = _offt._a_iid;
  overlap->g     = _gkp;

  if (_evalues)
    overlap->evalue(_evalues[_offt._overlapID++]);

  _offt._numOlaps--;


  return(1);
}



uint32
ovStore::numberOfOverlaps(void) {
  ovOverlap  *ovl  = NULL;
  uint32      novl = 0;

  return(readOverlaps(ovl, novl));
}



uint32
ovStore::readOverlaps(ovOverlap *&overlaps, uint32 &maxOverlaps, bool restrictToIID) {
  int    numOvl = 0;

  assert(_isOutput == FALSE);

  //  If we've finished reading overlaps for the current a_iid, get
  //  another a_iid.  If we hit EOF here, we're all done, no more
  //  overlaps.

  while (_offt._numOlaps == 0)
    if (0 == AS_UTL_safeRead(_offtFile, &_offt, "ovStore::readOverlaps::offset", sizeof(ovStoreOfft), 1))
      return(0);

  //  And if we've exited the range of overlaps requested, return.

  if (_offt._a_iid > _lastIIDrequested)
    return(0);

  //  Just a query?  Return the number of overlaps we'd want to read

  if ((overlaps == NULL) || (maxOverlaps == 0))
    return(_offt._numOlaps);

  //  Allocate more space, if needed

  if (maxOverlaps < _offt._numOlaps) {
    delete [] overlaps;

    while (maxOverlaps < _offt._numOlaps)
      maxOverlaps *= 2;

    overlaps = ovOverlap::allocateOverlaps(_gkp, maxOverlaps);
  }

  //  Read all the overlaps for this ID.

  while (((restrictToIID == true)  && (_offt._numOlaps > 0)) ||
         ((restrictToIID == false) && (_offt._numOlaps > 0) && (numOvl < maxOverlaps))) {

    //  Read an overlap.  If this fails, open the next partition and read from there.

    while ((_bof == NULL) ||
           (_bof->readOverlap(overlaps + numOvl) == false)) {
      char name[FILENAME_MAX];

      //  We read no overlap, open the next file and try again.

      delete _bof;
      _bof = NULL;

      _currentFileIndex++;

      if (_currentFileIndex > _info._highestFileIndex)
        //  No more files, stop trying to load an overlap.
        break;

      sprintf(name, "%s/%04d", _storePath, _currentFileIndex);
      _bof = new ovFile(name, ovFileNormal);
    }

    //  If the currentFileIndex is invalid, we ran out of overlaps to load.  Don't save that
    //  empty overlap to the list.

    if (_currentFileIndex <= _info._highestFileIndex) {
      overlaps[numOvl].a_iid = _offt._a_iid;
      overlaps[numOvl].g     = _gkp;

      if (_evalues)
        overlaps[numOvl].evalue(_evalues[_offt._overlapID++]);

      numOvl++;

      assert(_offt._numOlaps > 0);

      _offt._numOlaps--;
    }

    //  If restrictToIID == false, we're loading all overlaps up to the end of the store, or the
    //  request last IID.  If to the end of store, we never read a last 'offset' and so a_iid is
    //  still valid (below lastIIDrequested == infinity) but numOlaps is still zero, and the mail
    //  loop terminates.

    if (restrictToIID == false) {
      while (_offt._numOlaps == 0)
        if (0 == AS_UTL_safeRead(_offtFile, &_offt, "ovStore::readOverlap::offset", sizeof(ovStoreOfft), 1))
          break;
      if (_offt._a_iid > _lastIIDrequested)
        break;
    }
  }  //  while space for more overlaps, load overlaps

  assert(numOvl <= maxOverlaps);

  return(numOvl);
}











uint32
ovStore::readOverlaps(uint32         iid,
                      ovOverlap   *&ovl,
                      uint32        &ovlLen,
                      uint32        &ovlMax) {

  //  Allocate initial space if needed.

  if (ovl == NULL) {
    ovlLen = 0;
    ovlMax = 65 * 1024;
    ovl    = ovOverlap::allocateOverlaps(_gkp, ovlMax);
  }

  if (iid < ovl[0].a_iid)
    //  Overlaps loaded are for a future read.
    return(0);

  if (iid == ovl[0].a_iid)
    //  Overlaps loaded are for this read, YAY!  We assume that ALL overlaps are loaded
    //  for this iid.
    return(ovlLen);

  //  Until we load the correct overlap, repeat.

  do {
    //  Count the number of overlaps we would load
    ovlLen = numberOfOverlaps();

    if (ovlLen == 0)
      //  Quit now if there are no overlaps.  This simplifies the rest of the loop.
      return(0);

    //  Allocate space for these overlaps.
    while (ovlMax < ovlLen) {
      ovlMax *= 2;
      delete [] ovl;
      ovl = ovOverlap::allocateOverlaps(_gkp, ovlMax);
    }

    //  Load the overlaps
    ovlLen = readOverlaps(ovl, ovlMax);

    //  If we read overlaps for a fragment after 'iid', we're done.  The client will properly save
    //  these overlaps until the iid becomes active.
    //
    if (iid < ovl[0].a_iid)
      return(0);

    //  If we've found the overlaps, we're still done, we also return the number of overlaps.
    //
    if (iid == ovl[0].a_iid)
      return(ovlLen);

    //  On the otherhand, if we read overlaps for a fragment before 'iid', we can either keep
    //  reading until we find the overlaps for this fragment, or jump to the correct spot to read
    //  overlaps.
    //
    //  The rule is simple.  If we're within 50 of the correct IID, keep streaming.  Otherwise, make
    //  a jump.  setRange() seems to ALWAYS close and open a file, which is somewhat expensive,
    //  especially if the file doesn't actually change.
    //
    if (50 < iid - ovl[0].a_iid)
      setRange(iid, UINT32_MAX);

  } while (ovl[0].a_iid < iid);

  //  Code can't get here.
  //    If we ran out of overlaps, the first return is used.
  //    If we read past the iid we're looking for, the second return is used.
  //    If we found the overlaps we're looking for, the thrd return is used.
  //    If we didn't find the overlaps, we loop.

  assert(0);
  return(0);
}













void
ovStore::setRange(uint32 firstIID, uint32 lastIID) {
  char            name[FILENAME_MAX];

  //  make the index be one record per read iid, regardless, then we
  //  can quickly grab the correct record, and seek to the start of
  //  those overlaps

  if (firstIID > _info._largestIID)
    firstIID = _info._largestIID + 1;
  if (lastIID >= _info._largestIID)
    lastIID = _info._largestIID;

  //  If our range is invalid (firstIID > lastIID) we keep going, and
  //  let readOverlap() deal with it.

  AS_UTL_fseek(_offtFile, (size_t)firstIID * sizeof(ovStoreOfft), SEEK_SET);

  //  Unfortunately, we need to actually read the record to figure out
  //  where to position the overlap stream.  If the read fails, we
  //  silently return, letting readOverlap() deal with
  //  the problem.

  _offt.clear();

  //  Everything should notice that offsetFile is at EOF and not try
  //  to find overlaps, but, just in case, we set invalid first/last
  //  IIDs.
  //
  _firstIIDrequested = firstIID;
  _lastIIDrequested  = lastIID;

  if (0 == AS_UTL_safeRead(_offtFile, &_offt, "ovStore::setRange::offset", sizeof(ovStoreOfft), 1))
    return;

  _overlapsThisFile = 0;
  _currentFileIndex = _offt._fileno;

  delete _bof;

  sprintf(name, "%s/%04d", _storePath, _currentFileIndex);
  _bof = new ovFile(name, ovFileNormal);

  _bof->seekOverlap(_offt._offset);
}



void
ovStore::resetRange(void) {
  char            name[FILENAME_MAX];

  rewind(_offtFile);

  _offt.clear();

  _overlapsThisFile = 0;
  _currentFileIndex = 1;

  delete _bof;

  sprintf(name, "%s/%04d", _storePath, _currentFileIndex);
  _bof = new ovFile(name, ovFileNormal);

  _firstIIDrequested = _info._smallestIID;
  _lastIIDrequested  = _info._largestIID;
}





void
ovStore::writeOverlap(ovOverlap *overlap) {
  char            name[FILENAME_MAX];

  assert(_isOutput == TRUE);

  if (_offt._a_iid > overlap->a_iid) {
    //  Woah!  The last overlap we saw is bigger than the one we have now?!
    fprintf(stderr, "LAST:  a:"F_U32"\n", _offt._a_iid);
    fprintf(stderr, "THIS:  a:"F_U32" b:"F_U32"\n", overlap->a_iid, overlap->b_iid);
  }
  assert(_offt._a_iid <= overlap->a_iid);

  if (_info._smallestIID > overlap->a_iid)
    _info._smallestIID = overlap->a_iid;
  if (_info._largestIID < overlap->a_iid)
    _info._largestIID = overlap->a_iid;


  //  If we don't have an output file yet, or the current file is
  //  too big, open a new file.
  //
  if ((_bof) && (_overlapsThisFile >= 1024 * 1024 * 1024 / _bof->recordSize())) {
    delete _bof;

    _bof              = NULL;
    _overlapsThisFile = 0;
  }

  if (_bof == NULL) {
    char  name[FILENAME_MAX];

    _currentFileIndex++;

    sprintf(name, "%s/%04d", _storePath, _currentFileIndex);
    _bof = new ovFile(name, ovFileNormalWrite);
  }


  //  Put the index to disk, filling any gaps
  //
  if ((_offt._numOlaps != 0) &&
      (_offt._a_iid != overlap->a_iid)) {

    while (_offm._a_iid < _offt._a_iid) {
      _offm._fileno    = _offt._fileno;
      _offm._offset    = _offt._offset;
      _offm._overlapID = _offt._overlapID;  //  Not needed, but makes life easier

      AS_UTL_safeWrite(_offtFile, &_offm, "ovStore::writeOverlap::offset", sizeof(ovStoreOfft), 1);

      _offm._a_iid++;
    }

    _offm._a_iid++;  //  One more, since this iid is not missing -- we write it next!

    AS_UTL_safeWrite(_offtFile, &_offt, "AS_OVS_writeOverlapToStore offset", sizeof(ovStoreOfft), 1);

    _offt._numOlaps = 0;  //  Reset; this new id has no overlaps yet.
  }


  //  Update the index if this is the first overlap for this a_iid
  //
  if (_offt._numOlaps == 0) {
    _offt._a_iid     = overlap->a_iid;
    _offt._fileno    = _currentFileIndex;
    _offt._offset    = _overlapsThisFile;
    _offt._overlapID = _info._numOverlapsTotal;
  }

  //AS_OVS_accumulateStats(ovs, overlap);
  _bof->writeOverlap(overlap);

  _offt._numOlaps++;
  _info._numOverlapsTotal++;
  _overlapsThisFile++;
}



void
ovStore::writeOverlap(ovOverlap *overlap, uint32 maxOverlapsThisFile) {
  char            name[FILENAME_MAX];

  assert(_isOutput == TRUE);

  _currentFileIndex++;
  _overlapsThisFile = 0;

  for (uint64 i=0; i < maxOverlapsThisFile; i++ ) {
    //  All overlaps will be sorted by a_iid
    if (_offt._a_iid > overlap[i].a_iid) {
      fprintf(stderr, "LAST:  a:"F_U32"\n", _offt._a_iid);
      fprintf(stderr, "THIS:  a:"F_U32" b:"F_U32"\n", overlap[i].a_iid, overlap[i].b_iid);
    }

    assert(_offt._a_iid <= overlap[i].a_iid);

    if (_info._smallestIID > overlap[i].a_iid)
      _info._smallestIID = overlap[i].a_iid;
    if (_info._largestIID < overlap[i].a_iid)
      _info._largestIID = overlap[i].a_iid;


    //  Put the index to disk, filling any gaps
    if ((_offt._numOlaps != 0) && (_offt._a_iid != overlap[i].a_iid)) {

      while (_offm._a_iid < _offt._a_iid) {
        _offm._fileno    = _offt._fileno;
        _offm._offset    = _offt._offset;
        _offm._overlapID = _offt._overlapID;  //  Not needed, but makes life easier

        AS_UTL_safeWrite(_offtFile, &_offm, "AS_OVS_writeOverlapToStore offset", sizeof(ovStoreOfft), 1);

        _offm._a_iid++;
      }

      _offm._a_iid++;  //  One more, since this iid is not missing -- we write it next!

      AS_UTL_safeWrite(_offtFile, &_offt, "AS_OVS_writeOverlapToStore offset", sizeof(ovStoreOfft), 1);

      _offt._numOlaps  = 0;    //  Reset; this new id has no overlaps yet.
    }

    //  Update the index if this is the first overlap for this a_iid
    if (_offt._numOlaps == 0) {
      _offt._a_iid     = overlap[i].a_iid;
      _offt._fileno    = _currentFileIndex;
      _offt._offset    = _overlapsThisFile;
      _offt._overlapID = _info._numOverlapsTotal;
    }

    _offt._numOlaps++;
    _info._numOverlapsTotal++;
    _overlapsThisFile++;
  }

  fprintf(stderr,"Done building index for dumpfile %d.\n",_currentFileIndex);
}




uint64
ovStore::numOverlapsInRange(void) {
  size_t                     originalposition = 0;
  uint64                     i = 0;
  uint64                     len = 0;
  ovStoreOfft  *offsets = NULL;
  uint64                     numolap = 0;

  if (_firstIIDrequested > _lastIIDrequested)
    return(0);

  originalposition = AS_UTL_ftell(_offtFile);

  AS_UTL_fseek(_offtFile, (size_t)_firstIIDrequested * sizeof(ovStoreOfft), SEEK_SET);

  //  Even if we're doing a whole human-size store, this allocation is
  //  (a) temporary and (b) only 512MB.  The only current consumer of
  //  this code is FragCorrectOVL.c, which doesn't run on the whole
  //  human, it runs on ~24 pieces, which cuts this down to < 32MB.

  len = _lastIIDrequested - _firstIIDrequested + 1;
  offsets = new ovStoreOfft [len];

  if (len != AS_UTL_safeRead(_offtFile, offsets, "AS_OVS_numOverlapsInRange", sizeof(ovStoreOfft), len)) {
    fprintf(stderr, "AS_OVS_numOverlapsInRange()-- short read on offsets!\n");
    exit(1);
  }

  for (i=0; i<len; i++)
    numolap += offsets[i]._numOlaps;

  delete [] offsets;

  AS_UTL_fseek(_offtFile, originalposition, SEEK_SET);

  return(numolap);
}



uint32 *
ovStore::numOverlapsPerFrag(uint32 &firstFrag, uint32 &lastFrag) {

  if (_firstIIDrequested > _lastIIDrequested)
    return(NULL);

  firstFrag = _firstIIDrequested;
  lastFrag  = _lastIIDrequested;

  size_t originalPosition = AS_UTL_ftell(_offtFile);

  AS_UTL_fseek(_offtFile, (size_t)_firstIIDrequested * sizeof(ovStoreOfft), SEEK_SET);

  //  Even if we're doing a whole human-size store, this allocation is
  //  (a) temporary and (b) only 512MB.  The only current consumer of
  //  this code is FragCorrectOVL.c, which doesn't run on the whole
  //  human, it runs on ~24 pieces, which cuts this down to < 32MB.

  uint64 len = _lastIIDrequested - _firstIIDrequested + 1;

  ovStoreOfft  *offsets = new ovStoreOfft [len];
  uint32       *numolap = new uint32      [len];

  uint64 act = AS_UTL_safeRead(_offtFile, offsets, "ovStore::numOverlapsInRange::offsets", sizeof(ovStoreOfft), len);

  if (len != act)
    fprintf(stderr, "AS_OVS_numOverlapsPerFrag()-- short read on offsets!  Expected len="F_U64" read act="F_U64"\n", len, act), exit(1);

  for (uint64 i=0; i<len; i++)
    numolap[i] = offsets[i]._numOlaps;

  delete [] offsets;

  AS_UTL_fseek(_offtFile, originalPosition, SEEK_SET);

  return(numolap);
}










void
ovStore::addEvalues(uint32 bgnID, uint32 endID, uint16 *evalues, uint64 evaluesLen) {

  char  name[FILENAME_MAX];
  sprintf(name, "%s/evalues", _storePath);

  //  If we have an opened memory mapped file, and it isn't open for writing, close it.

  if ((_evaluesMap) && (_evaluesMap->type() == memoryMappedFile_readOnly)) {
    fprintf(stderr, "WARNING: closing read-only evalues file.\n");
    delete _evaluesMap;

    _evaluesMap = NULL;
    _evalues    = NULL;
  }

  //  Remove a bogus evalues file if one exists.

  if ((AS_UTL_fileExists(name) == true) &&
      (AS_UTL_sizeOfFile(name) != (sizeof(uint16) * _info._numOverlapsTotal))) {
    fprintf(stderr, "WARNING: existing evalues file is incorrect size: should be "F_U64" bytes, is "F_U64" bytes.  Removing.\n",
            (sizeof(uint16) * _info._numOverlapsTotal), AS_UTL_sizeOfFile(name));
    AS_UTL_unlink(name);
  }

  //  Make a new evalues file if one doesn't exist.

  if (AS_UTL_fileExists(name) == false) {
    fprintf(stderr, "Creating evalues file for "F_U64" overlaps.\r", _info._numOverlapsTotal);

    errno = 0;
    FILE *F = fopen(name, "w");
    if (errno)
      fprintf(stderr, "Failed to make evalues file '%s': %s\n", name, strerror(errno)), exit(1);

    uint16  *Z  = new uint16 [1048576];
    uint64   Zn = 0;

    memset(Z, 0, sizeof(uint16) * 1048576);

    while (Zn < _info._numOverlapsTotal) {
      uint64  S = (Zn + 1048576 < _info._numOverlapsTotal) ? 1048576 : _info._numOverlapsTotal - Zn;

      AS_UTL_safeWrite(F, Z, "zero evalues", sizeof(uint16), S);

      Zn += S;

      fprintf(stderr, "Creating evalues file for "F_U64" overlaps....%07.3f%%\r",
              _info._numOverlapsTotal, 100.0 * Zn / _info._numOverlapsTotal);
    }

    fprintf(stderr, "Creating evalues file for "F_U64" overlaps....%07.3f%%\n",
            _info._numOverlapsTotal, 100.0 * Zn / _info._numOverlapsTotal);

    fclose(F);
  }

  //  Open the evalues file if it isn't already opened

  if (_evalues == NULL) {
    _evaluesMap = new memoryMappedFile(name, memoryMappedFile_readWrite);
    _evalues    = (uint16 *)_evaluesMap->get(0);
  }

  //  Figure out the overlap ID for the first overlap associated with bgnID

  setRange(bgnID, endID);

  //  Load the evalues from 'evalues'

  for (uint64 ii=0; ii<evaluesLen; ii++)
    _evalues[_offt._overlapID + ii] = evalues[ii];

  //  That's it.  Deleting the ovStore object will close the memoryMappedFile.  It's left open
  //  for more updates.
}











//  For the parallel sort, write a block of sorted overlaps into a single file, with index and info.

void
writeOverlaps(char       *storePath,
              ovOverlap *ovls,
              uint64      ovlsLen,
              uint32      fileID) {

  char                        name[FILENAME_MAX];

  uint32                      currentFileIndex = fileID;
  uint64                      overlapsThisFile = 0;

  ovStoreInfo    info;

  info._ovsMagic              = 1;
  info._ovsVersion            = ovStoreVersion;
  info._UNUSED                = 0;
  info._smallestIID           = UINT64_MAX;
  info._largestIID            = 0;
  info._numOverlapsTotal      = 0;
  info._highestFileIndex      = 0;
  info._maxReadLenInBits      = AS_MAX_READLEN_BITS;

  ovStoreOfft    offt;
  ovStoreOfft    offm;

  offt._a_iid     = offm._a_iid     = ovls[0].a_iid;
  offt._fileno    = offm._fileno    = fileID;
  offt._offset    = offm._offset    = 0;
  offt._numOlaps  = offm._numOlaps  = 0;
  offt._overlapID = offm._overlapID = 0;

  //  Create the output file

  sprintf(name, "%s/%04d", storePath, fileID);
  ovFile *bof = new ovFile(name, ovFileNormalWrite);

  //  Create the index file

  sprintf(name,"%s/%04d.index", storePath, fileID);

  errno = 0;
  FILE *offtFile=fopen(name,"w");
  if (errno)
    fprintf(stderr, "ERROR: Failed to open '%s' for writing: %s\n", name, strerror(errno)), exit(1);

  //  Dump the overlaps

  fprintf(stderr, "Writing "F_U64" overlaps.\n", ovlsLen);

  for (uint64 i=0; i<ovlsLen; i++ ) {
    bof->writeOverlap(ovls + i);

    if (offt._a_iid > ovls[i].a_iid) {
      fprintf(stderr, "LAST:  a:"F_U32"\n", offt._a_iid);
      fprintf(stderr, "THIS:  a:"F_U32" b:"F_U32"\n", ovls[i].a_iid, ovls[i].b_iid);
    }
    assert(offt._a_iid <= ovls[i].a_iid);

    info._smallestIID = MIN(info._smallestIID, ovls[i].a_iid);
    info._largestIID  = MAX(info._largestIID,  ovls[i].a_iid);

    //  Put the index to disk, filling any gaps

    if ((offt._numOlaps != 0) && (offt._a_iid != ovls[i].a_iid)) {
      while (offm._a_iid < offt._a_iid) {
        offm._fileno     = offt._fileno;
        offm._offset     = offt._offset;
        offm._overlapID  = offt._overlapID;  //  Not needed, but makes life easier
        offm._numOlaps   = 0;

        AS_UTL_safeWrite(offtFile, &offm, "AS_OVS_writeOverlapToStore offt", sizeof(ovStoreOfft), 1);
        offm._a_iid++;
      }

      //  One more, since this iid is not offm -- we write it next!
      offm._a_iid++;

      AS_UTL_safeWrite(offtFile, &offt, "AS_OVS_writeOverlapToStore offt", sizeof(ovStoreOfft), 1);

      offt._overlapID += offt._numOlaps;  //  The next block of overlaps starts with this ID
      offt._numOlaps   = 0;               //  The next block has no overlaps yet.
    }

    //  Update the index if this is the first overlap for this a_iid

    if (offt._numOlaps == 0) {
      offt._a_iid   = ovls[i].a_iid;
      offt._fileno  = currentFileIndex;
      offt._offset  = overlapsThisFile;
    }

    offt._numOlaps++;

    info._numOverlapsTotal++;

    overlapsThisFile++;
  }

  //  Close the output file.

  delete bof;

  //  Write the final (empty) index entries.

  while (offm._a_iid < offt._a_iid) {
    offm._fileno     = offt._fileno;
    offm._offset     = offt._offset;
    offm._overlapID  = offt._overlapID;  //  Not needed, but makes life easier
    offm._numOlaps   = 0;

    AS_UTL_safeWrite(offtFile, &offm, "AS_OVS_writeOverlapToStore offt", sizeof(ovStoreOfft), 1);
    offm._a_iid++;
  }

  //  And the final (real) index entry.  We could, but don't need to, update overlapID with the
  //  number of overlaps in this block.

  AS_UTL_safeWrite(offtFile, &offt, "AS_OVS_writeOverlapToStore offt", sizeof(ovStoreOfft), 1);

  fclose(offtFile);

  //  In the nasty case that there were no overlaps in this slice, set meaningful smallest and
  //  largest.  Well, at least, set non-nonsense smallest and largest.

  if (overlapsThisFile == 0) {
    info._smallestIID = 0;
    info._largestIID  = 0;
  }

  //  Write the info, and some stats for the user.

  sprintf(name,"%s/%04d.info", storePath, fileID);

  errno = 0;
  FILE *F = fopen(name, "w");
  if (errno)
    fprintf(stderr, "ERROR: Failed to open '%s' for writing: %s\n", name, strerror(errno)), exit(1);

  AS_UTL_safeWrite(F, &info, "Partition ovs file", sizeof(ovStoreInfo), 1);

  fclose(F);

  fprintf(stderr, "Wrote "F_U64" overlaps into '%s'\n", info._numOverlapsTotal, name);
  fprintf(stderr, "  Smallest "F_U64"\n", info._smallestIID);
  fprintf(stderr, "  Largest  "F_U64"\n", info._largestIID);
}




//  For the parallel sort, but also generally applicable, test that the index is sane.

bool
testIndex(char *ovlName,
          bool  doFixes) {
  char name[FILENAME_MAX];
  FILE *I = NULL;
  FILE *F = NULL;

  sprintf(name, "%s/index", ovlName);

  errno = 0;
  I = fopen(name, "r");
  if (errno)
    fprintf(stderr, "ERROR: Failed to open '%s' for reading: %s\n", name, strerror(errno)), exit(1);

  //fprintf(stderr, "TESTING '%s'\n", name);

  if (doFixes) {
    sprintf(name, "%s/index.fixed", ovlName);

    errno = 0;
    F = fopen(name, "w");
    if (errno)
      fprintf(stderr, "ERROR: Failed to open '%s' for writing: %s\n", name, strerror(errno)), exit(1);

    //fprintf(stderr, "WITH FIXES TO '%s'\n", name);
  }

  ovStoreOfft  O;

  uint32  curIID = 0;
  uint32  minIID = UINT32_MAX;
  uint32  maxIID = 0;

  uint32  nErrs = 0;

  while (1 == AS_UTL_safeRead(I, &O, "offset", sizeof(ovStoreOfft), 1)) {
    bool  maxIncreases   = (maxIID < O._a_iid);
    bool  errorDecreased = ((O._a_iid < curIID));
    bool  errorGap       = ((O._a_iid > 0) && (curIID + 1 != O._a_iid));

    if (O._a_iid < minIID)
      minIID = O._a_iid;

    if (maxIncreases)
      maxIID = O._a_iid;

    if (errorDecreased)
      fprintf(stderr, "ERROR: index decreased from "F_U32" to "F_U32"\n", curIID, O._a_iid), nErrs++;
    else if (errorGap)
      fprintf(stderr, "ERROR: gap between "F_U32" and "F_U32"\n", curIID, O._a_iid), nErrs++;

    if ((maxIncreases == true) && (errorGap == false)) {
      if (doFixes)
        AS_UTL_safeWrite(F, &O, "offset", sizeof(ovStoreOfft), 1);

    } else if (O._numOlaps > 0) {
      fprintf(stderr, "ERROR: lost overlaps a_iid "F_U32" fileno "F_U32" offset "F_U32" numOlaps "F_U32"\n",
              O._a_iid, O._fileno, O._offset, O._numOlaps);
    }

    curIID = O._a_iid;
  }

  fclose(I);

  if (F)
    fclose(F);

  return(nErrs == 0);
}





//  For the parallel sort, merge index and info files into one, clean up the intermediates.

void
mergeInfoFiles(char       *storePath,
               uint32      nPieces) {
  ovStoreInfo    infopiece;
  ovStoreInfo    info;

  info._ovsMagic              = ovStoreMagic;
  info._ovsVersion            = ovStoreVersion;
  info._smallestIID           = UINT64_MAX;
  info._largestIID            = 0;
  info._numOverlapsTotal      = 0;
  info._highestFileIndex      = nPieces;
  info._maxReadLenInBits      = AS_MAX_READLEN_BITS;

  ovStoreOfft offm;

  offm._a_iid     = 0;
  offm._fileno    = 1;
  offm._offset    = 0;
  offm._numOlaps  = 0;
  offm._overlapID = 0;

  //  Open the new master index output file

  char            name[FILENAME_MAX];

  sprintf(name, "%s/index", storePath);

  errno = 0;
  FILE  *idx = fopen(name, "w");
  if (errno)
    fprintf(stderr, "ERROR: Failed to open '%s': %s\n", name, strerror(errno)), exit(1);

  //  Special case, we need an empty index for the zeroth fragment.

  AS_UTL_safeWrite(idx, &offm, "ovStore::mergeInfoFiles::offsetZero", sizeof(ovStoreOfft), 1);

  //  Sanity checking, compare the number of overlaps processed against the overlapID
  //  of each ovStoreOfft.

  uint64  totalOverlaps = 0;

  //  Process each

  for (uint32 i=1; i<=nPieces; i++) {
    sprintf(name, "%s/%04d.info", storePath, i);

    fprintf(stderr, "Processing '%s'\n", name);

    if (AS_UTL_fileExists(name, FALSE, FALSE) == false) {
      fprintf(stderr, "ERROR: file '%s' not found.\n", name);
      exit(1);
    }

    {
      errno = 0;
      FILE *F = fopen(name, "r");
      if (errno)
        fprintf(stderr, "ERROR: Failed to open '%s': %s\n", name, strerror(errno)), exit(1);
      AS_UTL_safeRead(F, &infopiece, "ovStore::mergeInfoFiles::infopiece", sizeof(ovStoreInfo), 1);
      fclose(F);
    }

    //  Add empty index elements for missing overlaps

    if (infopiece._numOverlapsTotal == 0) {
      fprintf(stderr, "  No overlaps found.\n");
      continue;
    }

    assert(infopiece._smallestIID <= infopiece._largestIID);

    if (info._largestIID + 1 < infopiece._smallestIID)
      fprintf(stderr, "  Adding empty records for fragments "F_U64" to "F_U64"\n",
              info._largestIID + 1, infopiece._smallestIID - 1);

    while (info._largestIID + 1 < infopiece._smallestIID) {
      offm._a_iid     = info._largestIID + 1;
      //offm._fileno    = set below, where the recs are written to the master file
      //offm._offset    = set below, where the recs are written to the master file

      AS_UTL_safeWrite(idx, &offm, "ovStore::mergeInfoFiles::offsets", sizeof(ovStoreOfft), 1);

      info._largestIID++;
    }

    //  Copy index elements for existing overlaps.  While copying, update the supposed position
    //  of any fragments with no overlaps.  Without doing this, accessing the store beginning
    //  or ending at such a fragment will fail.

    {
      sprintf(name, "%s/%04d.index", storePath, i);

      errno = 0;
      FILE  *F = fopen(name, "r");
      if (errno)
        fprintf(stderr, "ERROR: Failed to open '%s': %s\n", name, strerror(errno)), exit(1);

      uint32          recsLen = 0;
      uint32          recsMax = 1024 * 1024;
      ovStoreOfft    *recs    = new ovStoreOfft [recsMax];

      recsLen = AS_UTL_safeRead(F, recs, "ovStore::mergeInfoFiles::offsetsLoad", sizeof(ovStoreOfft), recsMax);

      if (recsLen > 0) {
        if (info._largestIID + 1 != recs[0]._a_iid)
          fprintf(stderr, "ERROR: '%s' starts with iid "F_U32", but store only up to "F_U64"\n",
                  name, recs[0]._a_iid, info._largestIID);
        assert(info._largestIID + 1 == recs[0]._a_iid);
      }

      while (recsLen > 0) {

        //  Update location of missing reads.

        offm._fileno     = recs[recsLen-1]._fileno;
        offm._offset     = recs[recsLen-1]._offset;

        //  Update overlapID for each record.

        for (uint32 rr=0; rr<recsLen; rr++) {
          recs[rr]._overlapID += info._numOverlapsTotal;

          if (recs[rr]._numOlaps > 0)
            assert(recs[rr]._overlapID == totalOverlaps);

          totalOverlaps += recs[rr]._numOlaps;
        }

        //  Write the records, read next batch

        AS_UTL_safeWrite(idx, recs, "ovStore::mergeInfoFiles::offsetsWrite", sizeof(ovStoreOfft), recsLen);

        recsLen = AS_UTL_safeRead(F, recs, "ovStore::mergeInfoFiles::offsetsReLoad", sizeof(ovStoreOfft), recsMax);
      }

      delete [] recs;

      fclose(F);
    }

    //  Update the info block to include the overlaps we just added

    info._smallestIID = MIN(info._smallestIID, infopiece._smallestIID);
    info._largestIID  = MAX(info._largestIID,  infopiece._largestIID);

    info._numOverlapsTotal += infopiece._numOverlapsTotal;

    fprintf(stderr, "  Now finished with fragments "F_U64" to "F_U64" -- "F_U64" overlaps.\n",
            info._smallestIID, info._largestIID, info._numOverlapsTotal);
  }

  fclose(idx);


  //  Dump the new store info file

  {
    sprintf(name, "%s/info", storePath);

    errno = 0;
    FILE  *F = fopen(name, "w");
    if (errno)
      fprintf(stderr, "ERROR: Failed to open '%s': %s\n", name, strerror(errno)), exit(1);

    AS_UTL_safeWrite(F, &info, "ovStore::mergeInfoFiles::finalInfo", sizeof(ovStoreInfo), 1);

    fclose(F);
  }

  fprintf(stderr, "\n");
  fprintf(stderr, "Index finalized for reads "F_U64" to "F_U64" with "F_U64" overlaps.\n",
          info._smallestIID,
          info._largestIID,
          info._numOverlapsTotal);
}







//
//
//  For overlap store building, both sequential and parallel.  Overlap filtering.
//
//



#define OBT_FAR5PRIME        (29)
#define OBT_MIN_LENGTH       (75)



//  Are the 5' end points very different?  If the overlap is flipped, then, yes, they are.
static
bool
isOverlapDifferent(ovOverlap &ol) {
  bool   isDiff = true;

  if (ol.flipped() == false) {
    if (ol.a_bgn() > ol.b_bgn())
      isDiff = ((ol.a_bgn() - ol.b_bgn()) > OBT_FAR5PRIME) ? (true) : (false);
    else
      isDiff = ((ol.b_bgn() - ol.a_bgn()) > OBT_FAR5PRIME) ? (true) : (false);
  }

  return(isDiff);
}


//  Is the overlap long?
static
bool
isOverlapLong(ovOverlap &ol) {
  int32 ab    = ol.a_bgn();
  int32 ae    = ol.a_end();
  int32 bb    = ol.b_bgn();
  int32 be    = ol.b_end();

  int32 Alength = ae - ab;
  int32 Blength = be - bb;

  if (be < bb)
    Blength = bb - be;

  return(((Alength > OBT_MIN_LENGTH) && (Blength > OBT_MIN_LENGTH)) ? (true) : (false));
}




void
ovStoreFilter::filterOverlap(ovOverlap       &foverlap,
                             ovOverlap       &roverlap) {

  //  Quick sanity check on IIDs.

  if ((foverlap.a_iid == 0) ||
      (foverlap.b_iid == 0) ||
      (foverlap.a_iid >= maxID) ||
      (foverlap.b_iid >= maxID)) {
    char ovlstr[256];

    fprintf(stderr, "Overlap has IDs out of range (maxID "F_U32"), possibly corrupt input data.\n", maxID);
    fprintf(stderr, "  coords -- %s\n", foverlap.toString(ovlstr, ovOverlapAsCoords, false));
    fprintf(stderr, "  hangs  -- %s\n", foverlap.toString(ovlstr, ovOverlapAsHangs, false));
    exit(1);
  }

  //  Make the reverse overlap (important, AFTER resetting the erate-based 'for' flags).

  roverlap.swapIDs(foverlap);


  //  Ignore high error overlaps

  if ((foverlap.evalue() > maxEvalue)) {
    foverlap.dat.ovl.forUTG = false;
    foverlap.dat.ovl.forOBT = false;
    foverlap.dat.ovl.forDUP = false;

    roverlap.dat.ovl.forUTG = false;
    roverlap.dat.ovl.forOBT = false;
    roverlap.dat.ovl.forDUP = false;

    skipERATE++;
    skipERATE++;
  }




  //  Don't OBT if not requested.

  if ((foverlap.dat.ovl.forOBT == false) && (skipReadOBT[foverlap.a_iid] == true)) {
    foverlap.dat.ovl.forOBT = false;
    skipOBT++;
  }

  if ((roverlap.dat.ovl.forOBT == false) && (skipReadOBT[roverlap.a_iid] == true)) {
    roverlap.dat.ovl.forOBT = false;
    skipOBT++;
  }

  //  If either overlap is good for either obt or dup, compute if it is different and long.  These
  //  are the same for both foverlap and roverlap.

  bool  isDiff = isOverlapDifferent(foverlap);
  bool  isLong = isOverlapLong(foverlap);

  //  Remove the bad-for-OBT overlaps.

  if ((isDiff == false) && (foverlap.dat.ovl.forOBT == true)) {
    foverlap.dat.ovl.forOBT = false;
    skipOBTbad++;
  }

  if ((isDiff == false) && (roverlap.dat.ovl.forOBT == true)) {
    roverlap.dat.ovl.forOBT = false;
    skipOBTbad++;
  }

  //  Remove the too-short-for-OBT overlaps.

  if ((isLong == false) && (foverlap.dat.ovl.forOBT == true)) {
    foverlap.dat.ovl.forOBT = false;
    skipOBTshort++;
  }

  if ((isLong == false) && (roverlap.dat.ovl.forOBT == true)) {
    roverlap.dat.ovl.forOBT = false;
    skipOBTshort++;
  }




  //  Don't dedupe if not requested.

  if ((foverlap.dat.ovl.forDUP == true) && (skipReadDUP[foverlap.a_iid] == true)) {
    foverlap.dat.ovl.forDUP = false;
    skipDUP++;
  }

  if ((roverlap.dat.ovl.forDUP == true) && (skipReadDUP[roverlap.b_iid] == true)) {
    roverlap.dat.ovl.forDUP = false;
    skipDUP++;
  }

  //  Remove the bad-for-DUP overlaps.

#if 0
  //  Nah, do this in dedupe, since parameters can change.
  if ((isDiff == true) && (foverlap.dat.ovl.forDUP == true)) {
    foverlap.dat.ovl.forDUP = false;
    skipDUPdiff++;
  }

  if ((isDiff == true) && (roverlap.dat.ovl.forDUP == true)) {
    roverlap.dat.ovl.forDUP = false;
    skipDUPdiff++;
  }
#endif

  //  Can't have duplicates between libraries.

  if (((foverlap.dat.ovl.forDUP == true) ||
       (roverlap.dat.ovl.forDUP == true)) &&
      (gkp->gkStore_getRead(foverlap.a_iid)->gkRead_libraryID() != gkp->gkStore_getRead(foverlap.b_iid)->gkRead_libraryID())) {

    if ((foverlap.dat.ovl.forDUP == true)) {
      foverlap.dat.ovl.forDUP = false;
      skipDUPlib++;
    }

    if ((roverlap.dat.ovl.forDUP == true)) {
      roverlap.dat.ovl.forDUP = false;
      skipDUPlib++;
    }
  }

  //  All done with the filtering, record some counts.

  if (foverlap.dat.ovl.forUTG == true)  saveUTG++;
  if (foverlap.dat.ovl.forOBT == true)  saveOBT++;
  if (foverlap.dat.ovl.forDUP == true)  saveDUP++;

  if (roverlap.dat.ovl.forUTG == true)  saveUTG++;
  if (roverlap.dat.ovl.forOBT == true)  saveOBT++;
  if (roverlap.dat.ovl.forDUP == true)  saveDUP++;
}




void
ovStoreFilter::reportFate(void) {
  fprintf(stderr, "overlap fate:\n");
  fprintf(stderr, "%16"F_U64P" SAVE  - overlaps output (for unitigging)\n", saveUTG);
  fprintf(stderr, "%16"F_U64P" SAVE  - overlaps output (for OBT)\n", saveOBT);
  fprintf(stderr, "%16"F_U64P" SAVE  - overlaps output (for dedupe)\n", saveDUP);
  fprintf(stderr, "\n");
  fprintf(stderr, "%16"F_U64P" ERATE - low quality, more than %.3f fraction error\n", skipERATE, AS_OVS_decodeEvalue(maxEvalue));
  fprintf(stderr, "\n");
  fprintf(stderr, "%16"F_U64P" OBT   - not requested\n", skipOBT);
  fprintf(stderr, "%16"F_U64P" OBT   - too similar\n", skipOBTbad);
  fprintf(stderr, "%16"F_U64P" OBT   - too short\n", skipOBTshort);
  fprintf(stderr, "\n");
  fprintf(stderr, "%16"F_U64P" DUP   - dedupe not requested\n", skipDUP);
  fprintf(stderr, "%16"F_U64P" DUP   - different library\n", skipDUPlib);
  fprintf(stderr, "%16"F_U64P" DUP   - obviously not duplicates\n", skipDUPdiff);
}


void
ovStoreFilter::resetCounters(void) {
  saveUTG         = 0;
  saveOBT         = 0;
  saveDUP         = 0;

  skipERATE       = 0;

  skipOBT         = 0;
  skipOBTbad      = 0;
  skipOBTshort    = 0;

  skipDUP         = 0;
  skipDUPdiff     = 0;
  skipDUPlib      = 0;
}
