
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

#include <vector>
#include <algorithm>

using namespace std;



ovStore::ovStore(const char *path, gkStore *gkp) {
  char  name[FILENAME_MAX];

  if (path == NULL)
    fprintf(stderr, "ovStore::ovStore()-- ERROR: no name supplied.\n"), exit(1);

  if ((path[0] == '-') &&
      (path[1] == 0))
    fprintf(stderr, "ovStore::ovStore()-- ERROR: name cannot be '-' (stdin).\n"), exit(1);

  memset(_storePath, 0, FILENAME_MAX);
  strncpy(_storePath, path, FILENAME_MAX-1);

  _info.clear();
  _gkp = gkp;

  _offtFile  = NULL;
  _offt.clear();
  _offm.clear();

  _evaluesMap = NULL;
  _evalues    = NULL;

  _overlapsThisFile  = 0;
  _currentFileIndex  = 0;
  _bof               = NULL;

  //  Now open the store

  if (_info.load(_storePath) == false)
    fprintf(stderr, "ERROR:  failed to intiialize ovStore '%s'.\n", path), exit(1);

  if (_info.checkIncomplete() == true)
    fprintf(stderr, "ERROR:  directory '%s' is an incomplete ovStore, remove and rebuild.\n", path), exit(1);

  if (_info.checkMagic() == false)
    fprintf(stderr, "ERROR:  directory '%s' is not an ovStore.\n", path), exit(1);

  if (_info.checkVersion() == false)
    fprintf(stderr, "ERROR:  directory '%s' is not a supported ovStore version (store version %u; supported version %u.\n",
            path, _info.getVersion(), _info.getCurrentVersion()), exit(1);

  if (_info.checkSize() == false)
    fprintf(stderr, "ERROR:  directory '%s' is not a supported read length (store is %u bits, AS_MAX_READLEN_BITS is %u).\n",
            path, _info.getSize(), AS_MAX_READLEN_BITS), exit(1);

  //  Open the index

  snprintf(name, FILENAME_MAX, "%s/index", _storePath);

  errno = 0;
  _offtFile = fopen(name, "r");
  if (errno)
    fprintf(stderr, "ERROR:  failed to open offset file '%s': %s\n", name, strerror(errno)), exit(1);

  //  Open and load erates

  snprintf(name, FILENAME_MAX, "%s/evalues", _storePath);

  if (AS_UTL_fileExists(name)) {
    _evaluesMap  = new memoryMappedFile(name, memoryMappedFile_readOnly);
    _evalues     = (uint16 *)_evaluesMap->get(0);
  }

  //  Set the initial range to everything.

  _firstIIDrequested      = _info.smallestID();
  _lastIIDrequested       = _info.largestID();
}




ovStore::~ovStore() {

  if (_evaluesMap) {
    delete _evaluesMap;

    _evaluesMap = NULL;
    _evalues    = NULL;
  }

  delete _bof;

  AS_UTL_closeFile(_offtFile);
}



uint32
ovStore::readOverlap(ovOverlap *overlap) {

  //  If we've finished reading overlaps for the current a_iid, get
  //  another a_iid.  If we hit EOF here, we're all done, no more
  //  overlaps.

  while (_offt._numOlaps == 0)
    if (0 == AS_UTL_safeRead(_offtFile, &_offt, "ovStore::readOverlap::offset", sizeof(ovStoreOfft), 1))
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

    snprintf(name, FILENAME_MAX, "%s/%04d", _storePath, _currentFileIndex);
    _bof = new ovFile(_gkp, name, ovFileNormal);
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

      if (_currentFileIndex > _info.lastFileIndex())
        //  No more files, stop trying to load an overlap.
        break;

      snprintf(name, FILENAME_MAX, "%s/%04d", _storePath, _currentFileIndex);
      _bof = new ovFile(_gkp, name, ovFileNormal);
    }

    //  If the currentFileIndex is invalid, we ran out of overlaps to load.  Don't save that
    //  empty overlap to the list.

    if (_currentFileIndex <= _info.lastFileIndex()) {
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

  if (firstIID > _info.largestID())
    firstIID = _info.largestID() + 1;
  if (lastIID >= _info.largestID())
    lastIID = _info.largestID();

  //  If our range is invalid (firstIID > lastIID) we keep going, and
  //  let readOverlap() deal with it.

  AS_UTL_fseek(_offtFile, (off_t)firstIID * sizeof(ovStoreOfft), SEEK_SET);

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

  snprintf(name, FILENAME_MAX, "%s/%04d", _storePath, _currentFileIndex);
  _bof = new ovFile(_gkp, name, ovFileNormal);

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

  snprintf(name, FILENAME_MAX, "%s/%04d", _storePath, _currentFileIndex);
  _bof = new ovFile(_gkp, name, ovFileNormal);

  _firstIIDrequested = _info.smallestID();
  _lastIIDrequested  = _info.largestID();
}



uint64
ovStore::numOverlapsInRange(void) {
  off_t                      originalposition = 0;
  uint64                     i = 0;
  uint64                     len = 0;
  ovStoreOfft  *offsets = NULL;
  uint64                     numolap = 0;

  if (_firstIIDrequested > _lastIIDrequested)
    return(0);

  originalposition = AS_UTL_ftell(_offtFile);

  AS_UTL_fseek(_offtFile, (off_t)_firstIIDrequested * sizeof(ovStoreOfft), SEEK_SET);

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



//  Returns array with the number of overlaps per read.  If numReads is more than zero,
//  only those reads will be loaded.
//
uint32 *
ovStore::numOverlapsPerRead(uint32 numReads) {

  if ((numReads == 0) && (_gkp != NULL))
    numReads = _gkp->gkStore_getNumReads();

  assert(numReads > 0);

  uint32       *olapsPerRead = new uint32      [numReads+1];
  ovStoreOfft  *offsets      = new ovStoreOfft [numReads+1];

  off_t  originalPosition = AS_UTL_ftell(_offtFile);

  AS_UTL_fseek(_offtFile, 0, SEEK_SET);

  uint32 act = AS_UTL_safeRead(_offtFile, offsets, "ovStore::numOverlapsPerRead", sizeof(ovStoreOfft), _info.largestID()+1);

  if (_info.largestID()+1 != act)
    fprintf(stderr, "ovStore::numOverlapsPerRead()-- short read on offsets!  Read %u entries, store has smallest %u largest %u\n",
            act, _info.smallestID(), _info.largestID()), exit(1);

  for (uint32 ii=0; ii<_info.largestID()+1; ii++)
    olapsPerRead[ii] = offsets[ii]._numOlaps;

  for (uint32 ii=_info.largestID()+1; ii<numReads+1; ii++)
    olapsPerRead[ii] = 0;

  delete [] offsets;

  AS_UTL_fseek(_offtFile, originalPosition, SEEK_SET);

  return(olapsPerRead);
}



class evalueFileMap {
public:
  evalueFileMap() {
    _name  = NULL;
    _bgnID = UINT32_MAX;
    _endID = 0;
    _Nolap = 0;
  }

  char     *_name;

  uint32    _bgnID;
  uint32    _endID;
  uint64    _Nolap;
};

bool
operator<(evalueFileMap const &a, evalueFileMap const &b) {
  return(a._bgnID < b._bgnID);
};



void
ovStore::addEvalues(vector<char *> &fileList) {
  char  evalueName[FILENAME_MAX+1];
  char  evalueTemp[FILENAME_MAX+1];

  //  Handy to have the name of the file we're working with.

  snprintf(evalueTemp, FILENAME_MAX, "%s/evalues.WORKING", _storePath);
  snprintf(evalueName, FILENAME_MAX, "%s/evalues",         _storePath);

  //  If we have an opened memory mapped file, close it.  There _shouldn't_ be one,
  //  as it would exist only if evalues were already added, but it might.  And if it
  //  does exist, nuke it from disk too (well, not quite yet).

  if (_evaluesMap) {
    delete _evaluesMap;

    _evaluesMap = NULL;
    _evalues    = NULL;
  }

  if (AS_UTL_fileExists(evalueName) == true) {
    fprintf(stderr, "WARNING:\n");
    fprintf(stderr, "WARNING: existing evalue file will be overwritten!\n");
    fprintf(stderr, "WARNING:\n");
  }

  //  Scan each file, finding the first overlap in it.

  fprintf(stderr, "\n");
  fprintf(stderr, "Scanning.\n");

  evalueFileMap   *emap = new evalueFileMap [fileList.size()];

  for (uint32 ii=0; ii<fileList.size(); ii++) {
    emap[ii]._name  = fileList[ii];
    emap[ii]._bgnID = 0;
    emap[ii]._endID = 0;
    emap[ii]._Nolap = 0;

    FILE *E = AS_UTL_openInputFile(fileList[ii]);

    AS_UTL_safeRead(E, &emap[ii]._bgnID, "bgnID", sizeof(uint32), 1);
    AS_UTL_safeRead(E, &emap[ii]._endID, "endID", sizeof(uint32), 1);
    AS_UTL_safeRead(E, &emap[ii]._Nolap, "Nolap", sizeof(uint64), 1);

    AS_UTL_closeFile(E);

    fprintf(stderr, "  '%s' covers reads %7" F_U32P "-%-7" F_U32P " with %10" F_U64P " overlaps.\n",
            emap[ii]._name, emap[ii]._bgnID, emap[ii]._endID, emap[ii]._Nolap);
  }

  //  Sort the emap by starting read.

  sort(emap, emap + fileList.size());

  //  Check that all IDs are present.

  for (uint32 ii=1; ii<fileList.size(); ii++) {
    if (emap[ii-1]._endID < emap[ii]._bgnID)
      fprintf(stderr, "Discontinuity between files '%s' and '%s'.\n", emap[ii-1]._name, emap[ii]._name);
  }

  //  Now just copy the new evalues to the real evalues file.  The inputs two 32-bit words and
  //  a 64-bit word at the start we need to ignore.  That's 8 16-bit words.

  fprintf(stderr, "\n");
  fprintf(stderr, "Merging.\n");

  FILE *EO = AS_UTL_openOutputFile(evalueTemp);

  for (uint32 ii=0; ii<fileList.size(); ii++) {
    uint16 *ev = new uint16 [emap[ii]._Nolap + 8];

    AS_UTL_loadFile(emap[ii]._name, ev, emap[ii]._Nolap + 8);

    AS_UTL_safeWrite(EO, ev + 8, "evalues", sizeof(uint16), emap[ii]._Nolap);

    delete [] ev;

    fprintf(stderr, "  '%s' covers reads %7" F_U32P "-%-7" F_U32P "; %10" F_U64P " with overlaps.\n",
            emap[ii]._name, emap[ii]._bgnID, emap[ii]._endID, emap[ii]._Nolap);
  }

  AS_UTL_closeFile(EO, evalueTemp);

  fprintf(stderr, "\n");
  fprintf(stderr, "Renaming.\n");

  AS_UTL_unlink(evalueName);
  AS_UTL_rename(evalueTemp, evalueName);

  fprintf(stderr, "\n");
  fprintf(stderr, "Success!\n");
  fprintf(stderr, "\n");
}
