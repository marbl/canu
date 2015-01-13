
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 2007, J. Craig Venter Institute. All rights reserved.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received (LICENSE.txt) a copy of the GNU General Public
 * License along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *************************************************************************/

static const char *rcsid = "$Id$";

#include "ovStore.H"

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>
#include <fcntl.h>
#include <assert.h>
#include <limits.h>


const uint64 ovStoreVersion         = 2;
const uint64 ovStoreMagic           = 0x4c564f3a67336163;   //  == "ca3g:OVS - store complete
const uint64 ovStoreMagicIncomplete = 0x50564f3a67336163;   //  == "ca3g:OVP - store under construction

//
//  Are backups used anymore??
//

void
ovStore::renameToBackup(char const *name, uint32 index) {
  char   orig[FILENAME_MAX];
  char   bkup[FILENAME_MAX];

  if (name == NULL) {
    sprintf(orig, "%s/%04"F_U32P,    _storePath, index);
    sprintf(bkup, "%s/%04"F_U32P"~", _storePath, index);
  } else {
    sprintf(orig, "%s/%s",  _storePath, name);
    sprintf(bkup, "%s/%s~", _storePath, name);
  }

  errno = 0;
  rename(orig, bkup);
  if (errno)
    fprintf(stderr, "ovStore::renameToBackup()- ERROR: failed to make backup of '%s' into '%s': %s\n", orig, bkup, strerror(errno)), exit(1);
}


void
ovStore::renameFromBackup(char const *name, uint32 index) {
  char   orig[FILENAME_MAX];
  char   bkup[FILENAME_MAX];

  if (name == NULL) {
    sprintf(orig, "%s/%04"F_U32P,    _storePath, index);
    sprintf(bkup, "%s/%04"F_U32P"~", _storePath, index);
  } else {
    sprintf(orig, "%s/%s",  _storePath, name);
    sprintf(bkup, "%s/%s~", _storePath, name);
  }

  errno = 0;
  rename(bkup, orig);
  if (errno)
    fprintf(stderr, "ovStore::renameFromBackup()-- ERROR: failed to restore backup of '%s' into '%s': %s\n", bkup, orig, strerror(errno)), exit(1);
}


void
ovStore::removeBackup(char const *name, uint32 index) {
  char   bkup[FILENAME_MAX];

  if (name == NULL) {
    sprintf(bkup, "%s/%04"F_U32P"~", _storePath, index);
  } else {
    sprintf(bkup, "%s/%s~",  _storePath, name);
  }

  errno = 0;
  unlink(bkup);
  if ((errno) && (errno != ENOENT))
    fprintf(stderr, "ovStore::removeBackup()-- WARNING: failed to remove backup '%s': %s\n", bkup, strerror(errno));
}


void
ovStore::createBackup(void) {

  if (_useBackup == false)
    return;

  renameToBackup("info");
  renameToBackup("index");

  for (uint32 i=1; i<=_info._highestFileIndex; i++)
    renameToBackup(i);
}


void
ovStore::restoreBackup(void) {

  if (_useBackup == false)
    return;

  renameFromBackup("info");
  renameFromBackup("index");

  for (uint32 i=1; i<=_info._highestFileIndex; i++)
    renameFromBackup(i);
}


void
ovStore::removeBackup(void) {

  if (_useBackup == false)
    return;

  removeBackup("info");
  removeBackup("index");

  for (uint32 i=1; i<=_info._highestFileIndex; i++)
    removeBackup(i);
}






ovStore::ovStore(const char *path, ovStoreType cType) {

  if (path == NULL)
    fprintf(stderr, "ovStore::ovStore()-- ERROR: no name supplied.\n"), exit(1);

  if ((path[0] == '-') &&
      (path[1] == 0))
    fprintf(stderr, "ovStore::ovStore()-- ERROR: name cannot be '-' (stdin).\n"), exit(1);

  strcpy(_storePath, path);

  //  If we're doing output, force useBackup and saveSpace to false.

  _isOutput  = (cType & ovStoreWrite)   ? true : false;
  _useBackup = (cType & ovStoreBackup)  ? true : false;  //  NOT USED?
  _saveSpace = (cType & ovStoreReclaim) ? true : false;

  _info._ovsMagic              = 0;
  _info._ovsVersion            = 0;
  _info._numOverlapsPerFile    = 0;  //  not used for reading
  _info._smallestIID           = UINT_MAX;
  _info._largestIID            = 0;
  _info._numOverlapsTotal      = 0;
  _info._highestFileIndex      = 0;
  _info._maxReadLenInBits      = AS_MAX_READLEN_BITS;

  //  If for output, create a new store.
  if (_isOutput == true) {
    AS_UTL_mkdir(path);

    char name[FILENAME_MAX];

    sprintf(name, "%s/info", path);

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
                path), exit(1);
    }

    //  Create a new incomplete info file.
#warning overlap files not exactly 1gb, dont know size of the overlap it is storing

    _info._ovsMagic              = ovStoreMagicIncomplete;
    _info._ovsVersion            = ovStoreVersion;
    _info._numOverlapsPerFile    = 1024 * 1024 * 1024 / (sizeof(ovsOverlap) * sizeof(uint32));
    _info._smallestIID           = UINT64_MAX;
    _info._largestIID            = 0;
    _info._numOverlapsTotal      = 0;
    _info._highestFileIndex      = 0;
    _info._maxReadLenInBits      = AS_MAX_READLEN_BITS;

    errno = 0;
    FILE *ovsinfo = fopen(name, "w");

    if (errno)
      fprintf(stderr, "failed to create overlap store '%s': %s\n", _storePath, strerror(errno)), exit(1);

    AS_UTL_safeWrite(ovsinfo, &_info, "ovStore::ovStore::saveinfo", sizeof(ovStoreInfo), 1);

    fclose(ovsinfo);

    sprintf(name, "%s/index", path);

    errno = 0;
    _offtFile = fopen(name, "w");
    if (errno)
      fprintf(stderr, "AS_OVS_createOverlapStore()-- failed to open offset file '%s': %s\n", name, strerror(errno)), exit(1);

    _overlapsThisFile = 0;
    _currentFileIndex = 0;
    _bof              = NULL;
  }

  //  Otherwise, load the info and maybe create a backup.
  else {
    char  name[FILENAME_MAX];

    sprintf(name, "%s/info", path);
    errno = 0;
    FILE *ovsinfo = fopen(name, "r");
    if (errno)
      fprintf(stderr, "ERROR: directory '%s' is not an ovelrapStore; failed to open info file '%s': %s\n",
              path, name, strerror(errno)), exit(1);

    AS_UTL_safeRead(ovsinfo, &_info, "ovStore::ovStore::info", sizeof(ovStoreInfo), 1);

    fclose(ovsinfo);

    if ((_info._ovsMagic != ovStoreMagic) && (_info._ovsMagic != ovStoreMagicIncomplete))
      fprintf(stderr, "ERROR:  directory '%s' is not an overlapStore; magic number 0x%016"F_X64P" incorrect.\n",
              path, _info._ovsMagic), exit(1);

    if ((_info._ovsMagic != ovStoreMagic) && (_info._ovsMagic != ovStoreMagicIncomplete))
      fprintf(stderr, "ERROR:  overlapStore '%s' is incomplate; creation crashed?\n",
              path), exit(1);

    if (_info._ovsVersion != ovStoreVersion)
      fprintf(stderr, "ERROR:  overlapStore '%s' is version "F_U64"; this code supports only version %d.\n",
              path, _info._ovsVersion, ovStoreVersion), exit(1);

    if (_info._maxReadLenInBits != AS_MAX_READLEN_BITS)
      fprintf(stderr, "ERROR:  overlapStore '%s' is for AS_MAX_READLEN_BITS="F_U64"; this code supports only %d bits.\n",
              path, _info._maxReadLenInBits, AS_MAX_READLEN_BITS), exit(1);

    //  If we're not supposed to be using the backup, load the stats.
    //
#if 0
    if (_useBackup == 0) {
      sprintf(name, "%s/statistics", _storePath);
      errno = 0;
      FILE *ost = fopen(name, "r");
      if (errno)
        fprintf(stderr, "failed to open the stats file '%s': %s\n", name, strerror(errno)), exit(1);
      AS_UTL_safeRead(ost, &_stats, "ovStore::ovStore::stats", sizeof(OverlapStoreStats), 1);
      fclose(ost);
    }
#endif

    createBackup();

    sprintf(name, "%s/index%c", path, _useBackup ? '~' : 0);

    errno = 0;
    _offtFile = fopen(name, "r");

    if (errno)
      fprintf(stderr, "AS_OVS_openOverlapStore()-- failed to open offset file '%s': %s\n", name, strerror(errno)), exit(1);
  }

  _offt._a_iid    = 0;
  _offt._fileno   = 0;
  _offt._offset   = 0;
  _offt._numOlaps = 0;

  _offm._a_iid    = 0;
  _offm._fileno   = 0;
  _offm._offset   = 0;
  _offm._numOlaps = 0;

  _firstIIDrequested = _info._smallestIID;
  _lastIIDrequested  = _info._largestIID;

  _overlapsThisFile  = 0;

  _currentFileIndex  = 0;
  _bof               = NULL;
}




ovStore::~ovStore() {

  removeBackup();

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
    fprintf(stderr, "  info._numOverlapsPerFile = "F_U64"\n", _info._numOverlapsPerFile);
    fprintf(stderr, "  info._smallestIID        = "F_U64"\n", _info._smallestIID);
    fprintf(stderr, "  info._largestIID         = "F_U64"\n", _info._largestIID);
    fprintf(stderr, "  info._numOverlapsTotal   = "F_U64"\n", _info._numOverlapsTotal);
    fprintf(stderr, "  info._highestFileIndex   = "F_U64"\n", _info._highestFileIndex);
    fprintf(stderr, "  info._maxReadLenInBits   = "F_U64"\n", _info._maxReadLenInBits);
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

  delete _gkp;
  delete _bof;

  fclose(_offtFile);
}















uint32
ovStore::readOverlap(ovsOverlap *overlap) {

  assert(_isOutput == FALSE);

  //  If we've finished reading overlaps for the current a_iid, get
  //  another a_iid.  If we hit EOF here, we're all done, no more
  //  overlaps.
  //
 again:

  while (_offt._numOlaps == 0)
    if (0 == AS_UTL_safeRead(_offtFile, &_offt, "ovStore::readOverlap::offset",
                             sizeof(ovStoreOfft), 1))
      return(0);

  //  And if we've exited the range of overlaps requested, return.
  //
  if (_offt._a_iid > _lastIIDrequested)
    return(0);

  while ((_bof == NULL) ||
         (_bof->readOverlap(overlap) == FALSE)) {
    char name[FILENAME_MAX];

    //  We read no overlap, open the next file and try again.

    if (_bof)
      delete _bof;

    if (_saveSpace)
      removeBackup(_storePath, _currentFileIndex);

    _currentFileIndex++;

    sprintf(name, "%s/%04d%c", _storePath, _currentFileIndex, _useBackup ? '~' : 0);
    _bof = new ovFile(name, ovFileNormal);
  }

  overlap->a_iid   = _offt._a_iid;

  _offt._numOlaps--;

  return(1);
}


uint32
ovStore::readOverlaps(ovsOverlap *overlaps, uint32 maxOverlaps, bool restrictToIID) {
  int    numOvl = 0;

  assert(_isOutput == FALSE);

  //  If we've finished reading overlaps for the current a_iid, get
  //  another a_iid.  If we hit EOF here, we're all done, no more
  //  overlaps.
  //
  while (_offt._numOlaps == 0)
    if (0 == AS_UTL_safeRead(_offtFile, &_offt, "ovStore::readOverlaps::offset", sizeof(ovStoreOfft), 1))
      return(0);

  //  And if we've exited the range of overlaps requested, return.
  //
  if (_offt._a_iid > _lastIIDrequested)
    return(0);

  //  Just a query?  Return the number of overlaps we'd want to read

  if ((overlaps == NULL) || (maxOverlaps == 0))
    return(_offt._numOlaps);

  //  Read all the overlaps for this ID.

  assert(_offt._numOlaps <= maxOverlaps);

  while (((restrictToIID == true)  && (_offt._numOlaps > 0)) ||
         ((restrictToIID == false) && (_offt._numOlaps > 0) && (numOvl < maxOverlaps))) {

    //  Read an overlap.  If this fails, open the next partition and read from there.

    while ((_bof == NULL) ||
           (_bof->readOverlap(overlaps + numOvl) == false)) {
      char name[FILENAME_MAX];

      //  We read no overlap, open the next file and try again.

      delete _bof;
      
      if (_saveSpace)
        removeBackup(_storePath, _currentFileIndex);

      _currentFileIndex++;

      if (_currentFileIndex > _info._highestFileIndex)
        //  No more files, stop trying to load an overlap.
        break;

      sprintf(name, "%s/%04d%c", _storePath, _currentFileIndex, _useBackup ? '~' : 0);
      _bof = new ovFile(name, ovFileNormal);
    }

    //  If the currentFileIndex is invalid, we ran out of overlaps to load.  Don't save that
    //  empty overlap to the list.

    if (_currentFileIndex <= _info._highestFileIndex) {
      overlaps[numOvl].a_iid = _offt._a_iid;

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

 returnOverlaps:
  assert(numOvl <= maxOverlaps);

  return(numOvl);
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

  _offt._a_iid    = 0;
  _offt._fileno   = 0;
  _offt._offset   = 0;
  _offt._numOlaps = 0;

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

  sprintf(name, "%s/%04d%c", _storePath, _currentFileIndex, _useBackup ? '~' : 0);
  _bof = new ovFile(name, ovFileNormal);

  _bof->seekOverlap(_offt._offset);
}



void
ovStore::resetRange(void) {
  char            name[FILENAME_MAX];

  rewind(_offtFile);

  _offt._a_iid    = 0;
  _offt._fileno   = 0;
  _offt._offset   = 0;
  _offt._numOlaps = 0;

  _overlapsThisFile = 0;
  _currentFileIndex = 1;

  delete _bof;

  sprintf(name, "%s/%04d%c", _storePath, _currentFileIndex, _useBackup ? '~' : 0);
  _bof = new ovFile(name, ovFileNormal);

  _firstIIDrequested = _info._smallestIID;
  _lastIIDrequested  = _info._largestIID;
}





void
ovStore::writeOverlap(ovsOverlap *overlap) {
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
  if (_overlapsThisFile >= _info._numOverlapsPerFile) {
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
      _offm._numOlaps  = 0;
      AS_UTL_safeWrite(_offtFile, &_offm, "ovStore::writeOverlap::offset", sizeof(ovStoreOfft), 1);
      _offm._a_iid++;
    }

    //  One more, since this iid is not missing -- we write it next!
    _offm._a_iid++;

    AS_UTL_safeWrite(_offtFile,
                     &_offt,
                     "AS_OVS_writeOverlapToStore offset",
                     sizeof(ovStoreOfft),
                     1);
    _offt._numOlaps  = 0;
  }


  //  Update the index if this is the first overlap for this a_iid
  //
  if (_offt._numOlaps == 0) {
    _offt._a_iid     = overlap->a_iid;
    _offt._fileno    = _currentFileIndex;
    _offt._offset    = _overlapsThisFile;
  }

  //AS_OVS_accumulateStats(ovs, overlap);
  _bof->writeOverlap(overlap);

  _offt._numOlaps++;
  _info._numOverlapsTotal++;
  _overlapsThisFile++;
}


// Create overlap dump index - Gregory E. Sims
// Assumes that overlaps are presorted by a_iid

void
ovStore::writeOverlap(ovsOverlap *overlap, uint32 maxOverlapsThisFile) {
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
				_offm._numOlaps  = 0;
				AS_UTL_safeWrite(_offtFile,
						&_offm,
						"AS_OVS_writeOverlapToStore offset",
						sizeof(ovStoreOfft),
						1);
				_offm._a_iid++;
			}

			//  One more, since this iid is not missing -- we write it next!
			_offm._a_iid++;
			AS_UTL_safeWrite(_offtFile,
					&_offt,
					"AS_OVS_writeOverlapToStore offset",
					sizeof(ovStoreOfft),
					1);
			_offt._numOlaps  = 0;
		}
		//  Update the index if this is the first overlap for this a_iid
		if (_offt._numOlaps == 0) {
			_offt._a_iid     = overlap[i].a_iid;
			_offt._fileno    = _currentFileIndex;
			_offt._offset    = _overlapsThisFile;
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
  offsets = (ovStoreOfft *)safe_malloc(sizeof(ovStoreOfft) * len);

  if (len != AS_UTL_safeRead(_offtFile, offsets, "AS_OVS_numOverlapsInRange", sizeof(ovStoreOfft), len)) {
    fprintf(stderr, "AS_OVS_numOverlapsInRange()-- short read on offsets!\n");
    exit(1);
  }

  for (i=0; i<len; i++)
    numolap += offsets[i]._numOlaps;

  safe_free(offsets);

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



//  For the parallel sort, write a block of sorted overlaps into a single file, with index and info.

void
writeOverlaps(char       *storePath,
              ovsOverlap *ovls,
              uint64      ovlsLen,
              uint32      fileID) {

	char                        name[FILENAME_MAX];

  uint32                      currentFileIndex = fileID;
  uint64                      overlapsThisFile = 0;

	ovStoreInfo    info;
 
	info._ovsMagic              = 1;
	info._ovsVersion            = ovStoreVersion;
  info._numOverlapsPerFile    = 1024 * 1024 * 1024;  //  arbitrary
  info._smallestIID           = UINT64_MAX;
  info._largestIID            = 0;
  info._numOverlapsTotal      = 0;
  info._highestFileIndex      = 0;
	info._maxReadLenInBits      = AS_MAX_READLEN_BITS;

	ovStoreOfft    offt; 
  ovStoreOfft    offm;

  offt._a_iid     = offm._a_iid    = ovls[0].a_iid;
	offt._fileno    = offm._fileno   = fileID;
  offt._offset    = offm._offset   = 0;
  offt._numOlaps  = offm._numOlaps = 0;

  //  Create the output file

  sprintf(name, "%s/%04d", storePath, fileID);
  ovFile *bof = new ovFile(name, ovFileFullWrite);

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
				offm._fileno   = offt._fileno;
				offm._offset   = offt._offset;
				offm._numOlaps = 0;

				AS_UTL_safeWrite(offtFile, &offm, "AS_OVS_writeOverlapToStore offt", sizeof(ovStoreOfft), 1);
				offm._a_iid++;
			}

			//  One more, since this iid is not offm -- we write it next!
			offm._a_iid++;

			AS_UTL_safeWrite(offtFile, &offt, "AS_OVS_writeOverlapToStore offt", sizeof(ovStoreOfft), 1);
			offt._numOlaps  = 0;
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

  //  Write the final index entries.

	while (offm._a_iid < offt._a_iid) {
		offm._fileno    = offt._fileno;
		offm._offset    = offt._offset;
		offm._numOlaps  = 0;

		AS_UTL_safeWrite(offtFile, &offm, "AS_OVS_writeOverlapToStore offt", sizeof(ovStoreOfft), 1);
		offm._a_iid++;
	}

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

	info._ovsMagic              = 1;
	info._ovsVersion            = ovStoreVersion;
  info._numOverlapsPerFile    = 1024 * 1024 * 1024;  //  arbitrary
  info._smallestIID           = UINT64_MAX;
  info._largestIID            = 0;
  info._numOverlapsTotal      = 0;
  info._highestFileIndex      = 0;
	info._maxReadLenInBits      = AS_MAX_READLEN_BITS;

  ovStoreOfft offm;

  offm._a_iid     = 0;
  offm._fileno    = 1;
  offm._offset    = 0;
  offm._numOlaps  = 0;

  //  Open the new master index output file

  char            name[FILENAME_MAX];

  sprintf(name, "%s/index", storePath);

  errno = 0;
  FILE  *idx = fopen(name, "w");
  if (errno)
    fprintf(stderr, "ERROR: Failed to open '%s': %s\n", name, strerror(errno)), exit(1);

  info._highestFileIndex = nPieces;

  //  Special case, we need an empty index for the zeroth fragment.

  AS_UTL_safeWrite(idx, &offm, "ovStore::mergeInfoFiles::offsetZero", sizeof(ovStoreOfft), 1);

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
      //offm._fileno    = set elsewhere
      //offm._offset    = set elsewhere
      //offm._numOlaps  = 0;

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
        offm._fileno = recs[recsLen-1]._fileno;  //  Update location of missing stuff.
        offm._offset = recs[recsLen-1]._offset;

				AS_UTL_safeWrite(idx, recs, "ovStore::mergeInfoFiles::offsetsWrite", sizeof(ovStoreOfft), recsLen);

        recsLen = AS_UTL_safeRead(F, recs, "ovStore::mergeInfoFiles::offsetsReLoad", sizeof(ovStoreOfft), recsMax);
      }

      delete [] recs;

      fclose(F);
    }

    //  Update

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


