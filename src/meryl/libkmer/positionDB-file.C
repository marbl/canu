
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
 *    Brian P. Walenz from 2003-AUG-14 to 2004-APR-01
 *      are Copyright 2003-2004 Applera Corporation, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz on 2004-APR-30
 *      are Copyright 2004 Brian P. Walenz, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2005-MAR-12 to 2014-APR-11
 *      are Copyright 2005,2007-2008,2014 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "positionDB.H"

static
char     magic[16] = { 'p', 'o', 's', 'i', 't', 'i', 'o', 'n', 'D', 'B', '.', 'v', '1', ' ', ' ', ' '  };
static
char     faild[16] = { 'p', 'o', 's', 'i', 't', 'i', 'o', 'n', 'D', 'B', 'f', 'a', 'i', 'l', 'e', 'd'  };





//  These came from kmer/libutil/file.c and need to be replaced with more standard functions


//  Split writes/reads into smaller pieces, check the result of each
//  piece.  Really needed by OSF1 (V5.1).
//
void
safeWrite(int filedes, const void *buffer, const char *desc, size_t nbytes) {
  size_t  position = 0;
  size_t  length   = 32 * 1024 * 1024;
  size_t  towrite  = 0;
  size_t  written  = 0;

  while (position < nbytes) {
    towrite = length;
    if (position + towrite > nbytes)
      towrite = nbytes - position;

    errno = 0;
    written = write(filedes, ((char *)buffer) + position, towrite);

    if ((errno) || (towrite != written)) {
      fprintf(stderr, "safeWrite()-- Write failure on %s: %s\n", desc, strerror(errno));
      fprintf(stderr, "safeWrite()-- Wanted to write "F_S64" bytes, wrote "F_S64".\n", (int64)towrite, (int64)written);
      exit(1);
    }

    position += written;
  }
}

int
safeRead(int filedes, const void *buffer, const char *desc, size_t nbytes) {
  size_t  position = 0;
  size_t  length   = 32 * 1024 * 1024;
  size_t  toread   = 0;
  size_t  written  = 0;  //  readen?
  int     failed   = 0;

  while (position < nbytes) {
    toread = length;
    if (position + toread > nbytes)
      toread = nbytes - position;

    errno = 0;
    written = read(filedes, ((char *)buffer) + position, toread);

    failed = errno;
#ifdef VERY_SAFE
    if (toread != written)
      failed = 1;
#endif

    if ((failed) && (errno != EINTR)) {
      fprintf(stderr, "safeRead()-- Read failure on %s: %s.\n", desc, strerror(errno));
      fprintf(stderr, "safeRead()-- Wanted to read "F_S64" bytes, read "F_S64".\n", (int64)toread, (int64)written);
      exit(1);
    }

    if (written == 0)
      break;

    position += written;
  }

  return(position);
}








void
positionDB::saveState(char const *filename) {

  fprintf(stderr, "Saving positionDB to '%s'\n", filename);

  errno = 0;
  int F = open(filename, O_RDWR | O_CREAT | O_LARGEFILE,
               S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH | S_IWOTH);
  if (errno) {
    fprintf(stderr, "Can't open '%s' for writing positionDB.\n%s\n", filename, strerror(errno));
    exit(1);
  }

  bool    magicFirst = false;

  //  Test if this is a pipe.  If so, we write the magic first,
  //  otherwise we write the magic last.
  //
  errno      = 0;
  lseek(F, 0, SEEK_SET);
  if (errno == ESPIPE)
    magicFirst = true;

  if (magicFirst)
    write(F, magic, sizeof(char) * 16);
  else
    write(F, faild, sizeof(char) * 16);

  if (errno) {
    fprintf(stderr, "positionDB::saveState()-- Write failure on magic first.\n%s\n", strerror(errno));
    exit(1);
  }

  //  If only to be completely annoying and anal, we clear the
  //  pointers before we write the data.  Sure, we could just write
  //  the stuff we care about, but this is easier.  This is easier.
  //  Before you go rip out this stuff, remember that you can now
  //  checksum the resulting files.  So don't do it.
  //
  uint32     *bs = _bucketSizes;
  uint64     *cb = _countingBuckets;
  uint64     *hp = _hashTable_BP;
  uint32     *hw = _hashTable_FW;
  uint64     *bu = _buckets;
  uint64     *ps = _positions;
  uint64     *he = _hashedErrors;

  _bucketSizes     = 0L;
  _countingBuckets = 0L;
  _hashTable_BP    = (uint64 *)((_hashTable_BP) ? uint64ONE : uint64ZERO);
  _hashTable_FW    = (uint32 *)((_hashTable_FW) ? uint64ONE : uint64ZERO);
  _buckets         = 0L;
  _positions       = 0L;
  _hashedErrors    = 0L;

  safeWrite(F, this,       "this",       sizeof(positionDB) * 1);

  _bucketSizes     = bs;
  _countingBuckets = cb;
  _hashTable_BP    = hp;
  _hashTable_FW    = hw;
  _buckets         = bu;
  _positions       = ps;
  _hashedErrors    = he;

  if (_hashTable_BP) {
    safeWrite(F, _hashTable_BP, "_hashTable_BP", sizeof(uint64) * (_tableSizeInEntries * _hashWidth / 64 + 1));
  } else {
    safeWrite(F, _hashTable_FW, "_hashTable_FW", sizeof(uint32) * (_tableSizeInEntries + 1));
  }

  safeWrite(F, _buckets,      "_buckets",      sizeof(uint64) * (_numberOfDistinct   * _wFin      / 64 + 1));
  safeWrite(F, _positions,    "_positions",    sizeof(uint64) * (_numberOfEntries    * _posnWidth / 64 + 1));
  safeWrite(F, _hashedErrors, "_hashedErrors", sizeof(uint64) * (_hashedErrorsLen));

  if (magicFirst == false) {
    lseek(F, 0, SEEK_SET);
    if (errno) {
      fprintf(stderr, "positionDB::saveState()-- Failed to seek to start of file -- write failed.\n%s\n", strerror(errno));
      exit(1);
    }

    write(F, magic, sizeof(char) * 16);
    if (errno) {
      fprintf(stderr, "positionDB::saveState()-- Write failure on magic last.\n%s\n", strerror(errno));
      exit(1);
    }
  }

  close(F);
}


bool
positionDB::loadState(char const *filename, bool beNoisy, bool loadData) {
  char   cigam[16] = { 0 };

  fprintf(stderr, "Loading positionDB from '%s'\n", filename);

  errno = 0;
  int F = open(filename, O_RDONLY | O_LARGEFILE, 0);
  if (errno) {
    fprintf(stderr, "Can't open '%s' for reading pre-built positionDB: %s\n", filename, strerror(errno));
    return(false);
  }

  safeRead(F, cigam, "Magic Number", sizeof(char) * 16);

  if        (strncmp(faild, cigam, 16) == 0) {
    if (beNoisy) {
      fprintf(stderr, "positionDB::loadState()-- Incomplete positionDB binary file.\n");
      fprintf(stderr, "positionDB::loadState()-- Read     '%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c'\n",
              cigam[0],  cigam[1],  cigam[2],  cigam[3],
              cigam[4],  cigam[5],  cigam[6],  cigam[7],
              cigam[8],  cigam[9],  cigam[10], cigam[11],
              cigam[12], cigam[13], cigam[14], cigam[15]);
      fprintf(stderr, "positionDB::loadState()-- Expected '%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c'\n",
              magic[0],  magic[1],  magic[2],  magic[3],
              magic[4],  magic[5],  magic[6],  magic[7],
              magic[8],  magic[9],  magic[10], magic[11],
              magic[12], magic[13], magic[14], magic[15]);
    }
    close(F);
    return(false);
  } else if (strncmp(magic, cigam, 16) != 0) {
    if (beNoisy) {
      fprintf(stderr, "positionDB::loadState()-- Not a positionDB binary file, maybe a sequence file?\n");
      fprintf(stderr, "positionDB::loadState()-- Read     '%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c'\n",
              cigam[0],  cigam[1],  cigam[2],  cigam[3],
              cigam[4],  cigam[5],  cigam[6],  cigam[7],
              cigam[8],  cigam[9],  cigam[10], cigam[11],
              cigam[12], cigam[13], cigam[14], cigam[15]);
      fprintf(stderr, "positionDB::loadState()-- Expected '%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c'\n",
              magic[0],  magic[1],  magic[2],  magic[3],
              magic[4],  magic[5],  magic[6],  magic[7],
              magic[8],  magic[9],  magic[10], magic[11],
              magic[12], magic[13], magic[14], magic[15]);
    }

    close(F);
    return(false);
  }

  safeRead(F, this, "positionDB", sizeof(positionDB) * 1);

  _bucketSizes     = 0L;
  _countingBuckets = 0L;
  _buckets         = 0L;
  _positions       = 0L;
  _hashedErrors    = 0L;

  if (loadData) {
    uint64  hs = _tableSizeInEntries * _hashWidth / 64 + 1;
    uint64  bs = _numberOfDistinct   * _wFin      / 64 + 1;
    uint64  ps = _numberOfEntries    * _posnWidth / 64 + 1;

    if (_hashTable_BP) {
      _hashTable_BP = new uint64 [hs];
      _hashTable_FW = 0L;
      safeRead(F, _hashTable_BP, "_hashTable_BP", sizeof(uint64) * hs);
    } else {
      _hashTable_BP = 0L;
      _hashTable_FW = new uint32 [_tableSizeInEntries + 1];
      safeRead(F, _hashTable_FW, "_hashTable_FW", sizeof(uint32) * (_tableSizeInEntries + 1));
    }

    _buckets      = new uint64 [bs];
    _positions    = new uint64 [ps];
    _hashedErrors = new uint64 [_hashedErrorsMax];

    safeRead(F, _buckets,      "_buckets",      sizeof(uint64) * bs);
    safeRead(F, _positions,    "_positions",    sizeof(uint64) * ps);
    safeRead(F, _hashedErrors, "_hashedErrors", sizeof(uint64) * _hashedErrorsLen);
  }

  close(F);

  return(true);
}



void
positionDB::printState(FILE *stream) {
  fprintf(stream, "merSizeInBases:       "F_U32"\n", _merSizeInBases);
  fprintf(stream, "merSkipInBases:       "F_U32"\n", _merSkipInBases);
  fprintf(stream, "tableSizeInBits:      "F_U32"\n", _tableSizeInBits);
  fprintf(stream, "tableSizeInEntries:   "F_U64"\n", _tableSizeInEntries);
  fprintf(stream, "hashWidth:            "F_U32"\n", _hashWidth);
  fprintf(stream, "chckWidth:            "F_U32"\n", _chckWidth);
  fprintf(stream, "posnWidth:            "F_U32"\n", _posnWidth);
  fprintf(stream, "numberOfMers:         "F_U64"\n", _numberOfMers);
  fprintf(stream, "numberOfPositions:    "F_U64"\n", _numberOfPositions);
  fprintf(stream, "numberOfDistinct:     "F_U64"\n", _numberOfDistinct);
  fprintf(stream, "numberOfUnique:       "F_U64"\n", _numberOfUnique);
  fprintf(stream, "numberOfEntries:      "F_U64"\n", _numberOfEntries);
  fprintf(stream, "maximumEntries:       "F_U64"\n", _maximumEntries);
}

