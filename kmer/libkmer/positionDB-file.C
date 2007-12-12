#include "positionDB.H"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>

static
char     magic[16] = { 'p', 'o', 's', 'i', 't', 'i', 'o', 'n', 'D', 'B', '.', 'v', '1', ' ', ' ', ' '  };
static
char     faild[16] = { 'p', 'o', 's', 'i', 't', 'i', 'o', 'n', 'D', 'B', 'f', 'a', 'i', 'l', 'e', 'd'  };

void
positionDB::saveState(char const *filename) {

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
  u32bit     *bs = _bucketSizes;
  u64bit     *cb = _countingBuckets;
  u64bit     *ht = _hashTable;
  u64bit     *bu = _buckets;
  u64bit     *ps = _positions;

  _bucketSizes     = 0L;
  _countingBuckets = 0L;
  _hashTable       = 0L;
  _buckets         = 0L;
  _positions       = 0L;

  safeWrite(F, this,       "this",       sizeof(positionDB) * 1);

  _bucketSizes     = bs;
  _countingBuckets = cb;
  _hashTable       = ht;
  _buckets         = bu;
  _positions       = ps;

#if 0
  fprintf(stderr, "hash: %d\n", sizeof(u64bit) * (_tableSizeInEntries * _hashWidth / 64 + 1));
  fprintf(stderr, "buck: %d\n", sizeof(u64bit) * (_numberOfDistinct   * _wFin      / 64 + 1));
  fprintf(stderr, "posn: %d\n", sizeof(u64bit) * (_numberOfEntries    * _posnWidth / 64 + 1));
#endif

  safeWrite(F, _hashTable, "_hashTable", sizeof(u64bit) * (_tableSizeInEntries * _hashWidth / 64 + 1));
  safeWrite(F, _buckets,   "_buckets",   sizeof(u64bit) * (_numberOfDistinct   * _wFin      / 64 + 1));
  safeWrite(F, _positions, "_positions", sizeof(u64bit) * (_numberOfEntries    * _posnWidth / 64 + 1));

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

  _hashTable = 0L;
  _buckets   = 0L;
  _positions = 0L;

  if (loadData) {
    u64bit  hs = _tableSizeInEntries * _hashWidth / 64 + 1;
    u64bit  bs = _numberOfDistinct   * _wFin      / 64 + 1;
    u64bit  ps = _numberOfEntries    * _posnWidth / 64 + 1;

    _hashTable = new u64bit [hs];
    _buckets   = new u64bit [bs];
    _positions = new u64bit [ps];

    safeRead(F, _hashTable, "_hashTable", sizeof(u64bit) * hs);
    safeRead(F, _buckets,   "_buckets",   sizeof(u64bit) * bs);
    safeRead(F, _positions, "_positions", sizeof(u64bit) * ps);
  }

  close(F);

  return(true);
}



void
positionDB::printState(FILE *stream) {
  fprintf(stream, "merSizeInBases:       "u32bitFMT"\n", _merSizeInBases);
  fprintf(stream, "merSkipInBases:       "u32bitFMT"\n", _merSkipInBases);
  fprintf(stream, "tableSizeInBits:      "u32bitFMT"\n", _tableSizeInBits);
  fprintf(stream, "tableSizeInEntries:   "u64bitFMT"\n", _tableSizeInEntries);
  fprintf(stream, "hashWidth:            "u32bitFMT"\n", _hashWidth);
  fprintf(stream, "chckWidth:            "u32bitFMT"\n", _chckWidth);
  fprintf(stream, "posnWidth:            "u32bitFMT"\n", _posnWidth);
  fprintf(stream, "numberOfMers:         "u64bitFMT"\n", _numberOfMers);
  fprintf(stream, "numberOfPositions:    "u64bitFMT"\n", _numberOfPositions);
  fprintf(stream, "numberOfDistinct:     "u64bitFMT"\n", _numberOfDistinct);
  fprintf(stream, "numberOfUnique:       "u64bitFMT"\n", _numberOfUnique);
  fprintf(stream, "numberOfEntries:      "u64bitFMT"\n", _numberOfEntries);
  fprintf(stream, "maximumEntries:       "u64bitFMT"\n", _maximumEntries);
}

