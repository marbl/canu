#include "positionDB.H"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

static
char     magic[16] = { 'p', 'o', 's', 'i', 't', 'i', 'o', 'n', 'D', 'B', '1', ' ', ' ', ' ', ' ', ' '  };

//  XXX: Includes ifdef'd support for writing the magic AFTER the data
//  has been written, ensuring that the whole write is complete.  Not
//  tested much, so not enabled.
//
//#define MAGIC_AFTER_DATA

void
positionDB::saveState(char const *filename) {

  errno = 0;
  FILE *F = fopen(filename, "wb");
  if (errno) {
    fprintf(stderr, "Can't open '%s' for writing positionDB.\n%s\n", filename, strerror(errno));
    exit(1);
  }

  u64bit  hs = _tableSizeInEntries * _hashWidth / 64 + 1;
  u64bit  bs = _numberOfDistinct   * _wFin      / 64 + 1;
  u64bit  ps = _numberOfEntries    * _posnWidth / 64 + 1;

#ifndef MAGIC_AFTER_DATA
  fwrite(magic,      sizeof(char),       16, F);
#else
  //  Test if this is a pipe.  If so, we write the magic first,
  //  otherwise we write the magic last.
  //
  char  cigam[16] = { 0 };
  bool  magicFirst = false;

  lseek(fileno(F), 0, SEEK_SET);
  if (errno == ESPIPE)
    magicFirst = true;

  if (magicFirst)
    fwrite(magic,    sizeof(char),       16, F);
  else
    fwrite(cigam,    sizeof(char),       16, F);
#endif
  fwrite(this,       sizeof(positionDB), 1,  F);
  fwrite(_hashTable, sizeof(u64bit),     hs, F);
  fwrite(_buckets,   sizeof(u64bit),     bs, F);
  fwrite(_positions, sizeof(u64bit),     ps, F);

#ifdef MAGIC_AFTER_DATA
  if (!magicFirst) {
    lseek(fileno(F), 0, SEEK_SET);
    if (errno) {
      fprintf(stderr, "positionDB::saveState()-- Failed to seek to start of file -- write failed.\n%s\n", strerror(errno));
      exit(1);
    }
    fwrite(magic, sizeof(char), 16, F);
  }
#endif

  fclose(F);

  if (errno) {
    fprintf(stderr, "positionDB::saveState()-- Write failure.\n%s\n", strerror(errno));
    exit(1);
  }
}


bool
positionDB::loadState(char const *filename, bool beNoisy, bool loadData) {
  char   cigam[16] = { 0 };

  errno = 0;
  FILE *F = fopen(filename, "rb");
  if (errno) {
    fprintf(stderr, "Can't open '%s' for reading pre-built positionDB.\n%s\n", filename, strerror(errno));
    return(false);
  }

  fread(cigam, sizeof(char), 16, F);

  if (strncmp(magic, cigam, 16) != 0) {
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

    fclose(F);
    return(false);
  }

  fread(this, sizeof(positionDB), 1, F);

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

    fread(_hashTable, sizeof(u64bit), hs, F);
    fread(_buckets,   sizeof(u64bit), bs, F);
    fread(_positions, sizeof(u64bit), ps, F);
  }

  fclose(F);

  if (errno) {
    fprintf(stderr, "positionDB::loadState()-- Read failure.\n%s\n", strerror(errno));
    exit(1);
  }

  return(true);
}



void
positionDB::printState(FILE *stream) {

#ifdef TRUE64BIT
  fprintf(stream, "merSizeInBases:       %u\n",  _merSizeInBases);
  fprintf(stream, "merSkipInBases:       %u\n",  _merSkipInBases);
  fprintf(stream, "tableSizeInBits:      %u\n",  _tableSizeInBits);
  fprintf(stream, "tableSizeInEntries:   %lu\n", _tableSizeInEntries);
  fprintf(stream, "hashWidth:            %u\n",  _hashWidth);
  fprintf(stream, "chckWidth:            %u\n",  _chckWidth);
  fprintf(stream, "posnWidth:            %u\n",  _posnWidth);
  fprintf(stream, "numberOfMers:         %lu\n", _numberOfMers);
  fprintf(stream, "numberOfPositions:    %lu\n", _numberOfPositions);
  fprintf(stream, "numberOfDistinct:     %lu\n", _numberOfDistinct);
  fprintf(stream, "numberOfUnique:       %lu\n", _numberOfUnique);
  fprintf(stream, "numberOfEntries:      %lu\n", _numberOfEntries);
  fprintf(stream, "maximumEntries:       %lu\n", _maximumEntries);
#else
  fprintf(stream, "merSizeInBases:       %lu\n",  _merSizeInBases);
  fprintf(stream, "merSkipInBases:       %lu\n",  _merSkipInBases);
  fprintf(stream, "tableSizeInBits:      %lu\n",  _tableSizeInBits);
  fprintf(stream, "tableSizeInEntries:   %llu\n", _tableSizeInEntries);
  fprintf(stream, "hashWidth:            %lu\n",  _hashWidth);
  fprintf(stream, "chckWidth:            %lu\n",  _chckWidth);
  fprintf(stream, "posnWidth:            %lu\n",  _posnWidth);
  fprintf(stream, "numberOfMers:         %llu\n", _numberOfMers);
  fprintf(stream, "numberOfPositions:    %llu\n", _numberOfPositions);
  fprintf(stream, "numberOfDistinct:     %llu\n", _numberOfDistinct);
  fprintf(stream, "numberOfUnique:       %llu\n", _numberOfUnique);
  fprintf(stream, "numberOfEntries:      %llu\n", _numberOfEntries);
  fprintf(stream, "maximumEntries:       %llu\n", _maximumEntries);
#endif


}

