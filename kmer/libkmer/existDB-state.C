#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include "existDB.H"
#include "positionDB.H"
#include "merstream.H"


char     magic[16] = { 'e', 'x', 'i', 's', 't', 'D', 'B', '1', 
                       ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' '  };


void
existDB::saveState(char const *filename) {

  errno = 0;
  FILE *F = fopen(filename, "wb");
  if (errno) {
    fprintf(stderr, "Can't open '%s' for writing\n%s\n", strerror(errno));
    exit(1);
  }

#ifdef COMPRESSED_HASH
  magic[8] = 'h';
#endif

#ifdef COMPRESSED_BUCKET
  magic[9] = 'b';
#endif

  fwrite(magic, sizeof(char), 16, F);

  fwrite(&_merSizeInBases, sizeof(u32bit), 1, F);
  fwrite(&_shift1, sizeof(u32bit), 1, F);
  fwrite(&_shift2, sizeof(u32bit), 1, F);
  fwrite(&_mask1, sizeof(u64bit), 1, F);
  fwrite(&_mask2, sizeof(u64bit), 1, F);

#ifdef COMPRESSED_HASH
  fwrite(&_hashWidth, sizeof(u32bit), 1, F);
#endif

#ifdef COMPRESSED_BUCKET
  fwrite(&_chckWidth, sizeof(u32bit), 1, F);
#endif

  fwrite(&_hashMask, sizeof(u64bit), 1, F);
  fwrite(&_chckMask, sizeof(u64bit), 1, F);

  fwrite(&_hashTableWords, sizeof(u64bit), 1, F);
  fwrite(&_bucketsWords,   sizeof(u64bit), 1, F);

  fwrite(_hashTable, sizeof(u64bit), _hashTableWords, F);
  fwrite(_buckets,   sizeof(u64bit), _bucketsWords,   F);

  fclose(F);

  if (errno) {
    fprintf(stderr, "existDB::saveState()-- Write failure.\n%s\n", strerror(errno));
    exit(1);
  }
}



bool
existDB::loadState(char const *filename,
                   bool        beNoisy,
                   bool        loadData) {
  char     cigam[16];

  errno = 0;
  FILE *F = fopen(filename, "rb");
  if (errno) {
    fprintf(stderr, "Can't open '%s' for reading pre-built existDB\n%s\n", strerror(errno));
    return(false);
  }

#ifdef COMPRESSED_HASH
  magic[8] = 'h';
#endif

#ifdef COMPRESSED_BUCKET
  magic[9] = 'b';
#endif

  fread(cigam, sizeof(char), 16, F);

  if (strncmp(magic, cigam, 16) != 0) {
    if (beNoisy) {
      fprintf(stderr, "existDB::loadState()-- Not an existDB binary file, maybe a sequence file?\n");
      fprintf(stderr, "existDB::loadState()-- Read     '%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c'\n",
              cigam[0],  cigam[1],  cigam[2],  cigam[3],
              cigam[4],  cigam[5],  cigam[6],  cigam[7],
              cigam[8],  cigam[9],  cigam[10], cigam[11],
              cigam[12], cigam[13], cigam[14], cigam[15]);
      fprintf(stderr, "existDB::loadState()-- Expected '%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c'\n",
              magic[0],  magic[1],  magic[2],  magic[3],
              magic[4],  magic[5],  magic[6],  magic[7],
              magic[8],  magic[9],  magic[10], magic[11],
              magic[12], magic[13], magic[14], magic[15]);
    }

    fclose(F);
    return(false);
  }

  fread(&_merSizeInBases, sizeof(u32bit), 1, F);
  fread(&_shift1, sizeof(u32bit), 1, F);
  fread(&_shift2, sizeof(u32bit), 1, F);
  fread(&_mask1, sizeof(u64bit), 1, F);
  fread(&_mask2, sizeof(u64bit), 1, F);

#ifdef COMPRESSED_HASH
  fread(&_hashWidth, sizeof(u32bit), 1, F);
#endif

#ifdef COMPRESSED_BUCKET
  fread(&_chckWidth, sizeof(u32bit), 1, F);
#endif

  fread(&_hashMask, sizeof(u64bit), 1, F);
  fread(&_chckMask, sizeof(u64bit), 1, F);

  fread(&_hashTableWords, sizeof(u64bit), 1, F);
  fread(&_bucketsWords,   sizeof(u64bit), 1, F);

  _hashTable = 0L;
  _buckets   = 0L;

  if (loadData) {
    _hashTable = new u64bit [_hashTableWords];
    _buckets   = new u64bit [_bucketsWords];

    fread(_hashTable, sizeof(u64bit), _hashTableWords, F);
    fread(_buckets,   sizeof(u64bit), _bucketsWords,   F);
  }

  fclose(F);

  if (errno) {
    fprintf(stderr, "existDB::loadState()-- Read failure.\n%s\n", strerror(errno));
    exit(1);
  }

  return(true);
}


void
existDB::printState(FILE *stream) {

  fprintf(stream, "merSizeInBases:   %u\n", _merSizeInBases);
  fprintf(stream, "tableBits         %u\n", 2 * _merSizeInBases - _shift1);
  fprintf(stream, "-----------------\n");
  fprintf(stream, "_hashTableWords   %lu (%lu KB)\n", _hashTableWords, _hashTableWords >> 7);
  fprintf(stream, "_bucketsWords     %lu (%lu KB)\n", _bucketsWords, _bucketsWords >> 7);
  fprintf(stream, "-----------------\n");
  fprintf(stream, "_shift1:          %u\n", _shift1);
  fprintf(stream, "_shift2           %u\n", _shift2);
  fprintf(stream, "_mask1            0x%016lx\n", _mask1);
  fprintf(stream, "_mask2            0x%016lx\n", _mask2);

#ifdef COMPRESSED_HASH
  fprintf(stream, "COMPRESSED_HASH   enabled\n");
  fprintf(stream, "_hashWidth        %u\n", _hashWidth);
#endif

#ifdef COMPRESSED_BUCKET
  fprintf(stream, "COMPRESSED_BUCKET enabled\n");
  fprintf(stream, "_chckWidth        %u\n", _chckWidth);
#endif

  fprintf(stream, "_hashMask         0x%016lx\n", _hashMask);
  fprintf(stream, "_chckMask         0x%016lx\n", _chckMask);
}

