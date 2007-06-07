#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include "existDB.H"
#include "bio++.H"


const char  magic[16] = { 'e', 'x', 'i', 's', 't', 'D', 'B', '1', 
                          ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' '  };


void
existDB::saveState(char const *filename) {
  char     cigam[16];

  errno = 0;
  FILE *F = fopen(filename, "wb");
  if (errno) {
    fprintf(stderr, "Can't open '%s' for writing\n%s\n", filename, strerror(errno));
    exit(1);
  }

  strncpy(cigam, magic, 16);

  if (_compressedHash)
    cigam[8] = 'h';
  if (_compressedBucket)
    cigam[9] = 'b';
  if (_isForward)
    cigam[10] = 'f';
  if (_isCanonical)
    cigam[10] = 'c';

  fwrite(cigam, sizeof(char), 16, F);

  fwrite(&_merSizeInBases, sizeof(u32bit), 1, F);
  fwrite(&_shift1, sizeof(u32bit), 1, F);
  fwrite(&_shift2, sizeof(u32bit), 1, F);
  fwrite(&_mask1, sizeof(u64bit), 1, F);
  fwrite(&_mask2, sizeof(u64bit), 1, F);
  fwrite(&_hshWidth, sizeof(u32bit), 1, F);  //  only valid if _compressedHash
  fwrite(&_chkWidth, sizeof(u32bit), 1, F);  //  only valid if _compressedBucket
  fwrite(&_cntWidth, sizeof(u32bit), 1, F);  //  only valid if _compressedCounts

  fwrite(&_hashTableWords, sizeof(u64bit), 1, F);
  fwrite(&_bucketsWords,   sizeof(u64bit), 1, F);
  fwrite(&_countsWords,    sizeof(u64bit), 1, F);

  fwrite(_hashTable, sizeof(u64bit), _hashTableWords, F);
  fwrite(_buckets,   sizeof(u64bit), _bucketsWords,   F);
  fwrite(_counts,    sizeof(u64bit), _countsWords,    F);

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
    //fprintf(stderr, "Can't open '%s' for reading pre-built existDB\n%s\n", strerror(errno));
    return(false);
  }

  fread(cigam, sizeof(char), 16, F);

  _compressedHash   = false;
  _compressedBucket = false;
  _isForward        = false;
  _isCanonical      = false;

  if (cigam[8] == 'h')
    _compressedHash = true;
  if (cigam[9] == 'b')
    _compressedBucket = true;
  if (cigam[10] == 'f')
    _isForward = true;
  if (cigam[10] == 'c')
    _isCanonical = true;

  cigam[ 8] = ' ';
  cigam[ 9] = ' ';
  cigam[10] = ' ';

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
  fread(&_hshWidth, sizeof(u32bit), 1, F);  //  only valid if _compressedHash
  fread(&_chkWidth, sizeof(u32bit), 1, F);  //  only valid if _compressedBucket
  fread(&_cntWidth, sizeof(u32bit), 1, F);  //  only valid if _compressedCounts

  fread(&_hashTableWords, sizeof(u64bit), 1, F);
  fread(&_bucketsWords,   sizeof(u64bit), 1, F);
  fread(&_countsWords,    sizeof(u64bit), 1, F);

  _hashTable = 0L;
  _buckets   = 0L;

  if (loadData) {
    _hashTable = new u64bit [_hashTableWords];
    _buckets   = new u64bit [_bucketsWords];
    _counts    = new u64bit [_countsWords];

    fread(_hashTable, sizeof(u64bit), _hashTableWords, F);
    fread(_buckets,   sizeof(u64bit), _bucketsWords,   F);
    fread(_counts,    sizeof(u64bit), _countsWords,    F);
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

  fprintf(stream, "merSizeInBases:   "u32bitFMT"\n", _merSizeInBases);
  fprintf(stream, "tableBits         "u32bitFMT"\n", 2 * _merSizeInBases - _shift1);
  fprintf(stream, "-----------------\n");
  fprintf(stream, "_hashTableWords   "u64bitFMT" ("u64bitFMT" KB)\n", _hashTableWords, _hashTableWords >> 7);
  fprintf(stream, "_bucketsWords     "u64bitFMT" ("u64bitFMT" KB)\n", _bucketsWords, _bucketsWords >> 7);
  fprintf(stream, "_countsWords      "u64bitFMT" ("u64bitFMT" KB)\n", _countsWords, _countsWords >> 7);
  fprintf(stream, "-----------------\n");
  fprintf(stream, "_shift1:          "u32bitFMT"\n", _shift1);
  fprintf(stream, "_shift2           "u32bitFMT"\n", _shift2);
  fprintf(stream, "_mask1            "u64bitHEX"\n", _mask1);
  fprintf(stream, "_mask2            "u64bitHEX"\n", _mask2);

  if (_compressedHash) {
    fprintf(stream, "_compressedHash   true\n");
    fprintf(stream, "_hshWidth         "u32bitFMT"\n", _hshWidth);
  } else {
    fprintf(stream, "_compressedHash   false\n");
    fprintf(stream, "_hshWidth         undefined\n");
  }

  if (_compressedBucket) {
    fprintf(stream, "_compressedBucket true\n");
    fprintf(stream, "_chkWidth         "u32bitFMT"\n", _chkWidth);
  } else {
    fprintf(stream, "_compressedBucket false\n");
    fprintf(stream, "_chkWidth         undefined\n");
  }

  if (_compressedCounts) {
    fprintf(stream, "_compressedCount  true\n");
    fprintf(stream, "_cntWidth         "u32bitFMT"\n", _cntWidth);
  } else {
    fprintf(stream, "_compressedCount  false\n");
    fprintf(stream, "_cntWidth         undefined\n");
  }
}

