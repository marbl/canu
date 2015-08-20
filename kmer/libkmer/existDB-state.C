
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
 *    Brian P. Walenz on 2003-SEP-08
 *      are Copyright 2003 Applera Corporation, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2004-APR-12 to 2004-OCT-10
 *      are Copyright 2004 Brian P. Walenz, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2005-MAR-20 to 2014-APR-11
 *      are Copyright 2005,2007,2013-2014 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include "existDB.H"
#include "bio++.H"


const char  magic[16] = { 'e', 'x', 'i', 's', 't', 'D', 'B', '2',
                          ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' '  };


void
existDB::saveState(char const *filename) {
  char     cigam[16] = { 0 };

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
  if (_compressedCounts)
    cigam[10] = 'c';

  if (_isForward)
    cigam[11] = 'F';
  if (_isCanonical)
    cigam[11] = 'C';

  fwrite(cigam, sizeof(char), 16, F);

  fwrite(&_merSizeInBases, sizeof(uint32), 1, F);
  fwrite(&_shift1, sizeof(uint32), 1, F);
  fwrite(&_shift2, sizeof(uint32), 1, F);
  fwrite(&_mask1, sizeof(uint64), 1, F);
  fwrite(&_mask2, sizeof(uint64), 1, F);
  fwrite(&_hshWidth, sizeof(uint32), 1, F);  //  only valid if _compressedHash
  fwrite(&_chkWidth, sizeof(uint32), 1, F);  //  only valid if _compressedBucket
  fwrite(&_cntWidth, sizeof(uint32), 1, F);  //  only valid if _compressedCounts

  fwrite(&_hashTableWords, sizeof(uint64), 1, F);
  fwrite(&_bucketsWords,   sizeof(uint64), 1, F);
  fwrite(&_countsWords,    sizeof(uint64), 1, F);

  fwrite(_hashTable, sizeof(uint64), _hashTableWords, F);
  fwrite(_buckets,   sizeof(uint64), _bucketsWords,   F);
  fwrite(_counts,    sizeof(uint64), _countsWords,    F);

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
  _compressedCounts = false;
  _isForward        = false;
  _isCanonical      = false;

  if (cigam[8] == 'h')
    _compressedHash = true;
  if (cigam[9] == 'b')
    _compressedBucket = true;
  if (cigam[10] == 'c')
    _compressedCounts = true;

  if (cigam[11] == 'F')
    _isForward = true;
  if (cigam[11] == 'C')
    _isCanonical = true;

  cigam[ 8] = ' ';
  cigam[ 9] = ' ';
  cigam[10] = ' ';
  cigam[11] = ' ';

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

  fread(&_merSizeInBases, sizeof(uint32), 1, F);
  fread(&_shift1, sizeof(uint32), 1, F);
  fread(&_shift2, sizeof(uint32), 1, F);
  fread(&_mask1, sizeof(uint64), 1, F);
  fread(&_mask2, sizeof(uint64), 1, F);
  fread(&_hshWidth, sizeof(uint32), 1, F);  //  only valid if _compressedHash
  fread(&_chkWidth, sizeof(uint32), 1, F);  //  only valid if _compressedBucket
  fread(&_cntWidth, sizeof(uint32), 1, F);  //  only valid if _compressedCounts

  fread(&_hashTableWords, sizeof(uint64), 1, F);
  fread(&_bucketsWords,   sizeof(uint64), 1, F);
  fread(&_countsWords,    sizeof(uint64), 1, F);

  _hashTable = 0L;
  _buckets   = 0L;
  _counts    = 0L;

  if (loadData) {
    _hashTable = new uint64 [_hashTableWords];
    _buckets   = new uint64 [_bucketsWords];

    if (_countsWords > 0)
      _counts  = new uint64 [_countsWords];

    fread(_hashTable, sizeof(uint64), _hashTableWords, F);
    fread(_buckets,   sizeof(uint64), _bucketsWords,   F);

    if (_countsWords > 0)
      fread(_counts,  sizeof(uint64), _countsWords,    F);
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

  fprintf(stream, "merSizeInBases:   "uint32FMT"\n", _merSizeInBases);
  fprintf(stream, "tableBits         "uint32FMT"\n", 2 * _merSizeInBases - _shift1);
  fprintf(stream, "-----------------\n");
  fprintf(stream, "_hashTableWords   "uint64FMT" ("uint64FMT" KB)\n", _hashTableWords, _hashTableWords >> 7);
  fprintf(stream, "_bucketsWords     "uint64FMT" ("uint64FMT" KB)\n", _bucketsWords, _bucketsWords >> 7);
  fprintf(stream, "_countsWords      "uint64FMT" ("uint64FMT" KB)\n", _countsWords, _countsWords >> 7);
  fprintf(stream, "-----------------\n");
  fprintf(stream, "_shift1:          "uint32FMT"\n", _shift1);
  fprintf(stream, "_shift2           "uint32FMT"\n", _shift2);
  fprintf(stream, "_mask1            "uint64HEX"\n", _mask1);
  fprintf(stream, "_mask2            "uint64HEX"\n", _mask2);

  if (_compressedHash) {
    fprintf(stream, "_compressedHash   true\n");
    fprintf(stream, "_hshWidth         "uint32FMT"\n", _hshWidth);
  } else {
    fprintf(stream, "_compressedHash   false\n");
    fprintf(stream, "_hshWidth         undefined\n");
  }

  if (_compressedBucket) {
    fprintf(stream, "_compressedBucket true\n");
    fprintf(stream, "_chkWidth         "uint32FMT"\n", _chkWidth);
  } else {
    fprintf(stream, "_compressedBucket false\n");
    fprintf(stream, "_chkWidth         undefined\n");
  }

  if (_compressedCounts) {
    fprintf(stream, "_compressedCount  true\n");
    fprintf(stream, "_cntWidth         "uint32FMT"\n", _cntWidth);
  } else {
    fprintf(stream, "_compressedCount  false\n");
    fprintf(stream, "_cntWidth         undefined\n");
  }
}

