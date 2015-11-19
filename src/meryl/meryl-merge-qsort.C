
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
 *    kmer/meryl/merge.qsort.C
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2008-JUN-20 to 2014-APR-11
 *      are Copyright 2008,2014 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz on 2014-DEC-05
 *      are Copyright 2014 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "meryl.H"
#include "libmeryl.H"


using namespace std;

#include <algorithm>

struct mMer {
  kMer     _mer;
  uint32   _cnt;
  uint32   _off;
};

static
int
mMerGreaterThan(void const *a, void const *b) {
  mMer const *A = (mMer const *)a;
  mMer const *B = (mMer const *)b;
  return(B->_mer.qsort_less(A->_mer));
}



class mMerList {
public:
  mMerList(uint32 maxSize) {
    _posLen = 0;
    _posMax = 2 * maxSize;
    _pos    = new uint32 [_posMax];

    _mmmLen = 0;
    _mmmMax = maxSize;
    _mmm    = new mMer [_mmmMax];
  };
  ~mMerList() {
    delete [] _pos;
    delete [] _mmm;
  };

  uint32  length(void) { return(_mmmLen);              };

  //  Until we sort, first() is the last thing loaded.
  //  After we sort, first() is the lowest mer in the set.

  kMer   &first(void)   { return(_mmm[_mmmLen-1]._mer); };
  //kMer   &last(void)    { return(_mmm[0]._mer);         };
  //kMer   &get(uint32 i) { return(_mmm[i]._mer); };

  //  Return the first (sorted order) thing in the list -- it's the last on the list.
  kMer   *pop(uint32 &cnt, uint32* &pos) {
    if (_mmmLen == 0)
      return(0L);

    _mmmLen--;

    assert(_sorted);

    cnt = _mmm[_mmmLen]._cnt;
    pos = 0L;

    if (_mmm[_mmmLen]._off != ~uint32ZERO)
      pos = _pos + _mmm[_mmmLen]._off;

    return(&_mmm[_mmmLen]._mer);
  }


  //  rebuild the position list, squeezes out empty items
  void    rebuild(void) {
    if (_posLen > 0) {
      uint32   *np = new uint32 [_posMax];

      _posLen = 0;

      for (uint32 i=0; i<_mmmLen; i++) {
        mMer *m = _mmm + i;

        if (m->_off != ~uint32ZERO) {
          _mmm[_mmmLen]._off = _posLen;

          for (uint32 p=0; p<m->_cnt; p++, _posLen++)
            np[_posLen] = _pos[p];
        }
      }

      delete [] _pos;
      _pos = np;
    }
  };



  //  Read more mers from the file
  void    read(merylStreamReader *R, uint32 num) {
    uint32 xxx = 0;

    if (_mmmLen + num >= _mmmMax) {
      fprintf(stderr, "Reallocate _mmm\n");
      _mmmMax = _mmmMax + 2 * num;
      mMer *tmp = new mMer [_mmmMax];
      memcpy(tmp, _mmm, sizeof(mMer) * _mmmLen);
      delete [] _mmm;
      _mmm = tmp;
    }

    _sorted = false;

    R->nextMer();

    for (xxx=0; (xxx < num) && (R->validMer()); xxx++) {
      if (_mmmMax <= _mmmLen) {
        fprintf(stderr, "Reallocate _mmm\n");
        _mmmMax *= 2;
        mMer *tmp = new mMer [_mmmMax];
        memcpy(tmp, _mmm, sizeof(mMer) * _mmmLen);
        delete [] _mmm;
        _mmm = tmp;
      }

      _mmm[_mmmLen]._mer = R->theFMer();
      _mmm[_mmmLen]._cnt = R->theCount();
      _mmm[_mmmLen]._off = ~uint32ZERO;

      uint32  *pos = R->thePositions();
      if (pos) {
        _mmm[_mmmLen]._off = _posLen;

        if (_posMax <= _posLen + _mmm[_mmmLen]._cnt) {
          fprintf(stderr, "Reallocate _pos\n");
          _posMax *= 2;
          uint32 *tmp = new uint32 [_posMax];
          memcpy(tmp, _pos, sizeof(uint32) * _posLen);
          delete [] _pos;
          _pos = tmp;
        }

        for (uint32 i=0; i<_mmm[_mmmLen]._cnt; i++, _posLen++)
          _pos[_posLen] = pos[i];
      }

      _mmmLen++;

      R->nextMer();
    }

    //fprintf(stderr, "read()-- now up to "uint32FMT" mers ("uint32FMT" pos); loaded "uint32FMT" out of "uint32FMT" requested.\n", _mmmLen, _posLen, xxx, num);
  };


  //  Sort our list of mers
  void    sort(void) {
    if (_sorted == false) {
      //fprintf(stderr, "SORT BEG\n");
      qsort_mt(_mmm, _mmmLen, sizeof(mMer), mMerGreaterThan, 8, 32 * 1024);
      _sorted = true;
      //fprintf(stderr, "SORT END\n");
    }
  };


private:
  bool    _sorted;

  uint32  _posLen;
  uint32  _posMax;
  uint32 *_pos;

  uint32  _mmmLen;
  uint32  _mmmMax;
  mMer   *_mmm;
};




void
multipleOperations(merylArgs *args) {

  char  debugstring[256];
  char  debugstring2[256];

  if (args->mergeFilesLen < 2) {
    fprintf(stderr, "ERROR - must have at least two databases (you gave "uint32FMT")!\n", args->mergeFilesLen);
    exit(1);
  }
  if (args->outputFile == 0L) {
    fprintf(stderr, "ERROR - no output file specified.\n");
    exit(1);
  }
  if ((args->personality != PERSONALITY_MERGE) &&
      (args->personality != PERSONALITY_MIN) &&
      (args->personality != PERSONALITY_MINEXIST) &&
      (args->personality != PERSONALITY_MAX) &&
      (args->personality != PERSONALITY_ADD) &&
      (args->personality != PERSONALITY_AND) &&
      (args->personality != PERSONALITY_NAND) &&
      (args->personality != PERSONALITY_OR) &&
      (args->personality != PERSONALITY_XOR)) {
    fprintf(stderr, "ERROR - only personalities min, minexist, max, add, and, nand, or, xor\n");
    fprintf(stderr, "ERROR - are supported in multipleOperations().\n");
    fprintf(stderr, "ERROR - this is a coding error, not a user error.\n");
    exit(1);
  }

  merylStreamReader  **R = new merylStreamReader* [args->mergeFilesLen];
  merylStreamWriter   *W = 0L;

  uint32               maxSize = 512 * 1024;

  mMerList            *M = new mMerList(maxSize + maxSize / 4);


  //  Open the input files and load some mers - we need to do this
  //  just so we can check the mersizes/compression next.
  //
  for (uint32 i=0; i<args->mergeFilesLen; i++) {
    R[i] = new merylStreamReader(args->mergeFiles[i]);
    M->read(R[i], 1 + i);
  }

  //  Verify that the mersizes are all the same
  //
  bool    fail       = false;
  uint32  merSize    = R[0]->merSize();
  uint32  merComp    = R[0]->merCompression();

  for (uint32 i=0; i<args->mergeFilesLen; i++) {
    fail |= (merSize != R[i]->merSize());
    fail |= (merComp != R[i]->merCompression());
  }

  if (fail)
    fprintf(stderr, "ERROR:  mer sizes (or compression level) differ.\n"), exit(1);

  //  Open the output file, using the largest prefix size found in the
  //  input/mask files.
  //
  uint32  prefixSize = 0;
  for (uint32 i=0; i<args->mergeFilesLen; i++)
    if (prefixSize < R[i]->prefixSize())
      prefixSize = R[i]->prefixSize();

  W = new merylStreamWriter(args->outputFile, merSize, merComp, prefixSize);


  kMer     lastLoaded;

  lastLoaded.setMerSize(merSize);
  lastLoaded.smallest();

  //  Load mers from all files, remember the largest mer we load.
  //
  for (uint32 i=0; i<args->mergeFilesLen; i++) {
    M->read(R[i], maxSize / args->mergeFilesLen);
    if (lastLoaded < M->first())
      lastLoaded = M->first();
  }

  //  Make sure all files have at least that largest mer loaded.
  //
  for (uint32 i=0; i<args->mergeFilesLen; i++)
    while (R[i]->validMer() && (R[i]->theFMer() <= lastLoaded))
      M->read(R[i], 2 * 1024);

  fprintf(stderr, "Initial load:  length="uint32FMT" lastLoaded=%s\n",
          M->length(), lastLoaded.merToString(debugstring));

  M->sort();

  bool     allLoaded = false;
  bool     moreStuff = true;

  kMer     currentMer;                      //  The current mer we're operating on
  uint32   currentCount     =  uint32ZERO;  //  The count (operation dependent) of this mer
  uint32   currentTimes     =  uint32ZERO;  //  Number of files it's in

  uint32   currentPositionsMax =  0;
  uint32  *currentPositions    =  0L;

  kMer    *thisMer;                         //  The mer we just read
  uint32   thisCount        =  uint32ZERO;  //  The count of the mer we just read
  uint32  *thisPositions    = 0L;

  speedCounter *C = new speedCounter("    %7.2f Mmers -- %5.2f Mmers/second\r", 1000000.0, 0x1fffff, args->beVerbose);

  currentMer.setMerSize(merSize);

  while (moreStuff) {

    //  Load more stuff if needed.  M is sorted, so first() is the
    //  smallest mer in the set - we're good up to and including
    //  lastLoaded.
    //
    if ((allLoaded == false) &&
        ((M->length() == 0) || (lastLoaded < M->first()))) {

#if 0
      if (M->length() > 0)
        fprintf(stderr, "LOADMORE length="uint32FMT" lastLoaded=%s first=%s\n",
                M->length(), lastLoaded.merToString(debugstring2), M->first().merToString(debugstring));
      else
        fprintf(stderr, "LOADMORE length="uint32FMT" lastLoaded=%s first=EMPTY\n",
                M->length(), lastLoaded.merToString(debugstring2));
#endif

      //  We need to copy all the mers currently loaded into fresh
      //  storage, so we can deallocate the position storage.  Yucky.
      //
      M->rebuild();

      allLoaded = true;

      //  Load more stuff to give us a large collection of mers
      //
      uint32 additionalLoading = 8192;

      if (maxSize / args->mergeFilesLen > M->length())
        additionalLoading = maxSize / args->mergeFilesLen - M->length();

      //fprintf(stderr, "LOADMORE adding "uint32FMT" from each file\n", additionalLoading);

      lastLoaded.setMerSize(merSize);
      lastLoaded.smallest();

      for (uint32 i=0; i<args->mergeFilesLen; i++) {
        if (R[i]->validMer()) {
          M->read(R[i], additionalLoading);
          if (lastLoaded < M->first())
            lastLoaded = M->first();
          allLoaded = false;
        }
      }

      //  Make sure all files have at least that largest mer loaded.
      //
      for (uint32 i=0; i<args->mergeFilesLen; i++)
        while (R[i]->validMer() && (R[i]->theFMer() <= lastLoaded))
          M->read(R[i], 2 * 1024);

      M->sort();
    }


    //  All done?  Exit.
    if (M->length() == 0)
      moreStuff = false;

    thisMer = M->pop(thisCount, thisPositions);

    //  If we've hit a different mer, write out the last one
    if ((M->length() == 0) || (*thisMer != currentMer)) {
      switch (args->personality) {
        case PERSONALITY_MIN:
          if (currentTimes == args->mergeFilesLen)
            W->addMer(currentMer, currentCount);
          break;
        case PERSONALITY_MERGE:
        case PERSONALITY_MINEXIST:
        case PERSONALITY_MAX:
        case PERSONALITY_ADD:
          W->addMer(currentMer, currentCount, currentPositions);
          break;
        case PERSONALITY_AND:
          if (currentTimes == args->mergeFilesLen)
            W->addMer(currentMer, currentCount);
          break;
        case PERSONALITY_NAND:
          if (currentTimes != args->mergeFilesLen)
            W->addMer(currentMer, currentCount);
          break;
        case PERSONALITY_OR:
          W->addMer(currentMer, currentCount);
          break;
        case PERSONALITY_XOR:
          if ((currentTimes % 2) == 1)
            W->addMer(currentMer, currentCount);
          break;
        default:
          fprintf(stderr, "ERROR - invalid personality in multipleOperations::write\n");
          fprintf(stderr, "ERROR - this is a coding error, not a user error.\n");
          exit(1);
          break;
      }

      currentMer = *thisMer;

      currentCount = uint32ZERO;
      currentTimes = uint32ZERO;

      C->tick();
    }


    if (moreStuff == false)
      break;


    //  Perform the operation
    switch (args->personality) {
      case PERSONALITY_MERGE:
        if (thisPositions) {

          if (currentPositionsMax == 0) {
            currentPositionsMax = 1048576;
            currentPositions    = new uint32 [currentPositionsMax];
          }

          if (currentPositionsMax < currentCount + thisCount) {
            while (currentPositionsMax < currentCount + thisCount)
              currentPositionsMax *= 2;

            uint32 *t = new uint32 [currentPositionsMax];
            memcpy(t, currentPositions, sizeof(uint32) * currentCount);
            delete [] currentPositions;
            currentPositions = t;
          }

          if (thisCount < 16) {
            for (uint32 i=0; i<thisCount; i++)
              currentPositions[currentCount + i] = thisPositions[i];
          } else {
            memcpy(currentPositions + currentCount, thisPositions, sizeof(uint32) * thisCount);
          }
        }
        //  Otherwise, we're the same as ADD.
        currentCount += thisCount;
        break;
      case PERSONALITY_MIN:
      case PERSONALITY_MINEXIST:
        if (currentTimes == 0) {
          currentCount = thisCount;
        } else {
          if (currentCount > thisCount)
            currentCount = thisCount;
        }
        break;
      case PERSONALITY_MAX:
        if (currentCount < thisCount)
          currentCount = thisCount;
        break;
      case PERSONALITY_ADD:
        currentCount += thisCount;
        break;
      case PERSONALITY_AND:
      case PERSONALITY_NAND:
      case PERSONALITY_OR:
      case PERSONALITY_XOR:
        currentCount = 1;
        break;
      default:
        fprintf(stderr, "ERROR - invalid personality in multipleOperations::operate\n");
        fprintf(stderr, "ERROR - this is a coding error, not a user error.\n");
        exit(1);
        break;
    }

    currentTimes++;
  }

  for (uint32 i=0; i<args->mergeFilesLen; i++)
    delete R[i];
  delete R;
  delete W;
  delete M;
  delete C;
}
