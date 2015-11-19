
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
 *    kmer/meryl/merge.listmerge.C
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
  uint32   _nxt;
  uint32   _stp;
};


class mMerList {
public:
  mMerList(uint32 maxSize) {
    _posLen = 0;
    _posMax = 2 * maxSize;
    _pos    = new uint32 [_posMax];

    _mmmLen = 0;
    _mmmMax = maxSize;
    _mmm    = new mMer [_mmmMax];

    _tip    = ~uint32ZERO;
    _fre    = 0;

    for (uint32 i=0; i<_mmmMax; i++) {
      _mmm[i]._cnt = 0;
      _mmm[i]._off = 0;
      _mmm[i]._nxt = i+1;
      _mmm[i]._stp = 0;
    }

    _mmm[_mmmMax-1]._nxt = ~uint32ZERO;
  };
  ~mMerList() {
    delete [] _pos;
    delete [] _mmm;
  };

  bool    loadMore(void)  { return((_mmmMax < _tip) || (_mmm[_tip]._stp == 1)); };

  uint32  length(void) { return(_mmmLen); };

  kMer   *pop(uint32 &cnt, uint32* &pos) {
    kMer  *ret = 0L;

    //fprintf(stderr, "POP tip="uint32FMT"\n", _tip);

    if (_tip < _mmmMax) {
      uint32 f = _tip;

      ret = &_mmm[f]._mer;
      cnt =  _mmm[f]._cnt;
      pos = (_mmm[f]._off != ~uint32ZERO) ? _pos + _mmm[f]._off : 0L;

      //  Move tip to the next thing
      _tip = _mmm[f]._nxt;

      //  And append this one to the free list.
      _mmm[f]._nxt = _fre;
      _fre = f;

      _mmmLen--;

      //fprintf(stderr, "POP f="uint32FMT" tip="uint32FMT" len="uint32FMT"\n", f, _tip, _mmmLen);
    }

    return(ret);
  };


  //  rebuild the position list, squeezes out empty items
  void    rebuild(void) {
    if (_posLen > 0) {
      assert(0);
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
  void    read(merylStreamReader *R, uint32 num, bool loadAll) {
    uint32 xxx  = 0;
    uint32 las  = ~uint32ZERO;
    uint32 pos  = _tip;
    bool   stop = false;

    //fprintf(stderr, "read()- loading "uint32FMT"\n", num);

    assert(_mmmLen + num < _mmmMax);

    //  Load until we hit the sentinal.
    if (loadAll == false)
      num = ~uint32ZERO;

    for (xxx=0; (xxx < num) && (stop == false) && (R->nextMer()); xxx++) {

      //  Insert into a free node
      uint32 fre = _fre;
      _fre = _mmm[fre]._nxt;

      _mmm[fre]._mer = R->theFMer();
      _mmm[fre]._cnt = R->theCount();
      _mmm[fre]._off = ~uint32ZERO;
      _mmm[fre]._stp = 0;

      uint32  *ppp = R->thePositions();
      if (ppp) {
        _mmm[fre]._off = _posLen;

        if (_posMax <= _posLen + _mmm[fre]._cnt) {
          fprintf(stderr, "Reallocate _pos\n");
          _posMax *= 2;
          uint32 *tmp = new uint32 [_posMax];
          memcpy(tmp, _pos, sizeof(uint32) * _posLen);
          delete [] _pos;
          _pos = tmp;
        }

        for (uint32 i=0; i<_mmm[fre]._cnt; i++, _posLen++)
          _pos[_posLen] = ppp[i];
      }

      //  Keep count
      _mmmLen++;

      //  Figure out where to put it in the list.  New duplicates must
      //  go AFTER the existing -- that's the job of <=.

      while ((pos < _mmmMax) && (_mmm[pos]._mer <= R->theFMer())) {
        las = pos;
        pos = _mmm[pos]._nxt;
      }

      if (_mmmMax < _tip) {
        //  No tip, make new list.
        _mmm[fre]._nxt = _tip;
        _tip           = fre;
        las = ~uint32ZERO;
        pos = _tip;
      } else if (_mmmMax < las) {
        //  Valid list, but we want to insert before the start
        _mmm[fre]._nxt = _tip;
        _tip           = fre;
        las = ~uint32ZERO;
        pos = _tip;
      } else if (pos < _mmmMax) {
        //  Valid pos, insert in the middle (after las, before pos)
        _mmm[fre]._nxt = _mmm[las]._nxt;
        _mmm[las]._nxt = fre;
        las = fre;
        //pos = _mmm[las]._nxt;
      } else {
        //  Have a list, but we ran off the end, append (after las)
        _mmm[fre]._nxt = ~uint32ZERO;
        _mmm[las]._nxt = fre;
        pos = fre;

        if (loadAll == false)
          stop = true;
      }
    }

    //  Set the sentinal.  This forces us to load more mers.
    //
    if (loadAll == true) {
      //fprintf(stderr, "read()-- stop on tip = "uint32FMT"\n", las);
      _mmm[las]._stp = 1;
    }

    //fprintf(stderr, "read()-- now up to "uint32FMT" mers ("uint32FMT" pos); loaded "uint32FMT" out of "uint32FMT" requested.\n", _mmmLen, _posLen, xxx, num);
  };

private:
  uint32  _posLen;
  uint32  _posMax;
  uint32 *_pos;

  uint32  _mmmLen;
  uint32  _mmmMax;
  mMer   *_mmm;

  uint32  _tip;
  uint32  _fre;
};




void
multipleOperations(merylArgs *args) {

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

  uint32               maxSize = 64 * 1024 * 1024;

  merylStreamReader  **R = new merylStreamReader* [args->mergeFilesLen];
  merylStreamWriter   *W = 0L;
  mMerList            *M = new mMerList(maxSize + maxSize / 4);

  for (uint32 i=0; i<args->mergeFilesLen; i++)
    R[i] = new merylStreamReader(args->mergeFiles[i]);

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
    fprintf(stderr, "ERROR:  mer size or compression level differ.\n"), exit(1);

  //  Open the output file, using the largest prefix size found in the
  //  input/mask files.
  //
  uint32  prefixSize = 0;
  for (uint32 i=0; i<args->mergeFilesLen; i++)
    if (prefixSize < R[i]->prefixSize())
      prefixSize = R[i]->prefixSize();

  W = new merylStreamWriter(args->outputFile, merSize, merComp, prefixSize);

  //  Load mers from all files, remember the largest mer we load.
  //
  bool     loadAll = true;
  for (uint32 i=0; i<args->mergeFilesLen; i++) {
    M->read(R[i], maxSize / args->mergeFilesLen, loadAll);
    loadAll = false;
  }

  fprintf(stderr, "Initial load:  length="uint32FMT"\n", M->length());

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

    //  Load more stuff if needed.
    //
    if (M->loadMore() == true) {
      M->rebuild();

      uint32 additionalLoading = 8192;

      if (maxSize / args->mergeFilesLen > M->length())
        additionalLoading = maxSize / args->mergeFilesLen - M->length();

      loadAll   = true;

      for (uint32 i=0; i<args->mergeFilesLen; i++) {
        if (R[i]->validMer()) {
          M->read(R[i], additionalLoading, loadAll);
          loadAll   = false;
        }
      }
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
