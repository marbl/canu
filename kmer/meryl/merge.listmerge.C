#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "meryl.H"
#include "libmeryl.H"


using namespace std;

#include <algorithm>

struct mMer {
  kMer     _mer;
  u32bit   _cnt;
  u32bit   _off;
  u32bit   _nxt;
  u32bit   _stp;
};


class mMerList {
public:
  mMerList(u32bit maxSize) {
    _posLen = 0;
    _posMax = 2 * maxSize;
    _pos    = new u32bit [_posMax];

    _mmmLen = 0;
    _mmmMax = maxSize;
    _mmm    = new mMer [_mmmMax];

    _tip    = ~u32bitZERO;
    _fre    = 0;

    for (u32bit i=0; i<_mmmMax; i++) {
      _mmm[i]._cnt = 0;
      _mmm[i]._off = 0;
      _mmm[i]._nxt = i+1;
      _mmm[i]._stp = 0;
    }

    _mmm[_mmmMax-1]._nxt = ~u32bitZERO;
  };
  ~mMerList() {
    delete [] _pos;
    delete [] _mmm;
  };

  bool    loadMore(void)  { return((_mmmMax < _tip) || (_mmm[_tip]._stp == 1)); };

  u32bit  length(void) { return(_mmmLen); };

  kMer   *pop(u32bit &cnt, u32bit* &pos) {
    kMer  *ret = 0L;

    //fprintf(stderr, "POP tip="u32bitFMT"\n", _tip);

    if (_tip < _mmmMax) {
      u32bit f = _tip;

      ret = &_mmm[f]._mer;
      cnt =  _mmm[f]._cnt;
      pos = (_mmm[f]._off != ~u32bitZERO) ? _pos + _mmm[f]._off : 0L;

      //  Move tip to the next thing
      _tip = _mmm[f]._nxt;

      //  And append this one to the free list.
      _mmm[f]._nxt = _fre;
      _fre = f;

      _mmmLen--;

      //fprintf(stderr, "POP f="u32bitFMT" tip="u32bitFMT" len="u32bitFMT"\n", f, _tip, _mmmLen);
    }

    return(ret);
  };


  //  rebuild the position list, squeezes out empty items
  void    rebuild(void) {
    if (_posLen > 0) {
      assert(0);
      u32bit   *np = new u32bit [_posMax];

      _posLen = 0;
    
      for (u32bit i=0; i<_mmmLen; i++) {
        mMer *m = _mmm + i;

        if (m->_off != ~u32bitZERO) {
          _mmm[_mmmLen]._off = _posLen;

          for (u32bit p=0; p<m->_cnt; p++, _posLen++)
            np[_posLen] = _pos[p];
        }
      }

      delete [] _pos;
      _pos = np;
    }
  };



  //  Read more mers from the file
  void    read(merylStreamReader *R, u32bit num, bool loadAll) {
    u32bit xxx  = 0;
    u32bit las  = ~u32bitZERO;
    u32bit pos  = _tip;
    bool   stop = false;

    //fprintf(stderr, "read()- loading "u32bitFMT"\n", num);

    assert(_mmmLen + num < _mmmMax);

    //  Load until we hit the sentinal.
    if (loadAll == false)
      num = ~u32bitZERO;

    for (xxx=0; (xxx < num) && (stop == false) && (R->nextMer()); xxx++) {

      //  Insert into a free node
      u32bit fre = _fre;
      _fre = _mmm[fre]._nxt;

      _mmm[fre]._mer = R->theFMer();
      _mmm[fre]._cnt = R->theCount();
      _mmm[fre]._off = ~u32bitZERO;
      _mmm[fre]._stp = 0;

      u32bit  *ppp = R->thePositions();
      if (ppp) {
        _mmm[fre]._off = _posLen;

        if (_posMax <= _posLen + _mmm[fre]._cnt) {
          fprintf(stderr, "Reallocate _pos\n");
          _posMax *= 2;
          u32bit *tmp = new u32bit [_posMax];
          memcpy(tmp, _pos, sizeof(u32bit) * _posLen);
          delete [] _pos;
          _pos = tmp;
        }

        for (u32bit i=0; i<_mmm[fre]._cnt; i++, _posLen++)
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
        las = ~u32bitZERO;
        pos = _tip;
      } else if (_mmmMax < las) {
        //  Valid list, but we want to insert before the start
        _mmm[fre]._nxt = _tip;
        _tip           = fre;
        las = ~u32bitZERO;
        pos = _tip;
      } else if (pos < _mmmMax) {
        //  Valid pos, insert in the middle (after las, before pos)
        _mmm[fre]._nxt = _mmm[las]._nxt;
        _mmm[las]._nxt = fre;
        las = fre;
        //pos = _mmm[las]._nxt;
      } else {
        //  Have a list, but we ran off the end, append (after las)
        _mmm[fre]._nxt = ~u32bitZERO;
        _mmm[las]._nxt = fre;
        pos = fre;

        if (loadAll == false)
          stop = true;
      }
    }

    //  Set the sentinal.  This forces us to load more mers.
    //
    if (loadAll == true) {
      //fprintf(stderr, "read()-- stop on tip = "u32bitFMT"\n", las);
      _mmm[las]._stp = 1;
    }

    //fprintf(stderr, "read()-- now up to "u32bitFMT" mers ("u32bitFMT" pos); loaded "u32bitFMT" out of "u32bitFMT" requested.\n", _mmmLen, _posLen, xxx, num);
  };

private:
  u32bit  _posLen;
  u32bit  _posMax;
  u32bit *_pos;

  u32bit  _mmmLen;
  u32bit  _mmmMax;
  mMer   *_mmm;

  u32bit  _tip;
  u32bit  _fre;
};




void
multipleOperations(merylArgs *args) {

  if (args->mergeFilesLen < 2) {
    fprintf(stderr, "ERROR - must have at least two databases (you gave "u32bitFMT")!\n", args->mergeFilesLen);
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

  u32bit               maxSize = 64 * 1024 * 1024;

  merylStreamReader  **R = new merylStreamReader* [args->mergeFilesLen];
  merylStreamWriter   *W = 0L;
  mMerList            *M = new mMerList(maxSize + maxSize / 4);

  for (u32bit i=0; i<args->mergeFilesLen; i++)
    R[i] = new merylStreamReader(args->mergeFiles[i]);

  //  Verify that the mersizes are all the same
  //
  bool    fail       = false;
  u32bit  merSize    = R[0]->merSize();
  u32bit  merComp    = R[0]->merCompression();

  for (u32bit i=0; i<args->mergeFilesLen; i++) {
    fail |= (merSize != R[i]->merSize());
    fail |= (merComp != R[i]->merCompression());
  }

  if (fail)
    fprintf(stderr, "ERROR:  mer size or compression level differ.\n"), exit(1);

  //  Open the output file, using the largest prefix size found in the
  //  input/mask files.
  //
  u32bit  prefixSize = 0;
  for (u32bit i=0; i<args->mergeFilesLen; i++)
    if (prefixSize < R[i]->prefixSize())
      prefixSize = R[i]->prefixSize();

  W = new merylStreamWriter(args->outputFile, merSize, merComp, prefixSize);

  //  Load mers from all files, remember the largest mer we load.
  //
  bool     loadAll = true;
  for (u32bit i=0; i<args->mergeFilesLen; i++) {
    M->read(R[i], maxSize / args->mergeFilesLen, loadAll);
    loadAll = false;
  }

  fprintf(stderr, "Initial load:  length="u32bitFMT"\n", M->length());

  bool     moreStuff = true;

  kMer     currentMer;                      //  The current mer we're operating on
  u32bit   currentCount     =  u32bitZERO;  //  The count (operation dependent) of this mer
  u32bit   currentTimes     =  u32bitZERO;  //  Number of files it's in

  u32bit   currentPositionsMax =  0;
  u32bit  *currentPositions    =  0L;

  kMer    *thisMer;                         //  The mer we just read
  u32bit   thisCount        =  u32bitZERO;  //  The count of the mer we just read
  u32bit  *thisPositions    = 0L;

  speedCounter *C = new speedCounter("    %7.2f Mmers -- %5.2f Mmers/second\r", 1000000.0, 0x1fffff, args->beVerbose);

  currentMer.setMerSize(merSize);

  while (moreStuff) {

    //  Load more stuff if needed.
    //
    if (M->loadMore() == true) {
      M->rebuild();

      u32bit additionalLoading = 8192;

      if (maxSize / args->mergeFilesLen > M->length())
        additionalLoading = maxSize / args->mergeFilesLen - M->length();

      loadAll   = true;

      for (u32bit i=0; i<args->mergeFilesLen; i++) {
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

      currentCount = u32bitZERO;
      currentTimes = u32bitZERO;

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
            currentPositions    = new u32bit [currentPositionsMax];
          }

          if (currentPositionsMax < currentCount + thisCount) {
            while (currentPositionsMax < currentCount + thisCount)
              currentPositionsMax *= 2;

            u32bit *t = new u32bit [currentPositionsMax];
            memcpy(t, currentPositions, sizeof(u32bit) * currentCount);
            delete [] currentPositions;
            currentPositions = t;
          }

          if (thisCount < 16) {
            for (u32bit i=0; i<thisCount; i++)
              currentPositions[currentCount + i] = thisPositions[i];
          } else {
            memcpy(currentPositions + currentCount, thisPositions, sizeof(u32bit) * thisCount);
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

  for (u32bit i=0; i<args->mergeFilesLen; i++)
    delete R[i];
  delete R;
  delete W;
  delete M;
  delete C;
}
