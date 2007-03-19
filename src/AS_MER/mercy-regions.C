
/**************************************************************************
 * This file is part of Celera Assembler, a software program that 
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 2006-2007, J. Craig Venter Institute
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

//  This reads the assembly frgctg, varctg and merQC badmers, computes
//  the number and location of bad-mer, bad-var regions, and their
//  depth, in contig space.
//
//  File paths are hardcoded.
//  This code ONLY works on 64-bit hardware, but it's easy to fix.

using namespace std;
#include <map>

//
//  Begin code from Bri's intervalList.H, intervalList.C, splitToWords.H
//
typedef unsigned long     u64bit;
typedef unsigned int      u32bit;

#define  u64bitFMT       "%lu"
#define  u32bitFMT       "%u"
#define  u32bitFMTW(X)   "%" #X "u"

#define  strtou32bit(N,O) (u32bit)strtoul(N, O, 10)
#define  strtou64bit(N,O) (u64bit)strtoul(N, O, 10)

class splitToWords {
public:
  splitToWords() {
    _argWords = 0;
    _maxWords = 0;
    _arg      = 0L;
    _maxChars = 0;
    _cmd      = 0L;
  };
  splitToWords(char *cmd) {
    _argWords = 0;
    _maxWords = 0;
    _arg      = 0L;
    _maxChars = 0;
    _cmd      = 0L;

    split(cmd);
  };
  ~splitToWords() {
    delete [] _cmd;
    delete [] _arg;
  };


  void   split(char *cmd) {

    //  Step Zero:
    //
    //  Count the length of the string, in words and in characters.
    //  For simplicity, we overcount words, by just counting white-space.
    //
    //  Then, allocate space for a temporary copy of the string, and a
    //  set of pointers into the temporary copy (much like argv).
    //
    u32bit   cmdChars = 1;  //  1 == Space for terminating 0
    u32bit   cmdWords = 2;  //  2 == Space for first word and terminating 0L

    for (char *tmp=cmd; *tmp; tmp++) {
      cmdWords += *tmp == ' ';
      cmdWords += *tmp == '\t';
      cmdChars++;
    }

    if (cmdChars > _maxChars) {
      delete [] _cmd;
      _cmd      = new char   [cmdChars];
      _maxChars = cmdChars;
    }
    if (cmdWords > _maxWords) {
      delete [] _arg;
      _arg      = new char * [cmdWords];
      _maxWords = cmdWords;
    }

    _argWords = 0;

    //  Step One:
    //
    //  Determine where the words are in the command string, copying the
    //  string to _cmd and storing words in _arg.
    //
    bool           isFirst  = true;
    char          *cmdI = cmd;
    char          *cmdO = _cmd;

    while (*cmdI) {

      //  If we are at a non-space character, we are in a word.  If
      //  this is the first character in the word, save the word in
      //  the args list.
      //
      //  Otherwise we are at a space and thus not in a word.  Make
      //  all spaces be string terminators, and declare that we are
      //  at the start of a word.
      //
      if ((*cmdI != ' ') && (*cmdI != '\t')) {
        *cmdO = *cmdI;

        if (isFirst) {
          _arg[_argWords++] = cmdO;
          isFirst           = false;
        }
      } else {
        *cmdO   = 0;
        isFirst = true;
      }

      cmdI++;
      cmdO++;
    }

    //  Finish off the list by terminating the last arg, and
    //  terminating the list of args.
    //
    *cmdO           = 0;
    _arg[_argWords] = 0L;
  };


  u32bit  numWords(void)        { return(_argWords); };
  char   *getWord(u32bit i)     { return(_arg[i]); };
  char   *operator[](u32bit i)  { return(_arg[i]); };
private:
  u32bit    _argWords;
  u32bit    _maxWords;
  char    **_arg;
  u32bit    _maxChars;
  char     *_cmd;
};




typedef u64bit  intervalNumber;

struct _intervalPair {
  intervalNumber    lo;
  intervalNumber    hi;
};

struct _intervalDepth {
  intervalNumber    lo;
  intervalNumber    hi;
  u32bit            de;
};


class intervalList {
public:
  intervalList();
  ~intervalList();

  intervalList &operator=(intervalList &src);

  //  Clear a list
  void        clear(void) {
    _isSorted = true;
    _isMerged = true;
    _listLen  = 0;
  }

  //  Insert a new interval into the list
  void        add(intervalNumber position, intervalNumber length);

  //  Sort the set of intervals by the lo value
  void        sort(void);

  //  Merge overlapping or adjacent intervals together.
  void        merge(void);

  void        invert(intervalNumber lo, intervalNumber hi);

  //  Returns the number of intervals
  u32bit      numberOfIntervals(void) {
    return(_listLen);
  };

  //  Returns the sum of the length of all intervals
  intervalNumber      sumOfLengths(void) {
    intervalNumber len = 0;
    u32bit         i   = numberOfIntervals();

    if (i > 0)
      while (i--)
        len += _list[i].hi - _list[i].lo;

    return(len);
  };

  //  Populates an array with the intervals that are within the
  //  supplied interval.  Return
  //
  u32bit      overlapping(intervalNumber    lo,
                          intervalNumber    hi,
                          u32bit          *&intervals,
                          u32bit           &intervalsLen,
                          u32bit           &intervalsMax);

  //  Populates this intervalList with the intersection of A and B.
  //  This intervalList is not cleared prior to adding new intervals.
  //
  //  Both A and B call merge().
  //
  void                intersect(intervalList &A,
                                intervalList &B);

  //  Populates this intervalList with regions in A that are completely
  //  contained in a region in B.
  //
  //  Both A and B call merge().
  //
  void                contained(intervalList &A,
                                intervalList &B);


  intervalNumber      lo(u32bit i) { return(_list[i].lo); };
  intervalNumber      hi(u32bit i) { return(_list[i].hi); };

private:
  bool                      _isSorted;
  bool                      _isMerged;
  u32bit                    _listLen;
  u32bit                    _listMax;
  _intervalPair            *_list;
};



//  Takes as input an intervalList, computes the number of intervals
//  covering every position in there, stores this as a new set of
//  intervals, annotated with the depth.
//
//  This is a static object, initialized once by the intervalList.
//
class intervalDepth {
public:
  intervalDepth(intervalList &IL);
  ~intervalDepth();

  //  Returns the number of intervals
  u32bit      numberOfIntervals(void) {
    return(_listLen);
  };

  intervalNumber      lo(u32bit i) { return(_list[i].lo); };
  intervalNumber      hi(u32bit i) { return(_list[i].hi); };
  u32bit              de(u32bit i) { return(_list[i].de); };

private:
  u32bit                    _listLen;
  u32bit                    _listMax;
  _intervalDepth           *_list;
};


intervalList::intervalList() {
  _isSorted = true;
  _isMerged = true;
  _listLen  = 0;
  _listMax  = 16;
  _list     = new _intervalPair [_listMax];
}


intervalList::~intervalList() {
  delete [] _list;
}


intervalList &
intervalList::operator=(intervalList &src) {
  _isSorted = src._isSorted;
  _isMerged = src._isMerged;
  _listLen = src._listLen;

  if (_listMax < src._listMax) {
    delete [] _list;
    _listMax = src._listMax;
    _list    = new _intervalPair [_listMax];
  }

  memcpy(_list, src._list, _listLen * sizeof(_intervalPair));

  return(*this);
}


void
intervalList::add(intervalNumber position, intervalNumber length) {

  if (_listLen >= _listMax) {
    _listMax *= 2;
    _intervalPair *l = new _intervalPair [_listMax];
    memcpy(l, _list, sizeof(_intervalPair) * _listLen);
    delete [] _list;
    _list = l;
  }

  _list[_listLen].lo   = position;
  _list[_listLen].hi   = position + length;

#if 0
  //  Aborted attempt to add a data field here.  Got stuck
  //  deciding how to handle merges lightweight

  _list[_listLen].data = 0L;

  if (data != ~u64bitZERO) {
    _list[_listLen].dataLen = 1;
    _list[_listLen].dataMax = 4;
    _list[_listLen].data    = new u64bit [_list[_listLen].dataMax];
    _list[_listLen].data[0] = data;
  }
#endif
    
  if ((_listLen > 0) &&
      (_list[_listLen-1].lo > _list[_listLen].lo)) {
    _isSorted = false;
    _isMerged = false;
  }

  _listLen++;
}


static
int
intervalList_sort_helper(const void *a, const void *b) {
  _intervalPair *A = (_intervalPair *)a;
  _intervalPair *B = (_intervalPair *)b;

  if (A->lo < B->lo) return(-1);
  if (A->lo > B->lo) return(1);
  if (A->hi < B->hi) return(-1);
  if (A->hi > B->hi) return(1);
  return(0);
}


void
intervalList::sort(void) {

  if (_isSorted)
    return;

  if (_listLen > 1)
    qsort(_list, _listLen, sizeof(_intervalPair), intervalList_sort_helper);

  _isSorted = true;
}


void
intervalList::merge(void) {
  u32bit  thisInterval  = 0;
  u32bit  nextInterval = 1;

  if (_listLen < 2)
    return;

  sort();

  while (nextInterval < _listLen) {

    if ((_list[thisInterval].lo == 0) &&
        (_list[thisInterval].hi == 0)) {

      //  Our interval is empty.  Copy in the interval we are
      //  examining and move to the next.

      //  XXX This is probably useless, thisInterval should always be
      //  valid.

      _list[thisInterval].lo = _list[nextInterval].lo;
      _list[thisInterval].hi = _list[nextInterval].hi;

      _list[nextInterval].lo = 0;
      _list[nextInterval].hi = 0;
      nextInterval++;
    } else {

      //  This interval is valid.  See if it overlaps with the next
      //  interval.

      if (_list[thisInterval].hi >= _list[nextInterval].lo) {

        //  Got an intersection.

        //  Merge nextInterval into thisInterval -- the hi range
        //  is extended if the nextInterval range is larger.
        //
        if (_list[thisInterval].hi < _list[nextInterval].hi)
          _list[thisInterval].hi = _list[nextInterval].hi;
        
        //  Clear the just merged nextInterval and move to the next one.
        //
        _list[nextInterval].lo = 0;
        _list[nextInterval].hi = 0;
        nextInterval++;
      } else {

        //  No intersection.  Move along.  Nothing to see here.

        //  If there is a gap between the target and the examine (we
        //  must have merged sometime in the past), copy examine to
        //  the next target.

        thisInterval++;

        if (thisInterval != nextInterval) {
          _list[thisInterval].lo = _list[nextInterval].lo;
          _list[thisInterval].hi = _list[nextInterval].hi;
        }

        nextInterval++;
      }
    }
  }

  if (thisInterval+1 < _listLen)
    _listLen = thisInterval + 1;

  _isMerged = true;
}


void
intervalList::invert(intervalNumber lo, intervalNumber hi) {

  if (!_isSorted || !_isMerged) {
    fprintf(stderr, "intervalList::invert()--  ERROR!  List is not sorted or not merged!\n");
    exit(1);
  }

  //  Create a new list to store the inversion
  //
  u32bit             invLen = 0;
  u32bit             invMax = _listLen + 2;
  _intervalPair     *inv    = new _intervalPair [invMax];

  //  Add the first
  //
  if (lo < _list[0].lo) {
    inv[invLen].lo = lo;
    inv[invLen].hi = _list[0].lo;
    invLen++;
  }

  //  Add the pieces
  for (u32bit i=1; i<_listLen; i++) {
    if (_list[i-1].hi < _list[i].lo) {
      inv[invLen].lo = _list[i-1].hi;
      inv[invLen].hi = _list[i].lo;
      invLen++;
    }
  }

  //  Add the last
  if (_list[_listLen-1].hi < hi) {
    inv[invLen].lo = _list[_listLen-1].hi;
    inv[invLen].hi = hi;
    invLen++;
  }

  //  Nuke the old list, swap in the new one
  delete [] _list;

  _list = inv;
  _listLen = invLen;
  _listMax = invMax;
}



u32bit
intervalList::overlapping(intervalNumber    rangelo,
                          intervalNumber    rangehi,
                          u32bit          *&intervals,
                          u32bit           &intervalsLen,
                          u32bit           &intervalsMax) {


  //  XXX: Naive implementation that is easy to verify (and that works
  //  on an unsorted list).

  if (intervals == 0L) {
    intervalsMax = 256;
    intervals    = new u32bit [intervalsMax];
  }

  intervalsLen = 0;

  for (u32bit i=0; i<_listLen; i++) {
    if ((rangelo <= _list[i].hi) &&
        (rangehi >= _list[i].lo)) {
      if (intervalsLen >= intervalsMax) {
        intervalsMax *= 2;
        u32bit *X = new u32bit [intervalsMax];
        memcpy(X, intervals, sizeof(u32bit) * intervalsLen);
        delete [] intervals;
        intervals = X;
      }

      intervals[intervalsLen++] = i;
    }
  }

  return(intervalsLen);
}



void
intervalList::intersect(intervalList &A,
                        intervalList &B) {
  A.merge();
  B.merge();

  u32bit  ai = 0;
  u32bit  bi = 0;

  while ((ai < A.numberOfIntervals()) &&
         (bi < B.numberOfIntervals())) {
    u32bit   al = A.lo(ai);
    u32bit   ah = A.hi(ai);
    u32bit   bl = B.lo(bi);
    u32bit   bh = B.hi(bi);
    u32bit   nl = 0;
    u32bit   nh = 0;

    //  If they intersect, make a new region
    //
    if ((al <= bl) && (bl < ah)) {
      nl = bl;
      nh = (ah < bh) ? ah : bh;
    }

    if ((bl <= al) && (al < bh)) {
      nl = al;
      nh = (ah < bh) ? ah : bh;
    }

    if (nl < nh)
      add(nl, nh - nl);

    //  Advance the list with the earlier region.
    //
    if        (ah < bh) {
      //  A ends before B
      ai++;
    } else if (ah > bh) {
      //  B ends before A
      bi++;
    } else {
      //  Exactly the same ending!
      ai++;
      bi++;
    }
  }
}

void
intervalList::contained(intervalList &A,
                        intervalList &B) {
  A.merge();
  B.merge();

  u32bit  ai = 0;
  u32bit  bi = 0;

  while ((ai < A.numberOfIntervals()) &&
         (bi < B.numberOfIntervals())) {
    u32bit   al = A.lo(ai);
    u32bit   ah = A.hi(ai);
    u32bit   bl = B.lo(bi);
    u32bit   bh = B.hi(bi);

    //  If A is contained in B, make a new region.
    //
    if ((bl <= al) && (ah <= bh))
      add(bl, bh - bl);

#if 0
    if ((al <= bl) && (bh <= ah))
      add(al, ah - al);
#endif

    //  Advance the list with the earlier region.
    //
    if        (ah < bh) {
      //  A ends before B
      ai++;
    } else if (ah > bh) {
      //  B ends before A
      bi++;
    } else {
      //  Exactly the same ending!
      ai++;
      bi++;
    }
  }
}






static
int
intervalDepth_sort_helper(const void *a, const void *b) {
  _intervalDepth *A = (_intervalDepth *)a;
  _intervalDepth *B = (_intervalDepth *)b;

  if (A->lo < B->lo) return(-1);
  if (A->lo > B->lo) return(1);
  return(0);
}


intervalDepth::intervalDepth(intervalList &IL) {

  u32bit           idlen = IL.numberOfIntervals() * 2;
  _intervalDepth  *id    = new _intervalDepth [idlen];

  for (u32bit i=0; i<IL.numberOfIntervals(); i++) {
    id[2*i  ].lo = IL.lo(i);
    id[2*i  ].hi = 0;
    id[2*i  ].de = 1;
    id[2*i+1].lo = IL.hi(i);
    id[2*i+1].hi = 0;
    id[2*i+1].de = 0;
  }

  qsort(id, idlen, sizeof(_intervalDepth), intervalDepth_sort_helper);

  //  Scan the list, counting how many times we change depth.
  //
  _listMax = 1;
  for (u32bit i=1; i<idlen; i++) {
    if (id[i-1].lo != id[i].lo)
      _listMax++;
  }

  //  Allocate the real depth of coverage intervals
  //
  _listLen = 0;
  _list    = new _intervalDepth [_listMax];

  //  Build new intervals
  //
  //  Initialize the first interval
  //
  _list[_listLen].lo = id[0].lo;
  _list[_listLen].hi = id[0].lo;
  _list[_listLen].de = 1;

  for (u32bit i=1; i<idlen; i++) {

    if (_list[_listLen].de == 0) {
      //  Update the start position if the current interval is at zero
      //  depth.
      //
      _list[_listLen].lo = id[i].lo;
    } else {

      //  If we are at a position different from the start, we need to
      //  close out the current interval and make a new one.
      //
      if (id[i-1].lo != id[i].lo) {
        _list[_listLen].hi = id[i].lo;

        _listLen++;

        _list[_listLen].lo = id[i].lo;
        _list[_listLen].hi = id[i].lo;
        _list[_listLen].de = _list[_listLen-1].de;
      }
    }

    //  Finally, update the depth of the current interval
    //
    if (id[i].de)
      _list[_listLen].de++;
    else
      _list[_listLen].de--;
  }

  //  Toss out the last one if it's zero length -- I think it's always
  //  zero length, just can convince myself.
  //
  if (_list[_listLen].lo == _list[_listLen].hi)
    _listLen--;

  delete [] id;
}

intervalDepth::~intervalDepth() {
  delete [] _list;
}

//
//  End code from Bri's libutil/util++.H
//



void
readDepth(char *depthname, map<u64bit,intervalDepth*> &lowCoverage) {
  char                         line[1024] = {0};
  map<u64bit,intervalList*>    ILs;

  fprintf(stderr, "Reading depth from '%s'\n", depthname);

  errno = 0;
  FILE *F = fopen(depthname, "r");
  if (errno)
    fprintf(stderr, "failed to open '%s': %s\n", depthname, strerror(errno)), exit(1);

  u32bit i=0;

  fgets(line, 1024, F);
  while (!feof(F)) {
    splitToWords   W(line);

    u64bit  uid = strtou64bit(W[1], 0L);
    u32bit  beg = strtou32bit(W[2], 0L);
    u32bit  end = strtou32bit(W[3], 0L);

    if (beg > end)
      fprintf(stderr, "ERROR: l="u32bitFMT" h="u32bitFMT"\n", beg, end);

    if (ILs[uid] == 0L)
      ILs[uid] = new intervalList();
    ILs[uid]->add(beg, end - beg);

    i++;

    fgets(line, 1024, F);
  }

  fclose(F);
  fprintf(stderr, " "u32bitFMT" lines.\n", i);

  map<u64bit,intervalList*>::iterator    it = ILs.begin();
  map<u64bit,intervalList*>::iterator    ed = ILs.end();

  while (it != ed) {
    lowCoverage[it->first] = new intervalDepth(*it->second);
    delete it->second;
    it->second = 0L;
    it++;
  }
}


void
readVariation(char *depthname, map<u64bit,intervalList*> &variation) {
  char                         line[1024 * 1024] = {0};

  fprintf(stderr, "Reading variation from '%s'\n", depthname);

  errno = 0;
  FILE *F = fopen(depthname, "r");
  if (errno)
    fprintf(stderr, "failed to open '%s': %s\n", depthname, strerror(errno)), exit(1);

  u32bit i=0;

  fgets(line, 1024 * 1024, F);
  while (!feof(F)) {
    splitToWords   W(line);

    u64bit  uid = strtou64bit(W[1], 0L);
    u32bit  beg = strtou32bit(W[2], 0L);
    u32bit  end = strtou32bit(W[3], 0L);

    if (variation[uid] == 0L)
      variation[uid] = new intervalList();
    variation[uid]->add(beg, end - beg);

    i++;

    fgets(line, 1024 * 1024, F);
  }

  fclose(F);
  fprintf(stderr, " "u32bitFMT" lines.\n", i);
}


void
readBadMers(char *depthname, map<u64bit,intervalList*> &badMers) {
  char                         line[1024] = {0};

  fprintf(stderr, "Reading badMers from '%s'\n", depthname);

  errno = 0;
  FILE *F = fopen(depthname, "r");
  if (errno)
    fprintf(stderr, "failed to open '%s': %s\n", depthname, strerror(errno)), exit(1);

  u32bit i=0;

  fgets(line, 1024, F);
  while (!feof(F)) {
    splitToWords   W(line);

    //  Change every non-digit to a space in the first word.
    for (u32bit z=strlen(W[0])-1; z--; )
      if (!isdigit(W[0][z]))
        W[0][z] = ' ';

    u64bit  uid = strtou64bit(W[0], 0L);
    u32bit  beg = strtou32bit(W[3], 0L);
    u32bit  end = strtou32bit(W[4], 0L);

    if (badMers[uid] == 0L)
      badMers[uid] = new intervalList();
    badMers[uid]->add(beg, end - beg);

    i++;

    fgets(line, 1024, F);
  }

  fclose(F);
  fprintf(stderr, " "u32bitFMT" lines.\n", i);
}



int
main(int argc, char **argv) {
  map<u64bit,intervalList*>    badMers;
  map<u64bit,intervalList*>    variation;
  map<u64bit,intervalDepth*>   lowCoverage;

  bool  showDepthIntersect    = false;
  bool  showVariantIntersect  = false;
  bool  showVarDepthIntersect = false;

  int arg=1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-D") == 0) {

    } else if (strcmp(argv[arg], "-pd") == 0) {
      showDepthIntersect = true;
    } else if (strcmp(argv[arg], "-pv") == 0) {
      showVariantIntersect = true;
    } else if (strcmp(argv[arg], "-pvd") == 0) {
      showVarDepthIntersect = true;
    } else {
      fprintf(stderr, "usage: %s [-D debugfile] [-pd] [-pv] [-pvd]\n", argv[0]);
      fprintf(stderr, " -pd    print bad mers regions isect depth\n");
      fprintf(stderr, " -pv    print bad mers regions isect variants\n");
      fprintf(stderr, " -pvd   print bad mers regions isect both variants and depth\n");
      exit(1);
    }
    arg++;
  }

#if 1
  //  HuRef6, in the assembly directory.
  //
  readDepth    ("/project/huref6/assembly/h6/9-terminator/h6.posmap.frgctg", lowCoverage);
  readVariation("/project/huref6/assembly/h6/9-terminator/h6.posmap.varctg", variation);
  readBadMers  ("/project/huref6/assembly/h6-mer-validation/h6-ms22-allfrags-normalcontigs.badmers.0.singlecontig.zerofrag.badmers", badMers);
#endif

#if 0
  //  HuRef6, ws=25, in the assembly directory.
  //
  readDepth    ("/project/huref6/assembly/h6/9-terminator-ws25/h6.posmap.frgctg", lowCoverage);
  readVariation("/project/huref6/assembly/h6/9-terminator-ws25/h6.posmap.varctg", variation);
  readBadMers  ("/project/huref6/assembly/h6-mer-validation/h6-version4-ws25/h6-ms22-allfrags-normalcontigs.badmers.0.singlecontig.zerofrag.badmers", badMers);
#endif

#if 0
  //  Our scratch huref
  //
  readDepth    ("/project/huref6/redo_consensus-gennady/mer-validation/h6tmp.posmap.frgctg", lowCoverage);
  readVariation("/project/huref6/redo_consensus-gennady/mer-validation/h6tmp.posmap.varctg", variation);
  readBadMers  ("/project/huref6/redo_consensus-gennady/mer-validation/h6tmp-ms22-allfrags-allcontigs.badmers.0.singlecontig.zerofrag.badmers", badMers);
#endif

  u32bit   badBegDepth[1024] = {0};
  u32bit   badEndDepth[1024] = {0};

  u32bit   badDepth[32][32];
  for (u32bit i=0; i<32; i++)
    for (u32bit j=0; j<32; j++)
      badDepth[i][j] = 0;

  map<u64bit,intervalList*>::iterator    it = badMers.begin();
  map<u64bit,intervalList*>::iterator    ed = badMers.end();
  while (it != ed) {
    u64bit         uid        = it->first;

    intervalList  *Iv = variation[uid];
    intervalList  *Ib = badMers[uid];
    intervalList  *Ii = 0L;
    intervalDepth *Id = lowCoverage[uid];

    if (Iv)
      Iv->merge();
    if (Ib)
      Ib->merge();

    if (Iv && Ib) {
      Ii = new intervalList();
      Ii->intersect(*Iv, *Ib);
    }


    if (Ii) {
      u32bit ii = 0;
      u32bit id = 0;

      while ((ii < Ii->numberOfIntervals()) &&
             (id < Id->numberOfIntervals())) {

        //  We want to count the number of times a badmer region
        //  begins/ends in some depth.

        //fprintf(stderr, "testing beg        "u32bitFMT" "u32bitFMT" -- "u32bitFMT" "u32bitFMT"\n",
        //        Ii->lo(ii), Ii->hi(ii), Id->lo(id), Id->hi(id));

        u32bit  beg = 0;
        u32bit  end = 0;

        //  Low points are not allowed to be equal to high points, skip to the next
        while ((id < Id->numberOfIntervals()) &&
               (Id->hi(id) <= Ii->lo(ii))) {
          id++;
          //fprintf(stderr, "testing beg (m)     "u32bitFMT" "u32bitFMT" -- "u32bitFMT" "u32bitFMT"\n",
          //        Ii->lo(ii), Ii->hi(ii), Id->lo(id), Id->hi(id));
        }
        if (id < Id->numberOfIntervals()) {
          u32bit lo = Id->lo(id);
          u32bit hi = Id->hi(id);

          //  Low points are not allowed to be equal to high points.
          if ((lo <= Ii->lo(ii)) && (Ii->lo(ii) < hi)) {
            beg = Id->de(id);
          } else {
            fprintf(stderr, "failed to find begin "u32bitFMT" "u32bitFMT" -- "u32bitFMT" "u32bitFMT" "u32bitFMT"\n",
                    Ii->lo(ii), Ii->hi(ii), Id->lo(id), Id->hi(id), Id->de(id));
            if (id > 0)
              fprintf(stderr, "                     "u32bitFMT" "u32bitFMT" -- "u32bitFMT" "u32bitFMT" "u32bitFMT"\n",
                      Ii->lo(ii), Ii->hi(ii), Id->lo(id-1), Id->hi(id-1), Id->de(id-1));
            //exit(1);
          }
        }

        //fprintf(stderr, "testing end        "u32bitFMT" "u32bitFMT" -- "u32bitFMT" "u32bitFMT"\n",
        //        Ii->lo(ii), Ii->hi(ii), Id->lo(id), Id->hi(id));

        //  High points can be equal.
        while ((id < Id->numberOfIntervals()) &&
               (Id->hi(id) < Ii->hi(ii))) {
          id++;
          //fprintf(stderr, "testing end (m)    "u32bitFMT" "u32bitFMT" -- "u32bitFMT" "u32bitFMT"\n",
          //        Ii->lo(ii), Ii->hi(ii), Id->lo(id), Id->hi(id));
        }
        if (id < Id->numberOfIntervals()) {
          u32bit lo = Id->lo(id);
          u32bit hi = Id->hi(id);

          //  High points aren't allowed to be equal to lo, but can be equal to hi.
          if ((lo < Ii->hi(ii)) && (Ii->hi(ii) <= hi)) {
            end = Id->de(id);
          } else {
            fprintf(stderr, "failed to find end "u32bitFMT" "u32bitFMT" -- "u32bitFMT" "u32bitFMT" "u32bitFMT"\n",
                    Ii->lo(ii), Ii->hi(ii), Id->lo(id), Id->hi(id), Id->de(id));
            if (id > 0)
              fprintf(stderr, "                     "u32bitFMT" "u32bitFMT" -- "u32bitFMT" "u32bitFMT" "u32bitFMT"\n",
                      Ii->lo(ii), Ii->hi(ii), Id->lo(id-1), Id->hi(id-1), Id->de(id-1));
            //exit(1);
          }
        }

        badBegDepth[beg]++;
        badEndDepth[end]++;

        fprintf(stdout, u64bitFMT"\t"u32bitFMT"\t"u32bitFMT"\tdepth="u32bitFMT","u32bitFMT"\n",
                uid, Ii->lo(ii), Ii->hi(ii), beg, end);

        if ((beg < 32) && (end < 32))
          badDepth[beg][end]++;

        ii++;
      }
    }

    it++;
  }

  u32bit bb = 0;
  u32bit be = 0;
  for (u32bit x=0; x<32; x++) {
    fprintf(stdout, u32bitFMT"\t"u32bitFMT"\t"u32bitFMT"\n", x, badBegDepth[x], badEndDepth[x]);
    bb += badBegDepth[x];
    be += badEndDepth[x];
  }
  fprintf(stdout, "total\t"u32bitFMT"\t"u32bitFMT"\n", bb, be);

  for (u32bit i=0; i<30; i++) {
    for (u32bit j=0; j<30; j++)
      fprintf(stdout, u32bitFMTW(5), badDepth[i][j]);
    fprintf(stdout, "\n");
  }

  return(0);
}
