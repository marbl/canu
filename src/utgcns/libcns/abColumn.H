
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
 *    src/AS_CNS/MultiAlignment_CNS.h
 *    src/AS_CNS/MultiAlignment_CNS_private.H
 *    src/AS_CNS/MultiAlignment_CNS_private.h
 *    src/utgcns/libcns/MultiAlignment_CNS_private.H
 *
 *  Modifications by:
 *
 *    Gennady Denisov from 2005-MAY-23 to 2007-OCT-25
 *      are Copyright 2005-2007 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2005-JUL-08 to 2013-AUG-01
 *      are Copyright 2005-2009,2011,2013 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Eli Venter from 2006-FEB-13 to 2008-FEB-13
 *      are Copyright 2006,2008 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Sergey Koren from 2008-JAN-28 to 2009-SEP-25
 *      are Copyright 2008-2009 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Jason Miller on 2011-SEP-21
 *      are Copyright 2011 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2014-NOV-17 to 2015-JUL-01
 *      are Copyright 2014-2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2015-DEC-18
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#ifndef ABCOLUMN_H
#define ABCOLUMN_H

#undef BASECOUNT  //  Used in abAbacus-refine.C, which is disabled.

#include "abBead.H"

class abAbacus;

class abColumn {
public:
  abColumn() {
    _columnPosition = INT32_MAX;
    _call           = '-';
    _qual           = 0;
    _prevColumn     = NULL;
    _nextColumn     = NULL;
    _beadsMax       = 0;
    _beadsLen       = 0;
    _beads          = NULL;
#if 0
    _beadReadIDs    = NULL;
#endif

#ifdef BASECOUNT
    for (int32 ii=0; ii<CNS_NUM_SYMBOLS; ii++)
      _baseCounts[ii]  = 0;
#endif
  };

  ~abColumn() {
    delete [] _beads;
#if 0
    delete [] _beadReadIDs;
#endif
  };


  int32     &position(void)      { return(_columnPosition); };

  uint32     depth(void)           { return(_beadsLen);      };

#ifdef BASECOUNT
  uint32     baseCount(char c)     { return(_baseCounts[baseToIndex[c]]); };  //  NO CHECKING.
  void       baseCountIncr(char c) { assert(baseToIndex[c] != UINT32_MAX);  _baseCounts[baseToIndex[c]]++; };
  void       baseCountDecr(char c) { assert(baseToIndex[c] != UINT32_MAX);  _baseCounts[baseToIndex[c]]--; };

  //  Return the character with the highest base count, breaking ties deterministically (previous
  //  versions tried to be random).  The last highest count is used - if the mapping is -acgt, then
  //  a tie with gap will never select gap.

  char       majorityBase(bool ignoreGap) {
    uint32 iMax = (ignoreGap == true) ? 1 : 0;  //  Start at entry 1 to ignore gap

    for (int32 i=iMax+1; i<CNS_NUM_SYMBOLS-1; i++)
      if (_baseCounts[i] >= _baseCounts[iMax])
        iMax = i;

    return(indexToBase[iMax]);
  };

  uint32     majorityCount(bool ignoreGap) {
    return(_baseCounts[majorityBase(ignoreGap)]);
  };

  uint32     mismatchCount(void) {
    return(depth() - majorityCount(false));
  }

#else
  uint32     baseCount(char )      { return(0); };  //  Just a stub so abAbacus-refine.C will compile
  char       majorityBase(bool )   { return(0); };
  uint32     mismatchCount(void)   { return(0); };
#endif



  char            baseCall(void)      { return(_call); };
  char            baseQual(void)      { return(_qual); };

  abColumn       *prev(void)          { return(_prevColumn); };
  abColumn       *next(void)          { return(_nextColumn); };


private:
  void            allocateInitialBeads(void);
  void            inferPrevNextBeadPointers(void);

public:
  uint16          insertAtBegin(abColumn *first, uint16 prevLink, char base, uint8 qual);
  uint16          insertAtEnd  (abColumn *prev,  uint16 prevLink, char base, uint8 qual);
  uint16          insertAfter  (abColumn *prev,  uint16 prevLink, char base, uint8 qual);

  uint16          alignBead(uint16 prevIndex, char base, uint8 qual);

  uint16          extendRead(abColumn *column, uint16 beadLink);
  bool            mergeWithNext(abAbacus *abacus, bool highQuality);

private:
  void            baseCallMajority(void);
  void            baseCallQuality(void);
public:
  char            baseCall(bool highQuality);

public:
  abColumn       *mergeNext(void);

private:
  int32            _columnPosition;  //  coordinate of the column in the multialign; needs to be signed for compat with tgTig bgn/end

  char             _call;            //  The base call for this column.
  uint8            _qual;            //  The quality of that base call.

  //  16 bytes of pointers.
  //  Alternate schemes:
  //    two 4 byte offsets into an arary of pointers, but that's also 16 bytes per column.
  //    two 4 byte offsets into an array of objects, but we'd then need to realloc sometime.
  //    two 4 byte offsets into a heap, but that adds a lot of complication.
  //
  abColumn        *_prevColumn;
  abColumn        *_nextColumn;

public:
  abBead           *bead(uint32 ii) { return(_beads + ii); };
private:
  uint16           _beadsMax;   //  Number of beads allocated
  uint16           _beadsLen;   //  Depth; number of reads that span this column
  abBead          *_beads;


  //  If allocated, the read idx (NOT gkpID) for each bead in the column.  This will
  //  be used to (efficiently) map arbitrary columns back to their reads when refining abacus
  //  layouts.  The data is never populated (needs to be done in applyAlignment()).
  //
#if 0
public:
  bool             beadReadIDsExist(void)   { return(_beadReadIDs != NULL); };
  uint32           beadReadID(uint32 xx)    { return(_beadReadIDs[xx]);     };
private:
  uint32          *_beadReadIDs;
#endif

#ifdef BASECOUNT
  uint16           _baseCounts[CNS_NUM_SYMBOLS];
#endif

public:
  void             checkLinks(void) {
    uint32  nErrors = 0;

    for (uint32 tt=0; tt<_beadsLen; tt++) {
      uint16 pl = _beads[tt].prevOffset();
      uint16 nl = _beads[tt].nextOffset();

      if ((prev() != NULL) && (pl != UINT16_MAX) && (prev()->_beads[pl].nextOffset() != tt))
        fprintf(stderr, "ERROR: prev bead pl=%d nextOffset=%d != tt=%d\n",
                pl, prev()->_beads[pl].nextOffset(), tt), nErrors++;

      if ((prev() == NULL) && (pl != UINT16_MAX))
        fprintf(stderr, "no prev column, yet defined link for column %p bead %d\n",
                this, tt), nErrors++;

      if ((next() != NULL) && (nl != UINT16_MAX) && (next()->_beads[nl].prevOffset() != tt))
        fprintf(stderr, "ERROR: next bead pl=%d prevOffset=%d != tt=%d\n",
                nl, next()->_beads[nl].prevOffset(), tt), nErrors++;

      if ((next() == NULL) && (nl != UINT16_MAX))
        fprintf(stderr, "no next column, yet defined link for column %p bead %d\n",
                this, tt), nErrors++;
    }

    if (nErrors == 0)
      return;

    fprintf(stderr, "column %d %p has %d error in links\n", position(), this, nErrors);
    showLinks();

    assert(nErrors == 0);
  };

  void             showLinks(uint16 ppos=UINT16_MAX, uint16 tpos=UINT16_MAX, uint16 npos=UINT16_MAX) {
    uint16   pmax = (_prevColumn == NULL) ? 0 : _prevColumn->_beadsLen;
    uint16   tmax =                                          _beadsLen;
    uint16   nmax = (_nextColumn == NULL) ? 0 : _nextColumn->_beadsLen;

    uint16   lim = max(max(pmax, tmax), nmax);

    fprintf(stderr, "\n");
    for (uint32 xx=0; xx<lim; xx++) {
      fprintf(stderr, "column %d %p link %3d: %5d-%2d-%-5d %c -- %5d-%2d-%-5d %c -- %5d-%2d-%-5d %c\n",
              position(), this, xx,
              (xx < pmax) ? _prevColumn->_beads[xx].prevOffset() : UINT16_MAX,
              xx,
              (xx < pmax) ? _prevColumn->_beads[xx].nextOffset() : UINT16_MAX,
              (xx == ppos) ? '*' : ' ',
              (xx < tmax) ?              _beads[xx].prevOffset() : UINT16_MAX,
              xx,
              (xx < tmax) ?              _beads[xx].nextOffset() : UINT16_MAX,
              (xx == tpos) ? '*' : ' ',
              (xx < nmax) ? _nextColumn->_beads[xx].prevOffset() : UINT16_MAX,
              xx,
              (xx < nmax) ? _nextColumn->_beads[xx].nextOffset() : UINT16_MAX,
              (xx == npos) ? '*' : ' ');
      }
  };


  friend class abAbacus;
  //  friend bool  mergeColumns(abColumn *lcolumn, abColumn *rcolumn);
};

#endif  //  ABCOLUMN_H
