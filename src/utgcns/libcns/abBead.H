
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
 *    Brian P. Walenz from 2014-NOV-17 to 2015-JAN-27
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

#ifndef ABBEAD_H
#define ABBEAD_H


class abBead {
public:
  abBead() {
    clear();
  };
  ~abBead() {
  };

  void         clear(void) {
    _unused     = 0;
    _isRead     = 0;
    _isUnitig   = 0;
    _base       = '-';
    _qual       = 0;
    _prevOffset = UINT16_MAX;
    _nextOffset = UINT16_MAX;
  };

  void         initialize(char   base,
                         uint8   qual,
                         uint16  prevOff,
                         uint16  nextOff) {
    _base       = base;
    _qual       = qual;
    _prevOffset = prevOff;
    _nextOffset = nextOff;
  };

public:  //  ACCESSORS
  char         base(void)          { return(_base);       };  //  Set only at init, or maybe
  char         qual(void)          { return(_qual);       };  //    reset with swap()

  uint16       prevOffset(void)    { return(_prevOffset); };  //  Generally set only at init,
  uint16       nextOffset(void)    { return(_nextOffset); };  //    columns are added

public:  //  METHODS

private:
  uint16       _unused:1;
  uint16       _isRead:1;     //  If set, base is from an actual read, use it for consensus.
  uint16       _isUnitig:1;   //  If set, base is from a unitig, don't use it for consensus (unless needed).
  uint16       _base:7;       //  Base at this position.  (eventually will be encoded to 3 bits)
  uint16       _qual:6;       //  Quality at this position.

  uint16       _prevOffset;  //  Position in the array of beads for the previous column
  uint16 unused_thisOffset;  //  Position in the arary of beads for this column
  uint16       _nextOffset;  //  Position in the array of beads for the next column

  friend class abAbacus;
  friend class abColumn;  //  inferPrevNextBeadPointers(), mergeWithNext()

  friend void swap(abBead &a, abBead &b);
};


inline
void
swap(abBead &a, abBead &b) {
  abBead c;

  c._isRead = b._isRead; c._isUnitig = b._isUnitig; c._base = b._base; c._qual = b._qual;
  b._isRead = a._isRead; b._isUnitig = a._isUnitig; b._base = a._base; b._qual = a._qual;
  a._isRead = c._isRead; a._isUnitig = c._isUnitig; a._base = c._base; a._qual = c._qual;
}


#endif  //  ABBEAD_H
