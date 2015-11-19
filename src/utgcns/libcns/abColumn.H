
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
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

class abBaseCount {
public:
  abBaseCount() {
    clear();
  };
  ~abBaseCount() {
  };

  void clear(void) {
    _depth    = 0;

    for (uint32 ii=0; ii<CNS_NUM_SYMBOLS; ii++)
      _count[ii] = 0;
  };

  uint32 IncBaseCount(char c) {
    uint32 r = baseToIndex[c];

    if (r > 5)
      fprintf(stderr, "IncBaseCount c=%d '%c' r=%u out of range.\n", c, c, r);
    assert(r <= 5);

    _depth++;

    return(_count[r]++);
  };

  uint32 DecBaseCount(char c) {
    uint32 r = baseToIndex[c];

    if (r > 5)
      fprintf(stderr, "DecBaseCount c=%d '%c' r=%u out of range.\n", c, c, r);
    assert(r <= 5);

    _depth--;

    return(_count[r]--);
  };

  uint32 GetBaseCount(char c) {
    uint32 r = baseToIndex[c];

    if (r > 5)
      fprintf(stderr, "GetBaseCount c=%d '%c' r=%u out of range.\n", c, c, r);
    assert(r <= 5);

    return(_count[r]);
  };



  uint32 GetDepth(void) {
    return(_depth);
  };


  //  Return the character with the highest base count, breaking ties randomly (excluding gaps).
  //  Start at 1 to disallow gaps.
  //
  //  BUG FIX:  a two-way tie was NOT broken, the lower letter always won
  //
  //  BUG:  This wants to return a letter, but doesn't know the mapping
  //
  uint32
  GetMaxBaseCountBase(bool skipGap) {
    uint32 iMax = (skipGap == false) ? 0 : 1;
    uint32 tied = 0;

    for (int32 i=iMax+1; i<CNS_NUM_SYMBOLS-1; i++) {
      if (_count[i] > _count[iMax]) {
        iMax = i;
        tied = 0;
      } else if (_count[i] == _count[iMax]) {
        tied++;
      }
    }

    if (tied > 0) {
      uint32   iTie = 0;

      for (uint32 i=1; i<CNS_NUM_SYMBOLS-1; i++) {
        if (_count[i] == _count[iMax]) {
          uint32 tieBreaker = random();

          if (tieBreaker > iTie) {
            iTie = tieBreaker;
            iMax = i;
          }
        }
      }
    }

    return(indexToBase[iMax]);
  };

private:
  uint32   _depth;
  uint32   _count[CNS_NUM_SYMBOLS];
};




//  A column of beads in the multialign
//
class abColumn {
public:
  abColumn() {
    ma_position = 0;
  };
  ~abColumn() {
  };

  uint32     GetDepth(void)                 { return(base_count.GetDepth());      };
  uint32     GetColumnBaseCount(char c)     { return(base_count.GetBaseCount(c)); };

  void       CheckBaseCounts(abAbacus *ab);

  //  Max over what?  Also needs abacus to get next beads!  Gah!  Two with the same name!
  //#warning GetMaxBaseCount confusion.
  //char       GetMaxBaseCount(abAbacus *ab, int32 start_index);
  char       GetMaxBaseCountBase(bool skipGap) {
    return(base_count.GetMaxBaseCountBase(skipGap));
  };

  abColID    const  ident(void)                    { return(lid);  };

  abBeadID   const  callID(void)                   { return(call); };

  abColID          &prevID(void)                   { return(prev); };
  abColID          &nextID(void)                   { return(next); };

  int32             mismatch(void) {
    return(base_count.GetDepth() - base_count.GetBaseCount( base_count.GetMaxBaseCountBase(false) ));
  };

  int32            &position(void)                 { return(ma_position); };

private:
  abColID          lid;           //  index in columnStore
  abBeadID         call;          //  global offset in beadStore;
  abColID          prev;          //  navigation in columnStore;
  abColID          next;          //
  abMultiAlignID   ma_id;         //  abMultiAlign this column belongs to
  int32            ma_position;   //  was ma_index; column ordinal, needs to be signed for compat with tgTig bgn/end
  abBaseCount      base_count;

  friend class abAbacus;
};




