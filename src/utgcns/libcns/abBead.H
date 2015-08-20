
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
 *    Brian P. Walenz beginning on 2014-NOV-17
 *      are Copyright 2014-2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

//  A single bead in the multialign
//
class abBead {
public:
  abBead() {
    clear();
  };
  ~abBead() {
  };

  //  ClearBead
  void        clear(void) {
    boffset      = abBeadID();
    soffset      = abBaseID();
    foffset      = UINT32_MAX;
    prev         = abBeadID();
    next         = abBeadID();
    up           = abBeadID();
    down         = abBeadID();
    frag_index   = abSeqID();
    column_index = abColID();
  };

  abBeadID   const ident(void)    const { return(boffset);      };

  abBeadID        &prevID(void)         { return(prev);         };
  abBeadID        &nextID(void)         { return(next);         };
  abBeadID        &upID(void)           { return(up);           };
  abBeadID        &downID(void)         { return(down);         };

  abBaseID   const baseIdx(void)  const { return(soffset);      };
  abSeqID    const seqIdx(void)   const { return(frag_index);   };
  abColID         &colIdx(void)         { return(column_index); };

private:
  abBeadID   boffset;       //  Location in abAbacus bead array
  abBaseID   soffset;       //  Location in abAbacus bases/quals array
  uint32     foffset;       //  Location in fragment sequence
  abBeadID   prev;          //  |
  abBeadID   next;          //  |
  abBeadID   up;            //  |- Navigation in abAbacus beads
  abBeadID   down;          //  |
  abSeqID    frag_index;    //  Location in abAbacus sequence array
  abColID    column_index;  //  Location in abAbacus column array

  friend class abAbacus;
};

