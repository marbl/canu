
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

//  Former Fragment, a sequence loaded into the multialign
//

class abSequence {
public:
  abSequence() {
  };
  ~abSequence() {
  };

  uint32                gkpIdent(void)          { return(_iid);        };
  abSeqID               ident(void)             { return(_lid);        };

  abBeadID              firstBead(void)         { return(_firstbead);  };
  abBaseID              firstBase(void)         { return(_firstbase);  };

  abBeadID              lastBead(void)          { return(_lastbead);   };
  abBaseID              lastBase(void)          { return(_lastbase);   };

  abBeadID              getBead(uint32 offset)  { return(abBeadID(_firstbead.get() + offset)); };

  bool                  isRead(void)            { return(_is_read);    };

  uint32                length(void)            { return(_length);     };

  bool                  isForward(void)         { return(_complement == false); };

  void                  initialize(uint32 readID,
                                   uint32 readIdx,
                                   uint32 length,
                                   uint32 complemented,
                                   uint32 firstBase,
                                   uint32 firstBead) {
    _iid           = readID;
    _lid.set(readIdx);                 // index in sequence/quality/fragment store

    _is_read       = true;

    _length        = length;

    _complement    = complemented;

    _firstbase.set(firstBase);
    _firstbead.set(firstBead);

    _lastbase.set(firstBase + _length - 1);
    _lastbead.set(firstBead + _length - 1);
  };

private:
  uint32           _iid;            //  external, aka gatekeeper, ID;
  abSeqID          _lid;            //  index in sequence/quality/fragment store

  bool             _is_read;

  uint32           _length;

  bool             _complement;

  abBaseID         _firstbase;      //  global index of first sequence character
  abBeadID         _firstbead;      //  global index of first "bead"

  abBaseID         _lastbase;       //  global index of first sequence character
  abBeadID         _lastbead;       //  global index of first "bead"
};

