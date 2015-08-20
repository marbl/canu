
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

//  Former Fragment, a sequence loaded into the multialign
//

class abSequence {
public:
  abSequence() {
  };
  ~abSequence() {
  };

  uint32                gkpIdent(void)  { return(_iid);        };
  abSeqID               ident(void)     { return(_lid);        };

  abBeadID              firstBead(void)   { return(_firstbead);  };
  abBaseID              firstBase(void)   { return(_firstbase);  };
  //abColID               firstColumn(void) { needs abacus to get columns };

  //abBeadID              lastBead(void)    { return(_firstbead + _length - 1);  };
  //abBaseID              lastBase(void)    { return(_firstbase + _length - 1);  };
  abBeadID              lastBead(void)    { return(_lastbead);  };
  abBaseID              lastBase(void)    { return(_lastbase);  };
  //abColID               lastColumn(void)  { needs avacus to get columns };

  bool             isRead(void)    { return(_is_read);    };

  uint32           length(void)    { return(_length);     };

  bool             isForward(void) { return(_complement == false); };

  //abMultiAlignID   multiAlignID(void) { return(_multiAlignID); };
  //void             addToMultiAlign(abMultiAlignID maid) {_multiAlignID = maid;};

  void             initialize(uint32 readID,
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
    //_container_iid = -1;
    //_is_contained  = false;
    //_deleted       = false;

    //_multiAlignID  = abMultiAlignID();

    _firstbase.set(firstBase);
    _firstbead.set(firstBead);

    _lastbase.set(firstBase + _length - 1);
    _lastbead.set(firstBead + _length - 1);

    //_n_components  = 0;
    //_components    = -1;
  };

private:
  uint32           _iid;            //  external, aka gatekeeper, ID;
  abSeqID          _lid;            //  index in sequence/quality/fragment store

  bool             _is_read;

  uint32           _length;

  bool             _complement;
  //uint32           _container_iid;  //  if non-zero, the iid of our container
  //bool             _is_contained;   //  if non-zero, consensus detected this fragment is contained
  //bool             _deleted;        //  deleted during consensus

  //  Never really used, was supposed to indicate if the read was present in this multialign.  The
  //  original code initialized this to -1, then reset it to NULL (ApplyAlignment.C:365).  On
  //  printing (PrintAlignment()), and maybe other places, reads with this value != mid (of the
  //  multialign we are printing) are skipped.  The multialign ID is _always_ 0.
  //
  //abMultiAlignID   _multiAlignID;

  abBaseID         _firstbase;      //  global index of first sequence character
  abBeadID         _firstbead;      //  global index of first "bead"

  abBaseID         _lastbase;       //  global index of first sequence character
  abBeadID         _lastbead;       //  global index of first "bead"

  //uint32           _n_components;   //  number of component frags (in case of "unitig" Fragments)
  //uint32           _components;     //  global index of first component frag
};

