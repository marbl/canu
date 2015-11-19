
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
 *    Brian P. Walenz from 2014-NOV-17 to 2015-JAN-31
 *      are Copyright 2014-2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

//  Iterate over all beads in a sequence
class abSeqBeadIterator {
public:
  abSeqBeadIterator(abAbacus *ab, abSequence *seq) {  //  used to take a fid, now needs _sequences + fid
    assert(seq != NULL);
    _abacus   = ab;
    _seq      = seq;
    _beadID   = _seq->firstBead();
    _isNull   = false;
  };
  ~abSeqBeadIterator() {
  };


  bool     isNull(void) {
    return(_isNull);
  };


  void     nullify(void) {
    _seq    = NULL;  //  used to be the 0th sequence
    _beadID = abBeadID();
    _isNull = true;
  };


  abBeadID  next(void) {
    assert(_isNull == false);

    abBeadID nid = _beadID;

    if (_beadID.isValid())
      _beadID = _abacus->getBead(_beadID)->nextID();

    return(nid);
  };


  abBeadID  prev(void) {
    assert(_isNull == false);

    abBeadID nid = _beadID;

    if (_beadID.isValid())
      _beadID      = _abacus->getBead(_beadID)->prevID();

    return(nid);
  };


private:
  abAbacus      *_abacus;
  abSequence    *_seq;
  abBeadID       _beadID;

  bool           _isNull;

  friend class abAbacus;
};




//  Iterate over all beads in a column
class abColBeadIterator {
public:
  abColBeadIterator(abAbacus *ab, abColumn *col) {
    assert(col != NULL);
    _abacus = ab;
    _col    = col;
    _beadID = _col->callID();
  };
  ~abColBeadIterator() {
  };

  abBeadID  next(void) {
    abBeadID nid;

    if (_beadID.isValid()) {
      abBead *bead = _abacus->getBead(_beadID);
      nid          = bead->downID();
      _beadID      = nid;
    }

    return(nid);
  }

private:
  abAbacus      *_abacus;
  abColumn      *_col;
  abBeadID       _beadID;

  friend class abAbacus;
};


