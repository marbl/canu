
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
 *    Brian P. Walenz beginning on 2015-NOV-23
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#ifndef ABSEQUENCE_H
#define ABSEQUENCE_H

#include "abColumn.H"

class abSequence {
public:
  abSequence() {
    _iid        = 0;

    _is_read    = false;
    _complement = false;

    _length     = 0;
    _bases      = NULL;
    _quals      = NULL;
  };

  abSequence(uint32  readID,
             uint32  length,
             char   *seq,
             char   *qlt,
             uint32  complemented);

  ~abSequence() {
    delete [] _bases;
    delete [] _quals;
  };


  uint32                gkpIdent(void)          { return(_iid);        };

  bool                  isRead(void)            { return(_is_read);    };
  bool                  isForward(void)         { return(_complement == false); };

  uint32                length(void)            { return(_length);     };
  char                  getBase(uint32 o)       { return(_bases[o]);   };
  uint8                 getQual(uint32 o)       { return(_quals[o]);   };

  char                 *getBases(void)          { return(_bases);      };

private:
  uint32           _iid;            //  external, aka gatekeeper, ID;

  bool             _is_read;
  bool             _complement;

  //  Until the read is added to the abacus, we need to store the correctly
  //  oriented bases/quals somewhere.

  uint32           _length;
  char            *_bases;
  uint8           *_quals;
};

#endif  //  ABSEQUENCE_H
