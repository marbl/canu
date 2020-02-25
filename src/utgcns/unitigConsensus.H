
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
 *    src/utgcns/libcns/unitigConsensus.H
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
 *    Brian P. Walenz from 2014-NOV-17 to 2015-AUG-05
 *      are Copyright 2014-2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2015-OCT-09
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *    Sergey Koren beginning on 2015-DEC-28
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#ifndef UNITIGCONSENSUS_H
#define UNITIGCONSENSUS_H

#include "runtime.H"

#include "tgStore.H"

class ALNoverlap;
class NDalign;


#define CNS_MIN_QV 0
#define CNS_MAX_QV 60


class abSequence {
public:
  abSequence() {
    _iid        = 0;

    _complement = false;

    _length     = 0;
    _bases      = NULL;
  };

  abSequence(uint32  readID,
             uint32  length,
             char   *seq,
             uint32  complemented);

  ~abSequence() {
    delete [] _bases;
  };


  uint32      seqIdent(void)          { return(_iid);        };

  bool        isForward(void)         { return(_complement == false); };

  uint32      length(void)            { return(_length);     };
  char        getBase(uint32 o)       { return(_bases[o]);   };

  char       *getBases(void)          { return(_bases);      };

private:
  uint32     _iid;            //  external, aka seqStore, ID;

  bool       _complement;

  uint32     _length;
  char      *_bases;
};



class unitigConsensus {
public:
  unitigConsensus(sqStore  *seqStore_,
                  double    errorRate_,
                  double    errorRateMax_,
                  uint32    minOverlap_);
  ~unitigConsensus();

private:
  void   addRead(uint32 readID,
                 uint32 askip, uint32 bskip,
                 bool complemented,
                 map<uint32, sqRead *>     *inPackageRead);

  bool   initialize(map<uint32, sqRead *>     *reads);

public:
  bool   generate(tgTig                     *tig_,
                  char                       algorithm_,
                  char                       aligner_,
                  map<uint32, sqRead *>     *reads_ = NULL);

private:
  void   switchToUncompressedCoordinates(void);

  char  *generateTemplateStitch(void);

  bool   initializeGenerate(tgTig                     *tig,
                            map<uint32, sqRead *>     *reads = NULL);


  bool   generatePBDAG(tgTig                     *tig,
                       char                       aligner,
                       map<uint32, sqRead *>     *reads = NULL);

  bool   generateQuick(tgTig                     *tig,
                       map<uint32, sqRead *>     *reads = NULL);

  bool   generateSingleton(tgTig                     *tig,
                           map<uint32, sqRead *>     *reads = NULL);

  void   findCoordinates(void);
  void   findRawAlignments(void);

public:
  bool   showProgress(void)         { return(_tig->_utgcns_verboseLevel >= 1); };  //  -V          displays which reads are processing
  bool   showAlgorithm(void)        { return(_tig->_utgcns_verboseLevel >= 2); };  //  -V -V       displays some details on the algorithm
  bool   showPlacement(void)        { return(_tig->_utgcns_verboseLevel >= 3); };  //  -V -V -V    displays aligns for placement
  bool   showAlignments(void)       { return(_tig->_utgcns_verboseLevel >= 4); };  //  -V -V -V -V displays aligns and multialigns

  abSequence *getSequence(uint32 id) {
    assert(id < _sequencesLen);
    return(_sequences[id]);
  };

private:
  sqStore        *_seqStore;

  tgTig          *_tig;
  uint32          _numReads;    //  == tig->numberOfChildren()

  uint32          _sequencesMax;
  uint32          _sequencesLen;
  abSequence    **_sequences;

  //  The two positions below are storing the low/high coords for the read.
  //  They do not encode the orientation in the coordinates.
  //
  tgPosition     *_utgpos;      //  Original unitigger location.
  tgPosition     *_cnspos;      //  Actual location in frankenstein.

  uint32          _minOverlap;
  double          _errorRate;
  double          _errorRateMax;
};


#endif
