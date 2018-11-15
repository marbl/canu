
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

#include "AS_global.H"

#include "tgStore.H"

class ALNoverlap;
class NDalign;


#define CNS_MIN_QV 0
#define CNS_MAX_QV 60


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
             uint8  *qlt,
             uint32  complemented);

  ~abSequence() {
    delete [] _bases;
    delete [] _quals;
  };


  uint32      seqIdent(void)          { return(_iid);        };

  bool        isRead(void)            { return(_is_read);    };
  bool        isForward(void)         { return(_complement == false); };

  uint32      length(void)            { return(_length);     };
  char        getBase(uint32 o)       { return(_bases[o]);   };
  uint8       getQual(uint32 o)       { return(_quals[o]);   };

  char       *getBases(void)          { return(_bases);      };

private:
  uint32     _iid;            //  external, aka seqStore, ID;

  bool       _is_read;
  bool       _complement;

  uint32     _length;
  char      *_bases;
  uint8     *_quals;
};



class unitigConsensus {
public:
  unitigConsensus(sqStore  *seqStore_,
                  double    errorRate_,
                  double    errorRateMax_,
                  uint32    minOverlap_);
  ~unitigConsensus();

  bool   generate(tgTig                     *tig,
                  char                       algorithm,
                  char                       aligner,
                  map<uint32, sqRead *>     *reads = NULL,
                  map<uint32, sqReadData *> *datas = NULL);

private:
  char *
  generateTemplateStitch(tgPosition  *utgpos,
                         uint32       numfrags,
                         double       errorRate,
                         bool         verbose);

  bool   generatePBDAG(tgTig                     *tig,
                       char                       aligner,
                       map<uint32, sqRead *>     *reads = NULL,
                       map<uint32, sqReadData *> *datas = NULL);

  bool   generateQuick(tgTig                     *tig,
                       map<uint32, sqRead *>     *reads = NULL,
                       map<uint32, sqReadData *> *datas = NULL);

  bool   generateSingleton(tgTig                     *tig,
                           map<uint32, sqRead *>     *reads = NULL,
                           map<uint32, sqReadData *> *datas = NULL);

private:
  void   addRead(sqStore *seqStore,
                 uint32 readID,
                 uint32 askip, uint32 bskip,
                 bool complemented,
                 map<uint32, sqRead *>     *inPackageRead,
                 map<uint32, sqReadData *> *inPackageReadData);

public:
  int32  initialize(map<uint32, sqRead *>     *reads,
                    map<uint32, sqReadData *> *datas);

#if 0
  void   setErrorRate(double errorRate_)   { errorRate  = errorRate_;  };
  void   setMinOverlap(uint32 minOverlap_) { minOverlap = minOverlap_; };
#endif

  //bool   showProgress(void)         { return(tig->_utgcns_verboseLevel >= 1); };  //  -V          displays which reads are processing
  bool   showAlgorithm(void)        { return(tig->_utgcns_verboseLevel >= 2); };  //  -V -V       displays some details on the algorithm
  //bool   showPlacementBefore(void)  { return(tig->_utgcns_verboseLevel >= 3); };  //  -V -V -V    displays placement info before each read
  //bool   showPlacement(void)        { return(tig->_utgcns_verboseLevel >= 3); };  //
  //bool   showAlignments(void)       { return(tig->_utgcns_verboseLevel >= 4); };  //  -V -V -V -V displays aligns and multialigns
  //bool   showMultiAlignments(void)  { return(tig->_utgcns_verboseLevel >= 4); };  //

#warning NEEDS GUARDS
  abSequence *getSequence(uint32 id) { return(_sequences[id]); };

  void   generateConsensus(tgTig *tig);

private:
  sqStore        *seqStore;

  tgTig          *tig;
  uint32          numfrags;    //  == tig->numberOfChildren()

  uint32          _sequencesMax;
  uint32          _sequencesLen;
  abSequence    **_sequences;

  //  The two positions below are storing the low/high coords for the read.
  //  They do not encode the orientation in the coordinates.
  //
  tgPosition     *utgpos;      //  Original unitigger location.
  tgPosition     *cnspos;      //  Actual location in frankenstein.

  uint32          minOverlap;
  double          errorRate;
  double          errorRateMax;
};


#endif
