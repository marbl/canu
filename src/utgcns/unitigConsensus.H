
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' r4587 (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' r1994 (http://kmer.sourceforge.net)
 *
 *  Except as indicated otherwise, this is a 'United States Government Work',
 *  and is released in the public domain.
 *
 *  File 'README.licenses' in the root directory of this distribution
 *  contains full conditions and disclaimers.
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
  void   trimCircular(void);

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
