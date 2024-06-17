
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


#include "tgStore.H"

#include <map>
#include <set>

class ALNoverlap;
class NDalign;


#define CNS_MIN_QV 0
#define CNS_MAX_QV 60


class abSequence {
public:
  abSequence(uint32  readID       = 0,
             uint32  length       = 0,
             char   *seq          = nullptr,
             bool    isReverse    = false);
  ~abSequence();

  uint32      seqIdent(void)          { return _iid;          }

  uint32      length(void)            { return _length;       }
  char        getBase(uint32 o)       { return _bases[o];     }

  char       *getBases(void)          { return _bases;        }

private:
  uint32     _iid = 0;            //  external, aka seqStore, ID;

  uint32     _length       = 0;
  char      *_bases        = nullptr;
};


//typedef std::map<uint32, sqRead *>  u32toRead;

class unitigConsensus {
public:
  unitigConsensus(sqStore  *seqStore_,
                  double    errorRate_,
                  double    errorRateMax_,
                  uint32    errorRateMaxID_,
                  uint32    minOverlap_,
                  uint32    minCoverage_);
  ~unitigConsensus();

private:
  void   addRead(uint32     readID,
                 uint32     askip,
                 uint32     bskip,
                 bool       complemented,
                 u32toRead &reads);

public:
  bool   generate(tgTig      *tig_,
                  char        algorithm_,
                  char        aligner_,
                  u32toRead  &reads_);

private:
  void   switchToUncompressedCoordinates(void);
  void   updateReadPositions(void);
                 
  char  *generateTemplateStitch(void);

  bool   generatePBDAG     (tgTig *tig, char aligner, u32toRead &reads);
  bool   generateQuick     (tgTig *tig,               u32toRead &reads);
  bool   generateSingleton (tgTig *tig,               u32toRead &reads);

  void   adjustPosition(tgPosition   utgpos,
                        tgPosition   cnspos,
                        tgPosition  &adjusted,
                        bool         isS);

  void   findCoordinates(char algorithm_, u32toRead &reads_);
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
  sqStore        *_seqStore = nullptr;

  tgTig          *_tig = nullptr;

  uint32          _numReads       = 0;   //  == tig->numberOfChildren()
  uint32          _numReadsUsable = 0;   //  == number not 'ont'

  uint32          _sequencesMax = 0;
  uint32          _sequencesLen = 0;
  abSequence    **_sequences    = nullptr;

  //  The two positions below are storing the low/high coords for the read.
  //  They do not encode the orientation in the coordinates.
  //
  tgPosition     *_utgpos = nullptr;     //  Original unitigger location.
  tgPosition     *_cnspos = nullptr;     //  Actual location in frankenstein.
  tgPosition     *_adjpos = nullptr;

  uint32         *_templateToCNS  = nullptr;
  uint32          _templateLength = 0;

  uint32          _minOverlap     = 0;
  uint32          _minCoverage    = 0;
  double          _errorRate      = 0;
  double          _errorRateMax   = 0;
  uint32          _errorRateMaxID = 0;
};


#endif
