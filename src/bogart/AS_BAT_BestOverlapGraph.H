
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
 *    src/AS_BAT/AS_BAT_BestOverlapGraph.H
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2010-NOV-23 to 2014-JAN-29
 *      are Copyright 2010-2014 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2014-OCT-21 to 2015-JUN-03
 *      are Copyright 2014-2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2016-JAN-11
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *    Sergey Koren beginning on 2016-MAR-11
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#ifndef INCLUDE_AS_BAT_BESTOVERLAPGRAPH
#define INCLUDE_AS_BAT_BESTOVERLAPGRAPH

#include "AS_global.H"
#include "AS_BAT_OverlapCache.H"

class ReadEnd {
public:
  ReadEnd() {
    _id  = 0;
    _e3p = false;
  };
  ReadEnd(uint32 id, bool e3p) {
    _id  = id;
    _e3p = e3p;
  };

  uint32  readId(void)  const { return(_id); };
  bool    read3p(void)  const { return(_e3p == true);  };
  bool    read5p(void)  const { return(_e3p == false); };

  bool operator==(ReadEnd const that) const {
    return((readId() == that.readId()) && (read3p() == that.read3p()));
  };

  bool operator!=(ReadEnd const that) const {
    return((readId() != that.readId()) || (read3p() != that.read3p()));
  };

  bool operator<(ReadEnd const that) const {
    if (readId() != that.readId())
      return readId() < that.readId();
    else
      return read3p() < that.read3p();
  };

private:
  uint32   _id:31;
  uint32   _e3p:1;
};



//  Stores an overlap from an 'a' read (implied by the index into the array of best edges) to a 'b'
//  read.  The hangs are relative to the 'a' read - just as a normal overlap would be.
//
class BestEdgeOverlap {
public:
  BestEdgeOverlap()                          { clear();  };
  BestEdgeOverlap(BAToverlap const &ovl)     { set(ovl); };
  ~BestEdgeOverlap()                         {           };

  void    clear(void) {
    _id     = 0;
    _e3p    = 0;
    _ahang  = 0;
    _bhang  = 0;
    _evalue = 0;
  };

  void    set(BAToverlap const &olap) {
    _id     = olap.b_iid;
    if (((olap.a_hang <= 0) && (olap.b_hang >= 0)) ||    //  If contained, _e3p just means
        ((olap.a_hang >= 0) && (olap.b_hang <= 0)))      //  the other read is flipped
      _e3p  = olap.flipped;
    else                                                 //  Otherwise, means the
      _e3p    = olap.BEndIs3prime();                     //    olap is to the 3' end

    _ahang  = olap.a_hang;
    _bhang  = olap.b_hang;
    _evalue = olap.evalue;
  };

  void    set(uint32 id, bool e3p, int32 ahang, int32 bhang, uint32 evalue) {
    _id     = id;
    _e3p    = e3p;
    _ahang  = ahang;
    _bhang  = bhang;
    _evalue = evalue;
  };


  uint32  readId(void)  const { return(_id); };
  bool    read3p(void)  const { return(_e3p == true);  };
  bool    read5p(void)  const { return(_e3p == false); };

  int32   ahang(void)   const { return(_ahang); };
  int32   bhang(void)   const { return(_bhang); };

  uint32  evalue(void)  const { return(_evalue); };
  double  erate(void)   const { return(AS_OVS_decodeEvalue(_evalue)); };

private:
  uint32            _id;
  uint64            _e3p    : 1;    //  Overlap with the 3' end of that read, or flipped contain
  int64             _ahang  : AS_MAX_READLEN_BITS+1;
  int64             _bhang  : AS_MAX_READLEN_BITS+1;
  uint64            _evalue : AS_MAX_EVALUE_BITS;
};

#if (1 + AS_MAX_READLEN_BITS + 1 + AS_MAX_READLEN_BITS + 1 + AS_MAX_EVALUE_BITS > 64)
#error not enough bits to store overlaps.  decrease AS_MAX_EVALUE_BITS or AS_MAX_READLEN_BITS.
#endif









class BestOverlaps {
public:
  BestEdgeOverlap     _best5;
  BestEdgeOverlap     _best3;
  uint32              _isC;
};


class BestScores {
public:
  BestScores() {
    _best5score = 0;
    _best3score = 0;
    _isC        = 0;
  };

  uint64  _best5score;
  uint64  _best3score;
  uint32  _isC;
};



class BestOverlapGraph {
private:
  void   removeSuspicious(const char *prefix);
  void   removeSpurs(const char *prefix);
  void   removeLopsidedEdges(const char *prefix);

  void   findEdges(void);

  void   removeHighErrorBestEdges(void);

  void   removeContainedDovetails(void);

public:
  BestOverlapGraph(double        erateGraph,
                   double        deviationGraph,
                   const char   *prefix,
                   bool          filterSuspicious,
                   bool          filterHighError,
                   bool          filterLopsided,
                   bool          filterSpur);

  ~BestOverlapGraph() {
    delete [] _bestA;
    delete [] _scorA;
  };

  //  Given a read UINT32 and which end, returns pointer to
  //  BestOverlap node.
  BestEdgeOverlap *getBestEdgeOverlap(uint32 readid, bool threePrime) {
    if (_bestA)
      return((threePrime) ? (&_bestA[readid]._best3) : (&_bestA[readid]._best5));
    return((threePrime) ? (&_bestM[readid]._best3) : (&_bestM[readid]._best5));
  };

  // given a ReadEnd sets it to the next ReadEnd after following the
  // best edge
  ReadEnd   followOverlap(ReadEnd end) {
    if (end.readId() == 0)
      return(ReadEnd());

    BestEdgeOverlap *edge = getBestEdgeOverlap(end.readId(), end.read3p());

    return(ReadEnd(edge->readId(), !edge->read3p()));
  };

  void setContained(const uint32 readid) {
    if (_bestA)
      _bestA[readid]._isC = true;
    else
      _bestM[readid]._isC = true;
  };

  bool isContained(const uint32 readid) {
    if (_bestA)
      return(_bestA[readid]._isC);
    return(_bestM[readid]._isC);
  };

  bool isSuspicious(const uint32 readid) {
    return(_suspicious.count(readid) > 0);
  };

  void      reportEdgeStatistics(const char *prefix, const char *label);
  void      reportBestEdges(const char *prefix, const char *label);

public:
  bool     isOverlapBadQuality(BAToverlap& olap);  //  Used in repeat detection
private:
  uint64   scoreOverlap(BAToverlap& olap);

private:
  void     scoreContainment(BAToverlap& olap);
  void     scoreEdge(BAToverlap& olap);

private:
  uint64  &best5score(uint32 id) {
    if (_restrictEnabled == false)
      return(_scorA[id]._best5score);
    return(_scorM[id]._best5score);
  };

  uint64  &best3score(uint32 id) {
    if (_restrictEnabled == false)
      return(_scorA[id]._best3score);
    return(_scorM[id]._best3score);
  };

private:
  BestOverlaps              *_bestA;
  BestScores                *_scorA;

  double                     _mean;
  double                     _stddev;

  double                     _median;
  double                     _mad;

  uint32                     _nSuspicious;          //  Stats for output
  uint32                     _n1EdgeFiltered;
  uint32                     _n2EdgeFiltered;
  uint32                     _n1EdgeIncompatible;
  uint32                     _n2EdgeIncompatible;

  set<uint32>                _suspicious;
  set<uint32>                _spur;

  map<uint32, BestOverlaps>  _bestM;
  map<uint32, BestScores>    _scorM;

  //  These restrict the best overlap graph to a set of reads, instead of all reads.
  //  Currently (Aug 2016) unused.  There used to be a constructor that would take
  //  a set(uint32) of reads we cared about, but it was quite stale and was removed.
private:
  bool     isOverlapRestricted(const BAToverlap &olap);
private:
  set<uint32>               *_restrict;
  bool                       _restrictEnabled;


public:
  double                     _erateGraph;
  double                     _deviationGraph;
private:
  double                     _errorLimit;
}; //BestOverlapGraph



extern BestOverlapGraph *OG;

#endif  //  INCLUDE_AS_BAT_BESTOVERLAPGRAPH
