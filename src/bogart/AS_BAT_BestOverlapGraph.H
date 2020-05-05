
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

#ifndef INCLUDE_AS_BAT_BESTOVERLAPGRAPH
#define INCLUDE_AS_BAT_BESTOVERLAPGRAPH

#include "runtime.H"

#include "AS_BAT_ReadInfo.H"
#include "AS_BAT_OverlapCache.H"

#include <set>
#include <map>
using namespace std;




class ReadEnd {
public:
  ReadEnd(uint32 id=0, bool e3p=false) : _id(id), _e3p(e3p) {};

  uint32  readId(void)  const { return(_id); };
  bool    read3p(void)  const { return(_e3p == true);  };
  bool    read5p(void)  const { return(_e3p == false); };

  bool operator==(ReadEnd const that) const;
  bool operator!=(ReadEnd const that) const;
  bool operator< (ReadEnd const that) const;

private:
  uint32   _id:31;
  uint32   _e3p:1;
};



//  Stores an overlap from an 'a' read (implied by the index into the array of best edges) to a 'b'
//  read.  The hangs are relative to the 'a' read - just as a normal overlap would be.
//
class BestEdgeOverlap {
public:
  BestEdgeOverlap() : _id(0), _e3p(false), _ahang(0), _bhang(0), _evalue(0)   {};

  BestEdgeOverlap(BAToverlap const &ovl) {
    set(ovl);
  };

  void    clear(void);

  void    set(BAToverlap const &olap);
  void    set(uint32 id, bool e3p, int32 ahang, int32 bhang, uint32 evalue);
  void    setReverse(BestEdgeOverlap *that, uint32 readId, bool read3p);

  bool    isValid(void) const { return(_id != 0); }
  bool    isUnset(void) const { return(_id == 0); }

  uint32  readId(void)  const { return(_id); };
  bool    read3p(void)  const { return(_e3p == true);  };
  bool    read5p(void)  const { return(_e3p == false); };

  int32   ahang(void)   const { return(_ahang); };
  int32   bhang(void)   const { return(_bhang); };

  uint32  evalue(void)  const { return(_evalue); };
  double  erate(void)   const { return(AS_OVS_decodeEvalue(_evalue)); };

  bool    operator==(BestEdgeOverlap const that) const { return((_id == that._id) && (_e3p == that._e3p)); };
  bool    operator!=(BestEdgeOverlap const that) const { return((_id != that._id) || (_e3p != that._e3p)); };

private:
  uint32            _id;
  uint64            _e3p    : 1;
  int64             _ahang  : AS_MAX_READLEN_BITS+1;
  int64             _bhang  : AS_MAX_READLEN_BITS+1;
  uint64            _evalue : AS_MAX_EVALUE_BITS;
};

#if (1 + AS_MAX_READLEN_BITS + 1 + AS_MAX_READLEN_BITS + 1 + AS_MAX_EVALUE_BITS > 64)
#error not enough bits to store overlaps.  decrease AS_MAX_EVALUE_BITS or AS_MAX_READLEN_BITS.
#endif



//  A node in the BestOverlapGraph is just a read
//  with two edges out of it, and some flags
//  indicating status of the read.
//
class BestEdgeRead {
public:
  BestEdgeRead() {
    _contained   = false;
    _ignored     = false;
    _coverageGap = false;
    _lopsided5   = false;
    _lopsided3   = false;

    _backbone    = false;
    _spur        = false;
    _bubble      = false;
    _orphan      = false;
    _delinquent  = false;
  };

private:
  BestEdgeOverlap   _best5;
  BestEdgeOverlap   _best3;

  //
  //////////
  //  Contained
  //   - the read has at least one overlap showing it is contained
  //     in some other read.
  //   - contained reads have chnk graph length of zero.
  //   - ignored during spur path detection
  //   - used all over the place to exclude useless reads from various
  //     bits.  should be converted to use backbone instead.
  //
  //////////
  //  Ignored
  //   - the read is flagged as an orphan or a bubble.  ONLY applies
  //     to a BestOverlapGraph constructed from the initial BOG.
  //   - used only to ignore reads when computing a second BOG.  This BOG
  //     is used only for generating a bubble-removed GFA output.
  //
  //////////
  //  CoverageGap
  //   - Probably a chimeric read, but could be a low coverage variant.
  //   - Generally excluded from the assembly, but edges from them are allowed
  //     so they can be possibly popped as bubbles.
  //   - Treated as terminal spur reads when finding spur paths.
  //   - Edges to these should not exist, they cannot seed unitigs, and
  //     will assert() if they're encountered when a unitig is constructed.
  //
  //////////
  //  Lopsided
  //   - Suspected bubble near the end of a read that disrupts all but
  //     short overlaps.  Could also be caused by repeats and low coverage.
  //   - Treated like a normal read, except they cannot seed unitigs.
  //   -
  //
  //////////
  //  Backbone   - Read was placed as part of the backbone of a contig.
  //  Orphan     - Read was placed into a contig as an orphan.
  //  Bubble     - Read can be placed into a contig as a bubble.
  //  Delinquint - Read cannot be plaed either as an orphan or a bubble.
  //

  uint32            _contained   : 1;
  uint32            _ignored     : 1;
  uint32            _coverageGap : 1;
  uint32            _lopsided5   : 1;
  uint32            _lopsided3   : 1;

  uint32            _backbone    : 1;   //  Read was placed as part of the backbone.
  uint32            _spur        : 1;   //  Read is part of a spur path.
  uint32            _bubble      : 1;   //  Read is part of a bubble.
  uint32            _orphan      : 1;   //  Read was placed as an orphan.
  uint32            _delinquent  : 1;   //  unable to plae the read as an orphan.

  friend class BestOverlapGraph;
};




class BestOverlapGraph {
private:
  void   removeReadsWithCoverageGap(const char *prefix, uint32 covGapOlap);
  void   removeLopsidedEdges(const char *prefix, const char *label, double lopsidedDiff);
  uint32 spurDistance(BestEdgeOverlap *edge, uint32 limit, uint32 distance=0);
  void   removeSpannedSpurs(const char *prefix, uint32 spurDepth);

  void   findContains(void);
  void   findEdges(bool redoAll);

  void   findErrorRateThreshold(void);

  void   checkForContainedDovetails(void) const;
  void   checkForCovGapEdges(void) const;

  void   emitGoodOverlaps(const char *prefix, const char *label);

public:
  BestOverlapGraph(double            erateGraph,
                   double            deviationGraph,    double  minOlapPercent,
                   const char       *prefix,
                   bool              filterCoverageGap, uint32  covGapOlap,
                   bool              filterHighError,
                   bool              filterLopsided,    double  lopsidedDiff,
                   bool              filterSpur,        uint32  spurDepth,
                   BestOverlapGraph *BOG = NULL);

  ~BestOverlapGraph() {
    delete [] _reads;
    delete [] _best5score;
    delete [] _best3score;
  };

  bool             bestEdgeExists(uint32 readid, bool threePrime) const {
    return((threePrime) ? (_reads[readid]._best3.readId() != 0) : (_reads[readid]._best5.readId() != 0));
  };

  BestEdgeOverlap *getBestEdgeOverlap(uint32 readid, bool threePrime) const {    //  const on this is cheating, because we can
    return((threePrime) ? (&_reads[readid]._best3) : (&_reads[readid]._best5));        //  still change the BestEdgeOverlap returned!
  };

  BestEdgeOverlap *getBestEdgeOverlap(ReadEnd re) const {
    return((re.read3p()) ? (&_reads[re.readId()]._best3) : (&_reads[re.readId()]._best5));
  };


  bool      isContained  (const uint32 r) const { return( _reads[r]._contained   == 1); };
  bool      isIgnored    (const uint32 r) const { return( _reads[r]._ignored     == 1); };
  bool      isCoverageGap(const uint32 r) const { return( _reads[r]._coverageGap == 1); };
  bool      isLopsided   (const uint32 r) const { return((_reads[r]._lopsided5   == 1) || (_reads[r]._lopsided3 == 1)); };
  bool      isLopsided2  (const uint32 r) const { return((_reads[r]._lopsided5   == 1) && (_reads[r]._lopsided3 == 1)); };
  bool      isBackbone   (const uint32 r) const { return( _reads[r]._backbone    == 1); };
  bool      isSpur       (const uint32 r) const { return( _reads[r]._spur        == 1); };
  bool      isBubble     (const uint32 r) const { return( _reads[r]._bubble      == 1); };
  bool      isOrphan     (const uint32 r) const { return( _reads[r]._orphan      == 1); };
  bool      isDelinquent (const uint32 r) const { return( _reads[r]._delinquent  == 1); };

  void      setContained  (const uint32 r, bool t=true) { _reads[r]._contained   = t; };
  void      setIgnored    (const uint32 r, bool t=true) { _reads[r]._ignored     = t; };
  void      setCoverageGap(const uint32 r, bool t=true) { _reads[r]._coverageGap = t; };
  void      setLopsided5  (const uint32 r, bool t=true) { _reads[r]._lopsided5   = t; };
  void      setLopsided3  (const uint32 r, bool t=true) { _reads[r]._lopsided3   = t; };
  void      setBackbone   (const uint32 r, bool t=true) { _reads[r]._backbone    = t; };
  void      setSpur       (const uint32 r, bool t=true) { _reads[r]._spur        = t; };
  void      setBubble     (const uint32 r, bool t=true) { _reads[r]._bubble      = t; };
  void      setOrphan     (const uint32 r, bool t=true) { _reads[r]._orphan      = t; };
  void      setDelinquent (const uint32 r, bool t=true) { _reads[r]._delinquent  = t; };

  uint32    numContained  (void) const { uint32 n=0;  for (uint32 fi=1; fi <= RI->numReads(); fi++)  if (isContained(fi))   n++;  return(n); };
  uint32    numIgnored    (void) const { uint32 n=0;  for (uint32 fi=1; fi <= RI->numReads(); fi++)  if (isIgnored(fi))     n++;  return(n); };
  uint32    numCoverageGap(void) const { uint32 n=0;  for (uint32 fi=1; fi <= RI->numReads(); fi++)  if (isCoverageGap(fi)) n++;  return(n); };
  uint32    numLopsided   (void) const { uint32 n=0;  for (uint32 fi=1; fi <= RI->numReads(); fi++)  if (isLopsided(fi))    n++;  return(n); };
  uint32    numLopsided2  (void) const { uint32 n=0;  for (uint32 fi=1; fi <= RI->numReads(); fi++)  if (isLopsided2(fi))   n++;  return(n); };
  uint32    numBackbone   (void) const { uint32 n=0;  for (uint32 fi=1; fi <= RI->numReads(); fi++)  if (isBackbone(fi))    n++;  return(n); };
  uint32    numSpur       (void) const { uint32 n=0;  for (uint32 fi=1; fi <= RI->numReads(); fi++)  if (isSpur(fi))        n++;  return(n); };
  uint32    numBubble     (void) const { uint32 n=0;  for (uint32 fi=1; fi <= RI->numReads(); fi++)  if (isBubble(fi))      n++;  return(n); };
  uint32    numOrphan     (void) const { uint32 n=0;  for (uint32 fi=1; fi <= RI->numReads(); fi++)  if (isOrphan(fi))      n++;  return(n); };
  uint32    numDelinquent (void) const { uint32 n=0;  for (uint32 fi=1; fi <= RI->numReads(); fi++)  if (isDelinquent(fi))  n++;  return(n); };

  void      reportEdgeStatistics(const char *prefix, const char *label);
  void      reportBestEdges(const char *prefix, const char *label);
  double    reportErrorLimit() const {return _errorLimit;};

public:
  bool      isOverlapBadQuality(BAToverlap& olap) const;  //  Used in repeat detection

private:
  void      scoreEdge(BAToverlap& olap, bool c5, bool c3);

private:
  BestEdgeRead              *_reads;        //  Nodes in the graph.

  uint64                    *_best5score;   //  Temporary data for computing
  uint64                    *_best3score;   //  best edges.

  double                     _mean;
  double                     _stddev;

  double                     _median;
  double                     _mad;

public:
  double                     _erateGraph;
  double                     _deviationGraph;
private:
  double                     _errorLimit;
  double                     _minOlapPercent;
};



#include "AS_BAT_BestOverlapGraph_implementation.H"



extern BestOverlapGraph *OG;

#endif  //  INCLUDE_AS_BAT_BESTOVERLAPGRAPH
