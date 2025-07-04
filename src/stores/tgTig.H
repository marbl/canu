
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

#ifndef TG_TIG_H
#define TG_TIG_H

#include "sqStore.H"
#include "bits.H"

#include <map>

//  Used in unitigConsensus and tgTig::dumpBAM.
using u32toRead = std::map<uint32, sqRead *>;

//  This stupid enum.  It used to be a legacy typedef, but gcc 9.2 in
//  holy-build-box started complaining about it (December 2024) not fitting
//  in the bitfield in tgTigRecord.  BPW switched it to the below
//  forced-unsigned but that didn't help; gcc 9.2 is just broken and needs
//  -fpermissive; see utility/scripts/build-release.sh.
//
enum tgTig_class : uint32 {
  tgTig_noclass      = 0x00,    //  0x02 was for bubbles, but that was moved
  tgTig_unassembled  = 0x01,    //  to be only a suggestion.
  tgTig_contig       = 0x03     //  tgTig_contig MUST be 3 for compatibility.
};

static
const
char *
toString(tgTig_class c) {
  switch (c) {
    case tgTig_noclass:      return("unsetc");  break;
    case tgTig_unassembled:  return("unassm");  break;
    case tgTig_contig:       return("contig");  break;
  }

  return("undefined-class");
}

static
tgTig_class
toClass(char const *c) {
  if (strcmp(c, "unsetc") == 0)
    return(tgTig_noclass);
  if (strcmp(c, "unassm") == 0)
    return(tgTig_unassembled);
  if (strcmp(c, "contig") == 0)
    return(tgTig_contig);

  fprintf(stderr, "WARNING:  tiTig_class toClass('%s') is not a valid class.\n", c);

  return(tgTig_noclass);
}




//  Info about the placement of an object in a tig.  For unitigs, this
//  will be just reads.  For contigs, this could be unitigs and reads.
//
class tgPosition {
public:
  void                 initialize(void);

  //  Accessors

  uint32               ident(void)     const { return(_objID); };

  bool                 isRead(void)    const { return(_isRead   == true); };
  bool                 isUnitig(void)  const { return(_isUnitig == true); };
  bool                 isContig(void)  const { return(_isContig == true); };

  //  These are not really used anywhere, just in tgTig::dumpLayout.
  //  See comments near the definition of _anchor below.
  uint32               anchor(void)          { return(_anchor); };
  int32                aHang(void)           { return(_ahang);  };
  int32                bHang(void)           { return(_bhang);  };

  bool                 isForward(void) const { return(_isReverse == false); };
  bool                 isReverse(void) const { return(_isReverse == true);  };

  //  We want to ignore some reads during consensus, like excess contained reads
  //  and ONT in verkko.  The flag is antisense for compatibility with existing
  //  Canu stores - those will have the _skipCNS bit set to 0, so if our
  //  flag was _useInCNS, no reads would be used!
  bool                 skipConsensus(void) const { return _skipCNS;      }
  void                 skipConsensus(bool v)     { _skipCNS = v;         }

  bool                 isLowQuality(void) const  { return _isLowQuality; }
  void                 isLowQuality(bool v)      { _isLowQuality = v;    }

  //  Position in the parent, both oriented (bgn/end) and unoriented (min/max).

  int32                bgn(void) const       { return((_isReverse == false) ? _min : _max); };
  int32                end(void) const       { return((_isReverse == false) ? _max : _min); };

  int32                min(void) const       { return(_min); };
  int32                max(void) const       { return(_max); };

  //  Amount of this object to ignore; e.g., trim from the ends.

  int32                askip(void) const       { return(_askip); };
  int32                bskip(void) const       { return(_bskip); };

  //  Delta values for the alignment to the parent.

  uint32               deltaOffset(void)     { return(_deltaOffset); };
  uint32               deltaLength(void)     { return(_deltaLen);    };

  uint32               getMAPQ(void)         {
    if (_numPlacement <= 1) {
        return 60;  // high cofidence
    }
    double p_err = 1.0 - (1.0 / _numPlacement);
    return int(floor(-10.0 * log10(p_err)));
  }

  //  Set everything.  This is to be used by unitigger.
  void                 set(uint32 id, uint32 anchor, int32 ahang, int32 bhang, int32  bgn, int32  end, uint32 askip = 0, uint32 bskip = 0, uint32 numPlacement = 1) {

    _objID         = id;

    _isRead        = true;
    _isUnitig      = false;
    _isContig      = false;
    //kipCNS       = _skipCNS;       //  Leave unchanged!
    //sLowQuality  = _isLowQuality;  //  Leave unchanged!

    _anchor        = anchor;
    _ahang         = ahang;
    _bhang         = bhang;

    _askip         = askip;
    _bskip         = bskip;

    _numPlacement = (numPlacement >= 32 ? 31 : numPlacement);

    _min           = std::min(bgn, end);
    _max           = std::max(bgn, end);
    _isReverse     = (bgn >= end);

    _deltaOffset   = 0;
    _deltaLen      = 0;
  };

  //  Reverse-complement the read.
  //  Anchor and hangs are invalidated.
  void                 reverseComplement(int32 len) {
    int32  min = len - _max;
    int32  max = len - _min;

    assert(_min < _max);
    assert( min <  max);

    _anchor    = 0;
    _ahang     = 0;    _askip = 0;    _deltaOffset = 0;
    _bhang     = 0;    _bskip = 0;    _deltaLen    = 0;

    _min       = min;
    _max       = max;

    _isReverse = !_isReverse;
  }

  //  Set the coords, ignoring orientation.
  void                 setMinMax(int32 min, int32 max) {
    _min = min;
    _max = max;
  }

  //  Sort by position, however, if:
  //   - I'm the anchor of the other tig, I must be before it.
  //   - My anchor is the other tig, I must be after it.
  //
  bool                 operator<(tgPosition const &that) const {
    if      (min()   != that.min())   return min()   < that.min();
    else if (_objID  == that._anchor) return true;
    else if (_anchor == that._objID)  return false;
    else if (max()   != that.max())   return max()   > that.max();
    else                              return ident() < that.ident();
  }

private:
public:
  uint32              _objID;             //  ID of this object

  uint32              _isRead       : 1;  //  A full length read alignment
  uint32              _isUnitig     : 1;  //
  uint32              _isContig     : 1;  //

  uint32              _isReverse    : 1;  //  Child is oriented forward relative to parent, used during consensus.

  uint32              _skipCNS      : 1;  //  Child should not be used for consensus.
  uint32              _isLowQuality : 1;  //  Child should not be used for consensus.

  uint32              _numPlacement : 5;  // Number of times this read is placed in the assembly, used to compute MAPQ in bam output, canu always sets this to 1
  uint32              _spare        : 21;

  //  These are not used any more; they were used by the early consensus
  //  algorithm; they are still set in bogart, and MUST remain in the structure
  //  for binary compatibility.
  uint32              _anchor;           //  ID of the like-type object we align to
  int32               _ahang;            //  Placement relative to anchor object
  int32               _bhang;            //

  int32               _askip;            //  Amount of sequence to not align on each end
  int32               _bskip;            //  Used in correction and verkko consensus.

  //  Must be signed, utgcns can push reads negative.
  int32               _min;
  int32               _max;

  uint32              _deltaOffset;      //  Our delta alignment to the parent object consensus
  uint32              _deltaLen;
};



//  On-disk tig, same as tgTig without the pointers.
struct tgTigRecordV1 {
  uint32              _tigID;

  uint32              _unused32;
  uint32              _trimBgn;
  uint32              _trimEnd;

  tgTig_class         _class           : 2;
  uint32              _suggestRepeat   : 1;
  uint32              _suggestCircular : 1;
  uint32              _suggestBubble   : 1;
  uint32              _circularLength  : 20;
  uint32              _suggestNoTrim   : 1;
  uint32              _spare           : 6;

  uint32              _layoutLen;
  uint32              _basesLen;
  uint32              _childrenLen;
  uint32              _childDeltaBitsLen;
};


//  On-disk tig, same as tgTig without the pointers.
//  Version 2 was introduced in Canu v2.3 to handle CIGAR strings
//  and BAM output.  Lengths were increased to 64-bit.
//
struct tgTigRecord {
  uint32              _tigID;

  tgTig_class         _class           : 2;
  uint32              _suggestRepeat   : 1;
  uint32              _suggestCircular : 1;
  uint32              _suggestBubble   : 1;
  uint32              _circularLength  : 20;
  uint32              _suggestNoTrim   : 1;
  uint32              _spare           : 6;

  uint32              _trimBgn;
  uint32              _trimEnd;

  uint64              _layoutLen;
  uint64              _basesLen;
  uint64              _childrenLen;
  uint64              _childDeltaBitsLen;
  uint64              _childCIGARLen;

  tgTigRecord&
  operator=(const tgTigRecordV1 &that) {
    _tigID             = that._tigID;

    _class             = that._class;
    _suggestRepeat     = that._suggestRepeat;
    _suggestCircular   = that._suggestCircular;
    _suggestBubble     = that._suggestBubble;
    _circularLength    = that._circularLength;
    _suggestNoTrim     = that._suggestNoTrim;
    _spare             = that._spare;
    _trimBgn           = that._trimBgn;
    _trimEnd           = that._trimEnd;

    _layoutLen         = that._layoutLen;
    _basesLen          = that._basesLen;
    _childrenLen       = that._childrenLen;
    _childDeltaBitsLen = that._childDeltaBitsLen;
    _childCIGARLen     = 0;

    return *this;
  }
};


//  In-core tig.
class tgTig {
public:
  tgTig()   {          }   //  All data unallocated, lengths set to zero
  ~tgTig()  { clear(); }   //  Releases memory

  void           clear(void);    //  Releases memory and clears data.

  //  Accessors

  uint32         tigID(void)                 { return(_tigID); };

  bool           consensusExists(void)       { return(_basesLen > 0); };
  double         computeCoverage(void);

  //  Should be private, but used in tgStoreLoad.C for testing cns files.
public:
  uint32         length(void)                { return(_layoutLen); };
  char          *bases(void)                 { return(_bases);     };
  uint8         *quals(void)                 { return(_quals);     };

public:
  uint32         numberOfChildren(void)      {                            return(_childrenLen); };
  tgPosition    *getChild(uint32 c)          { assert(c < _childrenLen);  return(_children + c);  };
  tgPosition    *addChild(void)              {
    increaseArray(_children, _childrenLen, _childrenMax, std::min(_childrenMax/4, 16u));

    _children[_childrenLen].initialize();

    return(_children + _childrenLen++);
  };

  void           allocateChildren(uint32 m) {
    resizeArray(_children, _childrenLen, _childrenMax, m, _raAct::doNothing);

    for (uint32 ii=0; ii<m; ii++)
      _children[ii].initialize();
  };

  //  Operators

  tgTig         &operator=(tgTig &tg);

  bool           loadFromStreamOrLayout(FILE *F);

private:
  void           sumCIGAR(void);
  void           writeCIGAR(writeBuffer *B, FILE *);
  void           loadCIGAR(readBuffer *B, FILE *F);
public:
  char const    *getChildCIGAR(uint32 child) {
    return (_childCIGAR == nullptr) ? nullptr : _childCIGAR[child];
  }

public:
  void           saveToBuffer(writeBuffer *B);
  void           saveToStream(FILE *F);

  bool           loadFromBuffer(readBuffer *B);
  bool           loadFromStream(FILE *F);

  void           dumpLayout(FILE *F, bool withSequence=true);
  bool           loadLayout(FILE *F);

  void           filterContains(double maxCov, bool enableStash);
  void           unstashContains(void);

  void           dropExcessCoverage(double maxCov) {
    filterContains(maxCov, true);
    clearStash();
  };

  void           cleanup(void);   //  Sort reads, make sure the start at zero, set length.

  //  Save and load a package of data needed to process this tig.

  void           exportData(writeBuffer                *exportDataFile,
                            sqStore                    *seqStore,
                            bool                        isForCorrection);

  void           exportData(writeBuffer                *exportDataFile,
                            sqCache                    *seqCache,
                            bool                        isForCorrection);

  bool           importData(readBuffer                 *importDataFile,
                            std::map<uint32, sqRead *> &reads,
                            FILE                       *layoutOutput,
                            FILE                       *sequenceOutput);


  void           reverseComplement(void);  //  Does NOT update childDeltas

  void           dumpFASTA(FILE *F);
  void           dumpFASTQ(FILE *F);
  void           dumpBAM(char const *prefix, sqStore *seqStore, u32toRead &seqReads);

  //  There are two multiAlign displays; this one, and one in abMultiAlign.
  void           display(FILE     *F,
                         sqStore  *seq,
                         uint32    displayWidth    = 100,    //  Width of display
                         uint32    displaySpacing  = 3,      //  Space between reads on the same line
                         bool      withDots        = true);

  //  Copy our important information to/from a clean struct we can then save to disk.
  void           saveToRecord(tgTigRecord &tr);
  void           restoreFromRecord(tgTigRecord &tr);

private:
public:
  uint32       _tigID = uint32max;             //  ID in the store, or UINT32_MAX if not set

  //  This is overloaded.
  //
  //  It's set in falconsense to track where we extracted the corrected
  //  read from.  But that information is never otherwise used, it was
  //  probably added for doing consensus with raw reads.
  //
  //  It's set during utgcns to indiacte what bases should be reported for a
  //  contig when it is trimmed back to circularize, or maybe later,
  //  bluntified.

  uint32       _trimBgn = 0;
  uint32       _trimEnd = 0;

  //  A variety of flags to suggest what type of unitig this is

  tgTig_class  _class           : 2  = tgTig_noclass;  //  Output classification: unassembled, bubble, contig
  uint32       _suggestRepeat   : 1  = false;          //  Bogart made this from detected repeat.
  uint32       _suggestCircular : 1  = false;          //  Bogart found overlaps making a circle.
  uint32       _suggestBubble   : 1  = false;          //  Bogart thinks this is a bubble.
  uint32       _circularLength  : 20 = 0;              //  Length of the end-to-end circularizing overlap.
  uint32       _suggestNoTrim   : 1 = false;           //  If set to true, ignore any trim value specified on the command line to consensus and do no trimming
  uint32       _spare           : 6  = 0;

  uint32       _layoutLen = 0;           //  The max coord in the layout.  Same as _basesLen if _bases exists.

  uint32       _basesLen = 0;            //  Doesn't include the NUL.
  uint32       _basesMax = 0;            //  Does    include the NUL.
  char        *_bases    = nullptr;      //  Linear consensus sequence.
  uint8       *_quals    = nullptr;      //

public:
  uint64       _childDeltaBitsLen = 0;         //  Created by utgcns or loaded from
  stuffedBits *_childDeltaBits    = nullptr;   //  tgStore output of consensus

  uint64       _childCIGARLen  = 0;
  char        *_childCIGARData = nullptr;
  char       **_childCIGAR     = nullptr;      //  Used by utgcns to output BAM format

private:
  tgPosition  *_children    = nullptr;   //  positions of objects that make up this tig
  uint32       _childrenLen = 0;
  uint32       _childrenMax = 0;

private:
  tgPosition  *_stashed    = nullptr;    //  children that we want to temporarily ignore
  uint32       _stashedLen = 0;
  uint32       _stashedMax = 0;

  uint64       _stashBack[2] = {0, 0};   //  Stash stats: read is backbone.
  uint64       _stashCont[2] = {0, 0};   //  Stash stats: read is contained and used for consensus.
  uint64       _stashStsh[2] = {0, 0};   //  Stash stats: read is stashed   NOT used for consensus.
  uint64       _stashIgnr[2] = {0, 0};   //  Stash stats: read is ignore    NOT used for consensus.

public:
  void                 clearStash(void) {
    delete [] _stashed;

    _stashed    = nullptr;
    _stashedLen = 0;
    _stashedMax = 0;

    _stashBack[0] = _stashBack[1] = 0;
    _stashCont[0] = _stashCont[1] = 0;
    _stashStsh[0] = _stashStsh[1] = 0;
    _stashIgnr[0] = _stashIgnr[1] = 0;
  }

  uint64        nStashBack(void) { return _stashBack[0]; }
  uint64        nStashCont(void) { return _stashCont[0]; }
  uint64        nStashStsh(void) { return _stashStsh[0]; }
  uint64        nStashIgnr(void) { return _stashIgnr[0]; }

  double        cStashBack(void) { return (double)_stashBack[1] / length(); }
  double        cStashCont(void) { return (double)_stashCont[1] / length(); }
  double        cStashStsh(void) { return (double)_stashStsh[1] / length(); }
  double        cStashIgnr(void) { return (double)_stashIgnr[1] / length(); }

  //  Flags for computing consensus/multialignments

  uint32       _utgcns_verboseLevel = 0;
};


#endif
