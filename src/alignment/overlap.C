


#include "AS_global.H"
#include "aligners.H"

//#include "overlapInCore.H"  //  For Extend_Alignment
#include "prefixEditDistance.H"

#include "AS_UTL_reverseComplement.H"

#include "kMer.H"
#include "merStream.H"

#include <map>
#include <vector>
#include <algorithm>

using namespace std;


//  Lots of code duplication with overlapPair.C

class exactMatch {
public:
  exactMatch(int32 a, int32 b, int32 l) {
    aBgn = a;
    bBgn = b;
    tLen = l;
  };

  int32  aBgn;  //  Signed to allow for easy compute of diagonal
  int32  bBgn;
  int32  tLen;

  bool operator<(exactMatch const that) const {
    if (tLen == that.tLen)
      return(aBgn < that.aBgn);

    return(tLen > that.tLen);
  };
};





ALNoverlap *
Overlap_forCNS(char  *aStr, uint32  aLen,
               char  *bStr, uint32  bLen,
               int32  minDiag,
               int32  maxDiag,
               double erate) {

  static ALNoverlap O;


  //  Find seeds - hash the kmer and position from the first read, then lookup
  //  each kmer in the second read.  For unique hits, save the diagonal.  Then what?
  //  If the diagonal is too far from the expected diagonal (based on the overlap),
  //  ignore the seed.

  map<uint64,int32>   aMap;  //  Signed, to allow for easy compute of diagonal
  map<uint64,int32>   bMap;

  uint32 merSize   = 18;
  bool   dupIgnore = false;  //  Ignore dups?  Otherwise, use the first occurrence

  //  Find mers in A
 findMersAgain:

  //fprintf(stderr, "Finding mers of size %u; dupIgnore=%s\n", merSize, dupIgnore ? "true" : "false");

  {
    kMerBuilder *kb = new kMerBuilder(merSize);
    seqStream   *ss = new seqStream(aStr, aLen);
    merStream   *ms = new merStream(kb, ss, true, true);

    while (ms->nextMer() == true) {
      uint64  kmer = ms->theFMer();

      if (aMap.find(kmer) != aMap.end()) {
        if (dupIgnore == true)
          aMap[kmer] = INT32_MAX;  //  Duplicate mer, now ignored!
      } else {
        aMap[kmer] = (int32)ms->thePositionInSequence();
      }
    }

    delete ms;
  }

  if (aMap.size() == 0) {
    aMap.clear();
    bMap.clear();

    merSize--;

    if ((merSize < 8) && (dupIgnore == true)) {
      merSize   = 20;
      dupIgnore =  false;
    }

    if (merSize >= 8)
      goto findMersAgain;
  }

  //  Find mers in B

  {
    kMerBuilder *kb = new kMerBuilder(merSize);
    seqStream   *ss = new seqStream(bStr, bLen);
    merStream   *ms = new merStream(kb, ss, true, true);

    while (ms->nextMer() == true) {
      uint64  kmer = ms->theFMer();

      if (aMap.find(kmer) == aMap.end())
        //  Doesn't exist in aSeq, don't care about it.
        continue;

      int32  apos = aMap[kmer];
      int32  bpos = (int32)ms->thePositionInSequence();

      if (apos == INT32_MAX)
        //  Exists too many times in aSeq, don't care about it.
        continue;

      if ((apos - bpos < minDiag) ||
          (apos - bpos > maxDiag))
        //  Too different.
        continue;

      if (bMap.find(kmer) != bMap.end()) {
        if (dupIgnore == true)
          bMap[kmer] = INT32_MAX;  //  Duplicate mer, now ignored!
      } else {
        bMap[kmer] = bpos;
      }
    }

    delete ms;
  }

  if (bMap.size() == 0) {
    aMap.clear();
    bMap.clear();

    merSize--;

    if ((merSize < 8) && (dupIgnore == true)) {
      merSize   = 20;
      dupIgnore =  false;
    }

    if (merSize >= 8)
      goto findMersAgain;
  }

  //  Still zero?  Didn't find any unique seeds anywhere.

  if (bMap.size() == 0)
    return(NULL);

  //  For unique mers in B, if the mer is also unique in A, add a hit.

  vector<exactMatch>   rawhits;
  vector<exactMatch>   hits;

  {
    for (map<uint64,int32>::iterator bit=bMap.begin(); bit != bMap.end(); bit++) {
      uint64  kmer = bit->first;
      int32   bpos = bit->second;

      if (bpos == INT32_MAX)
        //  Exists too many times in bSeq, don't care about it.
        continue;

      int32   apos = aMap[kmer];

      assert(apos != INT32_MAX);       //  Should never get a bMap for these

      assert(apos - bpos >= minDiag);  //  ...these too.
      assert(apos - bpos <= maxDiag);

      rawhits.push_back(exactMatch(apos, bpos, merSize));
    }
  }

  //  Sort by aPos (actually by length, then by aPos, but length is constant here).

  sort(rawhits.begin(), rawhits.end());

#if 0
  for (uint32 rr=0; rr<rawhits.size(); rr++)
    fprintf(stderr, "HIT: %d - %d diag %d\n", rawhits[rr].aBgn, rawhits[rr].bBgn, rawhits[rr].aBgn - rawhits[rr].bBgn);
#endif

  //  Chain the hits.

  hits.push_back(rawhits[0]);

  for (uint32 rr=1; rr<rawhits.size(); rr++) {
    uint32  hh = hits.size() - 1;

    int32   da = rawhits[rr].aBgn - hits[hh].aBgn;
    int32   db = rawhits[rr].bBgn - hits[hh].bBgn;

    if ((da > 0) && (da < 2 * merSize) && (da == db))
      hits[hh].tLen += da;
    else
      hits.push_back(rawhits[rr]);
  }

  //  Sort by longest

  sort(hits.begin(), hits.end());

#if 0
  for (uint32 hh=0; hh<hits.size(); hh++) {
    fprintf(stderr, "hit %02u %5d-%5d diag %d len %3u\n",
            hh,
            hits[hh].aBgn, hits[hh].bBgn,
            hits[hh].aBgn - hits[hh].bBgn,
            hits[hh].tLen);
  }
#endif

  //  Recompute.

  for (uint32 hh=0; hh<hits.size(); hh++) {
    Match_Node_t  match;

    match.Start  = hits[0].aBgn;   //  Begin position in a
    match.Offset = hits[0].bBgn;   //  Begin position in b
    match.Len    = merSize;        //  tLen can include mismatches
    match.Next   = 0;              //  Not used here


    int32 aLo=0, aHi=0;
    int32 bLo=0, bHi=0;

    prefixEditDistance   ed(false, 0.40);  // This will bomb.

    int32      errors  = 0;
    Overlap_t  ovltype = ed.Extend_Alignment(&match,       //  Initial exact match, relative to start of string
                                             aStr, aLen,
                                             bStr, bLen,
                                             aLo,  aHi,    //  Output: Regions which the match extends
                                             bLo,  bHi,
                                             errors,
                                             false);

    if (ovltype == NONE)
      //  Not a good overlap, keep searching for one.
      continue;

    int32  olapLen = 1 + min(aHi - aLo, bHi - bLo);
    double quality = (double)errors / olapLen;

    //  struct ALNoverlap {
    //  int begpos;
    //  int endpos;
    //  int length;
    //  int diffs;
    //  int comp;
    //  int *trace;
    //};  //  Former 'Overlap'

    //     Logic is same as in Optimal_Overlap_AS_forCNS

    O.begpos = (aLo        > 0) ? (aLo)        : -(bLo);
    O.endpos = (bLen - bHi > 0) ? (bLen - bHi) : -(aLen - aHi);
    O.length = olapLen;  //  NOT TRUE!  Want's the gapped length.
    O.diffs  = errors;
    O.comp   = 0;  //  0 if forward, 1 if complement
#warning NEED TO COPY DELTAS TO OUTPUT
    O.trace  = ed.Left_Delta;

    if (quality < erate)
      return(&O);
  }

  return(NULL);
}
