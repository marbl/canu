
#include <algorithm>

class tapperComputation {
public:
  tapperComputation(tapperTag *a, tapperTag *b) {
    u16bit  id[4];

    tag1f.clear();
    tag1r.clear();

    tag2f.clear();
    tag2r.clear();

    tag1size = 0;
    tag2size = 0;

    //  Process the tags.
    //
    //  It's not a trivial operation (probably not even possible) to
    //  reverse-complement a SOLiD read.  To reverse complement a
    //  read, we would need to construct a new reference base, but to
    //  construct that base, we need to decode the read from color to
    //  acgt.  Any errors in the read prevent precise decoding, and we
    //  end up building the new reverse-complemented read with all the
    //  errors at the start.  By adding the anchor base to the end,
    //  we're fixing all the bases in error.
    //  
    //  So, we need to handle reverse reads specially.  Reverse the
    //  read (color-space is self-complementing), and RC the anchor
    //  base.  Any downstream processes need to know that a read has
    //  the anchor at the start OR at the end.
    //
    //  When building the mers (the for loops), yes, i=2.  The first
    //  letter in the tag is the last in the adapter, and it's not in
    //  the tag.  We need to skip it.  The second letter (the first
    //  color) is biased by the adapter, and it will be an error 75%
    //  of the time.  Skip it too.
    
    if (a) {
      tag1id   = a->decode(id, tag1fseq, tag1fqlt);
      tag1size = strlen(tag1fseq);

      tag1f.setMerSize(tag1size-1); tag1f.setMerSpan(tag1size-1);
      tag1r.setMerSize(tag1size-1); tag1r.setMerSpan(tag1size-1);

      for (u32bit i=0, j=tag1size-1; i<tag1size; i++, j--) {
        tag1rseq[i] = tag1fseq[j];
        tag1rqlt[i] = tag1fqlt[j];
      }
      tag1rseq[tag1size-1] = complementSymbol[tag1rseq[tag1size-1]];
      tag1rseq[tag1size] = 0;
      tag1rqlt[tag1size] = 0;

      for (u32bit i=2; i<tag1size; i++)
        tag1f += letterToBits[tag1fseq[i]];

      for (u32bit i=0; i<tag1size-2; i++)
        tag1r += letterToBits[tag1rseq[i]];

      tag1size--;
    }

    if (b) {
      tag2id   = b->decode(id, tag2fseq, tag2fqlt);
      tag2size = strlen(tag2fseq);

      tag2f.setMerSize(tag2size-1); tag2f.setMerSpan(tag2size-1);
      tag2r.setMerSize(tag2size-1); tag2r.setMerSpan(tag2size-1);

      for (u32bit i=0, j=tag2size-1; i<tag2size; i++, j--) {
        tag2rseq[i] = tag2fseq[j];
        tag2rqlt[i] = tag2fqlt[j];
      }
      tag2rseq[tag2size-1] = complementSymbol[tag2rseq[tag2size-1]];
      tag2rseq[tag2size] = 0;
      tag2rqlt[tag2size] = 0;

      for (u32bit i=2; i<tag2size; i++)
        tag2f += letterToBits[tag2fseq[i]];

      for (u32bit i=0; i<tag2size-2; i++)
        tag2r += letterToBits[tag2rseq[i]];

      tag2size--;
    }

    //fprintf(stderr, "T1: %s/%s T2:%s/%s\n", tag1fseq, tag1rseq, tag2fseq, tag2rseq);

    tag1hitsLen = 0;
    tag1hitsMax = 16;
    tag1hits    = new tapperHit [tag1hitsMax];

    tag2hitsLen = 0;
    tag2hitsMax = 16;
    tag2hits    = new tapperHit [tag2hitsMax];

    mean   = 0;
    stddev = 0;

    alignQualHistogram           = 0L;
    alignQualHistogramLen        = 0;

    resultFragment               = 0L;
    resultFragmentLen            = 0L;

    resultSingleton              = 0L;
    resultSingletonLen           = 0L;

    resultTangledAlignment       = 0L;
    resultTangledAlignmentLen    = 0L;

    resultMated                  = 0L;
    resultMatedLen               = 0L;

    resultTangled                = 0L;
    resultTangledLen             = 0L;
  };

  ~tapperComputation() {
    delete [] tag1hits;
    delete [] tag2hits;

    delete [] resultFragment;
    delete [] resultSingleton;
    delete [] resultTangledAlignment;
    delete [] resultMated;
    delete [] resultTangled;
  };

  void addHit(tapperGlobalData *g, tapperHit& h, bool tag1) {

    u32bit  ii = g->TA->AQIindex(g->maxBaseError, g->maxColorError,
                                 h.numberOfBaseMismatches(),
                                 h.numberOfColorMismatches(),
                                 h.numberOfColorInconsistencies());
 
    if (alignQualHistogram == 0L) {
      alignQualHistogramLen = g->TA->AQIlength(g->maxBaseError, g->maxColorError);
      alignQualHistogram    = new u32bit [alignQualHistogramLen];

      memset(alignQualHistogram, 0, sizeof(u32bit) * alignQualHistogramLen);
    }

    alignQualHistogram[ii]++;

    if (tag1) {
      if (tag1hitsLen >= tag1hitsMax) {
        tag1hitsMax *= 2;
        tapperHit *nits = new tapperHit [tag1hitsMax];
        memcpy(nits, tag1hits, sizeof(tapperHit) * tag1hitsLen);
        delete [] tag1hits;
        tag1hits = nits;
      }
      tag1hits[tag1hitsLen++] = h;
    } else {
      if (tag2hitsLen >= tag2hitsMax) {
        tag2hitsMax *= 2;
        tapperHit *nits = new tapperHit [tag2hitsMax];
        memcpy(nits, tag2hits, sizeof(tapperHit) * tag2hitsLen);
        delete [] tag2hits;
        tag2hits = nits;
      }
      tag2hits[tag2hitsLen++] = h;
    }
  };

  void sortHitsByPosition(void) {
    tapperHitPositionCompare pc;
    std::sort(tag1hits, tag1hits+tag1hitsLen, pc);
    std::sort(tag2hits, tag2hits+tag2hitsLen, pc);
  };
  
public:
  kMer        tag1f, tag1r;
  kMer        tag2f, tag2r;

  u32bit      tag1size;
  u32bit      tag2size;

  u64bit      tag1id;
  u64bit      tag2id;

  char        tag1fseq[TAG_LEN_MAX], tag1rseq[TAG_LEN_MAX];
  char        tag2fseq[TAG_LEN_MAX], tag2rseq[TAG_LEN_MAX];

  u64bit      tag1fqlt[TAG_LEN_MAX], tag1rqlt[TAG_LEN_MAX];
  u64bit      tag2fqlt[TAG_LEN_MAX], tag2rqlt[TAG_LEN_MAX];

  u32bit      tag1hitsLen;
  u32bit      tag1hitsMax;
  tapperHit  *tag1hits;

  u32bit      tag2hitsLen;
  u32bit      tag2hitsMax;
  tapperHit  *tag2hits;

  u32bit                  mean;
  u32bit                  stddev;

  u32bit                 *alignQualHistogram;
  u32bit                  alignQualHistogramLen;

  tapperResultFragment   *resultFragment;
  u32bit                  resultFragmentLen;

  tapperResultFragment   *resultSingleton;
  u32bit                  resultSingletonLen;

  tapperResultFragment   *resultTangledAlignment;
  u32bit                  resultTangledAlignmentLen;

  tapperResultMated      *resultMated;
  u32bit                  resultMatedLen;

  tapperResultTangled    *resultTangled;
  u32bit                  resultTangledLen;
};

