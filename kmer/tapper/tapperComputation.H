

class tapperComputation {
public:
  tapperComputation(tapperTag *a, tapperTag *b) {

    tag1f.clear();
    tag1r.clear();

    tag2f.clear();
    tag2r.clear();

    //  Process the first tag.  It's the same for the second tag.
    //  When this gets stable, replicate.

    tag1id   = a->decode(tag1fseq, tag1fqlt);
    tag1size = strlen(tag1fseq);

    tag1f.setMerSize(tag1size-1); tag1f.setMerSpan(tag1size-1);
    tag1r.setMerSize(tag1size-1); tag1r.setMerSpan(tag1size-1);

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

    tag1rseq[0] = complementSymbol[tag1rseq[0]];

    reverseString(tag1rseq, tag1size + 1);

    //  Make sure things are terminated.

    tag1fseq[tag1size] = 0;
    tag1rseq[tag1size] = 0;

    //  Now just push the seqs into mers.

    //  Yes, i=2.  The first letter in the tag is the last in the
    //  adapter, and it's not in the tag.  We need to skip it.  The
    //  second letter (the first color) is biased by the adapter, and
    //  it will be an error 75% of the time.  Skip it too.
    //
    for (u32bit i=2; i<tag1size; i++) {
      tag1f += compressSymbol[tag1fseq[i]];
    }
    for (u32bit i=0; i<tag1size-2; i++) {
      tag1r += compressSymbol[tag1rseq[i]];
    }

    //  Now, strip out the first color.
    tag1size--;

    //fprintf(stderr, "T1: %s/%s T2:%s/%s\n", tag1fseq, tag1rseq, tag2fseq, tag2rseq);

    hitsLen = 0;
    hitsMax = 16;
    hits    = new tapperHit [hitsMax];
  };
  ~tapperComputation() {
    delete [] hits;
  };

  void addHit(tapperHit& h) {
    if (hitsLen >= hitsMax) {
      hitsMax *= 2;
      tapperHit *nits = new tapperHit [hitsMax];
      memcpy(nits, hits, sizeof(tapperHit) * hitsLen);
      delete [] hits;
      hits = nits;
    }
    hits[hitsLen++] = h;
  };

  void sortHits(void) {
    
  };

  void writeHits(FILE *OF) {
    for (u32bit i=0; i<hitsLen; i++)
      hits[i].writeHit(OF, tag1id);
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

  u32bit      hitsLen;
  u32bit      hitsMax;
  tapperHit  *hits;
};

