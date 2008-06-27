
#define TAGNAMEMAX  1024

class tapperComputation {
public:
  tapperComputation(seqInCore *a, seqInCore *b) {

    tag1f.clear();
    tag1r.clear();

    tag2f.clear();
    tag2r.clear();

    tag1size = a->sequenceLength() - 2;
    tag2size = b->sequenceLength() - 2;

    tag1f.setMerSize(tag1size); tag1f.setMerSpan(tag1size);
    tag1r.setMerSize(tag1size); tag1r.setMerSpan(tag1size);

    tag2f.setMerSize(tag1size); tag2f.setMerSpan(tag1size);
    tag2r.setMerSize(tag1size); tag2r.setMerSpan(tag2size);

    strncpy(tag1name, a->header() + 1, TAGNAMEMAX);
    tag1name[TAGNAMEMAX-1] = 0;

    strncpy(tag2name, b->header() + 1, TAGNAMEMAX);
    tag2name[TAGNAMEMAX-1] = 0;

    u32bit seqlen = a->sequenceLength();

    for (u32bit i=0; i<seqlen; i++) {
      tag1fseq[i] = a->sequence()[i];
      tag2fseq[i] = b->sequence()[i];
      tag1rseq[i] = tag1fseq[i];
      tag2rseq[i] = tag2fseq[i];
    }

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
    tag2rseq[0] = complementSymbol[tag2rseq[0]];

    reverseString(tag1rseq, seqlen);
    reverseString(tag2rseq, seqlen);

    //  Make sure things are terminated.

    tag1fseq[seqlen] = 0;
    tag1rseq[seqlen] = 0;
    tag2fseq[seqlen] = 0;
    tag2rseq[seqlen] = 0;

    //  Now just push the seqs into mers.

    //  Yes, i=2.  The first letter in the tag is the last in the
    //  adapter, and it's not in the tag.  We need to skip it.  The
    //  second letter (the first color) is biased by the adapter, and
    //  it will be an error 75% of the time.  Skip it too.
    //
    for (u32bit i=2; i<seqlen; i++) {
      tag1f += compressSymbol[tag1fseq[i]];
      tag2f += compressSymbol[tag2fseq[i]];
    }
    for (u32bit i=0; i<seqlen-2; i++) {
      tag1r += compressSymbol[tag1rseq[i]];
      tag2r += compressSymbol[tag2rseq[i]];
    }

    //fprintf(stderr, "T1: %s/%s T2:%s/%s\n", tag1fseq, tag1rseq, tag2fseq, tag2rseq);

    hitsLen = 0;
    hitsMax = 1024;
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

  void writeHits(FILE *OF) {
    for (u32bit i=0; i<hitsLen; i++)
      hits[i].writeHit(OF, tag1name, tag1f, tag1r);
  };

public:
  kMer        tag1f, tag1r;
  kMer        tag2f, tag2r;

  u32bit      tag1size;
  u32bit      tag2size;

  char        tag1name[1024];
  char        tag2name[1024];

  char        tag1fseq[64], tag1rseq[64];
  char        tag2fseq[64], tag2rseq[64];

  char        refseq[64];

  u32bit      hitsLen;
  u32bit      hitsMax;
  tapperHit  *hits;
};

