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

    strncpy(tag1name, a->header() + 1, 64);
    tag1name[63] = 0;

    strncpy(tag2name, b->header() + 1, 64);
    tag2name[63] = 0;

    u32bit seqlen = a->sequenceLength();

    //  Yes, i=2.  The first letter in the tag is the last in the
    //  adapter, and it's not in the tag.  We need to skip it.  The
    //  second letter (the first color) is biased by the adapter, and
    //  it will be an error 75% of the time.  Skip it too.
    //
    for (u32bit i=0; i<seqlen; i++) {
      tag1fseq[i] = a->sequence()[i];
      tag2fseq[i] = b->sequence()[i];
    }

    //  Rebuild the tag, for reverse-complement.  Add a new base to
    //  the end of the forward tag - this will become the new
    //  reference letter.  We'll add a 'T'.
    //
    //  Reverse the color encoding.  Drop the last color (the first
    //  color in forward).  Complement the T, append that to the
    //  start.
    //
    //  The forward tag:
    //    T 3 1 0 2 1 1
    //     A C C T G T
    //
    //  Add the T to the end
    //    T 3 1 0 2 1 1 0
    //     A C C T G T T
    //
    //  Reverse the colors, drop the last, add the complement of the
    //  new reference.
    //    A 0 1 1 2 0 1
    //     A C A G G T
    //

    char  l1 = tag1fseq[0];
    char  l2 = tag2fseq[0];

    tag1rseq[0] = tag1fseq[0];
    tag2rseq[0] = tag2fseq[0];

    for (u32bit i=1; i<seqlen; i++) {
      tag1rseq[i] = tag1fseq[i];
      tag2rseq[i] = tag2fseq[i];
      l1 = baseToColor[l1][tag1rseq[i]];
      l2 = baseToColor[l2][tag2rseq[i]];
    }

    tag1rseq[seqlen] = baseToColor[l1]['A'];
    tag2rseq[seqlen] = baseToColor[l2]['A'];

    tag1rseq[seqlen+1] = 'A';
    tag2rseq[seqlen+1] = 'A';

    tag1rseq[seqlen+2] = 0;
    tag2rseq[seqlen+2] = 0;

    reverseString(tag1rseq, seqlen+2);
    reverseString(tag2rseq, seqlen+2);

    //  Make sure things are terminated.

    tag1fseq[seqlen] = 0;
    tag1rseq[seqlen] = 0;
    tag2fseq[seqlen] = 0;
    tag2rseq[seqlen] = 0;

    //  Now just push the seqs into mers.

    for (u32bit i=2; i<seqlen; i++) {
      tag1f += compressSymbol[tag1fseq[i]];
      tag1r += compressSymbol[tag1rseq[i]];
      tag2f += compressSymbol[tag2fseq[i]];
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
      hits[i].writeHit(OF, tag1name, tag1f, tag1r, tag2name, tag2f, tag2r);
  };

public:
  kMer        tag1f, tag1r;
  kMer        tag2f, tag2r;

  u32bit      tag1size;
  u32bit      tag2size;

  char        tag1name[64];
  char        tag2name[64];

  char        tag1fseq[64], tag1rseq[64];
  char        tag2fseq[64], tag2rseq[64];

  char        refseq[64];

  u32bit      hitsLen;
  u32bit      hitsMax;
  tapperHit  *hits;
};

