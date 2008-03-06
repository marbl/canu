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
    for (u32bit i=2; i<seqlen; i++) {
      tag1f += compressSymbol[a->sequence()[i]];
      tag1r += compressSymbol[a->sequence()[seqlen-i+1]];
    }

    for (u32bit i=2; i<seqlen; i++) {
      tag2f += compressSymbol[b->sequence()[i]];
      tag2r += compressSymbol[b->sequence()[seqlen-i+1]];
    }

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
      hits[i].writeHit(OF, tag1name, tag1f, tag2r, tag2name, tag2f, tag2r);
  };

public:
  kMer        tag1f, tag1r;
  kMer        tag2f, tag2r;

  u32bit      tag1size;
  u32bit      tag2size;

  char        tag1name[64];
  char        tag2name[64];

  u32bit      hitsLen;
  u32bit      hitsMax;
  tapperHit  *hits;
};

