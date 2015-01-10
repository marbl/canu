
class abBaseCount {
public:
  abBaseCount() {
    clear();
  };
  ~abBaseCount() {
  };

  void clear(void) {
    _depth    = 0;

    for (uint32 ii=0; ii<CNS_NALPHABET; ii++)
      _count[ii] = 0;
  };

  uint32 IncBaseCount(char c) {
    uint32 r = RINDEX[c];

    if ((r < 0) || (r > 5))
      fprintf(stderr, "IncBaseCount c=%c r=%u out of range.\n", c, r);
    assert(0 <= r);
    assert(r <= 5);

    _depth++;

    return(_count[r]++);
  };

  uint32 DecBaseCount(char c) {
    uint32 r = RINDEX[c];

    if ((r < 0) || (r > 5))
      fprintf(stderr, "DecBaseCount c=%c r=%u out of range.\n", c, r);
    assert(0 <= r);
    assert(r <= 5);

    _depth--;

    return(_count[r]--);
  };

  uint32 GetBaseCount(char c) {
    uint32 r = RINDEX[c];

    if ((r < 0) || (r > 5))
      fprintf(stderr, "GetBaseCount c=%c r=%u out of range.\n", c, r);
    assert(0 <= r);
    assert(r <= 5);

    return(_count[r]);
  };



  uint32 GetDepth(void) {
    return(_depth);
  };


  //  Return the character with the highest base count, breaking ties randomly (excluding gaps).
  //  Start at 1 to disallow gaps.
  //
  //  BUG FIX:  a two-way tie was NOT broken, the lower letter always won
  //
  //  BUG:  This wants to return a letter, but doesn't know the mapping
  //
  uint32  GetMaxBaseCount(bool skipGap) {
    uint32 iMax = (skipGap == false) ? 0 : 1;
    uint32 tied = 0;

    for (int32 i=iMax+1; i<CNS_NALPHABET-1; i++) {
      if (_count[i] > _count[iMax]) {
        iMax = i;
        tied = 0;
      } else if (_count[i] == _count[iMax]) {
        tied++;
      }
    }

    if (tied > 0) {
      uint32   iTie = 0;
    
      for (uint32 i=1; i<CNS_NALPHABET-1; i++) {
        if (_count[i] == _count[iMax]) {
          uint32 tieBreaker = random();

          if (tieBreaker > iTie) {
            iTie = tieBreaker;
            iMax = i;
          }
        }
      }
    }

    return(iMax);
  };


#if 0
  void ShowBaseCount(void) {
    fprintf(stderr,"%d total\n", _depth);
    for (uint32 i=0;i<CNS_NALPHABET;i++)
      fprintf(stderr,"%c\t",ALPHABET[i]);

    fprintf(stderr,"\n");
    for (uint32 i=0;i<CNS_NALPHABET;i++)
      fprintf(stderr,"%d\t", _count[i]);
    fprintf(stderr,"\n");
  };

  void ShowBaseCountPlain(FILE *out) {
    fprintf(out,"%d\t",b->depth);
    for (uint32 i=0;i<CNS_NALPHABET;i++) {
      fprintf(out,"%d\t", _count[i]);
    }
  };
#endif

private:
  uint32   _depth;
  uint32   _count[CNS_NALPHABET];
};


//  A column of beads in the multialign
//
class abColumn {
public:
  abColumn();
  ~abColumn();

  uint32     GetDepth(void)                 { return(base_count.GetDepth());      };
  uint32     GetColumnBaseCount(char c)     { return(base_count.GetBaseCount(c)); };

  void       CheckBaseCounts(abAbacus *ab);

  //  Max over what?  Also needs abacus to get next beads!  Gah!  Two with the same name!
  //#warning GetMaxBaseCount confusion.
  //char       GetMaxBaseCount(abAbacus *ab, int32 start_index);
  char       GetMaxBaseCount(bool skipGap)  {  return(base_count.GetMaxBaseCount(skipGap)); };

  abColID    ident(void)    { return(lid); };

  abBeadID   callID(void)   { return(call);    };

  abColID    prevID(void)   { return(prev);    };
  abColID    nextID(void)   { return(next);    };

  int32      mismatch(void) {
    return(base_count.GetDepth() - base_count.GetBaseCount( base_count.GetMaxBaseCount(false) ));
  };

  uint32     position(void) { return(ma_position); };

private:
  abColID          lid;       //  index in columnStore
  abBeadID         call;      //  global offset in beadStore;
  abColID          prev;      //  navigation in columnStore;
  abColID          next;
  abMultiAlignID   ma_id;     //  external tigStore id ??
  uint32           ma_position;  //  was ma_index; index in MANode; // refreshed only periodically // seems to also be gapped position in the align
  abBaseCount      base_count;

  friend class abAbacus;
};




