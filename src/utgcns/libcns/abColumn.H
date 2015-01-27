
class abBaseCount {
public:
  abBaseCount() {
    clear();
  };
  ~abBaseCount() {
  };

  void clear(void) {
    _depth    = 0;

    for (uint32 ii=0; ii<CNS_NUM_SYMBOLS; ii++)
      _count[ii] = 0;
  };

  uint32 IncBaseCount(char c) {
    uint32 r = baseToIndex[c];

    if ((r < 0) || (r > 5))
      fprintf(stderr, "IncBaseCount c=%d '%c' r=%u out of range.\n", c, c, r);
    assert(0 <= r);
    assert(r <= 5);

    _depth++;

    return(_count[r]++);
  };

  uint32 DecBaseCount(char c) {
    uint32 r = baseToIndex[c];

    if ((r < 0) || (r > 5))
      fprintf(stderr, "DecBaseCount c=%d '%c' r=%u out of range.\n", c, c, r);
    assert(0 <= r);
    assert(r <= 5);

    _depth--;

    return(_count[r]--);
  };

  uint32 GetBaseCount(char c) {
    uint32 r = baseToIndex[c];

    if ((r < 0) || (r > 5))
      fprintf(stderr, "GetBaseCount c=%d '%c' r=%u out of range.\n", c, c, r);
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
  uint32
  GetMaxBaseCountBase(bool skipGap) {
    uint32 iMax = (skipGap == false) ? 0 : 1;
    uint32 tied = 0;

    for (int32 i=iMax+1; i<CNS_NUM_SYMBOLS-1; i++) {
      if (_count[i] > _count[iMax]) {
        iMax = i;
        tied = 0;
      } else if (_count[i] == _count[iMax]) {
        tied++;
      }
    }

    if (tied > 0) {
      uint32   iTie = 0;
    
      for (uint32 i=1; i<CNS_NUM_SYMBOLS-1; i++) {
        if (_count[i] == _count[iMax]) {
          uint32 tieBreaker = random();

          if (tieBreaker > iTie) {
            iTie = tieBreaker;
            iMax = i;
          }
        }
      }
    }

    return(indexToBase[iMax]);
  };

private:
  uint32   _depth;
  uint32   _count[CNS_NUM_SYMBOLS];
};




//  A column of beads in the multialign
//
class abColumn {
public:
  abColumn() {
    ma_position = 0;
  };
  ~abColumn() {
  };

  uint32     GetDepth(void)                 { return(base_count.GetDepth());      };
  uint32     GetColumnBaseCount(char c)     { return(base_count.GetBaseCount(c)); };

  void       CheckBaseCounts(abAbacus *ab);

  //  Max over what?  Also needs abacus to get next beads!  Gah!  Two with the same name!
  //#warning GetMaxBaseCount confusion.
  //char       GetMaxBaseCount(abAbacus *ab, int32 start_index);
  char       GetMaxBaseCountBase(bool skipGap) {
    return(base_count.GetMaxBaseCountBase(skipGap));
  };

  abColID    const  ident(void)                    { return(lid);  };

  abBeadID   const  callID(void)                   { return(call); };

  abColID          &prevID(void)                   { return(prev); };
  abColID          &nextID(void)                   { return(next); };

  int32             mismatch(void) {
    return(base_count.GetDepth() - base_count.GetBaseCount( base_count.GetMaxBaseCountBase(false) ));
  };

  int32            &position(void)                 { return(ma_position); };

private:
  abColID          lid;           //  index in columnStore
  abBeadID         call;          //  global offset in beadStore;
  abColID          prev;          //  navigation in columnStore;
  abColID          next;          //
  abMultiAlignID   ma_id;         //  abMultiAlign this column belongs to
  int32            ma_position;   //  was ma_index; column ordinal, needs to be signed for compat with tgTig bgn/end
  abBaseCount      base_count;

  friend class abAbacus;
};




