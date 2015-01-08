
class abBaseCount {
public:
  abBaseCount() {
    ResetBaseCount();
  };
  ~abBaseCount() {
  };

  int32 IncBaseCount(char c) {
    if ((RINDEX[c] < 0) || (RINDEX[c] > 5))
      fprintf(stderr, "IncBaseCount i out of range (possibly non ACGTN letter?)");
    assert(0 <= RINDEX[c]);
    assert(RINDEX[c] <= 5);

    _depth++;

    return(_count[RINDEX[c]]++);
  };

  int32 DecBaseCount(char c) {
    if ((RINDEX[c] < 0) || (RINDEX[c] > 5))
      fprintf(stderr, "IncBaseCount i out of range (possibly non ACGTN letter?)");
    assert(0 <= RINDEX[c]);
    assert(RINDEX[c] <= 5);

    _depth--;

    return(_count[RINDEX[c]]--);
  };

  int32 GetBaseCount(char c) {
    if ((RINDEX[c] < 0) || (RINDEX[c] > 5))
      fprintf(stderr, "IncBaseCount i out of range (possibly non ACGTN letter?)");
    assert(0 <= RINDEX[c]);
    assert(RINDEX[c] <= 5);

    return(_count[RINDEX[c]]);
  };

  int32 GetDepth(void) {
    return(_depth);
  };

  void ResetBaseCount(void) {
    _depth    = 0;

    for (uint32 ii=0; ii<CNS_NALPHABET; ii++)
      _count[ii] = 0;
  };

  //  Return the character with the highest base count, breaking ties randomly (excluding gaps).
  //  Start at 1 to disallow gaps.
  //
  //  BUG FIX:  a two-way tie was NOT broken, the lower letter always won
  //
  char  GetMaxBaseCount(int32 startCharacter) {
    int32 iMax = startCharacter;
    int32 tied = 0;

    for (int32 i=startCharacter+1; i<CNS_NALPHABET-1; i++) {
      if (count[i] > count[iMax]) {
        iMax = i;
        tied = 0;
      } else if (count[i] == count[iMax]) {
        tied++;
      }
    }

    if (tied > 0) {
      int32   iTie = 0;
    
      for (int32 i=1; i<CNS_NALPHABET-1; i++) {
        if (count[i] == count[iMax]) {
          int32 tieBreaker = random();

          if (tie_breaker > iTie) {
            iTie = tieBreaker;
            iMax = i;
          }
        }
      }
    }

    return toupper(ALPHABET[iMax]);
  };


  void ShowBaseCount(void) {
    fprintf(stderr,"%d total\n",b->depth);
    for (uint32 i=0;i<CNS_NALPHABET;i++)
      fprintf(stderr,"%c\t",ALPHABET[i]);

    fprintf(stderr,"\n");
    for (uint32 i=0;i<CNS_NALPHABET;i++)
      fprintf(stderr,"%d\t",b->count[i]);
    fprintf(stderr,"\n");
  };

  void ShowBaseCountPlain(FILE *out) {
    fprintf(out,"%d\t",b->depth);
    for (uint32 i=0;i<CNS_NALPHABET;i++) {
      fprintf(out,"%d\t",b->count[i]);
    }
  };

private:
  int32   _depth;
  int32   _count[CNS_NALPHABET];
};


