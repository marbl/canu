#include "alphabet.h"


class tapperHit {
public:
  //tapperHit() {};
  //~tapperHit();

  u32bit    numberOfBaseMismatches(void)        { return(_basesMismatch); };
  u32bit    numberOfColorMismatches(void)       { return(_colorMismatch); };
  u32bit    numberOfColorInconsistencies(void)  { return(_colorMismatch); };

  void writeHit(FILE *OF, u64bit  tagid) {

    fprintf(OF, "ID:"u64bitFMT"\t"u32bitFMT"\t"u32bitFMT"\t"u64bitFMT"/"u64bitFMT"/"u64bitFMT"\t%c\t%s/%s\t%s/%s\n",
            tagid,
            _seqIdx, _seqPos,
            _basesMismatch,
            _colorMismatch,
            _colorInconsistent,
            _rev ? '-' : '+',
            _tagCOLOR, _refCOLOR,
            _tagACGT, _refACGT);
  };

  void   alignToReference(tapperGlobalData *g,
                          u32bit so,
                          u32bit po,
                          char  *tag, u32bit len);

private:
  u32bit  _seqIdx;
  u32bit  _seqPos;

  u64bit  _tagIdx;               //  4e9 tags is only 34x of human

  u64bit  _pad:17;               //  Nothing
  u64bit  _len:6;                //  Length of tag
  u64bit  _rev:1;                //  Match is reversecomplement
  u64bit  _rank:16;              //  Rank in the sorted list of hits

  u64bit  _basesMismatch:6;      //  Number of mismatches in ACGT alignment

  u64bit  _colorMismatch:6;      //  Number of consistent color mismatches
  u64bit  _colorInconsistent:6;  //  Number of inconsistent color mismatches

  char    _tagCOLOR[TAG_LEN_MAX];
  char    _refCOLOR[TAG_LEN_MAX];

  char    _tagACGT[TAG_LEN_MAX];
  char    _refACGT[TAG_LEN_MAX];
};
