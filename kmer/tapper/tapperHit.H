



class tapperHit {
public:
  //tapperHit() {};
  //~tapperHit();

  u32bit    numberOfColorErrors(void)  { return(_colorErrors); };
  u32bit    numberOfBaseErrors(void)   { return(_basesErrors); };

  void writeHit(FILE *OF,
                char *tagname, kMer& tagf, kMer& tagr) {

    fprintf(OF, "%s\t"u32bitFMT"\t"u32bitFMT"\t"u32bitFMT"/"u32bitFMT"\t%s\t%s/%s\t%s/%s\n",
            tagname,
            _seqIdx, _seqPos,
            _colorErrors,
            _basesErrors,
            _rev ? "<-" : "->",
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

  u32bit  _rev:1;
  u32bit  _len:6;
  u32bit  _colorErrors:6;
  u32bit  _basesErrors:6;

  char    _tagCOLOR[64];
  char    _tagCOREC[64];  //  corrected
  char    _refCOLOR[64];

  char    _tagACGT[64];
  char    _refACGT[64];
};
