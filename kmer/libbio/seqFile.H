
class seqFile {
protected:
  seqFile(char *filename) {};

public:
  seqFile() {};
  virtual ~seqFile() {};

  virtual const char  *fileTypeName(void) {
    return("undefined sequence file type");
  };

public:
  virtual seqFile     *openFile(const char *name) = 0;

  virtual void         openIndex(u32bit indextypetoload = 0) = 0;

  virtual bool         rewind(void) = 0;
  virtual bool         eof(void) = 0;

  virtual bool         getSequence(u32bit &hLen, char *&h,
                                   u32bit &sLen, char *&s) {
    fprintf(stderr, "seqFile::getSequence()-- ERROR, default function not re-defined.\n");
    return(false);
  };

  virtual seqInCore   *getSequenceInCore(void) {
    char   *h = NULL, *s = NULL;
    u32bit  hLen=0,    sLen=0;

    if (getSequence(hLen, h, sLen, s) == false) {
      fprintf(stderr, "getSequenceInCore()-- Failed to getSequence()\n");
      return(0L);
    }

    return(new seqInCore(_curIID++, h, hLen, s, sLen));
  };

  virtual seqOnDisk   *getSequenceOnDisk(void) {
    char   *h = NULL, *s = NULL;
    u32bit  hLen=0,    sLen=0;

    if (getSequence(hLen, h, sLen, s) == false) {
      fprintf(stderr, "getSequenceOnDisk()-- Failed to getSequence()\n");
      return(0L);
    }

    return(new seqOnDisk(_curIID++, h, hLen, s, sLen));
  };

  virtual u64bit       timeStamp(void) = 0;

  virtual u32bit       currentIID(void) {
    return(_curIID);
  };

  virtual bool         find(seqIID  iid) = 0;
  virtual bool         find(char   *uid) = 0;

  virtual u32bit       sequenceLength(seqIID iid) = 0;
  virtual u32bit       nameLength(seqIID iid) = 0;

  virtual u32bit       getNumberOfSequences(void) = 0;
  virtual char const  *getSourceName(void) = 0;

  virtual void         printDescription(FILE *out, char *name) = 0;

  virtual bool         isIndexed(void) = 0;
  virtual bool         isSqueezed(void) = 0;

protected:
  seqIID   _curIID;
};
