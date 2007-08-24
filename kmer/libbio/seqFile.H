
class seqFile {
protected:
  seqFile(char *filename) {};

public:
  seqFile() {};
  virtual ~seqFile() {};

public:
  virtual seqFile     *openFile(const char *name) = 0;

  virtual void         openIndex(u32bit indextypetoload = 0) = 0;

  virtual bool         rewind(void) = 0;
  virtual bool         eof(void) = 0;

  virtual seqInCore   *getSequenceInCore(void) = 0;
  virtual seqOnDisk   *getSequenceOnDisk(void) = 0;

  virtual u64bit       timeStamp(void) = 0;

  virtual u32bit       currentIID(void) = 0;

  virtual bool         find(seqIID  iid) = 0;
  virtual bool         find(char   *uid) = 0;

  virtual u32bit       sequenceLength(seqIID iid) = 0;
  virtual u32bit       nameLength(seqIID iid) = 0;

  virtual u32bit       getNumberOfSequences(void) = 0;
  virtual char const  *getSourceName(void) = 0;

  virtual void         printDescription(FILE *out, char *name) = 0;

  virtual bool         isIndexed(void) = 0;
  virtual bool         isSqueezed(void) = 0;
};


