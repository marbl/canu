 

class seqStream {
private:
  void          clearGuts(void);
public:
  seqStream();
  seqStream(const char *filename, bool finishMe);
  seqStream(seqFile *S, bool finishMe);
  ~seqStream();

  u64bit        timeStamp(void) { return(_fileTimeStamp); };

  //  Adds sequences to the useList.  For parse(), sequences are
  //  specified in a character string, e.g., "10,20,30-40,50".  add()
  //  pushes one sequence iid onto the list.  Multiple calls to
  //  parse() and add() are allowed.  There is (currently) no way to
  //  remove sequences from the useList.  Once you're all done adding
  //  sequences, call finish().
  //
  void          setFile(const char *filename);
  void          setFile(seqFile *S);
  void          parse(char *line);
  void          add(u32bit v);
  void          finish(void);

  void          setSeparator(char sep, u32bit len);
  char          getSeparator(void) { return(_separator); };

  //  get() returns one letter per input letter -- a gap of size n
  //  will return n gap symbols.  seqPos() is the position we are at
  //  in the current sequence; seqIID() is the iid of that sequence.
  //
  unsigned char get(void);
  bool          eof(void)        { return(_eof); };

  bool          rewind(void);

  u64bit        seqPos(void)     { return(_currentPos - 1); };
  u64bit        seqIID(void)     { return(_useList[_currentSeq].iid); };
  u64bit        strPos(void)     { return(_useList[_currentSeq].start + _currentPos - 1); };

  //  Return the length of, position of (in the chain) and IID of the
  //  (s)th sequence in the chain.
  //
  u32bit        numberOfSequences(void) { return(_useListLen); };
  u64bit        lengthOfSequences(void) { return(_lengthOfSequences); };
  u64bit        lengthOf(u32bit s) { return((s >= _useListLen) ? ~u64bitZERO : _useList[s].length); };
  u64bit        startOf(u32bit s)  { return((s >= _useListLen) ? ~u64bitZERO : _useList[s].start);  };
  u64bit        IIDOf(u32bit s)    { return((s >= _useListLen) ? ~u64bitZERO : _useList[s].iid);    };

  //  For a chain position p, returns the s (above) for that position.
  u32bit        sequenceNumberOfPosition(u64bit p);


private:
  struct use_s {
    u32bit   iid;
    u32bit   length;
    u64bit   start;
  };

  char        *_filename;
  seqFile     *_file;
  seqFile     *_fileToDelete;
  u64bit       _fileTimeStamp;

  u32bit       _currentSeq;  //  index into useList of the current sequence
  u32bit       _currentPos;  //  position in the current sequence
  seqOnDisk   *_sequence;    //  the current sequence

  u32bit       _useListLen;
  u32bit       _useListMax;
  use_s       *_useList;

  u64bit       _lengthOfSequences;

  bool         _eof;

  char         _separator;
  bool         _separatorDone;
  u32bit       _separatorLength;
  u32bit       _separatorPosition;
};
