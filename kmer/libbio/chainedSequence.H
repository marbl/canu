#ifndef CHAINEDSEQUENCE_H
#define CHAINEDSEQUENCE_H


////////////////////////////////////////
//
//  chainedSequence
//
//  Build a chained sequence from a source file or source
//  FastABase.  In the case of a file, a FastABase (with iid
//  index) is opened.
//
//  It "assembles" sequences, separated by a fixed number of any
//  single character (say 10 N's) and allows streaming access, or
//  random block access.
//
class chainedSequence {
public:
  chainedSequence();
  ~chainedSequence();

  //  Sets the separator and the length of
  //
  void         setSeparator(char sep);
  void         setSeparatorLength(u32bit len);

  //  Sets the source file
  //
  void         setSource(char const *filename);
  void         setSource(FastABase *file);

  //  Adds sequences to the useList.  For parse(), sequences are
  //  specified in a character string, e.g., "10,20,30-40,50".  add()
  //  pushes one sequence iid onto the list.
  //
  //  Multiple calls to parse() and add() are allowed.  There is
  //  (currently) no way to remove sequences from the useList.
  //
  void         parse(char *line);
  void         add(u32bit v);
  void         finish(void);

  //  get() returns 0 at the end of the chained sequence.  It operates
  //  like getc() -- the stream is advanced after each call.
  //
  //  eof() returns true if the last get() returned eof()
  //
  char         get(void);
  bool         eof(void) { return(_eof); };

  //  These positions are pre-incremented, so they're off by one.

  u64bit       thePositionInSequence(void) { return(_positionInSequence - 1); };
  u64bit       thePositionInStream(void)   { return(_positionInStream - 1);   };

  //  This returns the sequence number in the chain, not the IID of the sequence.  See IIDOf().
  u64bit       theSequenceNumber(void) {
    return(_currentSeq);
  };

  //  Return the length of, position of (in the chain) and IID of the
  //  (s)th sequence in the chain.
  //
  u32bit       numberOfSequences(void) { return(_useListLen); };
  u64bit       lengthOfSequences(void) { return(_lengthOfSequences); };
  u64bit       lengthOf(u32bit s);
  u64bit       startOf(u32bit s);
  u64bit       IIDOf(u32bit s);

  u32bit       sequenceNumberOfPosition(u64bit p);

  //  Move the file pointer back to the start of the chained sequence
  //
  bool         rewind(void);

  //  Save everything, even the current position in the sequence.
  //
  void         saveState(char const *filename);
  bool         loadState(char const *filename, bool beNoisy=false, bool loadData=true);

private:
  struct use_s {
    u32bit   iid;
    u32bit   length;
    u32bit   start;
  };

  char                 *_filename;
  FastABase            *_file;
  FastABase            *_fileToDelete;
  FastASequenceOnDisk  *_sequence;

  u32bit                _useListLen;
  u32bit                _useListMax;
  use_s                *_useList;

  //  The position we are at in the stream -- each sequence knows
  //  where it is at, so we just need to remember which sequence we
  //  are at.
  //
  u32bit                _currentSeq;

  //  Various things the user will find useful
  //
  u64bit                _positionInSequence;
  u64bit                _positionInStream;
  u64bit                _lengthOfSequences;

  //  Geez, very obnoxiously, we can tell if the _next_ get()
  //  is going to return EOF, but we cannot tell if the _last_
  //  get() was an eof.  So we set a flag.
  //
  bool                  _eof;

  //  If _separatorPosition is positive, we are sending
  //  _separatorLength _separator's out.
  //
  bool                  _separatorDone;
  u32bit                _separatorLength;
  u32bit                _separatorPosition;

  char                  _separator;
};



inline
char
chainedSequence::get(void) {
  char ret = _sequence->get();

  _positionInStream++;

  if (ret == 0) {

    //  Are we at the end of the chained sequence?
    //
    if (_currentSeq + 1 >= _useListLen) {
      _eof = true;
      return(0);
    }

    //  Nope, are we doing the separator?
    //
    if (_separatorPosition > 0) {
      _separatorPosition--;
      return(_separator);
    }

    //  Nope, did we just finish the separator, or do we need to start
    //  it?  The flow being: the first time _sequence->get() returns
    //  0, _separatorDone is false, and we start the separator.  We do
    //  the separator.  It finishes, and we're back here, with
    //  _separatorDone set to true.
    //
    //  On _separatorDone, we initialize the next sequence, load the
    //  first letter and fall through to the normal "got a good read"
    //  exit.
    //
    if (_separatorDone) {
      delete _sequence;
      _currentSeq++;

      _file->find(_useList[_currentSeq].iid);

      _sequence           = _file->getSequenceOnDisk();
      _separatorDone      = false;
      _positionInSequence = u64bitZERO;
      ret = _sequence->get();
    } else {
      _separatorPosition  = _separatorLength;
      _separatorDone      = true;
      _positionInSequence = ~u64bitZERO;
      return(_separator);
    }
  }

  //  Got a good read, bump to the next character and return
  //
  _sequence->next();
  _positionInSequence++;
  return(ret);
}



#endif  //  CHAINEDSEQUENCE_H
