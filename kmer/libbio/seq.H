#ifndef SEQ_H
#define SEQ_H

#include "util++.H"
#include "alphabet.h"

class seqInCore;
class seqOnDisk;

class seqFile;
class seqStream;
class seqFactory;

typedef  u32bit  seqIID;

class seqInCore {
public:
  seqInCore(seqIID iid);
  seqInCore(seqIID iid, char *hdr, u32bit hdrlen, char *seq, u32bit seqlen);
  ~seqInCore();

  seqInCore *copy(void);

  char           *header(void) const         { return(_header); };
  u32bit          headerLength(void) const   { return(_headerLen); };

  char           *sequence(void) const       { return(_seq); };
  u32bit          sequenceLength(void) const { return(_seqLen); };

  u32bit          getIID(void) const         { return(_idx); };
private:
  u32bit         _headerLen;
  u32bit         _headerMax;
  char          *_header;

  u32bit         _seqLen;
  u32bit         _seqMax;
  char          *_seq;

  u32bit         _idx;
};


////////////////////////////////////////////////////////////////////////////////
//  Leave the sequence mostly on disk.  A small buffer is used to
//  cache reads.  The defline is read into core.
//
//  Two operations are currently supported:
//    char  getChar(void)
//    char *getChars(char *block, u32bit position, u32bit length)
//
//  When the end of the sequence is reached, getChar() returns 0.
//
//  Only supported for indexed and squeezed fasta file sources
//
//  Three types of sequence
//   -- fully squeezed: use a readBuffer directly, both stream and block copy
//   -- fixed width: stream and block copy
//   -- free form: only stream
//
class seqOnDisk {
public:
  seqOnDisk(char const *filename,
            u64bit hdrstart, u32bit hdrlen,
            u64bit seqstart, u32bit seqlen,
            seqIID iid,
            bool isSqueezed,
            bool isFixedWidth,
            u32bit sl,
            u32bit ss);
  seqOnDisk(seqIID iid,
            char *hdr, u32bit hdrlen,
            char *seq, u32bit seqlen);
  ~seqOnDisk();

  char           *header(void) const         { return(_header);       };
  u32bit          headerLength(void) const   { return(_headerLength); };

  char           *sequence(void) const       { return(_sequence);       };
  u32bit          sequenceLength(void) const { return(_sequenceLength); };

  u32bit          getIID(void) const         { return(_idx); };

  //  XXX:  Probably fails if the sequence ends with whitespace
  //#warning seqOnDisk might fail if the sequence ends with whitespace!

  //  next() returns true on eof.
  //  get() returns 0 on eof.
  //
  bool            next(void);
  char            get(void);
  char           *getChars(char *block, u32bit position, u32bit length);
private:

  readBuffer    *_readBuffer;

  //  The iid of this guy
  //
  u32bit         _idx;

  //  Information about the sequence on disk
  //
  u32bit         _headerLength;
  u32bit         _sequenceLength;
  u64bit         _headerStart;
  u64bit         _sequenceStart;

  //  A copy of the defline.  We get it for nearly free when we make
  //  the first read buffer.
  //
  char          *_header;

  //  If defined, this means we are not readBuffer backed, but
  //  we are faking it.
  //
  char          *_sequence;

  //  The position we are at in the sequence.  Used to decide when we
  //  have seen all the sequence.  If the source is not squeezed, this
  //  is not the same as the file position.
  //
  u32bit         _sequencePosition;

  //  Cache the type of the source
  //    0 - squeezed
  //    1 - fixed width
  //    2 - freeform
  //
  u32bit         _sourceType;

  //  If the source file is not squeezed, we'll do on-the-fly
  //  squeezing.  If the source file is fixed width, we know the
  //  length of each line, and the number of bytes between lines.
  //
  u32bit         _lineLength;
  u32bit         _lineSep;
};



//  Advances to the next non-space letter.
//  Returns false if there is no next non-space letter.
//
inline
bool
seqOnDisk::next(void) {

  //  Check if we're already at the end
  //
  if (_sequencePosition >= _sequenceLength)
    return(false);

  //  Nope, we know there is another non-space letter out there
  //  somewhere.

  if (_readBuffer) {
    _readBuffer->next();

    if (_sourceType != 0)
      while (whitespaceSymbol[_readBuffer->get()])
        _readBuffer->next();    
  }

  _sequencePosition++;

  return(true);
}


//  Returns 0 on eof (or end of sequence).  Does not advance the
//  file position.
//
inline
char
seqOnDisk::get(void) {

  if (_sequencePosition >= _sequenceLength)
    return(0);

  if (_readBuffer)
    return(_readBuffer->get());
  else
    return(_sequence[_sequencePosition]);
}

////////////////////////////////////////////////////////////////////////////////



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


class seqStream {
private:
  void          clearGuts(void);
public:
  seqStream();
  seqStream(const char *filename, bool finishMe);
  seqStream(seqFile *S, bool finishMe);
  ~seqStream();

  void          setFile(const char *filename);
  void          setFile(seqFile *S);

  void          setSeparator(char sep, u32bit len) {
    _separator = sep;
    if (len == 0)
      fprintf(stderr, "seqStream::setSeparator()-- ERROR: the separator length must be at least one.\n"), exit(1);
    _separatorLength = len - 1;    
  };
  char          getSeparator(void) {
    return(_separator);
  };

  //  Adds sequences to the useList.  For parse(), sequences are
  //  specified in a character string, e.g., "10,20,30-40,50".  add()
  //  pushes one sequence iid onto the list.  Multiple calls to
  //  parse() and add() are allowed.  There is (currently) no way to
  //  remove sequences from the useList.  Once you're all done adding
  //  sequences, call finish().
  //
  void          parse(char *line);
  void          add(u32bit v);
  void          finish(void);

  bool          rewind(void);
  bool          eof(void)                   { return(_eof); };
  unsigned char get(void);         //  unfiltered
  unsigned char nextSymbol(void);  //  filters out runs of invalid

  u64bit        thePositionInSequence(void) { return(_positionInSequence); };
  u64bit        thePositionInStream(void)   { return(_positionInStream);   };
  u64bit        theSequenceNumber(void)     { return(_theSequenceNumber);  };

  //  Return the length of, position of (in the chain) and IID of the
  //  (s)th sequence in the chain.
  //
  u32bit        numberOfSequences(void) { return(_useListLen); };
  u64bit        lengthOfSequences(void) { return(_lengthOfSequences); };
  u64bit        lengthOf(u32bit s) { return((s >= _useListLen) ? ~u32bitZERO : _useList[s].length); };
  u64bit        startOf(u32bit s)  { return((s >= _useListLen) ? ~u32bitZERO : _useList[s].start);  };
  u64bit        IIDOf(u32bit s)    { return((s >= _useListLen) ? ~u32bitZERO : _useList[s].iid);    };

  //  For a chain position p, returns the s (above) for that position.
  u32bit        sequenceNumberOfPosition(u64bit p);

private:

  struct use_s {
    u32bit   iid;
    u32bit   length;
    u32bit   start;
  };

  char                 *_filename;
  seqFile              *_file;
  seqFile              *_fileToDelete;
  seqOnDisk            *_sequence;

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
  u64bit                _theSequenceNumber;
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


class seqFactory {
protected:
  seqFactory();
  ~seqFactory();

public:
  static seqFactory *instance(void) {
    if (me == 0L)
      me = new seqFactory;
    return(me);
  };
  void           registerFile  (seqFile   *f);
  seqFile       *openFile  (const char *name);
private:
  static seqFactory  *me;

  u32bit         _filesNum;
  u32bit         _filesMax;
  seqFile      **_files;
};

#define openSeqFile(S)    seqFactory::instance()->openFile((S))

#endif  //  SEQ_H
