#ifndef BRI_PLUS_PLUS_H
#define BRI_PLUS_PLUS_H



#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>
#include <fcntl.h>
#include <string.h>
#include <ctype.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "bri.h"

//  C++ objects for libbri



//  These are all inlined, and C doesn't want to listen to that.
//
#include "bitOperations.h"
#include "bitPacking.h"

//  Various methods for encoding numbers into a bitstream.
//
//  Still missing:
//    minimal binary
//    golomb (actually rice, since power of two)
//    teuhola exponential golomb
//
//  And a nice way of getting parameters to those (and generalizedUnary)
//
#include "unaryEncoding.h"
#include "generalizedUnaryEncoding.h"
#include "eliasGammaEncoding.h"
#include "eliasDeltaEncoding.h"
#include "fibonacciEncoding.h"

#include "mers.h"

////////////////////////////////////////
//
//  britime
//

class speedCounter {
public:
  //  fmt specifies the status format.  And example:
  //    "    %8f [unit]things (%8.5f [unit]things/sec)\r"
  //
  speedCounter(char const   *fmt,
               double        unit,
               u64bit        freq,
               bool          enabled=true);
  ~speedCounter();

  void   tick(void) {
    if (_enabled && ((++_count & _freq) == u64bitZERO)) {
      double  v = _count / _unit;
      fprintf(stderr, _fmt, v, v / (getTime() - _startTime));
      fflush(stderr);
    }
  };

  void   finish(void) {
    if (_enabled && (_count >= _freq)) {
      double  v = _count / _unit;
      fprintf(stderr, _fmt, v, v / (getTime() - _startTime));
      fprintf(stderr, "\n");
      fflush(stderr);
    }
    _count = 0;
  };

private:
  u64bit          _count;
  double          _unit;
  u64bit          _freq;
  double          _startTime;
  char const     *_fmt;
  bool            _enabled;
};


////////////////////////////////////////
//
//  readBuffer
//
class readBuffer {
public:
  //  If bufferMax is zero, then the file is accessed using memory
  //  mapped I/O.  Otherwise, a small buffer is used.
  //
  readBuffer(const char *filename, u32bit bufferMax = 32 * 1024);

  //  This constructor always uses a small buffer; memory mapped I/O
  //  is not allowed.  This lets us pass in fileno(stdin).
  //
  readBuffer(int fileptr, const char *filename, u32bit bufferMax = 32 * 1024);

  ~readBuffer();

  bool      eof(void);
  bool      next(void);
  char      get(void);
  char      getnext(void);
  bool      seek(off_t pos);
  size_t    read(char *buf, size_t len);   //  read the next len bytes into the user buffer buf
  off_t     tell(void);

private:
  void      fillBuffer(void);
  void      init(int fileptr, const char *filename, u32bit bufferMax);

  char    *_filename;

  int      _file;
  int      _fileType;  //  0=stdin, 1=normal, 2=mmap
#ifdef __APPLE__
  int      _pad1;
#endif
  off_t    _filePos;

  bool     _eof;

  //  If bufferMax is zero, then we are using the mmapped interface, otherwise,
  //  we are using a open()/read() and a small buffer.

  size_t   _bufferPos;
  size_t   _bufferLen;
  size_t   _bufferMax;
  char    *_buffer;
#ifdef __APPLE__
  int      _pad2;
#endif
};


inline
bool
readBuffer::eof(void) {
  return(_eof);
}


inline
bool
readBuffer::next(void) {

  if (_eof)
    return(true);

  switch (_fileType) {
    case 0:
    case 1:
      _bufferPos++;
      _filePos++;

      if (_bufferPos >= _bufferLen)
        fillBuffer();
      break;
    case 2:
      _bufferPos++;
      if (_bufferPos >= _bufferLen)
        _eof = true;
      break;
  }

  return(_eof);
}


inline
char
readBuffer::get(void) {
  return(_buffer[_bufferPos]);
}


inline
char
readBuffer::getnext(void) {
  char  x = _buffer[_bufferPos];
  next();
  return(x);
}


inline
off_t
readBuffer::tell(void) {

  if (_fileType == 2)
    return((size_t)_bufferPos);

  return(_filePos);
}




////////////////////////////////////////
//
//  splittowords
//
class splitToWords {
public:
  splitToWords(char *cmd) {


    //  Step Zero:
    //
    //  Count the length of the string, in words and in characters.
    //  For simplicity, we overcount words, by just counting white-space.
    //
    //  Then, allocate space for a temporary copy of the string, and a
    //  set of pointers into the temporary copy (much like argv).
    //
    u32bit   cmdChars = 1;  //  1 == Space for terminating 0
    u32bit   cmdWords = 2;  //  2 == Space for first word and terminating 0L

    for (char *tmp=cmd; *tmp; tmp++) {
      cmdWords += *tmp == ' ';
      cmdWords += *tmp == '\t';
      cmdChars++;
    }

    _cmd      = new char   [cmdChars];
    _arg      = new char * [cmdWords];
    _argWords = 0;


    //  Step One:
    //
    //  Determine where the words are in the command string, copying the
    //  string to _cmd and storing words in _arg.
    //
    bool           isFirst  = true;
    char          *cmdI = cmd;
    char          *cmdO = _cmd;

    while (*cmdI) {

      //  If we are at a non-space character, we are in a word.  If
      //  this is the first character in the word, save the word in
      //  the args list.
      //
      //  Otherwise we are at a space and thus not in a word.  Make
      //  all spaces be string terminators, and declare that we are
      //  at the start of a word.
      //
      if ((*cmdI != ' ') && (*cmdI != '\t')) {
        *cmdO = *cmdI;

        if (isFirst) {
          _arg[_argWords++] = cmdO;
          isFirst           = false;
        }
      } else {
        *cmdO   = 0;
        isFirst = true;
      }

      cmdI++;
      cmdO++;
    }

    //  Finish off the list by terminating the last arg, and
    //  terminating the list of args.
    //
    *cmdO           = 0;
    _arg[_argWords] = 0L;
  };

  ~splitToWords() {
    delete [] _cmd;
    delete [] _arg;
  };

  u32bit  numWords(void)    { return(_argWords); };
  char   *getWord(u32bit i) { return(_arg[i]); };
private:
  char     *_cmd;
  char    **_arg;
  u32bit    _argWords;
};


////////////////////////////////////////
//
//  bit-packing
//







////////////////////////////////////////
//
//  bitPackedFile
//
class bitPackedFile {
public:
  bitPackedFile(char const *name);
  ~bitPackedFile();

  u64bit     getBits(u32bit size);
  u64bit     getNumber(void);

  void       putBits(u64bit bits, u32bit size);
  void       putNumber(u64bit val);

  void       seek(u64bit pos, bool forceLoad=false);

  void       showStats(FILE *f) {
    fprintf(f, "inside: "u64bitFMT"  outside: "u64bitFMT"\n", stat_seekInside, stat_seekOutside);
    fflush(f);
  };
private:

  //  Ensure that the buffer has enough space for any future
  //  operation.  If the buffer doesn't have 128 bits free (2 words
  //  currently) flush it (if dirty) and read in a new page.  128 was
  //  chosen because the fibonacci encoded numbers use up to 90-some
  //  bits.
  //
  void       sync(void) {

    //  This constant must be strictly less than the constant used in deciding
    //  if seek() is moving forward or backwards.
    //
    if (((_bit >> 6) + 31) >= _bfrmax)
      seek((_pos << 6) + _bit, true);
  };

  void       flushDirty(void);

  int       _file;

  u64bit   *_bfr;     //  A chunk of the bitPackedFile in core
  u64bit    _bfrmax;  //  Number of words in the buffer
  u64bit    _pos;     //  The location this chunk is from (in words)
  u64bit    _bit;     //  The bit position we are modifying relative to _pos

  bool      _bfrDirty;
#ifdef __APPLE__
  u32bit    _pad;
#endif

  //  For collecting statistics on our usage
  //
  u64bit  stat_seekInside;
  u64bit  stat_seekOutside;
  u64bit  stat_dirtyFlushes;
};




////////////////////////////////////////
//
//  fastastream
//

//
//  Reads sequence from a multiFastA file, letter by letter, saving
//  the defline.  Returns values are:
//
//  0    EOF
//  254  New sequence starting
//  253  Break in contiguous sequence
//
//  If the letter is a valid symbol the letter is returned.
//  Otherwise, one of code 253 or 254 is returned.
//

class FastAstream {
public:
  FastAstream(const char *filename) {
    if ((filename == 0L) ||
        ((filename[0] == '-') && (filename[1] == 0)))
      _rb = new readBuffer(fileno(stdin), "stdin");
    else
      _rb = new readBuffer(filename);

    _thePositionInSequence  = 0;
    _thePositionInStream    = 0;
    _theSequenceNumber      = 0;

    _theDefLineMax          = 2048;
    _theDefLinePos          = 0;
    _theDefLine             = new char [_theDefLineMax + 1];
  };

  ~FastAstream() {
    delete    _rb;
    delete [] _theDefLine;
  };

  bool              rewind(void) {
    _thePositionInSequence  = 0;
    _thePositionInStream    = 0;
    _theSequenceNumber      = 0;
    return(_rb->seek(0));
  };

  unsigned char     nextSymbol(void);

  u64bit            thePositionInSequence(void) {
    return(_thePositionInSequence);
  };
  u64bit            thePositionInStream(void) {
    //  See the comment below why there is a -1 here.
    return(_thePositionInStream - 1);
  };
  u64bit            theSequenceNumber(void) {
    return(_theSequenceNumber);
  };
  char             *theDefLine(void) {
    return(_theDefLine);
  };
private:
  u64bit       _thePositionInSequence;
  u64bit       _thePositionInStream;
  u64bit       _theSequenceNumber;

  readBuffer  *_rb;

  u32bit       _theDefLineMax;
  u32bit       _theDefLinePos;
  char        *_theDefLine;
};


inline
unsigned char
FastAstream::nextSymbol(void) {
  unsigned char  nc;

  //  Grab the next character
  //
  nc = (unsigned char)_rb->getnext();

  //  Skip any whitespace
  //
  while (isspace(nc) && !_rb->eof())
    nc = (unsigned char)_rb->getnext();

  //  Are we at the EOF?
  //
  if (_rb->eof())
    return(0);

  //  Looks like we are in the sequence.  Increment the position.  If
  //  this is the start of a defline, the position is reset to zero.
  //
  //  The stream position is off by one because we pre-increment.  It
  //  is adjusted on access.
  //  
  _thePositionInSequence++;
  _thePositionInStream++;

  //  If it's valid, just return it.
  //
  if (validSymbol[nc])
    return(nc);

  //  It's not a valid character.  If it's a defline, read until we
  //  hit a '\n' or '\r'.
  //
  if (nc == '>') {

    //  Reset thePosition and increment theSequenceNumber
    //
    _thePositionInSequence = 0;
    _theSequenceNumber++;

    _theDefLinePos = 0;

    nc = (unsigned char)_rb->getnext();
    while ((nc != '\n') && (nc != '\r') && !_rb->eof()) {

      //  More space?  Must be a monster defline!
      //
      if (_theDefLinePos >= _theDefLineMax) {
        _theDefLineMax <<= 1;
        char *d = new char [_theDefLineMax + 1];

        for (u32bit i=_theDefLinePos; --i; )
          d[i] = _theDefLine[i];

        delete [] _theDefLine;
        _theDefLine = d;
      }

      //  Push the letter onto the end of the defline
      //
      _theDefLine[_theDefLinePos++] = nc;

      nc = (unsigned char)_rb->getnext();
    }

    _theDefLine[_theDefLinePos] = 0;

    //  Return 254; we've just read in the defline.
    //
    return(254);
  }

  //  We've hit a break, return 253.
  //
  return(253);
}


////////////////////////////////////////
//
//  merstreamfile
//
//  Uses the merStream to build a compressed merStreamFile.  The
//  merStreamFile contains all the forward mers that the stream
//  returns.  It does not contain position information or deflines.
//  Just the mers, ma'am.
//
//  To get position, sequence number and defline:
//
//  Each mer block can store the position start, and the sequence
//  number.  If we then store the deflines, we get everything!


//  Returns true if the merStreamFile exists.
bool  merStreamFileExists(const char *i);

class merStream;

class merStreamFileBuilder {
public:
  merStreamFileBuilder(u32bit merSize, const char *inputFile, const char *outputFile);
  merStreamFileBuilder(merStream *MS, const char *outputFile);
  ~merStreamFileBuilder();

  //  Returns the number of mers written
  u64bit                 build(bool beVerbose=false);

private:
  u32bit                _merSize;
  merStream            *_merStream;
  const char           *_inputFile;
  const char           *_outputFile;
};


class merStreamFileReader {
public:
  merStreamFileReader(const char *filename);
  ~merStreamFileReader();

  u32bit                 merSize(void)      { return(_merSize); };
  u64bit                 numberOfMers(void) { return(_numMers); };

  u64bit                 theFMer(void) { return(_theFMer); };
  u64bit                 theRMer(void) { return(_theRMer); };

  char const            *theFMerString(void);
  char const            *theRMerString(void);

  u64bit                 thePositionInSequence(void) { return(_thePosInSeq); };
  u64bit                 thePositionInStream(void)   { return(_thePosInStr); };
  u64bit                 theSequenceNumber(void)     { return(_theSequence); };
  char const            *theDefLine(void)            { return(_theDefline);  };

  bool                   seekToMer(u64bit merNumber);
  bool                   seekToSequence(u64bit seqNumber);

  //  This will limit the number of mers returned to 'limit'.  The
  //  seek*() methods do not reset the count or the limit.  The
  //  current count is reset when this function is called.  Set to
  //  ~u64bitZERO for an infinite limit.  Zero really does mean zero,
  //  not infinity.
  //
  void                   setIterationLimit(u64bit limit) {
    _iterationLimit = limit;
    _iteration      = u64bitZERO;
  };

  bool                   nextMer(u32bit skip=0);
private:
  const char           *_inputFile;

  bitPackedFile        *_streamFile;

#ifdef __APPLE__
  u32bit                _pad;
#endif

  u32bit                _merSize;
  u64bit                _numMers;
  u64bit                _numBlocks;
  u64bit                _numDefs;
  u64bit                _defLength;
  off_t                 _blockFileSize;
  off_t                 _deflineFileSize;
  off_t                 _streamFileSize;
  off_t                 _blkStart;
  off_t                 _defStart;
  off_t                 _strStart;

  u32bit               *_blockSize;
  u32bit               *_blockSequence;
  u64bit               *_blockPosInSeq;
  u64bit               *_blockPosInStr;
  char                 *_deflineStorage;
  char                **_deflines;

  //  Things used to decode the stream
  //
  u64bit                _thisBlock;
  u64bit                _thisBlockSize;
  u64bit                _merMask;
  u64bit                _firstMer;

  //  Things used to stop the stream
  //
  u64bit                _iterationLimit;
  u64bit                _iteration;

  //  Things returned by the accessors.  These are set/updated by
  //  nextMer().
  //
  u64bit                _theFMer;
  u64bit                _theRMer;
  u64bit                _thePosInSeq;
  u64bit                _thePosInStr;
  u64bit                _theSequence;

  char                 *_theDefline;

  char                  _theMerString[36];
};




////////////////////////////////////////
//
//  This is the merStream.  Given a file of sequence and deflines, or
//  a string of sequence, it returns mers.
//
//  The file based merstream uses a fastastream to get sequence, and
//  that doesn't care if there is whitespace or deflines in the
//  stream.  It returns only valid sequence.
//
//  If the sequence/length method is used, we assume that 'seq'
//  contains no whitespace, and no extra characters (e.g., deflines),
//  just pure sequence.  It will work if you give it deflines and
//  whitespcae, but positions will be wrong.
//
class merStream {
public:
  merStream(merStreamFileReader *merstreamfile);
  merStream(u32bit merSize, const char *filename);
  merStream(u32bit merSize, const char *seq, u32bit len);
  ~merStream();


  //  The merstream has three modes of operation, depending
  //  on what the input format is.  We could:
  //  1)  subclass an interface template
  //  2)  use function pointers set at create time
  //  3)  use a set of if tests in each call that cares
  //
  //  1 and 2 are essentially the same, except 2 can use human
  //  intelligence to get rid of the function call on a few methods
  //  that aren't input specific.
  //
  //  One can claim that 3 is better because it doesn't need
  //  to mess with the stack, and if inlined, branch prediction
  //  will remove the overhead of the if tests.


  u64bit       theFMer(void)        { return(_theFMer); };
  u64bit       theRMer(void)        { return(_theRMer); };

  char const  *theFMerString(void);
  char const  *theRMerString(void);

  u32bit       theMerSize(void) { return(_merSize); };

  bool         rewind(void) {
    bool  ret = false;


    if (_theMers) {
      ret = _theMers->seekToMer(0);
      _theFMer = _theMers->theFMer();
      _theRMer = _theMers->theRMer();
    } else if (_theFile) {
      ret = _theFile->rewind();
      loadMer(_merSize - 1);
    } else {
      ret = true;
      _theStringPos = 0;
      _thePosInSeq  = 0;
      _thePosInStr  = 0;
      _theNum       = 0;
      loadMer(_merSize - 1);
      _thePosInStr--;
    }

    return(ret);
  };

  u64bit       thePositionInSequence(void) {
    if      (_theMers)
      return(_theMers->thePositionInSequence());
    else if (_theFile)
      return(_theFile->thePositionInSequence() - _merSize);
    else
      return(_thePosInSeq - _merSize);
  };

  u64bit       thePositionInStream(void) {
    if      (_theMers)
      return(_theMers->thePositionInStream());
    else if (_theFile)
      return(_theFile->thePositionInStream() - _merSize);
    else
      return(_thePosInStr - _merSize);
  };

  u64bit       theSequenceNumber(void) {
    if      (_theMers)
      return(_theMers->theSequenceNumber());
    else  if (_theFile)
      return(_theFile->theSequenceNumber());
    else
      return(_theNum);
  };

  char const  *theDefLine(void) {
    if      (_theMers)
      return(_theMers->theDefLine());
    else if (_theFile)
      return(_theFile->theDefLine());
    else
      return(">");
  };

  bool               nextMer(u32bit skip=0);
private:
  void               loadMer(u32bit s);

  char               nextSymbol(void) {
    if (_theFile) {
      return(_theFile->nextSymbol());
    } else {
      //  See also the constructor for character strings -- it skips
      //  white space and the first defline.

      //  Skip any whitespace
      //
      while (isspace(_theString[_theStringPos]) && (_theString[_theStringPos]) && (_theStringPos < _theStringLen))
        _theStringPos++;

      //  Are we at a defline?  If so, zip along until we hit a
      //  newline.  Then we return the code (fastastream) for a
      //  sequence break.
      //
      if (_theString[_theStringPos] == '>') {
        while ((_theString[_theStringPos] != '\n') && (_theString[_theStringPos]) && (_theStringPos < _theStringLen))
          _theStringPos++;

        _thePosInSeq = 0;
        _thePosInStr++;

        _theNum++;
        return(254);
      }

      //  Are we at the EOF?
      //
      if ((_theString[_theStringPos] == 0) || (_theStringPos >= _theStringLen))
        return(0);

      //  Update the position of the character in this sequence, and
      //  the positon in the stream.
      //
      _thePosInSeq++;
      _thePosInStr++;

      return(_theString[_theStringPos++]);
    }
  };


  //  Data specific to the merStreamFileReader input
  //
  merStreamFileReader  *_theMers;

  //  Data specific to the FastA file input
  //
  FastAstream          *_theFile;

  //  Data specific to the character string input
  //
  char const           *_theString;
  u64bit                _theStringPos;
  u64bit                _theStringLen;
  u64bit                _thePosInSeq;
  u64bit                _thePosInStr;
  u32bit                _theNum;

  //  Data specific to nobody
  //
  u32bit                _merSize;
  s32bit                _timeUntilValid;
  u32bit                _theRMerShift;

  u64bit                _theMerMask;

  u64bit                _theFMer;
  u64bit                _theRMer;

  char                  _theMerString[64];
};





//  We have a problem; at the start of a stream, we want to
//  initialize the mer with merSize-1 bases.  In the middle of the
//  stream, we need to load the mer with merSize bases.
//
//  loadMer will push s bases onto the mer, restarting if it hits a
//  mer break.
//
//  No masking of the mer is performed.
//
inline
void
merStream::loadMer(u32bit s) {
  u64bit   ch = 255;
  u64bit   cf = 0;
  u64bit   cr = 0;

  _timeUntilValid = s;

  //  While we are invalid, and still in the sequence
  //  push characters onto the mer.  The valid
  //  count is updated if we hit an invalid base.
  //
  while ((_timeUntilValid > 0) && (ch != 0)) {
    ch = nextSymbol();
    cf = validCompressedSymbol[ch];

    //  Rather than take the chance of generating a cache miss accessing
    //  another array, we just mask out the upper bits of validCompressedSymbol[];
    //  this is exactly the same as using compressSymbol[].
    //
    //  We need to mask the upper bits for reverse (but not for forward)
    //  because, in reverse, we shift to the right.  If we don't mask theseq
    //  out, we will have extra bits set in the mer.
    //
    //  Example: Consider placing 255 (== 11111111, the invalid symbol returned
    //  from validCompressedSymbol) at the fourth base in a four mer:
    //    000000**xxxxxx -- the ** are is the fourth base.
    //
    //  Without masking, we'd set all the 0's and all the *'s to one, in effect,
    //  preloading the next three bases of the mer.
    //
    cr = validCompressedSymbol[complementSymbol[ch]] & 0x03;

    _timeUntilValid--;

    if (cf & 0xfc)
      _timeUntilValid = s;

    //  If the ch is valid, we can obviously add it to the mer.
    //  If it is invalid, we don't care; by the time the mer is valid
    //  all bits of any invalid mers are removed.  Sure, we could
    //  put this into the else, but I suspect that it will be faster outside.
    //
    _theFMer <<= 2;
    _theRMer >>= 2;

    _theFMer |= cf;
    _theRMer |= cr << _theRMerShift;
  }
}




//  how to get thePositionInStream if we don't know the length of ....
//  we need to add one for each character we get from the source
//    does the fasta stream eat characters?
//    do we really need to add separators between sequences?
//


inline
bool
merStream::nextMer(u32bit skip) {

  if (_theMers) {
    bool ret = _theMers->nextMer(skip);

    _theFMer = _theMers->theFMer();
    _theRMer = _theMers->theRMer();

    return(ret);
  } else {
    do {
      u64bit  ch = nextSymbol();
      u64bit  cf = validCompressedSymbol[ch];
      u64bit  cr = validCompressedSymbol[complementSymbol[ch]] & 0x03;

      //  EOF?
      if (ch == 0)
        return(false);

      //  Push the ch onto the mer.
      //
      _theFMer <<= 2;
      _theRMer >>= 2;

      _theFMer |= cf;
      _theRMer |= cr << _theRMerShift;

      //  If the ch is invalid, we need to make a whole new mer.
      //  This also implies that we're done skipping.
      //
      if (cf & (unsigned char)0xfc) {
        loadMer(_merSize);
        skip = 0;
      }
    } while (skip--);

    _theFMer &= _theMerMask;
    _theRMer &= _theMerMask;

    //  Still need to check if we are valid -- we could
    //  have run off the end of the sequence before a valid
    //  mer was created.
    //
    return(_timeUntilValid <= 0);
  }
}








////////////////////////////////////////
//
//  intervalList
//

//
//  An interval list -- 
//

//  I don't want this to be public, but the sort function
//  in intervalList.C needs it so.
//

struct _intervalPair {
  u32bit    lo;
  u32bit    hi;
};



class intervalList {
public:
  intervalList();
  ~intervalList();

  intervalList &operator=(intervalList &src);


  //  Clear a list
  void        clear(void) {
    _isSorted = true;
    _listLen  = 0;
  }

  //  Insert a new interval into the list
  void        add(u32bit position, u32bit length);

  //  Sort the set of intervals by the lo value
  void        sort(void);

  //  Merge overlapping or adjacent intervals together.
  void        merge(void);

  //  Returns the number of intervals
  u32bit      numberOfIntervals(void) {
    return(_listLen);
  };

  //  Returns the sum of the length of all intervals
  u32bit      sumOfLengths(void) {
    u32bit len = 0;
    u32bit i   = numberOfIntervals();

    if (i > 0)
      while (i--)
        len += _list[i].hi - _list[i].lo;

    return(len);
  };

  u32bit      lo(u32bit i) {
    return(_list[i].lo);
  };
  u32bit      hi(u32bit i) {
    return(_list[i].hi);
  };

private:
  bool                      _isSorted;
  u32bit                    _listLen;
  u32bit                    _listMax;
  _intervalPair            *_list;
};


#include "merCovering.H"


#endif  //  BRI_PLUS_PLUS_H

