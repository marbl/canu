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
//
//  Things included here used to be in
//
//    #include "britime.H"
//    #include "readBuffer.H"
//    #include "splittowords.H"
//    #include "bit-operations.H"
//    #include "bit-packing.H"
//    #include "bitPackedFile.H"
//    #include "fastastream.H"
//    #include "merstream.H"
//



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

  void   tick(void) {
    if (_enabled && ((++_count & _freq) == u64bitZERO)) {
      double  v = _count / _unit;
      fprintf(stderr, _fmt, v, v / (getTime() - _startTime));
      fflush(stderr);
    }
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
  void      seek(off_t pos);
  size_t    read(char *buf, size_t len);   //  read the next len bytes into the user buffer buf
  off_t     tell(void);

private:
  void      fillBuffer(void);
  void      init(int fileptr, const char *filename, u32bit bufferMax);

  char    *_filename;

  int      _file;
  int      _fileType;  //  0=stdin, 1=normal, 2=mmap
  off_t    _filePos;

  bool     _eof;

  //  If bufferMax is zero, then we are using the mmapped interface, otherwise,
  //  we are using a open()/read() and a small buffer.

  size_t   _bufferPos;
  size_t   _bufferLen;
  size_t   _bufferMax;
  char    *_buffer;
};


inline
void
readBuffer::fillBuffer(void) {

  _bufferPos = 0;
  errno = 0;
  _bufferLen = (u32bit)::read(_file, _buffer, _bufferMax * sizeof(char));
  if (errno) {
    fprintf(stderr, "readBuffer()-- couldn't read %d bytes from '%s'\n%s\n",
            (int)(_bufferMax * sizeof(char)), _filename, strerror(errno));
    exit(1);
  }
  if (_bufferLen == 0)
    _eof = true;
}


inline
void
readBuffer::init(int fileptr, const char *filename, u32bit bufferMax) {

  _filename  = new char [strlen(filename) + 1];
  strcpy(_filename, filename);

  if (bufferMax == 0) {
    _file        = 0;
    _fileType    = 2;
    _filePos     = 0;
    _eof         = false;
    _bufferPos   = 0;
    _bufferLen   = 0;
    _bufferMax   = 0;
    _buffer      = (char *)mapFile(_filename, &_bufferLen, 'r');
  } else {
    if (fileptr == -1) {
      _fileType = 1;
      errno = 0;
      fileptr = open(filename, O_RDONLY | O_LARGEFILE);
      if (errno) {
        fprintf(stderr, "readBuffer()-- ERROR: couldn't open the file '%s'.\n%s\n",
                filename, strerror(errno));
        exit(1);
      }
    } else {
      _fileType    = 0;
    }
    _file        = fileptr;
    _filePos     = 0;
    _eof         = false;
    _bufferPos   = 0;          //  Position we are at in the buffer
    _bufferLen   = 0;          //  Valid length of the buffer
    _bufferMax   = bufferMax;  //  Maximum size of the buffer
    _buffer      = new char [_bufferMax];

    fillBuffer();
  }
}


inline
readBuffer::readBuffer(const char *filename, u32bit bufferMax) {
  init(-1, filename, bufferMax);
}


inline
readBuffer::readBuffer(int fileptr, const char *filename, u32bit bufferMax) {
  if (bufferMax == 0) {
    fprintf(stderr, "readBuffer::readBuffer()-- file '%s' requested memory mapped I/O on previously\n", filename);
    fprintf(stderr, "readBuffer::readBuffer()-- opened file, but that's not allowed.  bufferMax reset to 32KB.\n");
    bufferMax = 32 * 1024;
  }
  init(fileptr, filename, bufferMax);
}


inline
readBuffer::~readBuffer() {

  switch (_fileType) {
    case 0:
      delete [] _filename;
      delete [] _buffer;
      break;
    case 1:
      delete [] _filename;
      delete [] _buffer;
      errno = 0;
      close(_file);
      if (errno) {
        fprintf(stderr, "readBuffer()-- WARNING: couldn't close the file '%s'\n%s\n",
                _filename, strerror(errno));
      }
      break;
    case 2:
      delete [] _filename;
      unmapFile(_buffer, _bufferLen);
      break;
  }
}


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
void
readBuffer::seek(off_t pos) {

  switch(_fileType) {
    case 0:
      fprintf(stderr, "readBuffer()-- seek() not available for file '%s'.\n", _filename);
      break;
    case 1:
      errno = 0;
      lseek(_file, pos, SEEK_SET);
      if (errno) {
        fprintf(stderr, "readBuffer()-- '%s' couldn't seek to position "s64bitFMT".\n%s\n",
                _filename, pos, strerror(errno));
        exit(1);
      }

      _filePos = pos;
      _eof     = false;

      //  We really do need to do a fillBuffer() here.  If the next
      //  operation is a get() we should return something valid.
      //
      //  Of course, we could put a check into get() to see if bufferPos
      //  is valid, but I don't want the overhead on EVERY get().
      //
      fillBuffer();
      break;
    case 2:
      if (pos < _bufferLen) {
        _bufferPos = pos;
        _eof       = false;
      } else {
        _eof       = true;
      }
      break;
  }
}


inline
size_t
readBuffer::read(char *buf, size_t len) {

  if (_fileType == 2) {
    size_t c = 0;

    while ((_bufferPos < _bufferLen) && (c < len))
      buf[c++] = _buffer[_bufferPos++];

    return(c);
  } else {
    //  The trick here is to use the existing buffered input first,
    //  then do a direct read to get the rest.
    //
    //  We fill the buffer again if it is empty.
    //
    //  The number of bytes actually put into buf is returned.
      
    size_t   bCopied = 0;   //  Number of bytes copied into the buffer
    size_t   bRead   = 0;   //  Number of bytes read into the buffer
    size_t   bAct    = 0;   //  Number of bytes actually read from disk

    //  Easy case; the next len bytes are already in the buffer; just
    //  copy and move the position.
    //
    //  XXX:  Check the zero-left-in-buffer case
    //
    if (_bufferLen - _bufferPos > len) {
      bCopied = len;
      bRead   = 0;

      memcpy(buf, _buffer + _bufferPos, sizeof(char) * len);
      _bufferPos += (u32bit)len;
    } else {

      //  Existing buffer not big enough.  Copy what's there, then finish
      //  with a read.
      //
      memcpy(buf, _buffer + _bufferPos, (_bufferLen - _bufferPos) * sizeof(char));
      bCopied    = _bufferLen - _bufferPos;
      _bufferPos = _bufferLen;

      while (bCopied + bRead < len) {
        errno = 0;
        bAct = (u32bit)::read(_file, buf + bCopied + bRead, (len - bCopied - bRead) * sizeof(char));
        if (errno) {
          fprintf(stderr, "readBuffer()-- couldn't read %d bytes from '%s'\n%s\n",
                  len * sizeof(char), _filename, strerror(errno));
          exit(1);
        }
        if (bAct == 0) {
          fprintf(stderr, "readBuffer()-- EOF hit while reading %d bytes from '%s'\n%s\n",
                  len * sizeof(char), _filename, strerror(errno));
          exit(1);
        }
        bRead += bAct;
      }
    }

    if (_bufferPos == _bufferLen)
      fillBuffer();

    //fprintf(stdout, "bCopied = %d  bRead = %d\n", bCopied, bRead);

    return(bCopied + bRead);
  }
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

#ifdef DEBUG
    fprintf(stderr, "splitToWords()-- found %u characters.\n", cmdChars);
    fprintf(stderr, "splitToWords()-- found %u potential words.\n", cmdWords);
#endif

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

#ifdef DEBUG
    fprintf(stderr, "splitToWords()-- found %u words.\n", _argWords);
#endif
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
      _buffer = new readBuffer(fileno(stdin), "stdin");
    else
      _buffer = new readBuffer(filename);

    _thePosition       = 0;
    _theSequenceNumber = 0;

    _theDefLineMax     = 2048;
    _theDefLinePos     = 0;
    _theDefLine        = new char [_theDefLineMax + 1];
  };

  ~FastAstream() {
    delete    _buffer;
    delete [] _theDefLine;
  };

  unsigned char     nextSymbol(void);

  u64bit            thePosition(void) {
    return(_thePosition);
  };
  u64bit            theSequenceNumber(void) {
    return(_theSequenceNumber);
  };
  char             *theDefLine(void) {
    return(_theDefLine);
  };
private:
  readBuffer  *_buffer;

  u64bit       _thePosition;
  u64bit       _theSequenceNumber;

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
  nc = (unsigned char)_buffer->getnext();

  //  Skip any whitespace
  //
  while (isspace(nc) && !_buffer->eof())
    nc = (unsigned char)_buffer->getnext();

  //  Are we at the EOF?
  //
  if (_buffer->eof())
    return(0);

  //  Looks like we are in the sequence.  Increment the position.  If
  //  this is the start of a defline, the position is reset to zero.
  //
  _thePosition++;

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
    _thePosition = 0;
    _theSequenceNumber++;

    _theDefLinePos = 0;

    nc = (unsigned char)_buffer->getnext();
    while ((nc != '\n') && (nc != '\r') && !_buffer->eof()) {

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

      nc = (unsigned char)_buffer->getnext();
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
//  merstream
//


//
//  This is the string-based merStreamer.  It gets sequence from
//  a character string.
//

class merStream {
public:
  merStream(u32bit merSize, const char *filename);
  merStream(u32bit merSize, char const *seq, u32bit len);
  ~merStream();

  u64bit       theFMer(void)        { return(_theFMer); };
  u64bit       theRMer(void)        { return(_theRMer); };

  char const  *theFMerString(void);
  char const  *theRMerString(void);

  u64bit       thePosition(void) {
    if (_theFile)
      return(_theFile->thePosition() - _merSize);
    else
      return(_thePos - _merSize);
  };
  u64bit       theSequenceNumber(void) {
    if (_theFile)
      return(_theFile->theSequenceNumber());
    else
      return(u64bitZERO);
  };
  char const  *theDefLine(void) {
    if (_theFile)
      return(_theFile->theDefLine());
    else
      return(">");
  };

  bool               nextMer(u32bit skip=0);
private:
  void               loadMer(u32bit s);
  char               nextSymbol(void) {
    if (_theFile)
      return(_theFile->nextSymbol());
    else
      return(_theSeq[_thePos++]);
  };

  FastAstream          *_theFile;
  char const           *_theSeq;
  u32bit                _theLen;
  u32bit                _thePos;

  u32bit                _merSize;
  s32bit                _timeUntilValid;
  u64bit                _theMerMask;

  u64bit                _theFMer;
  u64bit                _theRMer;

  u32bit                _theRMerShift;

  char                  _theMerString[33];
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
    //  because, in reverse, we shift to the right.  If we don't mask these
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

inline
bool
merStream::nextMer(u32bit skip) {

 again:
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

  if (skip--)
    goto again;

  _theFMer &= _theMerMask;
  _theRMer &= _theMerMask;

  //  Still need to check if we are valid -- we could
  //  have run off the end of the sequence before a valid
  //  mer was created.
  //
  return(_timeUntilValid <= 0);
}





////////////////////////////////////////
//
//  merstreamfile
//


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


class merStreamFileBuilder {
public:
  merStreamFileBuilder(u32bit merSize, const char *inputFile, const char *outputFile);
  ~merStreamFileBuilder();

  //  Returns the number of mers written
  u64bit                 build(bool beVerbose=false);

private:
  u32bit                _merSize;
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

  u64bit                 thePosition(void)       { return(_thePosition); };
  u64bit                 theSequenceNumber(void) { return(_theSequence); };
  char const            *theDefLine(void)        { return(_theDefline);  };

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
  u64bit               *_blockPosition;
  char                 *_deflineStorage;
  char                **_deflines;

  bitPackedFile        *_streamFile;

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
  u64bit                _thePosition;
  u64bit                _theSequence;
  char                 *_theDefline;

  char                  _theMerString[33];
};


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

