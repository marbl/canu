#ifndef BIO_PLUS_PLUS_H
#define BIO_PLUS_PLUS_H



#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>
#include <fcntl.h>
#include <string.h>
#include <ctype.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "bio.h"
#include "util++.H"

#include "mers.h"
#include "fasta.H"
#include "fasta-cache.H"

////////////////////////////////////////
//
//  chainedSequence
//
//  Build a chained sequence from a source file or source
//  FastAWrapper.  In the case of a file, a FastAWrapper (with iid
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
  void         setSource(char *filename);
  void         setSource(FastAWrapper *file);

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
  bool         eof(void);

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


  //  Move the file pointer back to the start of the chained sequence
  //
  bool         rewind(void);

private:
  struct use_s {
    u32bit   iid;
    u32bit   length;
    u32bit   start;
  };

  char                 *_filename;
  FastAWrapper         *_file;
  FastAWrapper         *_fileToDelete;
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
  u32bit                _separator;
  bool                  _separatorDone;
  u32bit                _separatorLength;
  u32bit                _separatorPosition;

#ifdef __APPLE__
  u32bit                _pad;
#endif
};




inline
bool
chainedSequence::eof(void) {

  //  It's easy to tell if the _next_ get() is going to return
  //  eof, but that's not what we need to return.
  //return((_sequence->get() == 0) && (_currentSeq + 1 >= _useListLen));
  //
  return(_eof);
}



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





////////////////////////////////////////
//
//  fastastream
//

#include "fastastream.H"



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

//
//  N.B. because we create a merStream internally, and that uses
//  the FastAstream to get sequence from the file, the sequence
//  file MUST BE SQUEEZED to get the correct positionInStream.
//

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
  FastAstream          *_fastaStream;
  merStream            *_merStream;
  const char           *_outputFile;
};


class merStreamFileReader {
public:
  merStreamFileReader(const char *filename);
  ~merStreamFileReader();

  u32bit                 merSize(void)      { return(_merSize); };
  u64bit                 numberOfMers(void) { return(_numMers); };

  u64bit                 theFMer(void) { return(_fMer); };
  u64bit                 theRMer(void) { return(_rMer); };

  char const            *theFMerString(void);
  char const            *theRMerString(void);

  u64bit                 thePositionInSequence(void) { return(_posInSeq); };
  u64bit                 thePositionInStream(void)   { return(_posInStr); };
  u64bit                 theSequenceNumber(void)     { return(_sequence); };
  char const            *theDefLine(void)            { return(_defline);  };

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
  u64bit                _fMer;
  u64bit                _rMer;
  u64bit                _posInSeq;
  u64bit                _posInStr;
  u64bit                _sequence;

  char                 *_defline;

  char                  _merString[36];
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
//  **********
//
//  N.B. white-space in the source (applies to FastAstream and char*
//  only) result in an incorrect positionInStream().  If this troubles
//  you, use a chainedSequence instead of a FastAstream, or squeeze
//  the sequences.
//

class merStream {
public:
  merStream(merStreamFileReader *merstreamfile);
  merStream(u32bit merSize, FastAstream *str);
  merStream(u32bit merSize, chainedSequence *cs);
  merStream(u32bit merSize, const char *seq, u32bit len);
  ~merStream();


  //  The merstream has many modes of operation, depending
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


  u64bit       theFMer(void)        { return(_fMer); };
  u64bit       theRMer(void)        { return(_rMer); };

  char const  *theFMerString(void);
  char const  *theRMerString(void);

  u32bit       theMerSize(void) { return(_merSize); };

  bool         rewind(void);

  //  This returns the position relative to the start of the sequence.
  //  The sequences are squeezed.
  //
  //  Why the +1?  FastAstream->thePositionInSequence() returns the
  //  position that the last character retrieved is at.  We need to go
  //  back _merSize-1 positions to get to the start of the mer.
  //
  //  ?
  //
  u64bit       thePositionInSequence(void) {
    if      (_ms_mers)
      return(_ms_mers->thePositionInSequence());
    else if (_fs_file)
      return(_fs_file->thePositionInSequence() - _merSize + 1);
    else if (_cs_chainedSeq)
      return(_cs_chainedSeq->thePositionInSequence() - _merSize + 1);
    else
      return(_st_posInSeq - _merSize);
  };

  //  This returns the position relative to the start of the source.
  //  It is a _byte_ position.  No squeezing.  If you give it a fasta
  //  file, it includes any deflines and white space.
  //
  u64bit       thePositionInStream(void) {
    if      (_ms_mers)
      return(_ms_mers->thePositionInStream());
    else if (_fs_file)
      return(_fs_file->thePositionInStream() - _merSize + 1);
    else if (_cs_chainedSeq)
      return(_cs_chainedSeq->thePositionInStream() - _merSize + 1);
    else
      return(_st_stringPos - _merSize);  //  was _st_posInStr
  };

  u64bit       theSequenceNumber(void) {
    if      (_ms_mers)
      return(_ms_mers->theSequenceNumber());
    else  if (_fs_file)
      return(_fs_file->theSequenceNumber());
    else if (_cs_chainedSeq)
      return(_cs_chainedSeq->theSequenceNumber());
    else
      return(_st_num);
  };

  char const  *theDefLine(void) {
    if      (_ms_mers)
      return(_ms_mers->theDefLine());
    else if (_fs_file)
      return(_fs_file->theDefLine());
    else
      return(">");
  };

  bool               nextMer(u32bit skip=0);
private:
  void               loadMer(u32bit s);
  u64bit             nextSymbol(void);


  //  Data specific to the merStreamFileReader input
  //
  merStreamFileReader  *_ms_mers;


  //  Data specific to the FastA file input
  //
  FastAstream          *_fs_file;


  //  Data specific to the chainedSequence input
  //
  chainedSequence      *_cs_chainedSeq;


  //  Data specific to the character string input
  //
  char const           *_st_string;
  u64bit                _st_stringPos;
  u64bit                _st_stringLen;
  u64bit                _st_posInSeq;
  u64bit                _st_posInStr;
  u32bit                _st_num;


  //  Data specific to nobody or everybody
  //
  u32bit                _merSize;
  s64bit                _timeUntilValid;
  u32bit                _rMerShift;

  u64bit                _merMask;

  u64bit                _fMer;
  u64bit                _rMer;

  //  Used for returning the mer as a character string.
  //
  char                  _merString[64];
};



inline
u64bit
merStream::nextSymbol(void) {

  //  This is called only for FastAstreams, chainedSequences and char* input.

  if (_fs_file)
    return(_fs_file->nextSymbol());

  if (_cs_chainedSeq)
    return(_cs_chainedSeq->get());

  //  See also the constructor for character strings -- it skips
  //  white space and the first defline.
  
  //  Skip any whitespace
  //
  while (isspace(_st_string[_st_stringPos]) && (_st_string[_st_stringPos]) && (_st_stringPos < _st_stringLen))
    _st_stringPos++;

  //  Are we at a defline?  If so, zip along until we hit a
  //  newline.  Then we return the code (fastastream) for a
  //  sequence break.
  //
  if (_st_string[_st_stringPos] == '>') {
    while ((_st_string[_st_stringPos] != '\n') && (_st_string[_st_stringPos]) && (_st_stringPos < _st_stringLen))
      _st_stringPos++;

    _st_posInSeq = 0;
    _st_posInStr++;

    _st_num++;
    return(254);
  }

  //  Are we at the EOF?
  //
  if ((_st_string[_st_stringPos] == 0) || (_st_stringPos >= _st_stringLen))
    return(0);

  //  Update the position of the character in this sequence, and
  //  the positon in the stream.
  //
  _st_posInSeq++;
  _st_posInStr++;

  return(_st_string[_st_stringPos++]);
}






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
    _fMer <<= 2;
    _rMer >>= 2;

    _fMer |= cf;
    _rMer |= cr << _rMerShift;
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

  if (_ms_mers) {
    bool ret = _ms_mers->nextMer(skip);

    _fMer = _ms_mers->theFMer();
    _rMer = _ms_mers->theRMer();

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
      _fMer <<= 2;
      _rMer >>= 2;

      _fMer |= cf;
      _rMer |= cr << _rMerShift;

      //  If the ch is invalid, we need to make a whole new mer.
      //  This also implies that we're done skipping.
      //
      if (cf & (unsigned char)0xfc) {
        loadMer(_merSize);
        skip = 0;
      }
    } while (skip--);

    _fMer &= _merMask;
    _rMer &= _merMask;

    //  Still need to check if we are valid -- we could
    //  have run off the end of the sequence before a valid
    //  mer was created.
    //
    return(_timeUntilValid <= 0);
  }
}




#include "merCovering.H"
#include "merList.H"


#endif  //  BIO_PLUS_PLUS_H

