#ifndef UTIL_PLUS_PLUS_H
#define UTIL_PLUS_PLUS_H



#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>
#include <fcntl.h>
#include <string.h>
#include <ctype.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "util.h"

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

////////////////////////////////////////
//
//  time
//

class speedCounter {
public:
  //  fmt specifies the status format.  An example:
  //    "    %8f [unit]things (%8.5f [unit]things/sec)\r"
  //
  speedCounter(char const   *fmt,
               double        unit,
               u64bit        freq,
               bool          enabled=true);
  ~speedCounter();

  void   enableSpinner(void) { _spin = true; };
  void   enableLiner(void)   { _line = true; };

  void   tick(void) {
    if (_enabled && ((++_count & _freq) == u64bitZERO)) {
      double  v = _count / _unit;
      if (_spin) fputs(_spinr[_draws %  4], stderr);
      if (_line) fputs(_liner[_draws % 19], stderr);
      _draws++;
      fprintf(stderr, _fmt, v, v / (getTime() - _startTime));
      fflush(stderr);
    }
  };

  void   finish(void) {
    if (_enabled && (_count >= _freq)) {
      double  v = _count / _unit;
      if (_spin) fputs(_spinr[_draws %  4], stderr);
      if (_line) fputs(_liner[_draws % 19], stderr);
      fprintf(stderr, _fmt, v, v / (getTime() - _startTime));
      fprintf(stderr, "\n");
      fflush(stderr);
    }
    _count = 0;
  };

private:
  static char    *_spinr[4];
  static char    *_liner[19];
  u64bit          _count;
  u64bit          _draws;
  double          _unit;
  u64bit          _freq;
  double          _startTime;
  char const     *_fmt;
  bool            _spin;
  bool            _line;
  bool            _enabled;
};


////////////////////////////////////////
//
//  readBuffer
//
#include "readBuffer.H"



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
//  bitPackedArray
//
//  implements an integer array using bit-widths less than word-sizes,
//  e.g., a memory efficient way to store 23 bit numbers.  Numbers may
//  be up to 64 bits wide.
//
//  The array is variable length, and it is implemented as an array,
//  not a list or tree -- accessing element 1,000,000 will allocate
//  elements 0 through 999,999.
//
class bitPackedArray {
public:

  //  Create a bitpacked array with elements of width 'width' using
  //  'segmentSize' KB per segment.  If you know your array is going
  //  to be much bigger or smaller, crank this value.
  //
  bitPackedArray(u32bit valueWidth, u32bit segmentSize = 1024);
  ~bitPackedArray();

  //  No array operator is provided, because we cannot return a
  //  reference to a value that is split across two words (or even a
  //  reference to a value that is not bit aligned in the word).
  //
  u64bit   get(u64bit idx);
  void     set(u64bit idx, u64bit val);

  //  Clear the array.  Since the array is variable sized, you must add
  //  things to a new array before clearing it.
  void     clear(void);

private:
  u32bit   _valueWidth;
  u32bit   _segmentSize;
  u32bit   _nextElement;  //  the first invalid element
  u32bit   _valuesPerSegment;

  u32bit   _numSegments;
  u32bit   _maxSegments;
  u64bit **_segments;
};




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




#endif  //  UTIL_PLUS_PLUS_H

