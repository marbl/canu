#ifndef UTIL_PLUS_PLUS_H
#define UTIL_PLUS_PLUS_H

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <ctype.h>
#include <errno.h>
#include <string.h>

//
//  This came from Brian Walenz's ESTmapper/ATAC/meryl/snapper utility
//  library.
//


//  Useful types.
//
//  *MASK(x) is only defined for unsigned types, with x != 0 and less
//  than the datawidth.

#if defined(__alpha) || defined(_AIX) || defined(__LP64__) || defined(_LP64)
#define TRUE64BIT
#endif

#ifdef TRUE64BIT

//  Generic 64-bit platform

typedef unsigned long       u64bit;
typedef unsigned int        u32bit;
typedef signed long         s64bit;
typedef signed int          s32bit;
typedef unsigned short      u16bit;
typedef signed short        s16bit;
typedef unsigned char       u8bit;
typedef signed char         s8bit;

#define  u64bitNUMBER(X) X ## LU
#define  u64bitZERO      (0x0000000000000000LU)
#define  u64bitONE       (0x0000000000000001LU)
#define  u64bitMAX       (0xffffffffffffffffLU)
#define  u64bitMASK(X)   ((~u64bitZERO) >> (64 - (X)))
#define  u64bitFMTW(X)   "%" #X "lu"
#define  u64bitFMT       "%lu"
#define  u64bitHEX       "0x%016lx"
#define  s64bitFMTW(X)   "%" #X "ld"
#define  s64bitFMT       "%ld"

#define  u32bitNUMBER(X) X ## U
#define  u32bitZERO      (0x00000000U)
#define  u32bitONE       (0x00000001U)
#define  u32bitMAX       (0xffffffffU)
#define  u32bitMASK(X)   ((~u32bitZERO) >> (32 - (X)))
#define  u32bitFMTW(X)   "%" #X "u"
#define  u32bitFMT       "%u"
#define  u32bitHEX       "0x%08x"
#define  s32bitFMTW(X)   "%" #X "d"
#define  s32bitFMT       "%d"

#define  u16bitZERO      (0x0000)
#define  u16bitONE       (0x0001)
#define  u16bitMAX       (0xffff)
#define  u16bitMASK(X)   ((~u16bitZERO) >> (16 - (X)))
#define  u16bitFMT       "%hd"

#define  u8bitZERO       (0x00)
#define  u8bitONE        (0x01)
#define  u8bitMAX        (0xff)
#define  u8bitMASK(X)    ((~u8bitZERO) >> (8 - (X)))

#define  strtou32bit(N,O) (u32bit)strtoul(N, O, 0)
#define  strtou64bit(N,O) (u64bit)strtoul(N, O, 0)

#else

//  Generic 32-bit platform

typedef unsigned long long  u64bit;
typedef unsigned long       u32bit;
typedef signed long long    s64bit;
typedef signed long         s32bit;
typedef unsigned short      u16bit;
typedef signed short        s16bit;
typedef unsigned char       u8bit;
typedef signed char         s8bit;

#define  u64bitNUMBER(X) X ## LLU
#define  u64bitZERO      (0x0000000000000000LLU)
#define  u64bitONE       (0x0000000000000001LLU)
#define  u64bitMAX       (0xffffffffffffffffLLU)
#define  u64bitMASK(X)   ((~u64bitZERO) >> (64 - (X)))
#define  u64bitFMTW(X)   "%" #X "llu"
#define  u64bitFMT       "%llu"
#define  u64bitHEX       "0x%016llx"
#define  s64bitFMTW(X)   "%" #X "lld"
#define  s64bitFMT       "%lld"

#define  u32bitNUMBER(X) X ## LU
#define  u32bitZERO      (0x00000000LU)
#define  u32bitONE       (0x00000001LU)
#define  u32bitMAX       (0xffffffffLU)
#define  u32bitMASK(X)   ((~u32bitZERO) >> (32 - (X)))
#define  u32bitFMTW(X)   "%" #X "lu"
#define  u32bitFMT       "%lu"
#define  u32bitHEX       "0x%08lx"
#define  s32bitFMTW(X)   "%" #X "ld"
#define  s32bitFMT       "%ld"

#define  u16bitZERO      (0x0000)
#define  u16bitONE       (0x0001)
#define  u16bitMAX       (0xffff)
#define  u16bitMASK(X)   ((~u16bitZERO) >> (16 - (X)))
#define  u16bitFMT       "%hd"

#define  u8bitZERO       (0x00)
#define  u8bitONE        (0x01)
#define  u8bitMAX        (0xff)
#define  u8bitMASK(X)    ((~u8bitZERO) >> (8 - (X)))

#define  strtou32bit(N,O) (u32bit)strtoul(N, O, 0)
#define  strtou64bit(N,O) (u64bit)strtoull(N, O, 0)

#endif  //  TRUE64BIT




double
getTime(void);


inline
off_t
sizeOfFile(const char *path) {
  struct stat s;

  if (stat(path, &s) != 0) {
    fprintf(stderr, "Couldn't stat() '%s'\n%s\n", path, strerror(errno));
    exit(1);
  }

  return(s.st_size);
}







////////////////////////////////////////
//
//  perl's chomp is pretty nice
//
#define chomp(S) { char *t=S; while (*t) t++; t--; while (isspace(*t)) *t--=0; }






#ifndef SPLITTOWORDS_H
#define SPLITTOWORDS_H


class splitToWords {
public:
  splitToWords(char *cmd) {
    _argWords = 0;
    _maxWords = 0;
    _arg      = 0L;
    _maxChars = 0;
    _cmd      = 0L;

    split(cmd);
  };
  ~splitToWords() {
    delete [] _cmd;
    delete [] _arg;
  };


  void   split(char *cmd) {

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

    if (cmdChars > _maxChars) {
      delete [] _cmd;
      _cmd      = new char   [cmdChars];
      _maxChars = cmdChars;
    }
    if (cmdWords > _maxWords) {
      delete [] _arg;
      _arg      = new char * [cmdWords];
      _maxWords = cmdWords;
    }

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


  u32bit  numWords(void)        { return(_argWords); };
  char   *getWord(u32bit i)     { return(_arg[i]); };
  char   *operator[](u32bit i)  { return(_arg[i]); };
private:
  u32bit    _argWords;
  u32bit    _maxWords;
  char    **_arg;
  u32bit    _maxChars;
  char     *_cmd;
};


#endif  //  SPLITTOWORDS_H







#ifndef INTERVALLIST_H
#define INTERVALLIST_H


//  I don't want this to be public, but the sort function
//  in intervalList.C needs it so.
//
typedef u64bit  intervalNumber;

struct _intervalPair {
  intervalNumber    lo;
  intervalNumber    hi;
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
  void        add(intervalNumber position, intervalNumber length);

  //  Sort the set of intervals by the lo value
  void        sort(void);

  //  Merge overlapping or adjacent intervals together.
  void        merge(void);

  //  Returns the number of intervals
  u32bit      numberOfIntervals(void) {
    return(_listLen);
  };

  //  Returns the sum of the length of all intervals
  intervalNumber      sumOfLengths(void) {
    intervalNumber len = 0;
    u32bit         i   = numberOfIntervals();

    if (i > 0)
      while (i--)
        len += _list[i].hi - _list[i].lo;

    return(len);
  };

  //  Populates an array with the intervals that are within the
  //  supplied interval.  Return
  //
  u32bit      overlapping(intervalNumber    lo,
                          intervalNumber    hi,
                          u32bit          *&intervals,
                          u32bit           &intervalsLen,
                          u32bit           &intervalsMax);

  intervalNumber      lo(u32bit i) { return(_list[i].lo); };
  intervalNumber      hi(u32bit i) { return(_list[i].hi); };

private:
  bool                      _isSorted;
  u32bit                    _listLen;
  u32bit                    _listMax;
  _intervalPair            *_list;
};


#endif  //  INTERVALLIST_H






#ifndef SPEEDCOUNTER_H
#define SPEEDCOUNTER_H

#include <stdio.h>

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
  static const char  *_spinr[4];
  static const char  *_liner[19];
  u64bit              _count;
  u64bit              _draws;
  double              _unit;
  u64bit              _freq;
  double              _startTime;
  char const         *_fmt;
  bool                _spin;
  bool                _line;
  bool                _enabled;
};


#endif  //  SPEEDCOUNTER_H









#endif  //  UTIL_PLUS_PLUS_H

