#ifndef UTIL_PLUS_PLUS_H
#define UTIL_PLUS_PLUS_H

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <ctype.h>
#include <errno.h>
#include <string.h>

extern "C" {
#include "AS_global.h"
}

//  perl's chomp is pretty nice
//
#define chomp(S) { char *t=S; while (*t) t++; t--; while (isspace(*t)) *t--=0; }


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
    uint32   cmdChars = 1;  //  1 == Space for terminating 0
    uint32   cmdWords = 2;  //  2 == Space for first word and terminating 0L

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


  uint32  numWords(void)        { return(_argWords); };
  char   *getWord(uint32 i)     { return(_arg[i]); };
  char   *operator[](uint32 i)  { return(_arg[i]); };
private:
  uint32    _argWords;
  uint32    _maxWords;
  char    **_arg;
  uint32    _maxChars;
  char     *_cmd;
};






//  I don't want this to be public, but the sort function
//  in intervalList.C needs it so.
//
typedef uint64  intervalNumber;

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
  uint32      numberOfIntervals(void) {
    return(_listLen);
  };

  //  Returns the sum of the length of all intervals
  intervalNumber      sumOfLengths(void) {
    intervalNumber len = 0;
    uint32         i   = numberOfIntervals();

    if (i > 0)
      while (i--)
        len += _list[i].hi - _list[i].lo;

    return(len);
  };

  //  Populates an array with the intervals that are within the
  //  supplied interval.  Return
  //
  uint32      overlapping(intervalNumber    lo,
                          intervalNumber    hi,
                          uint32          *&intervals,
                          uint32           &intervalsLen,
                          uint32           &intervalsMax);

  intervalNumber      lo(uint32 i) { return(_list[i].lo); };
  intervalNumber      hi(uint32 i) { return(_list[i].hi); };

private:
  bool                      _isSorted;
  uint32                    _listLen;
  uint32                    _listMax;
  _intervalPair            *_list;
};




#endif  //  UTIL_PLUS_PLUS_H

