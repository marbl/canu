
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' (http://kmer.sourceforge.net)
 *  both originally distributed by Applera Corporation under the GNU General
 *  Public License, version 2.
 *
 *  Canu branched from Celera Assembler at its revision 4587.
 *  Canu branched from the kmer project at its revision 1994.
 *
 *  This file is derived from:
 *
 *    src/utility/AS_UTL_decodeRange.H
 *    src/utility/findKeyAndValue.H
 *    src/utility/splitToWords.H
 *
 *  Modifications by:
 *
 *    Brian P. Walenz beginning on 2018-JUL-20
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#ifndef STRINGS_H
#define STRINGS_H

#include "AS_global.H"

#include <set>
#include <vector>

using namespace std;


//  perl's chomp is pretty nice
//  Not a great place to put this, but it's getting used all over.
#ifndef chomp
#define chomp(S)  { char *t=(S); while (*t) t++; t--; while (t >= S && isspace(*t)) *t--=0; }
#endif

#ifndef munch
#define munch(S)  { while (*(S) &&  isspace(*(S))) (S)++; }
#endif

#ifndef crunch
#define crunch(S) { while (*(S) && !isspace(*(S))) (S)++; }
#endif



//  For pretty-printing numbers.  Converts integers to, e.g., 123 k, 123 M, G, T, P.
uint64      scaledNumber(uint64 n, uint32 div=1024);
char        scaledUnit  (uint64 n, uint32 div=1024);
const char *scaledName  (uint64 n, uint32 div=1024);





template<typename T>
char *
decodeRange(char *range, T &lo, T &hi) {
  char    *ap = range;

  strtonumber(lo, ap, &ap);     //  Grab the first number.

  hi = lo;                      //  Set the second to that.

  if ((*ap == '-') ||           //  If this is a range,
      (*ap == '/')) {           //  or a one-of-many selection,
    ap++;                       //  grab the second number
    strtonumber(hi, ap, &ap);
  }

  if (*ap == ',')               //  If the next letter continues
    return(ap + 1);             //  move past that and return.

  if (*ap == 0)                 //  If the next letter is the end
    return(NULL);               //  of the string, return NULL.

  //  Otherwise, we can't decode this range.

  fprintf(stderr, "ERROR: invalid range '%s'\n", range);
  exit(1);

  return(NULL);
}



template<typename T>
void
decodeRange(char *range, vector<T> &bgn, vector<T> &end) {
  char    *ap = range;
  T        av = 0;
  T        bv = 0;

  while ((ap != NULL) && (*ap != 0)) {
    ap = decodeRange(ap, av, bv);

    bgn.push_back(av);
    end.push_back(bv);
  }
}



template<typename T>
void
decodeRange(char *range, set<T> &ranges) {
  char    *ap = range;
  T        av = 0;
  T        bv = 0;

  while ((ap != NULL) && (*ap != 0)) {
    ap = decodeRange(ap, av, bv);

    for (T xx=av; xx<=bv; xx++)
      ranges.insert(xx);
  }
}



//  Decodes a string with 0/1, false/true, no/yes into an integer flag.
bool   decodeBoolean(char *value);



class KeyAndValue {
public:
  KeyAndValue(char *line = NULL)   {  find(line);  };
  ~KeyAndValue()                   {               };

  bool     find(char *line);

  char    *key(void)           { return(key_); };
  char    *value(void)         { return(val_); };

  bool     value_bool(void)    { return(decodeBoolean(val_)); };

  int32    value_int32(void)   { return(strtol (val_, NULL, 10)); };
  int64    value_int64(void)   { return(strtoll(val_, NULL, 10)); };

  uint32   value_uint32(void)  { return(strtoul (val_, NULL, 10)); };
  uint64   value_uint64(void)  { return(strtoull(val_, NULL, 10)); };

  float    value_float(void)   { return(strtof(val_, NULL)); };
  double   value_double(void)  { return(strtod(val_, NULL)); };

public:
  bool    iscomment(char c)    { return((c == '!') || (c == '#') || (c == 0)); };
  bool    isdelimiter(char c)  { return((c == ':') || (c == '=') || isspace(c)); };

  char   *key_;
  char   *val_;
};



enum splitType {
  splitWords  = 0,
  splitPaths  = 1
};

class splitToWords {
public:
  splitToWords(const char *string=NULL, splitType type=splitWords);
  ~splitToWords();

private:
  bool   isPath(char c) {
    return(c == '/');
  };

  bool   isSpace(char c) {
    return((c == ' ')  ||
           (c == '\t') ||
           (c == '\n') ||
           (c == '\r'));
  };

  bool   isSeparator(char c, splitType type) {
    return(((type == splitWords) && (isSpace(c))) ||
           ((type == splitPaths) && (isPath(c))));
  };

public:
  void    split(const char *line, splitType type=splitWords);

  uint32  numWords(void)        { return(_wordsLen); };

  char   *operator[](uint32 i)  { return(first(i)); };

  char   *first(uint32 i=0)     { return((_wordsLen <= i) ? NULL : _words[i]);  };
  char   *shift(void)           {
    if (_wordsLen == 0)                     //  If no words, nothing to return.
      return(NULL);

    for (uint32 ii=1; ii<_wordsLen; ii++)   //  Shift all words down one place, moving
      swap(_words[ii-1], _words[ii]);       //  the word to shift off to the end.

    return(_words[--_wordsLen]);            //  Return the word we shifted out.
  };

  char   *last(uint32 i=0)      { return((_wordsLen == 0) ? NULL : _words[_wordsLen - i - 1]); };
  char   *pop(void)             { return((_wordsLen == 0) ? NULL : _words[--_wordsLen]);       };

  int32   toint32(uint32 i)     { return(strtoint32 (_words[i])); };
  uint32  touint32(uint32 i)    { return(strtouint32(_words[i])); };
  int64   toint64(uint32 i)     { return(strtoint64 (_words[i])); };
  uint64  touint64(uint32 i)    { return(strtouint64(_words[i])); };
  double  todouble(uint32 i)    { return(strtodouble(_words[i])); };

private:
  uint32    _wordsLen;
  uint32    _wordsMax;
  char    **_words;

  uint32    _charsLen;
  uint32    _charsMax;
  char     *_chars;
};







#endif  //  STRINGS_H
