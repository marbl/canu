
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
 *    src/utility/splitToWords.H
 *
 *  Modifications by:
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#ifndef SPLITTOWORDS_H
#define SPLITTOWORDS_H

#include "AS_global.H"



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




enum splitType {
  splitWords  = 0,
  splitPaths  = 1
};

class splitToWords {
public:
  splitToWords(const char *string=NULL, splitType type=splitWords) {
    _wordsLen  = 0;
    _wordsMax  = 0;
    _words     = NULL;

    _charsLen = 0;
    _charsMax = 0;
    _chars    = NULL;

    if (string)
      split(string, type);
  };

  ~splitToWords() {
    delete [] _chars;
    delete [] _words;
  };

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
  void   split(const char *line, splitType type=splitWords) {

    _wordsLen = 0;        //  Initialize to no words
    _charsLen = 0;        //  and no characters.

    if (line == NULL)     //  Bail if there isn't a line to process.
      return;

    //  Count the number of words and chars in the input line, then make
    //  sure there is space for us to store them.

    while (line[_charsLen] != 0)
      if (isSeparator(line[_charsLen++], type))
        _wordsLen++;

    resizeArray(_words, 0, _wordsMax, _wordsLen + 1, resizeArray_doNothing);
    resizeArray(_chars, 0, _charsMax, _charsLen + 1, resizeArray_doNothing);

    //  Clear all the words pointers, and copy the input line to our storage.
    //  This greatly simplifies the loop, as we don't need to worry about
    //  terminating the final word.

    memset(_words, 0,    sizeof(char *) * (_wordsLen + 1));
    memcpy(_chars, line, sizeof(char)   * (_charsLen + 1));

    //  Scan the line copy, converting word separators to NUL bytes.
    //  counting and saving the start of each word in _words.

    _wordsLen = 0;

    for (uint32 st=1, ii=0; ii < _charsLen; ii++) {
      if (isSeparator(line[ii], type)) {      //  If the character is a word
        _chars[ii] = 0;                       //  separator, convert to NUL,
        st         = true;                    //  and flag the next character
      }                                       //  as the start of a new word.

      else if (st) {                          //  Otherwise, if this is the
        _words[_wordsLen++] = _chars + ii;    //  start of a word, make
        st                  = false;          //  a new word.
      }
    }
  };

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

#endif  //  SPLITTOWORDS_H
