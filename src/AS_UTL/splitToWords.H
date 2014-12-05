#ifndef SPLITTOWORDS_H
#define SPLITTOWORDS_H


class splitToWords {
public:
  splitToWords() {
    _argWords = 0;
    _maxWords = 0;
    _arg      = 0L;
    _maxChars = 0;
    _cmd      = 0L;
  };
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
  int64  operator()(uint32 i)  { return(strtoull(_arg[i], NULL, 10)); };
private:
  uint32    _argWords;
  uint32    _maxWords;
  char    **_arg;
  uint32    _maxChars;
  char     *_cmd;
};


#endif  //  SPLITTOWORDS_H
