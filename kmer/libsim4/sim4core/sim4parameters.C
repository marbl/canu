#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "fasta.H"
#include "sim4parameters.H"

//
//  if -d dbFile, read in that file here.  Sim4 ALWAYS uses
//  (dbIdx, dbLo, dbHi) interface.
//

//  Converts a character string into a sim4parameter structure
//
sim4parameters::sim4parameters(char *cmd, SequenceManager *ESTs, SequenceManager *DBs) {

  //
  //  Clear the stuff
  //
  _numESTs = 0;
  _maxESTs = 16;
  _ESTlist = new unsigned int [_maxESTs];

  _ESTs         = ESTs;
  _DBs          = DBs;
  _DBsIsPrivate = false;

  _dbIdx  = 0;
  _dbLo   = 0;
  _dbHi   = 0;

  _doForward = true;
  _doReverse = true;

  _strandIndicator = 0L;


  //  Step Zero:
  //
  //  Count the length of the string, in words and in characters.
  //  For simplicity, we overcount words, by just counting white-space.
  //
  unsigned int   cmdChars = 1;  //  1 == Space for terminating 0
  unsigned int   cmdWords = 2;  //  2 == Space for first word and terminating 0L

  for (_cmd=cmd; *_cmd; _cmd++) {

    //  This is *significantly* faster than using isspace(), and
    //  somewhat faster than the usual if-based expression.
    //
    cmdWords += *_cmd == ' ';
    cmdWords += *_cmd == '\t';

    cmdChars++;
  }


  //  Now we can allocate space for a temporary copy of the string,
  //  and a set of pointers into the temporary copy that will act like
  //  argv.
  //
  _cmd = new char   [cmdChars];
  _arg = new char * [cmdWords];


  //  Step One:
  //
  //  Determine where the words are in the command string, copying the
  //  string to _cmd and storing words in _arg.
  //
  bool           isFirst  = true;
  unsigned int   argWords = 0;
  char          *cmdI     = cmd;
  char          *cmdO     = _cmd;

  while (*cmdI) {
    *cmdO = *cmdI;

    if ((*cmdI != ' ') && (*cmdI != '\t')) {
      //  We are at a non-space character, a word.  If this is the first
      //  character in the word, save the word in the args list.
      //
      if (isFirst) {
        _arg[argWords++] = cmdO;
        isFirst          = false;
      }
    } else {
      //  We are at a space and thus not in a word.  Make all spaces
      //  be string terminators, and declare that we are at the start
      //  of a word.
      //
      *cmdO   = 0;
      isFirst = true;
    }

    cmdI++;
    cmdO++;
  }

  //  Terminate the last arg
  //
  *cmdO   = 0;

  //  Terminate the args list
  //
  _arg[argWords] = 0L;



  //  Step Two:
  //
  //  Parse the arguments
  //
  //  Valid arguments:
  //    -f  Forward only
  //    -r  Reverse only
  //    -s  strand indicator string
  //    -d  dbFile
  //    -D  dbSeq dbLo dbHi
  //    -e  est-list
  //
  argWords = 0;

  while (_arg[argWords]) {
    switch (_arg[argWords][1]) {
    case 'f':
      _doForward = true;
      _doReverse = false;
      break;
    case 'r':
      _doForward = false;
      _doReverse = true;
      break;
    case 'S':
      argWords++;
      _strandIndicator = new char [strlen(_arg[argWords]) + 1];
      strcpy(_strandIndicator, _arg[argWords]);
      break;
    case 'd':
      if (_DBs != 0L) {
        fprintf(stderr, "ERROR:  Was given _DBs, but told to make a new one.\n");
        exit(1);
      }

      argWords++;
      _DBs          = new SequenceManager(_arg[argWords], true);
      _DBsIsPrivate = true;
      break;
    case 'D':
      argWords++;  _dbIdx = atoi(_arg[argWords]);
      argWords++;  _dbLo  = atoi(_arg[argWords]);
      argWords++;  _dbHi  = atoi(_arg[argWords]);
      break;
    case 'e':
      argWords++;

      while ((_arg[argWords]) && (_arg[argWords][0] != '-')) {
        if (_numESTs >= _maxESTs) {

          _maxESTs *= 2;
          unsigned int *newEST = new unsigned int [_maxESTs];

          for (unsigned int i=0; i<_numESTs; i++)
            newEST[i] = _ESTlist[i];

          delete [] _ESTlist;
          _ESTlist = newEST;
        }

        _ESTlist[_numESTs] = atoi(_arg[argWords]);

        _numESTs++;
        argWords++;
      }
      argWords--;
      break;
    default:
      //fprintf(stderr, "Unknown option '%s'\n", _arg[argWords]);
      break;
    }
    argWords++;
  }


  //  Make sure that the ranges of the database sequence are correct.
  //
  if (_dbHi <= _dbLo) {
    _dbLo = 0;
    _dbHi = getDBlength(_dbIdx);
  }
}


sim4parameters::sim4parameters(unsigned int ESTid,
                               unsigned int DBid,
                               SequenceManager *ESTs, SequenceManager *DBs) {
  _cmd = 0L;
  _arg = 0L;

  _numESTs = 1;
  _maxESTs = 1;
  _ESTlist = new unsigned int [1];
  _ESTlist[0] = ESTid;

  _ESTs         = ESTs;
  _DBs          = DBs;
  _DBsIsPrivate = false;

  _dbIdx  = DBid;
  _dbLo   = 0;
  _dbHi   = getDBlength(_dbIdx);

  _doForward = true;
  _doReverse = true;

  _strandIndicator = 0L;
}


sim4parameters::~sim4parameters() {

  delete [] _cmd;
  delete [] _arg;
  delete [] _ESTlist;

  if (_strandIndicator)
    delete [] _strandIndicator;

  if (_DBsIsPrivate)
    delete [] _DBs;
}
