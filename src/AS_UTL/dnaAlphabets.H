#include <stdio.h>
#include <ctype.h>



class dnaAlphabets {
public:
  dnaAlphabets() {
    initTablesForACGTSpace();
  };
  ~dnaAlphabets() {
  };

  void             initTablesForACGTSpace(void);
  void             initTablesForColorSpace(void);

  bool             isWhitespace(unsigned char x)      { return(_whitespaceSymbol[x]); };
  unsigned char    toLower(unsigned char x)           { return(_toLower[x]);          };
  unsigned char    toUpper(unsigned char x)           { return(_toUpper[x]);          };

  unsigned char    letterToBits(unsigned char x)      { return(_letterToBits[x]);     };
  unsigned char    bitsToLetter(unsigned char x)      { return(_bitsToLetter[x]);     };

  unsigned char    complementSymbol(unsigned char x)  { return(_complementSymbol[x]); };

  bool             validCompressedSymbol(unsigned char x) {
    return(_validCompressedSymbol[x]);
  };

private:
  unsigned char   _whitespaceSymbol[256];
  unsigned char   _toLower[256];
  unsigned char   _toUpper[256];

  unsigned char   _letterToBits[256];
  unsigned char   _bitsToLetter[256];
  unsigned char   _bitsToColor[256];

  unsigned char   _complementSymbol[256];
  unsigned char   _validCompressedSymbol[256];

  unsigned char   _IUPACidentity[128][128];
  unsigned char   _baseToColor[128][128];
};


extern dnaAlphabets  alphabet;
