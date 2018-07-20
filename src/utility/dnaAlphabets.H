
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
 *  Modifications by:
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_global.H"

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
