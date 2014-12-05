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

  bool             iswhitespace(unsigned char x)    { return(whitespaceSymbol[x]); };
  unsigned char    tolower(unsigned char x)         { return(toLower[x]);          };
  unsigned char    toupper(unsigned char x)         { return(toUpper[x]);          };

  unsigned char    tobits(unsigned char x)          { return(letterToBits[x]);     };
  unsigned char    toletter(unsigned char x)        { return(bitsToLetter[x]);     };

  unsigned char    complement(unsigned char x)      { return(complementSymbol[x]); };

  bool             validcompressedsymbol(unsigned char x) { return(validCompressedSymbol[x]); };

private:
  unsigned char   whitespaceSymbol[256];
  unsigned char   toLower[256];
  unsigned char   toUpper[256];

  unsigned char   letterToBits[256];
  unsigned char   bitsToLetter[256];
  unsigned char   bitsToColor[256];

  unsigned char   complementSymbol[256];
  unsigned char   validCompressedSymbol[256];

  unsigned char   IUPACidentity[128][128];
  unsigned char   baseToColor[128][128];
};



void
dnaAlphabets::initTablesForACGTSpace(void) {
  int i, j;

  for (i=0; i<256; i++) {
    whitespaceSymbol[i]      = isspace(i) ? 1 : 0;
    toLower[i]               = tolower(i);
    toUpper[i]               = toupper(i);
    letterToBits[i]          = (unsigned char)0xff;
    bitsToLetter[i]          = (unsigned char)'?';
    bitsToColor[i]           = (unsigned char)'?';
    complementSymbol[i]      = (unsigned char)'?';
  }

  for (i=0; i<128; i++)
    for (j=0; j<128; j++)
      IUPACidentity[i][j] = 0;

  letterToBits['a'] = letterToBits['A'] = (unsigned char)0x00;
  letterToBits['c'] = letterToBits['C'] = (unsigned char)0x01;
  letterToBits['g'] = letterToBits['G'] = (unsigned char)0x02;
  letterToBits['t'] = letterToBits['T'] = (unsigned char)0x03;

  letterToBits['0'] = (unsigned char)0x00;
  letterToBits['1'] = (unsigned char)0x01;
  letterToBits['2'] = (unsigned char)0x02;
  letterToBits['3'] = (unsigned char)0x03;

  bitsToLetter[0x00] = 'A';
  bitsToLetter[0x01] = 'C';
  bitsToLetter[0x02] = 'G';
  bitsToLetter[0x03] = 'T';

  bitsToColor[0x00] = '0';
  bitsToColor[0x01] = '1';
  bitsToColor[0x02] = '2';
  bitsToColor[0x03] = '3';

  complementSymbol['a'] = 't';  //  a
  complementSymbol['t'] = 'a';  //  t
  complementSymbol['u'] = 'a';  //  u, Really, only for RNA
  complementSymbol['g'] = 'c';  //  g
  complementSymbol['c'] = 'g';  //  c
  complementSymbol['y'] = 'r';  //  c t
  complementSymbol['r'] = 'y';  //  a g
  complementSymbol['s'] = 'w';  //  g c
  complementSymbol['w'] = 's';  //  a t
  complementSymbol['k'] = 'm';  //  t/u g
  complementSymbol['m'] = 'k';  //  a c
  complementSymbol['b'] = 'v';  //  c g t
  complementSymbol['d'] = 'h';  //  a g t
  complementSymbol['h'] = 'd';  //  a c t
  complementSymbol['v'] = 'b';  //  a c g
  complementSymbol['n'] = 'n';  //  a c g t

  complementSymbol['A'] = 'T';  //  a
  complementSymbol['T'] = 'A';  //  t
  complementSymbol['U'] = 'A';  //  u, Really, only for RNA
  complementSymbol['G'] = 'C';  //  g
  complementSymbol['C'] = 'G';  //  c
  complementSymbol['Y'] = 'R';  //  c t
  complementSymbol['R'] = 'Y';  //  a g
  complementSymbol['S'] = 'W';  //  g c
  complementSymbol['W'] = 'S';  //  a t
  complementSymbol['K'] = 'M';  //  t/u g
  complementSymbol['M'] = 'K';  //  a c
  complementSymbol['B'] = 'V';  //  c g t
  complementSymbol['D'] = 'H';  //  a g t
  complementSymbol['H'] = 'D';  //  a c t
  complementSymbol['V'] = 'B';  //  a c g
  complementSymbol['N'] = 'N';  //  a c g t

  complementSymbol['0'] = '0';  //  ColorSpace is self-complementing
  complementSymbol['1'] = '1';
  complementSymbol['2'] = '2';
  complementSymbol['3'] = '3';

  IUPACidentity['A']['A'] = 1;
  IUPACidentity['C']['C'] = 1;
  IUPACidentity['G']['G'] = 1;
  IUPACidentity['T']['T'] = 1;
  IUPACidentity['M']['A'] = 1;
  IUPACidentity['M']['C'] = 1;
  IUPACidentity['R']['A'] = 1;
  IUPACidentity['R']['G'] = 1;
  IUPACidentity['W']['A'] = 1;
  IUPACidentity['W']['T'] = 1;
  IUPACidentity['S']['C'] = 1;
  IUPACidentity['S']['G'] = 1;
  IUPACidentity['Y']['C'] = 1;
  IUPACidentity['Y']['T'] = 1;
  IUPACidentity['K']['G'] = 1;
  IUPACidentity['K']['T'] = 1;
  IUPACidentity['V']['A'] = 1;
  IUPACidentity['V']['C'] = 1;
  IUPACidentity['V']['G'] = 1;
  IUPACidentity['H']['A'] = 1;
  IUPACidentity['H']['C'] = 1;
  IUPACidentity['H']['T'] = 1;
  IUPACidentity['D']['A'] = 1;
  IUPACidentity['D']['G'] = 1;
  IUPACidentity['D']['T'] = 1;
  IUPACidentity['B']['C'] = 1;
  IUPACidentity['B']['G'] = 1;
  IUPACidentity['B']['T'] = 1;

  IUPACidentity['N']['A'] = 1;
  IUPACidentity['N']['C'] = 1;
  IUPACidentity['N']['G'] = 1;
  IUPACidentity['N']['T'] = 1;
  
  IUPACidentity['M']['M'] = 1;
  IUPACidentity['R']['R'] = 1;
  IUPACidentity['W']['W'] = 1;
  IUPACidentity['S']['S'] = 1;
  IUPACidentity['Y']['Y'] = 1;
  IUPACidentity['K']['K'] = 1;
  IUPACidentity['V']['V'] = 1;
  IUPACidentity['H']['W'] = 1;
  IUPACidentity['D']['D'] = 1;
  IUPACidentity['B']['B'] = 1;
  IUPACidentity['N']['N'] = 1;

  //  Order isn't important
  //
  for (i='A'; i<'Z'; i++)
    for (j='A'; j<'Z'; j++) {
      if (IUPACidentity[j][i])
        IUPACidentity[i][j] = 1;
    }

  //  Case isn't important
  //
  for (i='A'; i<'Z'; i++)
    for (j='A'; j<'Z'; j++) {
      if (IUPACidentity[j][i]) {
        IUPACidentity[tolower(i)][tolower(j)] = 1;
        IUPACidentity[tolower(i)][j         ] = 1;
        IUPACidentity[i         ][tolower(j)] = 1;
      }
    }
}



void
dnaAlphabets::initTablesForColorSpace(void) {
  int i, j;

  for (i=0; i<128; i++)
    for (j=0; j<128; j++)
      baseToColor[i][j] = '.';  //  Invalid

  //  Supports transforming a base sequence to a color sequence.

  //  Not sure how valid this is; treat every letter like it's a gap.
  //  We then override ACGT to be the correct encoding.
  for (i='a'; i<='z'; i++) {
    baseToColor['a'][i] = '4';
    baseToColor['c'][i] = '4';
    baseToColor['g'][i] = '4';
    baseToColor['t'][i] = '4';
    baseToColor['n'][i] = '4';
  }
  for (i='a'; i<='z'; i++) {
    baseToColor[i]['a'] = '0';
    baseToColor[i]['c'] = '1';
    baseToColor[i]['g'] = '2';
    baseToColor[i]['t'] = '3';
    baseToColor[i]['n'] = '4';
  }

  baseToColor['a']['a'] = '0';
  baseToColor['a']['c'] = '1';
  baseToColor['a']['g'] = '2';
  baseToColor['a']['t'] = '3';
  baseToColor['a']['n'] = '4';

  baseToColor['c']['a'] = '1';
  baseToColor['c']['c'] = '0';
  baseToColor['c']['g'] = '3';
  baseToColor['c']['t'] = '2';
  baseToColor['c']['n'] = '4';

  baseToColor['g']['a'] = '2';
  baseToColor['g']['c'] = '3';
  baseToColor['g']['g'] = '0';
  baseToColor['g']['t'] = '1';
  baseToColor['g']['n'] = '4';

  baseToColor['t']['a'] = '3';
  baseToColor['t']['c'] = '2';
  baseToColor['t']['g'] = '1';
  baseToColor['t']['t'] = '0';
  baseToColor['t']['n'] = '4';

  for (i='a'; i<='z'; i++)
    for (j='a'; j<='z'; j++) {
      baseToColor[toupper(i)][toupper(j)] = baseToColor[i][j];
      baseToColor[tolower(i)][toupper(j)] = baseToColor[i][j];
      baseToColor[toupper(i)][tolower(j)] = baseToColor[i][j];
      baseToColor[tolower(i)][tolower(j)] = baseToColor[i][j];
    }

  //  Supports composing colors

  baseToColor['0']['0'] = '0';
  baseToColor['0']['1'] = '1';
  baseToColor['0']['2'] = '2';
  baseToColor['0']['3'] = '3';
  baseToColor['0']['4'] = '4';

  baseToColor['1']['0'] = '1';
  baseToColor['1']['1'] = '0';
  baseToColor['1']['2'] = '3';
  baseToColor['1']['3'] = '2';
  baseToColor['1']['4'] = '4';

  baseToColor['2']['0'] = '2';
  baseToColor['2']['1'] = '3';
  baseToColor['2']['2'] = '0';
  baseToColor['2']['3'] = '1';
  baseToColor['2']['4'] = '4';

  baseToColor['3']['0'] = '3';
  baseToColor['3']['1'] = '2';
  baseToColor['3']['2'] = '1';
  baseToColor['3']['3'] = '0';
  baseToColor['3']['4'] = '4';

  //  Supports transforming color sequence to base sequence.

  baseToColor['a']['0'] = baseToColor['A']['0'] = 'a';
  baseToColor['a']['1'] = baseToColor['A']['1'] = 'c';
  baseToColor['a']['2'] = baseToColor['A']['2'] = 'g';
  baseToColor['a']['3'] = baseToColor['A']['3'] = 't';
  baseToColor['a']['4'] = baseToColor['A']['4'] = 'n';

  baseToColor['c']['0'] = baseToColor['C']['0'] = 'c';
  baseToColor['c']['1'] = baseToColor['C']['1'] = 'a';
  baseToColor['c']['2'] = baseToColor['C']['2'] = 't';
  baseToColor['c']['3'] = baseToColor['C']['3'] = 'g';
  baseToColor['c']['4'] = baseToColor['C']['4'] = 'n';

  baseToColor['g']['0'] = baseToColor['G']['0'] = 'g';
  baseToColor['g']['1'] = baseToColor['G']['1'] = 't';
  baseToColor['g']['2'] = baseToColor['G']['2'] = 'a';
  baseToColor['g']['3'] = baseToColor['G']['3'] = 'c';
  baseToColor['g']['4'] = baseToColor['G']['4'] = 'n';

  baseToColor['t']['0'] = baseToColor['T']['0'] = 't';
  baseToColor['t']['1'] = baseToColor['T']['1'] = 'g';
  baseToColor['t']['2'] = baseToColor['T']['2'] = 'c';
  baseToColor['t']['3'] = baseToColor['T']['3'] = 'a';
  baseToColor['t']['4'] = baseToColor['T']['4'] = 'n';

  baseToColor['n']['0'] = baseToColor['N']['0'] = 'a';
  baseToColor['n']['1'] = baseToColor['N']['1'] = 'c';
  baseToColor['n']['2'] = baseToColor['N']['2'] = 'g';
  baseToColor['n']['3'] = baseToColor['N']['3'] = 't';
  baseToColor['n']['4'] = baseToColor['N']['4'] = 'n';
}

