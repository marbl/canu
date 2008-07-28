#include <stdio.h>
#include <ctype.h>
#include "alphabet.h"

void
initCompressionTablesForACGTSpace(void) {
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



