#include <stdio.h>
#include <ctype.h>
#include "alphabet.h"

void
initCompressionTablesForACGTSpace(void) {
  int i, j;

  for (i=0; i<256; i++) {
    whitespaceSymbol[i]      = (unsigned char)0x00;
    toLower[i]               = (unsigned char)0x00;
    toUpper[i]               = (unsigned char)0x00;
    compressSymbol[i]        = (unsigned char)0x00;
    validSymbol[i]           = (unsigned char)0x00;
    decompressSymbol[i]      = (unsigned char)0x00;
    complementSymbol[i]      = (unsigned char)i;
    validCompressedSymbol[i] = (unsigned char)0xff;
  }

  for (i=0; i<256; i++) {
    whitespaceSymbol[i] = isspace(i) ? 1 : 0;
    toLower[i]          = tolower(i);
    toUpper[i]          = toupper(i);
  }

  for (i=0; i<128; i++)
    for (j=0; j<128; j++)
      IUPACidentity[i][j] = 0;

  compressSymbol['a'] = compressSymbol['A'] = (unsigned char)0x00;
  compressSymbol['c'] = compressSymbol['C'] = (unsigned char)0x01;
  compressSymbol['g'] = compressSymbol['G'] = (unsigned char)0x02;
  compressSymbol['t'] = compressSymbol['T'] = (unsigned char)0x03;

  validSymbol['a'] = validSymbol['A'] = 1;
  validSymbol['c'] = validSymbol['C'] = 1;
  validSymbol['g'] = validSymbol['G'] = 1;
  validSymbol['t'] = validSymbol['T'] = 1;

  validCompressedSymbol['a'] = validCompressedSymbol['A'] = (unsigned char)0x00;
  validCompressedSymbol['c'] = validCompressedSymbol['C'] = (unsigned char)0x01;
  validCompressedSymbol['g'] = validCompressedSymbol['G'] = (unsigned char)0x02;
  validCompressedSymbol['t'] = validCompressedSymbol['T'] = (unsigned char)0x03;

  decompressSymbol[0x00] = 'A';
  decompressSymbol[0x01] = 'C';
  decompressSymbol[0x02] = 'G';
  decompressSymbol[0x03] = 'T';

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



