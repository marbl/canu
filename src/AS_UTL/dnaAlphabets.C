
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
 *    Brian P. Walenz on 2014-DEC-08
 *      are Copyright 2014 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "dnaAlphabets.H"


dnaAlphabets  alphabet;


void
dnaAlphabets::initTablesForACGTSpace(void) {
  int i, j;

  for (i=0; i<256; i++) {
    _whitespaceSymbol[i]      = isspace(i) ? 1 : 0;
    _toLower[i]               = tolower(i);
    _toUpper[i]               = toupper(i);
    _letterToBits[i]          = (unsigned char)0xff;
    _bitsToLetter[i]          = (unsigned char)'?';
    _bitsToColor[i]           = (unsigned char)'?';
    _complementSymbol[i]      = (unsigned char)'?';
  }

  for (i=0; i<128; i++)
    for (j=0; j<128; j++)
      _IUPACidentity[i][j] = 0;

  _letterToBits['a'] = _letterToBits['A'] = (unsigned char)0x00;
  _letterToBits['c'] = _letterToBits['C'] = (unsigned char)0x01;
  _letterToBits['g'] = _letterToBits['G'] = (unsigned char)0x02;
  _letterToBits['t'] = _letterToBits['T'] = (unsigned char)0x03;

  _letterToBits['0'] = (unsigned char)0x00;
  _letterToBits['1'] = (unsigned char)0x01;
  _letterToBits['2'] = (unsigned char)0x02;
  _letterToBits['3'] = (unsigned char)0x03;

  _bitsToLetter[0x00] = 'A';
  _bitsToLetter[0x01] = 'C';
  _bitsToLetter[0x02] = 'G';
  _bitsToLetter[0x03] = 'T';

  _bitsToColor[0x00] = '0';
  _bitsToColor[0x01] = '1';
  _bitsToColor[0x02] = '2';
  _bitsToColor[0x03] = '3';

  _complementSymbol['a'] = 't';  //  a
  _complementSymbol['t'] = 'a';  //  t
  _complementSymbol['u'] = 'a';  //  u, Really, only for RNA
  _complementSymbol['g'] = 'c';  //  g
  _complementSymbol['c'] = 'g';  //  c
  _complementSymbol['y'] = 'r';  //  c t
  _complementSymbol['r'] = 'y';  //  a g
  _complementSymbol['s'] = 'w';  //  g c
  _complementSymbol['w'] = 's';  //  a t
  _complementSymbol['k'] = 'm';  //  t/u g
  _complementSymbol['m'] = 'k';  //  a c
  _complementSymbol['b'] = 'v';  //  c g t
  _complementSymbol['d'] = 'h';  //  a g t
  _complementSymbol['h'] = 'd';  //  a c t
  _complementSymbol['v'] = 'b';  //  a c g
  _complementSymbol['n'] = 'n';  //  a c g t

  _complementSymbol['A'] = 'T';  //  a
  _complementSymbol['T'] = 'A';  //  t
  _complementSymbol['U'] = 'A';  //  u, Really, only for RNA
  _complementSymbol['G'] = 'C';  //  g
  _complementSymbol['C'] = 'G';  //  c
  _complementSymbol['Y'] = 'R';  //  c t
  _complementSymbol['R'] = 'Y';  //  a g
  _complementSymbol['S'] = 'W';  //  g c
  _complementSymbol['W'] = 'S';  //  a t
  _complementSymbol['K'] = 'M';  //  t/u g
  _complementSymbol['M'] = 'K';  //  a c
  _complementSymbol['B'] = 'V';  //  c g t
  _complementSymbol['D'] = 'H';  //  a g t
  _complementSymbol['H'] = 'D';  //  a c t
  _complementSymbol['V'] = 'B';  //  a c g
  _complementSymbol['N'] = 'N';  //  a c g t

  _complementSymbol['0'] = '0';  //  ColorSpace is self-complementing
  _complementSymbol['1'] = '1';
  _complementSymbol['2'] = '2';
  _complementSymbol['3'] = '3';

  _IUPACidentity['A']['A'] = 1;
  _IUPACidentity['C']['C'] = 1;
  _IUPACidentity['G']['G'] = 1;
  _IUPACidentity['T']['T'] = 1;
  _IUPACidentity['M']['A'] = 1;
  _IUPACidentity['M']['C'] = 1;
  _IUPACidentity['R']['A'] = 1;
  _IUPACidentity['R']['G'] = 1;
  _IUPACidentity['W']['A'] = 1;
  _IUPACidentity['W']['T'] = 1;
  _IUPACidentity['S']['C'] = 1;
  _IUPACidentity['S']['G'] = 1;
  _IUPACidentity['Y']['C'] = 1;
  _IUPACidentity['Y']['T'] = 1;
  _IUPACidentity['K']['G'] = 1;
  _IUPACidentity['K']['T'] = 1;
  _IUPACidentity['V']['A'] = 1;
  _IUPACidentity['V']['C'] = 1;
  _IUPACidentity['V']['G'] = 1;
  _IUPACidentity['H']['A'] = 1;
  _IUPACidentity['H']['C'] = 1;
  _IUPACidentity['H']['T'] = 1;
  _IUPACidentity['D']['A'] = 1;
  _IUPACidentity['D']['G'] = 1;
  _IUPACidentity['D']['T'] = 1;
  _IUPACidentity['B']['C'] = 1;
  _IUPACidentity['B']['G'] = 1;
  _IUPACidentity['B']['T'] = 1;

  _IUPACidentity['N']['A'] = 1;
  _IUPACidentity['N']['C'] = 1;
  _IUPACidentity['N']['G'] = 1;
  _IUPACidentity['N']['T'] = 1;

  _IUPACidentity['M']['M'] = 1;
  _IUPACidentity['R']['R'] = 1;
  _IUPACidentity['W']['W'] = 1;
  _IUPACidentity['S']['S'] = 1;
  _IUPACidentity['Y']['Y'] = 1;
  _IUPACidentity['K']['K'] = 1;
  _IUPACidentity['V']['V'] = 1;
  _IUPACidentity['H']['W'] = 1;
  _IUPACidentity['D']['D'] = 1;
  _IUPACidentity['B']['B'] = 1;
  _IUPACidentity['N']['N'] = 1;

  //  Order isn't important
  //
  for (i='A'; i<'Z'; i++)
    for (j='A'; j<'Z'; j++) {
      if (_IUPACidentity[j][i])
        _IUPACidentity[i][j] = 1;
    }

  //  Case isn't important
  //
  for (i='A'; i<'Z'; i++)
    for (j='A'; j<'Z'; j++) {
      if (_IUPACidentity[j][i]) {
        _IUPACidentity[tolower(i)][tolower(j)] = 1;
        _IUPACidentity[tolower(i)][j         ] = 1;
        _IUPACidentity[i         ][tolower(j)] = 1;
      }
    }
}



void
dnaAlphabets::initTablesForColorSpace(void) {
  int i, j;

  for (i=0; i<128; i++)
    for (j=0; j<128; j++)
      _baseToColor[i][j] = '.';  //  Invalid

  //  Supports transforming a base sequence to a color sequence.

  //  Not sure how valid this is; treat every letter like it's a gap.
  //  We then override ACGT to be the correct encoding.
  for (i='a'; i<='z'; i++) {
    _baseToColor['a'][i] = '4';
    _baseToColor['c'][i] = '4';
    _baseToColor['g'][i] = '4';
    _baseToColor['t'][i] = '4';
    _baseToColor['n'][i] = '4';
  }
  for (i='a'; i<='z'; i++) {
    _baseToColor[i]['a'] = '0';
    _baseToColor[i]['c'] = '1';
    _baseToColor[i]['g'] = '2';
    _baseToColor[i]['t'] = '3';
    _baseToColor[i]['n'] = '4';
  }

  _baseToColor['a']['a'] = '0';
  _baseToColor['a']['c'] = '1';
  _baseToColor['a']['g'] = '2';
  _baseToColor['a']['t'] = '3';
  _baseToColor['a']['n'] = '4';

  _baseToColor['c']['a'] = '1';
  _baseToColor['c']['c'] = '0';
  _baseToColor['c']['g'] = '3';
  _baseToColor['c']['t'] = '2';
  _baseToColor['c']['n'] = '4';

  _baseToColor['g']['a'] = '2';
  _baseToColor['g']['c'] = '3';
  _baseToColor['g']['g'] = '0';
  _baseToColor['g']['t'] = '1';
  _baseToColor['g']['n'] = '4';

  _baseToColor['t']['a'] = '3';
  _baseToColor['t']['c'] = '2';
  _baseToColor['t']['g'] = '1';
  _baseToColor['t']['t'] = '0';
  _baseToColor['t']['n'] = '4';

  for (i='a'; i<='z'; i++)
    for (j='a'; j<='z'; j++) {
      _baseToColor[toupper(i)][toupper(j)] = _baseToColor[i][j];
      _baseToColor[tolower(i)][toupper(j)] = _baseToColor[i][j];
      _baseToColor[toupper(i)][tolower(j)] = _baseToColor[i][j];
      _baseToColor[tolower(i)][tolower(j)] = _baseToColor[i][j];
    }

  //  Supports composing colors

  _baseToColor['0']['0'] = '0';
  _baseToColor['0']['1'] = '1';
  _baseToColor['0']['2'] = '2';
  _baseToColor['0']['3'] = '3';
  _baseToColor['0']['4'] = '4';

  _baseToColor['1']['0'] = '1';
  _baseToColor['1']['1'] = '0';
  _baseToColor['1']['2'] = '3';
  _baseToColor['1']['3'] = '2';
  _baseToColor['1']['4'] = '4';

  _baseToColor['2']['0'] = '2';
  _baseToColor['2']['1'] = '3';
  _baseToColor['2']['2'] = '0';
  _baseToColor['2']['3'] = '1';
  _baseToColor['2']['4'] = '4';

  _baseToColor['3']['0'] = '3';
  _baseToColor['3']['1'] = '2';
  _baseToColor['3']['2'] = '1';
  _baseToColor['3']['3'] = '0';
  _baseToColor['3']['4'] = '4';

  //  Supports transforming color sequence to base sequence.

  _baseToColor['a']['0'] = _baseToColor['A']['0'] = 'a';
  _baseToColor['a']['1'] = _baseToColor['A']['1'] = 'c';
  _baseToColor['a']['2'] = _baseToColor['A']['2'] = 'g';
  _baseToColor['a']['3'] = _baseToColor['A']['3'] = 't';
  _baseToColor['a']['4'] = _baseToColor['A']['4'] = 'n';

  _baseToColor['c']['0'] = _baseToColor['C']['0'] = 'c';
  _baseToColor['c']['1'] = _baseToColor['C']['1'] = 'a';
  _baseToColor['c']['2'] = _baseToColor['C']['2'] = 't';
  _baseToColor['c']['3'] = _baseToColor['C']['3'] = 'g';
  _baseToColor['c']['4'] = _baseToColor['C']['4'] = 'n';

  _baseToColor['g']['0'] = _baseToColor['G']['0'] = 'g';
  _baseToColor['g']['1'] = _baseToColor['G']['1'] = 't';
  _baseToColor['g']['2'] = _baseToColor['G']['2'] = 'a';
  _baseToColor['g']['3'] = _baseToColor['G']['3'] = 'c';
  _baseToColor['g']['4'] = _baseToColor['G']['4'] = 'n';

  _baseToColor['t']['0'] = _baseToColor['T']['0'] = 't';
  _baseToColor['t']['1'] = _baseToColor['T']['1'] = 'g';
  _baseToColor['t']['2'] = _baseToColor['T']['2'] = 'c';
  _baseToColor['t']['3'] = _baseToColor['T']['3'] = 'a';
  _baseToColor['t']['4'] = _baseToColor['T']['4'] = 'n';

  _baseToColor['n']['0'] = _baseToColor['N']['0'] = 'a';
  _baseToColor['n']['1'] = _baseToColor['N']['1'] = 'c';
  _baseToColor['n']['2'] = _baseToColor['N']['2'] = 'g';
  _baseToColor['n']['3'] = _baseToColor['N']['3'] = 't';
  _baseToColor['n']['4'] = _baseToColor['N']['4'] = 'n';
}

