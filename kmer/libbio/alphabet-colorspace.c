#include <stdio.h>
#include <ctype.h>
#include "alphabet.h"

void
initCompressionTablesForColorSpace(void) {
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
