#include <stdio.h>
#include <ctype.h>
#include "alphabet.h"

void
initCompressionTablesForColorSpace(void) {
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

  compressSymbol['1'] = (unsigned char)0x00;
  compressSymbol['2'] = (unsigned char)0x01;
  compressSymbol['3'] = (unsigned char)0x02;
  compressSymbol['4'] = (unsigned char)0x03;

  validSymbol['1'] = (unsigned char)0x00;
  validSymbol['2'] = (unsigned char)0x01;
  validSymbol['3'] = (unsigned char)0x02;
  validSymbol['4'] = (unsigned char)0x03;

  validCompressedSymbol['1'] = (unsigned char)0x00;
  validCompressedSymbol['2'] = (unsigned char)0x01;
  validCompressedSymbol['3'] = (unsigned char)0x02;
  validCompressedSymbol['4'] = (unsigned char)0x03;

  decompressSymbol[0x00] = '1';
  decompressSymbol[0x01] = '2';
  decompressSymbol[0x02] = '3';
  decompressSymbol[0x03] = '4';

  complementSymbol['1'] = '1';  //  Self complementing
  complementSymbol['2'] = '2';
  complementSymbol['3'] = '3';
  complementSymbol['4'] = '4';
}
