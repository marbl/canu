#include <stdio.h>

//
//  Instead of forcing client applications to explicitly call
//  initCompressionTables(), static tables are now generated.
//
//  This file can be made much prettier, but it has not so that we can
//  still use the initCompressionTables() method (though I can't think
//  of a reason to do that).
//
//  Unless you *REALLY* *REALLY* know what you're doing, DO NOT
//  EVER change these constants.  Various bits depend on:
//    ACGT being two bits
//    A=00, C=01, G=10 and T=11
//

unsigned char  compressSymbol[256];
unsigned char  validSymbol[256];
unsigned char  decompressSymbol[256];
unsigned char  complementSymbol[256];
unsigned char  validCompressedSymbol[256];

void
initCompressionTables(void) {

  //
  //  Non-bug fix: complementSymbol['x'] would return 0, which would
  //  terminate a reverse complemented sequence abnormally early.
  //  Make it return 'x' instead.
  //
  //  bpw, Tue Oct 22 15:31:03 EDT 2002
  //

  unsigned int i=256;
  while (i) {
    i--;
    compressSymbol[i]        = (unsigned char)0x00;
    validSymbol[i]           = (unsigned char)0x00;
    decompressSymbol[i]      = (unsigned char)0x00;
    complementSymbol[i]      = (unsigned char)i;
    validCompressedSymbol[i] = (unsigned char)0xff;
  }

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
}

#ifdef MAIN
int
main(int argc, char **argv) {
  unsigned int i;

  FILE *C = fopen("alphabet.c", "w");
  FILE *H = fopen("alphabet.h", "w");

  initCompressionTables();

  fprintf(H, "//\n");
  fprintf(H, "//  Automagically generated -- DO NOT EDIT!\n");
  fprintf(H, "//  See %s for details.\n", __FILE__);
  fprintf(H, "//\n");
  fprintf(H, "\n");
  fprintf(H, "#ifdef __cplusplus\n");
  fprintf(H, "extern \"C\" {\n");
  fprintf(H, "#endif\n");
  fprintf(H, "\n");
  fprintf(H, "extern unsigned char const   compressSymbol[256];\n");

  fprintf(C, "//\n");
  fprintf(C, "//  Automagically generated -- DO NOT EDIT!\n");
  fprintf(C, "//  See %s for details.\n", __FILE__);
  fprintf(C, "//\n");
  fprintf(C, "\n");
  fprintf(C, "unsigned char const   compressSymbol[256] = { %d", compressSymbol[0]);
  for (i=1; i<256; i++)
    fprintf(C, ", %d", compressSymbol[i]);
  fprintf(C, " };\n");

  fprintf(H, "extern const unsigned char   validSymbol[256];\n");
  fprintf(C, "const unsigned char   validSymbol[256] = { %d", validSymbol[0]);
  for (i=1; i<256; i++)
    fprintf(C, ", %d", validSymbol[i]);
  fprintf(C, " };\n");

  fprintf(H, "extern const unsigned char   decompressSymbol[256];\n");
  fprintf(C, "const unsigned char   decompressSymbol[256] = { %d", decompressSymbol[0]);
  for (i=1; i<256; i++)
    fprintf(C, ", %d", decompressSymbol[i]);
  fprintf(C, " };\n");

  fprintf(H, "extern const unsigned char   complementSymbol[256];\n");
  fprintf(C, "const unsigned char   complementSymbol[256] = { %d", complementSymbol[0]);
  for (i=1; i<256; i++)
    fprintf(C, ", %d", complementSymbol[i]);
  fprintf(C, " };\n");

  fprintf(H, "extern const unsigned char   validCompressedSymbol[256];\n");
  fprintf(C, "const unsigned char   validCompressedSymbol[256] = { %d", validCompressedSymbol[0]);
  for (i=1; i<256; i++)
    fprintf(C, ", %d", validCompressedSymbol[i]);
  fprintf(C, " };\n");

  fprintf(H, "#ifdef __cplusplus\n");
  fprintf(H, "}\n");
  fprintf(H, "#endif\n");

  return(0);
}
#endif
