#include <stdio.h>
#include <ctype.h>

//  Instead of forcing client applications to explicitly call
//  initCompressionTables(), static tables are now generated.

#include "alphabet.h"  //  auto ganerated
#include "alphabet.c"  //  auto generated
#include "alphabet-acgtspace.c"
#include "alphabet-colorspace.c"

int
main(int argc, char **argv) {
  int i, j;

  FILE *C = fopen("alphabet.c", "w");
  FILE *H = fopen("alphabet.h", "w");

  initCompressionTablesForACGTSpace();

  fprintf(H, "//\n");
  fprintf(H, "//  Automagically generated -- DO NOT EDIT!\n");
  fprintf(H, "//  See libbri/alphabet-generate.c for details.\n");
  fprintf(H, "//\n");
  fprintf(H, "\n");
  fprintf(H, "#ifdef __cplusplus\n");
  fprintf(H, "extern \"C\" {\n");
  fprintf(H, "#endif\n");
  fprintf(H, "\n");

  fprintf(C, "//\n");
  fprintf(C, "//  Automagically generated -- DO NOT EDIT!\n");
  fprintf(C, "//  See %s for details.\n", __FILE__);
  fprintf(C, "//\n");

  fprintf(H, "extern unsigned char   whitespaceSymbol[256];\n");
  fprintf(C, "unsigned char   whitespaceSymbol[256] = { %d", whitespaceSymbol[0]);
  for (i=1; i<256; i++)
    fprintf(C, ", %d", whitespaceSymbol[i]);
  fprintf(C, " };\n");

  fprintf(H, "extern unsigned char   toLower[256];\n");
  fprintf(C, "unsigned char   toLower[256] = { %d", toLower[0]);
  for (i=1; i<256; i++)
    fprintf(C, ", %d", toLower[i]);
  fprintf(C, " };\n");

  fprintf(H, "extern unsigned char   toUpper[256];\n");
  fprintf(C, "unsigned char   toUpper[256] = { %d", toUpper[0]);
  for (i=1; i<256; i++)
    fprintf(C, ", %d", toUpper[i]);
  fprintf(C, " };\n");

  fprintf(H, "extern unsigned char   compressSymbol[256];\n");
  fprintf(C, "unsigned char   compressSymbol[256] = { %d", compressSymbol[0]);
  for (i=1; i<256; i++)
    fprintf(C, ", %d", compressSymbol[i]);
  fprintf(C, " };\n");

  fprintf(H, "extern unsigned char   validSymbol[256];\n");
  fprintf(C, "unsigned char   validSymbol[256] = { %d", validSymbol[0]);
  for (i=1; i<256; i++)
    fprintf(C, ", %d", validSymbol[i]);
  fprintf(C, " };\n");

  fprintf(H, "extern unsigned char   decompressSymbol[256];\n");
  fprintf(C, "unsigned char   decompressSymbol[256] = { %d", decompressSymbol[0]);
  for (i=1; i<256; i++)
    fprintf(C, ", %d", decompressSymbol[i]);
  fprintf(C, " };\n");

  fprintf(H, "extern unsigned char   complementSymbol[256];\n");
  fprintf(C, "unsigned char   complementSymbol[256] = { %d", complementSymbol[0]);
  for (i=1; i<256; i++)
    fprintf(C, ", %d", complementSymbol[i]);
  fprintf(C, " };\n");

  fprintf(H, "extern unsigned char   validCompressedSymbol[256];\n");
  fprintf(C, "unsigned char   validCompressedSymbol[256] = { %d", validCompressedSymbol[0]);
  for (i=1; i<256; i++)
    fprintf(C, ", %d", validCompressedSymbol[i]);
  fprintf(C, " };\n");

  fprintf(H, "extern unsigned char   IUPACidentity[128][128];\n");
  fprintf(C, "unsigned char   IUPACidentity[128][128] = {\n");
  for (i=0; i<128; i++) {
    fprintf(C, " {");
    if (IUPACidentity[i][0])
      fprintf(C, "1");
    else
      fprintf(C, "0");
    for (j=1;j<128; j++) {
      if (IUPACidentity[i][j])
        fprintf(C, ",1");
      else
        fprintf(C, ",0");
    }
    fprintf(C, "},\n");
  }
  fprintf(C, "};\n");

  fprintf(H, "\n");
  fprintf(H, "void initCompressionTablesForACGTSpace(void);\n");
  fprintf(H, "void initCompressionTablesForColorSpace(void);\n");

  fprintf(H, "\n");
  fprintf(H, "#ifdef __cplusplus\n");
  fprintf(H, "}\n");
  fprintf(H, "#endif\n");

  return(0);
}
