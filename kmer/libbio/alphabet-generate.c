#include <stdio.h>
#include <ctype.h>

//  Instead of forcing client applications to explicitly call
//  initCompressionTables(), static tables are now generated.

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

void initCompressionTablesForACGTSpace(void);
void initCompressionTablesForColorSpace(void);

#include "alphabet-acgtspace.c"
#include "alphabet-colorspace.c"

int
main(int argc, char **argv) {
  int i, j;

  FILE *C = fopen("alphabet.c", "w");
  FILE *H = fopen("alphabet.h", "w");

  initCompressionTablesForACGTSpace();
  initCompressionTablesForColorSpace();

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
    fprintf(C, ",%d", whitespaceSymbol[i]);
  fprintf(C, " };\n");

  fprintf(H, "extern unsigned char   toLower[256];\n");
  fprintf(C, "unsigned char   toLower[256] = { %d", toLower[0]);
  for (i=1; i<256; i++)
    fprintf(C, ",%d", toLower[i]);
  fprintf(C, " };\n");

  fprintf(H, "extern unsigned char   toUpper[256];\n");
  fprintf(C, "unsigned char   toUpper[256] = { %d", toUpper[0]);
  for (i=1; i<256; i++)
    fprintf(C, ",%d", toUpper[i]);
  fprintf(C, " };\n");

  fprintf(H, "extern unsigned char   letterToBits[256];\n");
  fprintf(C, "unsigned char   letterToBits[256] = { %d", letterToBits[0]);
  for (i=1; i<256; i++)
    fprintf(C, ",%d", letterToBits[i]);
  fprintf(C, " };\n");

  fprintf(H, "extern unsigned char   bitsToLetter[256];\n");
  fprintf(C, "unsigned char   bitsToLetter[256] = { %d", bitsToLetter[0]);
  for (i=1; i<256; i++)
    fprintf(C, ",%d", bitsToLetter[i]);
  fprintf(C, " };\n");

  fprintf(H, "extern unsigned char   bitsToColor[256];\n");
  fprintf(C, "unsigned char   bitsToColor[256] = { %d", bitsToColor[0]);
  for (i=1; i<256; i++)
    fprintf(C, ",%d", bitsToColor[i]);
  fprintf(C, " };\n");

  fprintf(H, "extern unsigned char   complementSymbol[256];\n");
  fprintf(C, "unsigned char   complementSymbol[256] = { %d", complementSymbol[0]);
  for (i=1; i<256; i++)
    fprintf(C, ",%d", complementSymbol[i]);
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


  fprintf(H, "extern unsigned char   baseToColor[128][128];\n");
  fprintf(C, "unsigned char   baseToColor[128][128] = {\n");
  for (i=0; i<128; i++) {
    fprintf(C, " {%d", baseToColor[i][0]);
    for (j=1;j<128; j++)
      fprintf(C, ",%d", baseToColor[i][j]);
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
