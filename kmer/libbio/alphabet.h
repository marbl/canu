//
//  Automagically generated -- DO NOT EDIT!
//  See libbri/alphabet-generate.c for details.
//

#ifdef __cplusplus
extern "C" {
#endif

extern unsigned char   whitespaceSymbol[256];
extern unsigned char   toLower[256];
extern unsigned char   toUpper[256];
extern unsigned char   letterToBits[256];
extern unsigned char   bitsToLetter[256];
extern unsigned char   bitsToColor[256];
extern unsigned char   complementSymbol[256];
extern unsigned char   IUPACidentity[128][128];
extern unsigned char   baseToColor[128][128];

void initCompressionTablesForACGTSpace(void);
void initCompressionTablesForColorSpace(void);

#ifdef __cplusplus
}
#endif
