#ifndef FASTA_C_H
#define FASTA_C_H

#include <stdio.h>

//
//  A limited FastA reader, callable from C
//

#ifdef __cplusplus
extern "C" {
#endif

void *createFastA(char *file);

unsigned char  *getFastAsequence(void *c, int idx);
unsigned char  *getFastAsequenceReverseComplement(void *c, int idx);
int             getFastAsequenceLength(void *c, int idx);
unsigned char  *getFastAheader(void *c, int idx);

void  destroyFastA(void *c);

#ifdef __cplusplus
}
#endif

#endif  //  FASTA_C_H
