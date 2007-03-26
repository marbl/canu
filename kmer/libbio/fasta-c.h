#ifndef FASTA_C_H
#define FASTA_C_H

#include <stdio.h>

//
//  A limited FastABase wrapper that is callable from C
//

#ifdef __cplusplus
extern "C" {
#endif

void           *createFastA(char *file);
void            destroyFastA(void *c);

char           *getFastAsequence(void *c, int idx);
int             getFastAsequenceLength(void *c, int idx);
char           *getFastAheader(void *c, int idx);

#ifdef __cplusplus
}
#endif

#endif  // FASTA_C_H
