#ifndef FASTA_SIMPLE_H
#define FASTA_SIMPLE_H

//
//  A very simple C based fasta reader.  No random access, just open a file,
//  get sequence, get sequence, get sequence, close file.
//

#include "libbritypes.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
  char     *header;
  u32bit    headerLen;
  char     *seq;
  u32bit    seqLen;
} fastaSequence;


typedef struct {
  FILE     *fileptr;
  char     *filename;
} fastaReader;


fastaReader    *FastAopen(char *filename);
fastaSequence  *FastAget(fastaReader *r);
void            FastAfree(fastaSequence *f);
void            FastAclose(fastaReader *r);


#ifdef __cplusplus
}
#endif

#endif  //  FASTA_SIMPLE_H
