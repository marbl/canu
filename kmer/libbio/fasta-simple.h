#ifndef FASTA_SIMPLE_H
#define FASTA_SIMPLE_H

//  A very simple C-callable multi-fasta file reader.  No random
//  access or other fancy features, just a simple sequence-by-sequence
//  reader.
//
//  Example usage:
//
//  fastaReader   *r = FastAopen("file.fasta");
//  fastaSequence *s = FastAget(r);
//
//  while (s) {
//    fprintf(stdout, "%s\n%s\n", s->header, s->seq);
//    FastAfree(s);
//    s = FastAget(r);
//  }
//
//  FastAclose(r);
//

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
  char     *header;
  int       headerLen;
  char     *seq;
  int       seqLen;
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
