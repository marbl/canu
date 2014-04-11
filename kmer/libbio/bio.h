#ifndef BIO_H
#define BIO_H

#include "util.h"


#ifdef __cplusplus
extern "C" {
#endif


////////////////////////////////////////
//
//  alphabet
//
#include "alphabet.h"


////////////////////////////////////////
//
//  reversecomplement.c
//
char *reverseComplementSequence(char *seq, uint32 seqlen);
char *reverseString(char *seq, uint32 seqlen);


//  halign
//
//  N.B. align() (aka halign) was switched over to palloc() -- this
//  fixed any memory leaks, and gives a 30%-ish speed increase.  This
//  is thread safe (unless someone breaks palloc2()).
//
void
halign(const char *string1,
       const char *string2,
       const int len1,
       const int len2,
       char *alnline1,
       char *alnline2);


#ifdef __cplusplus
}
#endif


#endif  //  BIO_H
