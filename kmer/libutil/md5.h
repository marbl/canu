#ifndef MD5_H
#define MD5_H

#include "libbritypes.h"

#ifdef __cplusplus
extern "C" {
#endif

//  XXX:  Notably missing is somethig that will do a checksum piece by piece


typedef struct {
  u64bit  a;
  u64bit  b;
  u32bit  i;
} md5_s;

//  Returns -1, 0, 1 depending on if a <, ==, > b.  Suitable for
//  qsort().
//
int     md5_compare(void const *a, void const *b);

//  Converts an md5_s into a character string.  s must be at least
//  33 bytes long.
char   *md5_toascii(md5_s *m, char *s);

//  Computes the md5 checksum on the string s.
md5_s  *md5_string(md5_s *m, char *s, u32bit l);


#ifdef __cplusplus
}
#endif

#endif  //  MD5_H
