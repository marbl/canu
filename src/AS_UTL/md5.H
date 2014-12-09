#include "AS_global.H"



typedef struct {
  uint64  a;
  uint64  b;
  uint32  i;    //  the iid, used in leaff
  uint32  pad;  //  keep us size compatible between 32- and 64-bit machines.
} md5_s;

#define MD5_BUFFER_SIZE   32*1024

typedef struct {
  uint64           a;
  uint64           b;
  void            *context;
  int              bufferPos;
  unsigned char    buffer[MD5_BUFFER_SIZE];
} md5_increment_s;


//  Returns -1, 0, 1 depending on if a <, ==, > b.  Suitable for
//  qsort().
//
int     md5_compare(void const *a, void const *b);


//  Converts an md5_s into a character string.  s must be at least
//  33 bytes long.
//
char   *md5_toascii(md5_s *m, char *s);


//  Computes the md5 checksum on the string s.
//
md5_s  *md5_string(md5_s *m, char *s, uint32 l);


//  Computes an md5 checksum piece by piece.
//
//  If m is NULL, a new md5_increment_s is allocated and returned.
//
md5_increment_s  *md5_increment_char(md5_increment_s *m, char s);
md5_increment_s  *md5_increment_block(md5_increment_s *m, char *s, uint32 l);
void              md5_increment_finalize(md5_increment_s *m);
void              md5_increment_destroy(md5_increment_s *m);
