#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>

#include "bri.h"

#include "../external/md5lib/global.h"
#include "../external/md5lib/md5.h"


int
md5_compare(void const *a, void const *b) {
  md5_s *A = (md5_s *)a;
  md5_s *B = (md5_s *)b;

  if (A->a < B->a) return(-1);
  if (A->a > B->a) return(1);
  if (A->b < B->b) return(-1);
  if (A->b > B->b) return(1);
  return(0);
}


static const char *md5_letters = "0123456789abcdef";

char*
md5_toascii(md5_s *m, char *s) {
  int i;
  for (i=0; i<16; i++) {
    s[15-i   ] = md5_letters[(m->a >> 4*i) & 0x0f];
    s[15-i+16] = md5_letters[(m->b >> 4*i) & 0x0f];
  }
  s[32] = 0;

  return(s);
}



md5_s*
md5_string(md5_s *m, char *s, u32bit l) {
  MD5_CTX         ctx;
  unsigned char   dig[16];
  int             i = 0;

  if (m == NULL) {
    errno = 0;
    m = (md5_s *)malloc(sizeof(md5_s *));
    if (errno) {
      fprintf(stderr, "md5_string()-- Can't allocate a md5_s.\n%s\n", strerror(errno));
      exit(1);
    }
  }

  MD5Init(&ctx);
  MD5Update(&ctx, (unsigned char*)s, l);
  MD5Final(dig, &ctx);

  m->a = dig[0];
  while (i<8) {
    m->a <<= 8;
    m->a |= dig[i++];
  }

  m->b  = dig[i++];
  while (i<16) {
    m->b <<= 8;
    m->b |= dig[i++];
  }

  return(m);
}



static
md5_increment_s*
md5_increment_initialize(void) {
  md5_increment_s *m;

  errno = 0;
  m = (md5_increment_s *)malloc(sizeof(md5_increment_s));
  if (errno) {
    fprintf(stderr, "md5_increment_*()-- Can't allocate a md5_increment_s.\n%s\n", strerror(errno));
    exit(1);
  }

  m->context = (MD5_CTX *)malloc(sizeof(MD5_CTX));
  if (errno) {
    fprintf(stderr, "md5_increment_*()-- Can't allocate a md5 context.\n%s\n", strerror(errno));
    exit(1);
  }
  MD5Init((MD5_CTX *)m->context);

  m->bufferPos = 0;

  return(m);
}


md5_increment_s*
md5_increment_char(md5_increment_s *m, char s) {

  if (m == NULL)
    m = md5_increment_initialize();

  m->buffer[m->bufferPos++] = s;

  if (m->bufferPos == MD5_BUFFER_SIZE) {
    MD5Update((MD5_CTX *)m->context, m->buffer, m->bufferPos);
    m->bufferPos = 0;
  }

  return(m);
}


md5_increment_s*
md5_increment_block(md5_increment_s *m, char *s, u32bit l) {

  if (m == NULL)
    m = md5_increment_initialize();

  MD5Update((MD5_CTX *)m->context, (unsigned char*)s, l);

  return(m);
}


void
md5_increment_finalize(md5_increment_s *m) {
  MD5_CTX        *ctx = (MD5_CTX *)m->context;
  unsigned char   dig[16];
  int             i = 0;

  if (m->bufferPos > 0) {
    MD5Update((MD5_CTX *)m->context, m->buffer, m->bufferPos);
    m->bufferPos = 0;
  }

  MD5Final(dig, ctx);

  m->a = dig[0];
  while (i<8) {
    m->a <<= 8;
    m->a |= dig[i++];
  }

  m->b  = dig[i++];
  while (i<16) {
    m->b <<= 8;
    m->b |= dig[i++];
  }

  m->context = 0L;

  free(ctx);
}

void
md5_increment_destroy(md5_increment_s *m) {
  free(m);
}












#ifdef MAIN

//
//  Performs the md5 test suite.
//

int
testit(char *str, char *ans) {
  md5_s  m;
  char   r[33];
  int    ret = 0;

  md5_toascii(md5_string(&m, str, strlen(str)), r);
  ret = strcmp(r, ans);
  printf("%s should equal %s -- %s\n",
         r, ans,
         ret ? "FAILED" : "PASSED");
  return(ret == 0);
}

int
main(int argc, char **argv) {
  int ret = 7;
  
  ret -= testit("", "d41d8cd98f00b204e9800998ecf8427e");
  ret -= testit("a", "0cc175b9c0f1b6a831c399e269772661");
  ret -= testit("abc", "900150983cd24fb0d6963f7d28e17f72");
  ret -= testit("message digest", "f96b697d7cb7938d525a2f31aaf161d0");
  ret -= testit("abcdefghijklmnopqrstuvwxyz", "c3fcd3d76192e4007dfb496cca67e13b");
  ret -= testit("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789", "d174ab98d277d9f5a5611c2c9f419d9f");
  ret -= testit("12345678901234567890123456789012345678901234567890123456789012345678901234567890", "57edf4a22be3c955ac49da2e2107b67a");
  if (ret == 0) {
    printf("All md5 tests PASSED.\n");
  } else {
    printf("Some md5 Tests FAILED.\n");
  }
  exit(ret);
}
#endif
