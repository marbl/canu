#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>

#include "util.h"

#include "md5lib/global.h"
#include "md5lib/md5.h"


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
    m = (md5_s *)malloc(sizeof(md5_s));
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

