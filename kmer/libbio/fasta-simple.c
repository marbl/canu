#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <ctype.h>
#include "fasta-simple.h"


fastaReader*
FastAopen(char *filename) {
  fastaReader *r;

  errno = 0;
  r = (fastaReader *)malloc(sizeof(fastaReader));
  if (r == 0L) {
    fprintf(stderr, "FastAopen()-- Couldn't allocate fastaReader structure for '%s'\n%s\n", filename, strerror(errno));
    exit(1);
  }

  r->fileptr = fopen(filename, "r");
  if (errno) {
    fprintf(stderr, "FastAopen()-- Couldn't open '%s'\n%s\n", filename, strerror(errno));
    exit(1);
  }

  r->filename = filename;

  return(r);
}


fastaSequence*
FastAget(fastaReader *r) {
  fastaSequence *fs;
  char           ch;

  int            hLen = 0;
  int            hMax = 128;
  char          *h    = (char *)malloc(sizeof(char) * (hMax + 1));
  int            sLen = 0;
  int            sMax = 16 * 1024 * 1024;
  char          *s    = (char *)malloc(sizeof(char) * (sMax + 1));

  //  Return end if we are at the end of the file
  //
  if (feof(r->fileptr))
    return(NULL);

  //  Skip whitespace at the start of the sequence.
  //
  ch = getc(r->fileptr);
  while ((!feof(r->fileptr)) && isspace(ch))
    ch = getc(r->fileptr);

  //  We should be at a '>' character now.  Fail if not.
  //
  if (ch != '>') {
    fprintf(stderr, "FastAget()-- In %s, expected '>' at beginning of defline, got '%c' instead.\n",
            r->filename, ch);
    exit(1);
  }

  //  Read the header
  //
  while ((!feof(r->fileptr)) && (ch != '\r') && (ch != '\n')) {
    if (hLen >= hMax) {
      hMax += 128;
      char *htmp = (char *)malloc(sizeof(char) * (hMax + 1));
      memcpy(htmp, h, sizeof(char) * hLen);
      free(h);
      h = htmp;
    }
    h[hLen++] = ch;

    ch = getc(r->fileptr);
  }
  h[hLen] = 0;

  //  Read the sequence, ignoring whitespace
  //
  while ((!feof(r->fileptr)) && (ch != '>')) {
    if (!isspace(ch)) {
      if (sLen >= sMax) {
        sMax += 32 * 1024 * 1024;
        char *stmp = (char *)malloc(sizeof(char) * (sMax + 1));
        memcpy(stmp, s, sizeof(char) * sLen);
        free(s);
        s = stmp;
      }
      s[sLen++] = ch;
    }

    ch = getc(r->fileptr);
  }
  s[sLen] = 0;

  //  Skip whitespace at the end of the sequence (gets us to eof at
  //  the very end).
  //
  while ((!feof(r->fileptr)) && isspace(ch))
    ch = getc(r->fileptr);

  //  If not eof, unget the last character.  This pushes back the '>'
  //  of the next sequence.
  //
  if (!feof(r->fileptr))
    ungetc(ch, r->fileptr);

  fs = (fastaSequence *)malloc(sizeof(fastaSequence));
  fs->header    = h;
  fs->headerLen = hLen;
  fs->seq       = s;
  fs->seqLen    = sLen; 

  return(fs);
}


void
FastAfree(fastaSequence *f) {
  if (f) {
    free(f->header);
    free(f->seq);
    free(f);
  }
}


void
FastAclose(fastaReader *r) {
  fclose(r->fileptr);
  free(r);
}
