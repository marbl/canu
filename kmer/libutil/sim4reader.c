#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include "sim4reader.h"


int
estIDcompare(const void *a, const void *b) {
  sim4polish *A = (*(sim4polish **)a);
  sim4polish *B = (*(sim4polish **)b);

  if (A == 0L)  return(1);
  if (B == 0L)  return(-1);

  if (A->estID < B->estID) return(-1);
  if (A->estID > B->estID) return(1);
  if (A->genID < B->genID) return(-1);
  if (A->genID > B->genID) return(1);
  if (A->genLo < B->genLo) return(-1);
  if (A->genLo > B->genLo) return(1);
  if (A->exons[0].genFrom < B->exons[0].genFrom) return(-1);
  if (A->exons[0].genFrom > B->exons[0].genFrom) return(1);
#ifdef DEVELOPMENT_CODE
  fprintf(stderr, "WARN: Sort not stable!\n");
#endif
  return(0);
}


int
genIDcompare(const void *a, const void *b) {
  sim4polish *A = (*(sim4polish **)a);
  sim4polish *B = (*(sim4polish **)b);

  if (A == 0L)  return(1);
  if (B == 0L)  return(-1);

  if (A->genID < B->genID) return(-1);
  if (A->genID > B->genID) return(1);
  if (A->estID < B->estID) return(-1);
  if (A->estID > B->estID) return(1);
  if (A->genLo < B->genLo) return(-1);
  if (A->genLo > B->genLo) return(1);
  if (A->exons[0].genFrom < B->exons[0].genFrom) return(-1);
  if (A->exons[0].genFrom > B->exons[0].genFrom) return(1);
#ifdef DEVELOPMENT_CODE
  fprintf(stderr, "WARN: Sort not stable!\n");
#endif
  return(0);
}



_line *newLine(void) {
  _line *L;
  errno = 0;
  L = (_line *)malloc(sizeof(_line));
  if (errno) {
    fprintf(stderr, "malloc() error in sim4reader::newLine\n%s\n", strerror(errno));
    exit(1);
  }
  L->l = 0;
  L->a = 10;
  L->lineNumber = 0;
  errno = 0;
  L->s = (char *)malloc(sizeof(char) * L->a);
  if (errno) {
    fprintf(stderr, "malloc() error in sim4reader::newLine\n%s\n", strerror(errno));
    exit(1);
  }
  return(L);
}

void
deleteLine(_line *L) {
  free(L->s);
  free(L);
}

void
readLine(FILE *F, _line *L) {

  L->l = 0;

  if (fgets(L->s, L->a, F)) {
    L->lineNumber++;

    L->l = strlen(L->s);

    while ((L->l + 1 >= L->a) &&
           (L->s[L->l - 1] != '\n')) {
      L->a *= 2;
      errno = 0;
      L->s  = (char *)realloc(L->s, sizeof(char) * L->a);
      if (errno) {
        fprintf(stderr, "realloc() error in sim4reader::readLine -- can't allocate more char\n%s\n", strerror(errno));
        exit(1);
      }
      fgets(L->s + L->l, L->a - L->l, F);
      L->l += strlen(L->s + L->l);
    }

    //  Trim newlines from the end of the string
    //
    while ((L->l > 0) &&
           ((L->s[L->l - 1] == '\r') ||
            (L->s[L->l - 1] == '\n')))
      L->s[--L->l] = 0;
  }
}

//
//  XXX:  How does this handle truncated files?
//

sim4polish*
readPolish(FILE *F) {
  sim4polish      *p = 0L;
  _line           *l = newLine();
  int              r;
  int              ef, et, gf, gt, nm, nn, id;
  int              el = 0;
  int              em = 256;
  sim4polishExon  *ex = 0L;
  char             mOri[65];
  char             sOri[65];

  readLine(F, l);
  while(!feof(F) && strcmp(l->s, "sim4begin")) {
    fprintf(stderr, "sim4reader: Got %d:'%s', expecting 'sim4begin'\n", l->lineNumber, l->s);
    readLine(F, l);
  }

  if (feof(F))
    return(0L);

  errno = 0;
  p = (sim4polish *)malloc(sizeof(sim4polish));
  if (errno) {
    fprintf(stderr, "malloc() error in sim4reader::readPolish -- can't allocate sim4polish\n%s\n", strerror(errno));
    exit(1);
  }

  //  Read the description line
  //
  readLine(F, l);

  //  Convert '-' into ' ', on the assumption that this is the description
  //  line.  This allows us to use scanf properly.
  //
  for (r=l->l; r--;)
    if (l->s[r] == '-')
      l->s[r] = ' ';

  mOri[0] = 0;
  sOri[0] = 0;
  r = sscanf(l->s, "%d[%d %d %d] %d[%d %d] <%d %d %d %s %s>",
             &p->estID,
             &p->estLen,
             &p->estPolyA,
             &p->estPolyT,
             &p->genID,
             &p->genLo,
             &p->genHi,
             &p->numMatches,
             &p->numMatchesN,
             &p->percentIdentity,
             mOri, sOri);
  if (r != 12) {
    fprintf(stderr, "sim4reader: Got %d:'%s'\n", l->lineNumber, l->s);
    fprintf(stderr, "sim4reader: Expecting description line, found %d tokens instead of 12.\n", r);
    fprintf(stderr, "sim4reader: mOri='%s'\n", mOri);
    fprintf(stderr, "sim4reader: sOri='%s'\n", sOri);
  }

  switch (mOri[0]) {
    case 'f':
      p->matchOrientation = MATCH_FORWARD;
      break;
    case 'c':
      p->matchOrientation = MATCH_COMPLEMENT;
      break;
    case 'r':
      //  BUG FIX -- old version of sim4 used "reverse-intractable"
      //  instead of "complement-intractable"
      p->matchOrientation = MATCH_COMPLEMENT;
      break;
    default:
      fprintf(stderr, "sim4reader: Got %d:'%s'\n", l->lineNumber, l->s);
      fprintf(stderr, "sim4reader: Unknown match orientation.\n");
      break;
  }

  switch (sOri[2]) {
    case 'r':
      p->strandOrientation = STRAND_POSITIVE;
      break;
    case 'v':
      p->strandOrientation = STRAND_NEGATIVE;
      break;
    case 'k':
      p->strandOrientation = STRAND_UNKNOWN;
      break;
    case 't':
      p->strandOrientation = STRAND_INTRACTABLE;
      break;
    case 'i':
      p->strandOrientation = STRAND_FAILED;
      break;
    default:
      fprintf(stderr, "sim4reader: Got %d:'%s'\n", l->lineNumber, l->s);
      fprintf(stderr, "sim4reader: Unknown strand orientation.\n");
      break;
  }

  p->querySeqIdentity = (int)floor(100 * (double)(p->numMatches) / (double)(p->estLen - p->estPolyA - p->estPolyT));

  readLine(F, l);

  p->comment = 0L;
  if (strncmp(l->s, "comment", 7) == 0) {
    errno = 0;
    p->comment = (char *)malloc(sizeof(char) * l->l - 7);
    if (errno) {
      fprintf(stderr, "malloc() error in sim4reader::readPolish -- can't allocate comment\n'%s'\n%s\n", l->s, strerror(errno));
      exit(1);
    }
    strcpy(p->comment, l->s + 8);

    //fprintf(stdout, "sim4reader: Got comment='%s'\n", p->comment);

    readLine(F, l);
  }

  p->estDefLine = 0L;
  if (strncmp(l->s, "edef", 4) == 0) {
    errno = 0;
    p->estDefLine = (char *)malloc(sizeof(char) * l->l - 4);
    if (errno) {
      fprintf(stderr, "malloc() error in sim4reader::readPolish -- can't allocate estDefLine\n'%s'\n%s\n", l->s, strerror(errno));
      exit(1);
    }
    strcpy(p->estDefLine, l->s + 5);

    //fprintf(stdout, "sim4reader: Got edef='%s'\n", p->estDefLine);

    readLine(F, l);
  }

  p->genDefLine = 0L;
  if (strncmp(l->s, "ddef", 4) == 0) {
    errno = 0;
    p->genDefLine = (char *)malloc(sizeof(char) * l->l - 4);
    if (errno) {
      fprintf(stderr, "malloc() error in sim4reader::readPolish -- can't allocate genDefLine\n'%s'\n%s\n", l->s, strerror(errno));
      exit(1);
    }
    strcpy(p->genDefLine, l->s + 5);

    //fprintf(stdout, "sim4reader: Got ddef='%s'\n", p->genDefLine);

    readLine(F, l);
  }



  //
  //  While we get exons, make exons.
  //
  errno = 0;
  ex = (sim4polishExon *)malloc(sizeof(sim4polishExon) * em);
  if (errno) {
    fprintf(stderr, "malloc() error in sim4reader::readPolish -- can't allocate sim4polishExon\n%s\n", strerror(errno));
    exit(1);
  }

  while (sscanf(l->s, "%d-%d (%d-%d) <%d-%d-%d>",
                &ef, &et, &gf, &gt, &nm, &nn, &id) == 7) {
    if (el >= em) {
      errno = 0;
      em *= 2;
      ex  = (sim4polishExon *)realloc(ex, sizeof(sim4polishExon) * em);
      if (errno) {
        fprintf(stderr, "realloc() error in sim4reader::readPolish -- can't allocate more sim4polishExon\n%s\n", strerror(errno));
        exit(1);
      }
    }

    ex[el].estFrom           = ef;
    ex[el].estTo             = et;
    ex[el].genFrom           = gf;
    ex[el].genTo             = gt;
    ex[el].numMatches        = nm;
    ex[el].numMatchesN       = nn;
    ex[el].percentIdentity   = id;
    ex[el].intronOrientation = '.';
    ex[el].estAlignment      = 0L;
    ex[el].genAlignment      = 0L;

    if ((l->s[l->l-2] == '-') && (l->s[l->l-1] == '>'))
      ex[el].intronOrientation = '>';
    if ((l->s[l->l-2] == '<') && (l->s[l->l-1] == '-'))
      ex[el].intronOrientation = '<';
    if ((l->s[l->l-2] == '-') && (l->s[l->l-1] == '-'))
      ex[el].intronOrientation = '-';
    if ((l->s[l->l-2] == '=') && (l->s[l->l-1] == '='))
      ex[el].intronOrientation = '=';

    el++;

    readLine(F, l);
  }

  if (el == 0) {
    fprintf(stderr, "WARNING: sim4reader::readPolish -- found ZERO exons?\n");
  }

  //  All done.  Save the exons to the sim4polish.
  //
  errno = 0;
  p->numExons = el;
  p->exons    = (sim4polishExon *)malloc(sizeof(sim4polishExon) * el);
  if (errno) {
    fprintf(stderr, "malloc() error in sim4reader::readPolish -- can't allocate %d exons\n%s\n", el, strerror(errno));
    exit(1);
  }
  memcpy(p->exons, ex, sizeof(sim4polishExon) * el);

  free(ex);

  //  Now, if we are not at 'sim4end', assume that there
  //  are alignment lines for each exon.
  //
  if (strcmp(l->s, "sim4end") != 0) {
    for (el=0; el<p->numExons; el++) {
      errno = 0;
      p->exons[el].estAlignment = (char *)malloc(sizeof(char) * (l->l + 1));
      if (errno) {
        fprintf(stderr, "malloc() error in sim4reader::readPolish -- can't allocate estAlignment\n%s\n", strerror(errno));
        exit(1);
      }
      strcpy(p->exons[el].estAlignment, l->s);
      readLine(F, l);

      errno = 0;
      p->exons[el].genAlignment = (char *)malloc(sizeof(char) * (l->l + 1));
      if (errno) {
        fprintf(stderr, "malloc() error in sim4reader::readPolish -- can't allocate genAlignment\n%s\n", strerror(errno));
        exit(1);
      }
      strcpy(p->exons[el].genAlignment, l->s);
      readLine(F, l);

      //  XXX: Sanity.  If we get a sim4end (and we have not just read
      //  the last exon's alignment), we need to stop right now!  To
      //  make the match consistent, we remove all alignment lines.
      //
      if ((el+1 < p->numExons) && (strcmp(l->s, "sim4end") == 0)) {
        int i;
        for (i=0; i<el; i++) {
          free(p->exons[i].estAlignment);
          free(p->exons[i].genAlignment);
          p->exons[i].estAlignment = 0L;
          p->exons[i].genAlignment = 0L;
        }
        el = p->numExons;

        fprintf(stderr, "Full alignment not found for EST=%d GEN=%d; all alignment removed.\n",
                p->estID, p->genID);
      }
    }
  }

  deleteLine(l);

  return(p);
}












void
readString(FILE *F, _line *L) {

  L->l = 0;

  if (fgets(L->s, L->a, F)) {
    L->l = strlen(L->s);

    while ((L->l + 1 >= L->a) &&
           (L->s[L->l - 1] != '\n')) {
      L->a *= 2;
      errno = 0;
      L->s  = (char *)realloc(L->s, sizeof(char) * L->a);
      if (errno) {
        fprintf(stderr, "realloc() error in sim4reader::readString -- can't allocate more char\n%s\n", strerror(errno));
        exit(1);
      }
      fgets(L->s + L->l, L->a - L->l, F);
      L->l += strlen(L->s + L->l);
    }

    //  Trim newlines from the end of the string
    //
    while ((L->l > 0) &&
           ((L->s[L->l - 1] == '\r') ||
            (L->s[L->l - 1] == '\n')))
      L->s[--L->l] = 0;
  }
}



//
//  DESTRUCTIVE to the string s!
//

sim4polish*
s4p_stringToPolish(char *s) {
  sim4polish      *p = 0L;
  int              r;
  int              ef, et, gf, gt, nm, nn, id;
  int              el = 0;
  int              em = 256;
  sim4polishExon  *ex = 0L;
  char             mOri[65];
  char             sOri[65];

  int      maxLine = 0;
  int      curLine = 0;
  char   **lines   = 0L;
  int     *lengths = 0L;

  if ((s == 0L) || (s[0] == 0))
    return(0L);

  errno = 0;
  ex = (sim4polishExon *)malloc(sizeof(sim4polishExon) * em);
  if (errno) {
    fprintf(stderr, "malloc() error in sim4reader::s4p_stringToPolish -- can't allocate sim4polishExon\n%s\n", strerror(errno));
    exit(1);
  }


  //  Convert s into an array of lines
  //
  //  1)  Count the number of lines in the string
  //
  maxLine = 0;
  for (r=0; s[r] != 0; r++)
    if (s[r] == '\n')
      maxLine++;

#if 0
  fprintf(stderr, "----------------------------------------\n");
  fprintf(stderr, "%s", s);
  fprintf(stderr, "----------------------------------------\n");
  fprintf(stderr, "Found %d lines\n", maxLine);
#endif

  if (maxLine == 0)
    return(0L);

  //  2)  Allocate space for the lines
  //
  errno = 0;
  lines = (char **)malloc(sizeof(char *) * maxLine);
  if (errno) {
    fprintf(stderr, "malloc() error in sim4reader::s4p_stringToPolish -- can't allocate space for line pointers\n%s\n", strerror(errno));
    exit(1);
  }
  lengths = (int *)malloc(sizeof(int) * maxLine);
  if (errno) {
    fprintf(stderr, "malloc() error in sim4reader::s4p_stringToPolish -- can't allocate space for line lengths\n%s\n", strerror(errno));
    exit(1);
  }

  //  3)  Set pointers to the character after each \n, change \n's to 0's.
  //
  curLine = 0;

  lines[curLine++] = s;

  for (r=0; s[r] != 0; r++) {
    if (s[r] == '\n') {
      s[r] = 0;

      //  We don't want to create a line for the last '\n'
      //
      if (curLine < maxLine) {
        lines[curLine]   = s + r + 1;
        curLine++;
      }
    }
  }

  for (r=0; r<maxLine; r++)
    lengths[r] = strlen(lines[r]);

  curLine = 0;



  ///////////////////////////////////////
  //
  //  Begin of normal readPolish() code
  //
  ///////////////////////////////////////




  while((curLine < maxLine) && strcmp(lines[curLine], "sim4begin")) {
    fprintf(stderr, "sim4reader: Got '%s', expecting 'sim4begin'\n", lines[curLine]);
    curLine++;
  }

  if (curLine >= maxLine)
    return(0L);

  errno = 0;
  p = (sim4polish *)malloc(sizeof(sim4polish));
  if (errno) {
    fprintf(stderr, "malloc() error in sim4reader::readPolish -- can't allocate sim4polish\n%s\n", strerror(errno));
    exit(1);
  }

  //  Read the description line
  //
  curLine++;

  //  Convert '-' into ' ', on the assumption that this is the description
  //  line.  This allows us to use scanf properly.
  //
  for (r=0; lines[curLine][r] != 0; r++)
    if (lines[curLine][r] == '-')
      lines[curLine][r] = ' ';

  mOri[0] = 0;
  sOri[0] = 0;
  r = sscanf(lines[curLine], "%d[%d %d %d] %d[%d %d] <%d %d %d %s %s>",
             &p->estID,
             &p->estLen,
             &p->estPolyA,
             &p->estPolyT,
             &p->genID,
             &p->genLo,
             &p->genHi,
             &p->numMatches,
             &p->numMatchesN,
             &p->percentIdentity,
             mOri, sOri);
  if (r != 12) {
    fprintf(stderr, "sim4reader: Got '%s'\n", lines[curLine]);
    fprintf(stderr, "sim4reader: Expecting description line (%d).\n", r);
    fprintf(stderr, "sim4reader: mOri='%s'\n", mOri);
    fprintf(stderr, "sim4reader: sOri='%s'\n", sOri);
  }

  switch (mOri[0]) {
    case 'f':
      p->matchOrientation = MATCH_FORWARD;
      break;
    case 'c':
      p->matchOrientation = MATCH_COMPLEMENT;
      break;
    case 'r':
      //  BUG FIX -- old version of sim4 used "reverse-intractable"
      //  instead of "complement-intractable"
      p->matchOrientation = MATCH_COMPLEMENT;
      break;
    default:
      fprintf(stderr, "sim4reader: Got '%s'\n", lines[curLine]);
      fprintf(stderr, "sim4reader: Unknown match orientation.\n");
      break;
  }

  switch (sOri[2]) {
    case 'r':
      p->strandOrientation = STRAND_POSITIVE;
      break;
    case 'v':
      p->strandOrientation = STRAND_NEGATIVE;
      break;
    case 'k':
      p->strandOrientation = STRAND_UNKNOWN;
      break;
    case 't':
      p->strandOrientation = STRAND_INTRACTABLE;
      break;
    case 'i':
      p->strandOrientation = STRAND_FAILED;
      break;
    default:
      fprintf(stderr, "sim4reader: Got '%s'\n", lines[curLine]);
      fprintf(stderr, "sim4reader: Unknown strand orientation.\n");
      break;
  }

  p->querySeqIdentity = (int)floor(100 * (double)(p->numMatches) / (double)(p->estLen - p->estPolyA - p->estPolyT));

  curLine++;

  p->comment = 0L;
  if (strncmp(lines[curLine], "comment", 7) == 0) {
    errno = 0;
    p->comment = (char *)malloc(sizeof(char) * lengths[curLine] - 7);
    if (errno) {
      fprintf(stderr, "malloc() error in sim4reader::readPolish -- can't allocate comment\n'%s'\n%s\n", lines[curLine], strerror(errno));
      exit(1);
    }
    strcpy(p->comment, lines[curLine] + 8);

    //fprintf(stdout, "sim4reader: Got comment='%s'\n", p->comment);

    curLine++;
  }

  p->estDefLine = 0L;
  if (strncmp(lines[curLine], "edef", 4) == 0) {
    errno = 0;
    p->estDefLine = (char *)malloc(sizeof(char) * lengths[curLine] - 4);
    if (errno) {
      fprintf(stderr, "malloc() error in sim4reader::readPolish -- can't allocate estDefLine\n'%s'\n%s\n", lines[curLine], strerror(errno));
      exit(1);
    }
    strcpy(p->estDefLine, lines[curLine] + 5);

    //fprintf(stdout, "sim4reader: Got edef='%s'\n", p->estDefLine);

    curLine++;
  }

  p->genDefLine = 0L;
  if (strncmp(lines[curLine], "ddef", 4) == 0) {
    errno = 0;
    p->genDefLine = (char *)malloc(sizeof(char) * lengths[curLine] - 4);
    if (errno) {
      fprintf(stderr, "malloc() error in sim4reader::readPolish -- can't allocate genDefLine\n'%s'\n%s\n", lines[curLine], strerror(errno));
      exit(1);
    }
    strcpy(p->genDefLine, lines[curLine] + 5);

    //fprintf(stdout, "sim4reader: Got ddef='%s'\n", p->genDefLine);

    curLine++;
  }



  //
  //  While we get exons, make exons.
  //

  while (sscanf(lines[curLine], "%d-%d (%d-%d) <%d-%d-%d>",
                &ef, &et, &gf, &gt, &nm, &nn, &id) == 7) {
    if (el >= em) {
      em *= 2;
      errno = 0;
      ex  = (sim4polishExon *)realloc(ex, sizeof(sim4polishExon) * em);
      if (errno) {
        fprintf(stderr, "realloc() error in sim4reader::readPolish -- can't allocate more sim4polishExon\n%s\n", strerror(errno));
        exit(1);
      }
    }

    ex[el].estFrom           = ef;
    ex[el].estTo             = et;
    ex[el].genFrom           = gf;
    ex[el].genTo             = gt;
    ex[el].numMatches        = nm;
    ex[el].numMatchesN       = nn;
    ex[el].percentIdentity   = id;
    ex[el].intronOrientation = '.';
    ex[el].estAlignment      = 0L;
    ex[el].genAlignment      = 0L;
 
    if ((lines[curLine][lengths[curLine]-2] == '-') && (lines[curLine][lengths[curLine]-1] == '>'))
      ex[el].intronOrientation = '>';
    if ((lines[curLine][lengths[curLine]-2] == '<') && (lines[curLine][lengths[curLine]-1] == '-'))
      ex[el].intronOrientation = '<';
    if ((lines[curLine][lengths[curLine]-2] == '-') && (lines[curLine][lengths[curLine]-1] == '-'))
      ex[el].intronOrientation = '-';
    if ((lines[curLine][lengths[curLine]-2] == '=') && (lines[curLine][lengths[curLine]-1] == '='))
      ex[el].intronOrientation = '=';

    el++;

    curLine++;
  }

  if (el == 0) {
    fprintf(stderr, "WARNING: sim4reader::readPolish -- found ZERO exons?\n");
  }

  //  All done.  Save the exons to the sim4polish.
  //
  errno = 0;
  p->numExons = el;
  p->exons    = (sim4polishExon *)malloc(sizeof(sim4polishExon) * el);
  if (errno) {
    fprintf(stderr, "malloc() error in sim4reader::readPolish -- can't allocate %d exons\n%s\n", el, strerror(errno));
    exit(1);
  }
  memcpy(p->exons, ex, sizeof(sim4polishExon) * el);

  free(ex);

  //  Now, if we are not at 'sim4end', assume that there
  //  are alignment lines for each exon.
  //
  if (strcmp(lines[curLine], "sim4end") != 0) {
    for (el=0; el<p->numExons; el++) {
      errno = 0;
      p->exons[el].estAlignment = (char *)malloc(sizeof(char) * (lengths[curLine] + 1));
      if (errno) {
        fprintf(stderr, "malloc() error in sim4reader::readPolish -- can't allocate estAlignment\n%s\n", strerror(errno));
        exit(1);
      }
      strcpy(p->exons[el].estAlignment, lines[curLine]);
      curLine++;

      errno = 0;
      p->exons[el].genAlignment = (char *)malloc(sizeof(char) * (lengths[curLine] + 1));
      if (errno) {
        fprintf(stderr, "malloc() error in sim4reader::readPolish -- can't allocate genAlignment\n%s\n", strerror(errno));
        exit(1);
      }
      strcpy(p->exons[el].genAlignment, lines[curLine]);
      curLine++;

      //  XXX: Sanity.  If we get a sim4end (and we have not just read
      //  the last exon's alignment), we need to stop right now!  To
      //  make the match consistent, we remove all alignment lines.
      //
      if ((el+1 < p->numExons) && (strcmp(lines[curLine], "sim4end") == 0)) {
        int i;
        for (i=0; i<el; i++) {
          free(p->exons[i].estAlignment);
          free(p->exons[i].genAlignment);
          p->exons[i].estAlignment = 0L;
          p->exons[i].genAlignment = 0L;
        }
        el = p->numExons;

        fprintf(stderr, "Full alignment not found for EST=%d GEN=%d; all alignment removed.\n",
                p->estID, p->genID);
      }
    }
  }


  ///////////////////////////////////////
  //
  //  Begin of normal readPolish() code
  //
  ///////////////////////////////////////

  free(lines);

  return(p);
}























void
destroyPolish(sim4polish *p) {
  int i;

  if (p) {
    for (i=0; i<p->numExons; i++) {
      free(p->exons[i].estAlignment);
      free(p->exons[i].genAlignment);
    }
    free(p->exons);
    free(p->comment);
    free(p->estDefLine);
    free(p->genDefLine);
    free(p);
  }
}


void
printPolish(FILE *O, sim4polish *p) {
  int          i;
  char const  *mOri;
  char const  *sOri;

  if (p) {
    fprintf(O, "sim4begin\n");

    switch (p->matchOrientation) {
      case MATCH_FORWARD:
        mOri = "forward";
        break;
      case MATCH_COMPLEMENT:
        mOri = "complement";
        break;
      default:
        fprintf(stderr, "sim4reader: Unknown matchOrientation '%d' in printPolish()\n", p->matchOrientation);
        mOri = "UNKNOWN";
        break;
    }

    switch (p->strandOrientation) {
      case STRAND_POSITIVE:
        sOri = "forward";
        break;
      case STRAND_NEGATIVE:
        sOri = "reverse";
        break;
      case STRAND_UNKNOWN:
        sOri = "unknown";
        break;
      case STRAND_INTRACTABLE:
        sOri = "intractable";
        break;
      case STRAND_FAILED:
        sOri = "aborted";
        break;
      default:
        fprintf(stderr, "sim4reader: Unknown strandOrientation '%d' in printPolish()\n", p->matchOrientation);
        sOri = "UNKNOWN";
        break;
    }

    fprintf(O, "%d[%d-%d-%d] %d[%d-%d] <%d-%d-%d-%s-%s>\n",
            p->estID,
            p->estLen,
            p->estPolyA,
            p->estPolyT,
            p->genID,
            p->genLo,
            p->genHi,
            p->numMatches,
            p->numMatchesN,
            p->percentIdentity,
            mOri, sOri);

    if (p->comment)
      fprintf(O, "comment=%s\n", p->comment);
    if (p->estDefLine)
      fprintf(O, "edef=%s\n", p->estDefLine);
    if (p->genDefLine)
      fprintf(O, "ddef=%s\n", p->genDefLine);
    for (i=0; i<p->numExons; i++) {
      fprintf(O, "%d-%d (%d-%d) <%d-%d-%d>",
              p->exons[i].estFrom,
              p->exons[i].estTo,
              p->exons[i].genFrom,
              p->exons[i].genTo,
              p->exons[i].numMatches,
              p->exons[i].numMatchesN,
              p->exons[i].percentIdentity);
      switch (p->exons[i].intronOrientation) {
        case INTRON_POSITIVE:
          fprintf(O, " ->\n");
          break;
        case INTRON_NEGATIVE:
          fprintf(O, " <-\n");
          break;
        case INTRON_AMBIGUOUS:
          fprintf(O, " --\n");
          break;
        case INTRON_GAP:
          fprintf(O, " ==\n");
          break;
        case INTRON_NONE:
          fprintf(O, "\n");
          break;
        default:
          fprintf(O, "\n");
          break;
      }
    }
    for (i=0; i<p->numExons; i++) {
      if (p->exons[i].estAlignment)
        fprintf(O, "%s\n", p->exons[i].estAlignment);
      if (p->exons[i].genAlignment)
        fprintf(O, "%s\n", p->exons[i].genAlignment);
    }
    fprintf(O, "sim4end\n");
  }
}

void
printPolishColumn(FILE *O, sim4polish *p) {
  int             i;

  //  XXX:  Original version printed strings for estID and genID.
  //  estID was the first word of the est defline
  //  genID was the ga_uid

  fprintf(O, "%d %d %c %d %d %d",
          p->estID,
          p->genID, 
          p->matchOrientation,
          p->percentIdentity,
          p->numMatches + p->numMatchesN,
          p->numExons);

  for (i=0; i<p->numExons; i++)
    fprintf(O, " %d %d %d %d %d %c",
            p->exons[i].estFrom,
            p->exons[i].estTo,
            p->exons[i].genFrom,
            p->exons[i].genTo,
            p->exons[i].percentIdentity,
            p->exons[i].intronOrientation);

  fprintf(O,"\n");
}

void
printPolishNormalized(FILE *O, sim4polish *p) {
  int  i;
  int  genLo = p->genLo;

  p->genLo = 0;

  for (i=0; i<p->numExons; i++) {
    p->exons[i].genFrom += genLo;
    p->exons[i].genTo   += genLo;
  }

  printPolish(O, p);

  for (i=0; i<p->numExons; i++) {
    p->exons[i].genFrom -= genLo;
    p->exons[i].genTo   -= genLo;
  }

  p->genLo = genLo;
}


void *
copyPolishDupHelper(void *orig, size_t size) {
  void *rslt;
  errno = 0;
  rslt = (char *)malloc(size);
  if (errno) {
    fprintf(stderr, "Out of memory in copyPolish.\n%s\n", strerror(errno));
    exit(1);
  }
  memcpy(rslt, orig, size);
  return(rslt);
}

sim4polish *
copyPolish(sim4polish *orig) {
  int         i;
  sim4polish *copy;

  errno = 0;
  copy = (sim4polish *)malloc(sizeof(sim4polish));
  if (errno) {
    fprintf(stderr, "malloc() error in sim4reader::copyPolish\n%s\n", strerror(errno));
    exit(1);
  }

  memcpy(copy, orig, sizeof(sim4polish));

  //  Deflines are optional -- if they exist, make copies.
  //
  if (orig->estDefLine)
    copy->estDefLine = (char *)copyPolishDupHelper(orig->estDefLine, sizeof(char) * (strlen(orig->estDefLine) + 1));

  if (orig->genDefLine)
    copy->genDefLine = (char *)copyPolishDupHelper(orig->genDefLine, sizeof(char) * (strlen(orig->genDefLine) + 1));

  //  Exons are not optional.  Copy them.
  //
  copy->exons = (sim4polishExon *)copyPolishDupHelper(orig->exons, sizeof(sim4polishExon) * copy->numExons);

  //  Finally, copy all the exon alignments
  //
  for (i=0; i<copy->numExons; i++) {
    if (copy->exons[i].estAlignment)
      copy->exons[i].estAlignment = (char *)copyPolishDupHelper(orig->exons[i].estAlignment, sizeof(char) * (strlen(copy->exons[i].estAlignment) + 1));

    if (copy->exons[i].genAlignment)
      copy->exons[i].genAlignment = (char *)copyPolishDupHelper(orig->exons[i].genAlignment, sizeof(char) * (strlen(copy->exons[i].genAlignment) + 1));
  }

  return(copy);
}



void
s4p_swapExons(sim4polish *p, int a, int b) {
  sim4polishExon  copyofa;

  //  Save the a exon into copyofa; copy b into a; copy copyofa into b
  //
  memcpy(&copyofa,   p->exons+a, sizeof(sim4polishExon));
  memcpy(p->exons+a, p->exons+b, sizeof(sim4polishExon));
  memcpy(p->exons+b, &copyofa,   sizeof(sim4polishExon));
}


void
s4p_deleteExon(sim4polish *p, int a) {
  char  *ed, *gd;
  int    i;
  int    editDistance    = 0;
  int    alignmentLength = 0;

  //  Warn if we don't have alignments -- this is now done by the
  //  driver (e.g., cleanPolishes.C)
  //
#if 0
  if ((p->exons[0].estAlignment == 0L) || (p->exons[0].genAlignment == 0L))
    fprintf(stderr, "s4p_deleteExon()-- Need alignments to recompute scores correctly!\n");
#endif

  //  Set the intron orientation for the exon before the one we are
  //  deleting:
  //    If we are deleting the first exon, there is no previous exon
  //    If we are deleting the last exon, set the previous to INTRON_NONE
  //    Otherwise, set the previous to INTRON_GAP
  //
  if (a == p->numExons-1) {
    p->exons[a-1].intronOrientation = INTRON_NONE;
  } else {
    if (a > 0)
      p->exons[a-1].intronOrientation = INTRON_GAP;
  }

  //  Update the match scores
  //
  p->numMatches  -= p->exons[a].numMatches;
  p->numMatchesN -= p->exons[a].numMatchesN;

  //  Delete any alignment in the soon to be deleted exon
  //
  free(p->exons[a].estAlignment);
  free(p->exons[a].genAlignment);

  //  Shift all the exons down by one, and decrement the number of
  //  exons present in the list.
  //
  for (i=a+1; i<p->numExons; i++)
    memcpy(p->exons+i-1, p->exons+i, sizeof(sim4polishExon));

  p->numExons--;

  //  To be safe, zero out the space for the (still allocated, but
  //  unused) last exon
  //
  memset(p->exons+p->numExons, 0, sizeof(sim4polishExon));

  //  The strand orientation doesn't change if we are deleting the
  //  first or last exon....unless we now have only one exon.
  //
  if (((0 < a) && (a < p->numExons-1)) ||
      (p->numExons == 1))
    p->strandOrientation = STRAND_UNKNOWN;



  //  Compute the alignment length and the number of edits.
  //
  alignmentLength = 0;
  editDistance    = 0;

  for (i=0; i<p->numExons; i++) {
    ed = p->exons[i].estAlignment;
    gd = p->exons[i].genAlignment;

    if (ed && gd) {
      alignmentLength += 2 * strlen(ed);
      for (; *ed && *gd; ed++, gd++) {
        if (*ed != *gd)
          editDistance++;
      }
    } else {
      int len = p->exons[i].estTo - p->exons[i].estFrom + 1 + p->exons[i].estTo - p->exons[i].estFrom + 1;

      alignmentLength += len;
      editDistance    += len / 2 - p->exons[i].numMatches - p->exons[i].numMatchesN;
    }
  }

#if 0
  fprintf(stdout, "Found (new)alignLen = %d\n", alignmentLength);
  fprintf(stdout, "Found (new)editDist = %d\n", editDistance);
#endif

  //  Fix the scores for the match.  Special case; if there is only
  //  one exon left, the score for the exon is the score for the
  //  match.
  //
  if (p->numExons == 1)
    p->percentIdentity = p->exons[0].percentIdentity;
  else
    if (alignmentLength > 0)
      p->percentIdentity = (int)floor(100 * (1 - 2 * (double)editDistance / (double)alignmentLength));
    else
      p->percentIdentity = 0;

  //  Update the query sequence identity
  //
  p->querySeqIdentity = (int)floor(100 * (double)(p->numMatches) / (double)(p->estLen - p->estPolyA - p->estPolyT));
}


// p[a] --> e
void
s4p_copyExon(sim4polish *p, int a, sim4polishExon *e) {
  memcpy(e, p->exons+a, sizeof(sim4polishExon));
}


// e --> p[a]
void
s4p_overwriteExon(sim4polish *p, sim4polishExon *e, int a) {
  memcpy(p->exons+a, e, sizeof(sim4polishExon));
}


void
s4p_insertExon(sim4polish *p, int a, sim4polishExon *e) {
  fprintf(stderr, "s4p_insertExon() not implemented!\n");
}
