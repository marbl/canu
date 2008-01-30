#include "sim4polish.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <assert.h>


//  Utility for reading a whole line, safely, from a file.
//
typedef struct {
  u32bit   l;
  u32bit   a;
  char    *s;
  u32bit   lineNumber;
} _line;

_line *newLine(void) {
  _line *L;
  errno = 0;
  L = (_line *)malloc(sizeof(_line));
  if (errno) {
    fprintf(stderr, "malloc() error in sim4reader::newLine\n%s\n", strerror(errno));
    assert(0), exit(1);
  }
  L->l = 0;
  L->a = 10;
  L->lineNumber = 0;
  errno = 0;
  L->s = (char *)malloc(sizeof(char) * L->a);
  if (errno) {
    fprintf(stderr, "malloc() error in sim4reader::newLine\n%s\n", strerror(errno));
    assert(0), exit(1);
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
        fprintf(stderr, "realloc() error in sim4reader::readLine("u32bitFMT") -- can't allocate more char\n%s\n", L->lineNumber, strerror(errno));
        assert(0), exit(1);
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
s4p_readPolish(FILE *F) {
  sim4polish      *p = 0L;
  _line           *l = newLine();
  int              r;
  u32bit           ef, et, gf, gt, nm, nn, id;
  int              el = 0;
  int              em = 256;
  sim4polishExon  *ex = 0L;
  char             mOri[65];
  char             sOri[65];

  readLine(F, l);
  while(!feof(F) && strcmp(l->s, "sim4begin")) {
    fprintf(stderr, "sim4reader: Got '%s', expecting 'sim4begin' at line "u32bitFMT"\n",
            l->s, l->lineNumber);
    readLine(F, l);
  }

  if (feof(F))
    return(0L);

  errno = 0;
  p = (sim4polish *)malloc(sizeof(sim4polish));
  if (errno) {
    fprintf(stderr, "malloc() error in sim4reader::readPolish("u32bitFMT") -- can't allocate sim4polish\n%s\n", l->lineNumber, strerror(errno));
    assert(0), exit(1);
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
  r = sscanf(l->s, ""u32bitFMT"["u32bitFMT" "u32bitFMT" "u32bitFMT"] "u32bitFMT"["u32bitFMT" "u32bitFMT"] <"u32bitFMT" "u32bitFMT" "u32bitFMT" %s %s>",
             &p->estID,
             &p->estLen,
             &p->estPolyA,
             &p->estPolyT,
             &p->genID,
             &p->genRegionOffset,
             &p->genRegionLength,
             &p->numMatches,
             &p->numMatchesN,
             &p->percentIdentity,
             mOri, sOri);
  if (r != 12) {
    fprintf(stderr, "sim4reader: line "u32bitFMT": '%s'\n", l->lineNumber, l->s);
    fprintf(stderr, "sim4reader: Expecting description line, found %d tokens instead of 12.\n", r);
  }

  switch (mOri[0]) {
    case 'f':
      p->matchOrientation = SIM4_MATCH_FORWARD;
      break;
    case 'c':
      p->matchOrientation = SIM4_MATCH_COMPLEMENT;
      break;
    case 'r':
      //  BUG FIX -- old version of sim4 used "reverse-intractable"
      //  instead of "complement-intractable"
      p->matchOrientation = SIM4_MATCH_COMPLEMENT;
      break;
    default:
      fprintf(stderr, "sim4reader: line "u32bitFMT": '%s'\n", l->lineNumber, l->s);
      fprintf(stderr, "sim4reader: unknown match orientation\n");
      break;
  }

  switch (sOri[2]) {
    case 'r':
      p->strandOrientation = SIM4_STRAND_POSITIVE;
      break;
    case 'v':
      p->strandOrientation = SIM4_STRAND_NEGATIVE;
      break;
    case 'k':
      p->strandOrientation = SIM4_STRAND_UNKNOWN;
      break;
    case 't':
      p->strandOrientation = SIM4_STRAND_INTRACTABLE;
      break;
    case 'i':
      p->strandOrientation = SIM4_STRAND_FAILED;
      break;
    default:
      fprintf(stderr, "sim4reader: line "u32bitFMT": '%s'\n", l->lineNumber, l->s);
      fprintf(stderr, "sim4reader: unknown strand orientation\n");
      break;
  }

  readLine(F, l);

  p->comment = 0L;
  if (strncmp(l->s, "comment", 7) == 0) {
    errno = 0;
    p->comment = (char *)malloc(sizeof(char) * l->l - 7);
    if (errno) {
      fprintf(stderr, "malloc() error in sim4reader::readPolish("u32bitFMT") -- can't allocate comment\n'%s'\n%s\n", l->lineNumber, l->s, strerror(errno));
      assert(0), exit(1);
    }
    strcpy(p->comment, l->s + 8);

    readLine(F, l);
  }

  p->estDefLine = 0L;
  if (strncmp(l->s, "edef", 4) == 0) {
    errno = 0;
    p->estDefLine = (char *)malloc(sizeof(char) * l->l - 4);
    if (errno) {
      fprintf(stderr, "malloc() error in sim4reader::readPolish("u32bitFMT") -- can't allocate estDefLine\n'%s'\n%s\n", l->lineNumber, l->s, strerror(errno));
      assert(0), exit(1);
    }
    strcpy(p->estDefLine, l->s + 5);

    readLine(F, l);
  }

  p->genDefLine = 0L;
  if (strncmp(l->s, "ddef", 4) == 0) {
    errno = 0;
    p->genDefLine = (char *)malloc(sizeof(char) * l->l - 4);
    if (errno) {
      fprintf(stderr, "malloc() error in sim4reader::readPolish("u32bitFMT") -- can't allocate genDefLine\n'%s'\n%s\n", l->lineNumber, l->s, strerror(errno));
      assert(0), exit(1);
    }
    strcpy(p->genDefLine, l->s + 5);

    readLine(F, l);
  }



  //
  //  While we get exons, make exons.
  //
  errno = 0;
  ex = (sim4polishExon *)malloc(sizeof(sim4polishExon) * em);
  if (errno) {
    fprintf(stderr, "malloc() error in sim4reader::readPolish("u32bitFMT") -- can't allocate sim4polishExon\n%s\n", l->lineNumber, strerror(errno));
    assert(0), exit(1);
  }

  p->numCovered = 0;

  while (sscanf(l->s, ""u32bitFMT"-"u32bitFMT" ("u32bitFMT"-"u32bitFMT") <"u32bitFMT"-"u32bitFMT"-"u32bitFMT">",
                &ef, &et, &gf, &gt, &nm, &nn, &id) == 7) {
    if (el >= em) {
      errno = 0;
      em *= 2;
      ex  = (sim4polishExon *)realloc(ex, sizeof(sim4polishExon) * em);
      if (errno) {
        fprintf(stderr, "realloc() error in sim4reader::readPolish("u32bitFMT") -- can't allocate more sim4polishExon\n%s\n", l->lineNumber, strerror(errno));
        assert(0), exit(1);
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

    p->numCovered += ex[el].estTo - ex[el].estFrom + 1;

    el++;

    readLine(F, l);
  }

  if (el == 0) {
    fprintf(stderr, "WARNING: sim4reader::readPolish("u32bitFMT") -- found ZERO exons?\n", l->lineNumber);
  }


  p->querySeqIdentity = s4p_percentCoverageApprox(p);


  //  All done.  Save the exons to the sim4polish.
  //
  errno = 0;
  p->numExons = el;
  p->exons    = (sim4polishExon *)malloc(sizeof(sim4polishExon) * el);
  if (errno) {
    fprintf(stderr, "malloc() error in sim4reader::readPolish("u32bitFMT") -- can't allocate %d exons\n%s\n", l->lineNumber, el, strerror(errno));
    assert(0), exit(1);
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
        fprintf(stderr, "malloc() error in sim4reader::readPolish("u32bitFMT") -- can't allocate estAlignment\n%s\n", l->lineNumber, strerror(errno));
        assert(0), exit(1);
      }
      strcpy(p->exons[el].estAlignment, l->s);
      readLine(F, l);

      errno = 0;
      p->exons[el].genAlignment = (char *)malloc(sizeof(char) * (l->l + 1));
      if (errno) {
        fprintf(stderr, "malloc() error in sim4reader::readPolish("u32bitFMT") -- can't allocate genAlignment\n%s\n", l->lineNumber, strerror(errno));
        assert(0), exit(1);
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

        fprintf(stderr, "WARNING: sim4reader::readPolish("u32bitFMT"): Full alignment not found for EST="u32bitFMT" GEN="u32bitFMT"; all alignment removed.\n",
                l->lineNumber, p->estID, p->genID);
      }
    }
  }

  deleteLine(l);

  return(p);
}




