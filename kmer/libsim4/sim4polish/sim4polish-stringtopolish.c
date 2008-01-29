#include "sim4polish.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>


//
//  DESTRUCTIVE to the string s!
//

sim4polish*
s4p_stringToPolish(char *s) {
  sim4polish      *p = 0L;
  int              r;
  u32bit           ef, et, gf, gt, nm, nn, id;
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
  r = sscanf(lines[curLine], ""u32bitFMT"["u32bitFMT" "u32bitFMT" "u32bitFMT"] "u32bitFMT"["u32bitFMT" "u32bitFMT"] <"u32bitFMT" "u32bitFMT" "u32bitFMT" %s %s>",
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
    fprintf(stderr, "sim4reader: Got '%s'\n", lines[curLine]);
    fprintf(stderr, "sim4reader: Expecting description line (%d).\n", r);
    fprintf(stderr, "sim4reader: mOri='%s'\n", mOri);
    fprintf(stderr, "sim4reader: sOri='%s'\n", sOri);
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
      fprintf(stderr, "sim4reader: Got '%s'\n", lines[curLine]);
      fprintf(stderr, "sim4reader: Unknown match orientation.\n");
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
      fprintf(stderr, "sim4reader: Got '%s'\n", lines[curLine]);
      fprintf(stderr, "sim4reader: Unknown strand orientation.\n");
      break;
  }

  p->querySeqIdentity = s4p_percentCoverageApprox(p);

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

  while (sscanf(lines[curLine], ""u32bitFMT"-"u32bitFMT" ("u32bitFMT"-"u32bitFMT") <"u32bitFMT"-"u32bitFMT"-"u32bitFMT">",
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

        fprintf(stderr, "Full alignment not found for EST="u32bitFMT" GEN="u32bitFMT"; all alignment removed.\n",
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

