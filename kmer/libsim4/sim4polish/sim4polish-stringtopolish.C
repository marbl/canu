#include "sim4polish.H"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>

#include <assert.h>

void
sim4polish::s4p_linesToPolish(u32bit lineNumber, u32bit maxLines, char **lines, u32bit *lengths) {
  char             mOri[65];
  char             sOri[65];

  assert(_comment    == 0L);
  assert(_estDefLine == 0L);
  assert(_genDefLine == 0L);
  assert(_exons      == 0L);
  assert(_numExons   == 0);

  if (strcmp(lines[0], "sim4begin")) {
    fprintf(stderr, "sim4polish::s4p_linesToPolish()-- Invalid sim4db format.  Cannot convert.\n");
    return;
  }

  u32bit cl = 1;

  //  Convert '-' into ' ', on the assumption that this is the description line.  This allows us to
  //  use scanf properly.
  //
  for (u32bit i=0; i<lengths[cl]; i++)
    if (lines[cl][i] == '-')
      lines[cl][i] = ' ';

  mOri[0] = 0;
  sOri[0] = 0;
  u32bit r = sscanf(lines[cl], ""u32bitFMT"["u32bitFMT" "u32bitFMT" "u32bitFMT"] "u32bitFMT"["u32bitFMT" "u32bitFMT"] <"u32bitFMT" "u32bitFMT" "u32bitFMT" %s %s>",
                    &_estID,
                    &_estLen,
                    &_estPolyA,
                    &_estPolyT,
                    &_genID,
                    &_genRegionOffset,
                    &_genRegionLength,
                    &_numMatches,
                    &_numMatchesN,
                    &_percentIdentity,
                    mOri, sOri);
  if (r != 12) {
    fprintf(stderr, "sim4polish::s4p_linesToPolish()--  line "u32bitFMT": '%s'\n", lineNumber, lines[cl]);
    fprintf(stderr, "sim4polish::s4p_linesToPolish()--  Expecting description line, found %d tokens instead of 12.\n", r);
  }

  switch (mOri[0]) {
    case 'f':
      _matchOrientation = SIM4_MATCH_FORWARD;
      break;
    case 'c':
      _matchOrientation = SIM4_MATCH_COMPLEMENT;
      break;
    case 'r':
      //  BUG FIX -- old version of sim4 used "reverse-intractable"
      //  instead of "complement-intractable"
      _matchOrientation = SIM4_MATCH_COMPLEMENT;
      break;
    default:
      fprintf(stderr, "sim4polish::s4p_linesToPolish()--  line "u32bitFMT": '%s'\n", lineNumber, lines[cl]);
      fprintf(stderr, "sim4polish::s4p_linesToPolish()--  unknown match orientation\n");
      break;
  }

  switch (sOri[2]) {
    case 'r':
      _strandOrientation = SIM4_STRAND_POSITIVE;
      break;
    case 'v':
      _strandOrientation = SIM4_STRAND_NEGATIVE;
      break;
    case 'k':
      _strandOrientation = SIM4_STRAND_UNKNOWN;
      break;
    case 't':
      _strandOrientation = SIM4_STRAND_INTRACTABLE;
      break;
    case 'i':
      _strandOrientation = SIM4_STRAND_FAILED;
      break;
    default:
      fprintf(stderr, "sim4polish::s4p_linesToPolish()--  line "u32bitFMT": '%s'\n", lineNumber, lines[cl]);
      fprintf(stderr, "sim4polish::s4p_linesToPolish()--  unknown strand orientation\n");
      break;
  }

  cl++;

  _comment = 0L;
  if (strncmp(lines[cl], "comment", 7) == 0) {
    _comment = new char [lengths[cl] - 7];
    strcpy(_comment, lines[cl] + 8);

    cl++;
  }

  _estDefLine = 0L;
  if (strncmp(lines[cl], "edef", 4) == 0) {
    _estDefLine = new char [lengths[cl] - 4];
    strcpy(_estDefLine, lines[cl] + 5);

    cl++;
  }

  _genDefLine = 0L;
  if (strncmp(lines[cl], "ddef", 4) == 0) {
    _genDefLine = new char [lengths[cl] - 4];
    strcpy(_genDefLine, lines[cl] + 5);

    cl++;
  }

  //
  //  While we get exons, make exons.
  //

  sim4polishExon    exon;
  u32bit            maxExons = 1024;

  _numExons = 0;
  _exons    = new sim4polishExon [maxExons];

  _numCovered = 0;

  while (sscanf(lines[cl], ""u32bitFMT"-"u32bitFMT" ("u32bitFMT"-"u32bitFMT") <"u32bitFMT"-"u32bitFMT"-"u32bitFMT">",
                &exon._estFrom, &exon._estTo,
                &exon._genFrom, &exon._genTo,
                &exon._numMatches,
                &exon._numMatchesN,
                &exon._percentIdentity) == 7) {

    //  Dang, out of space!  This would be a chore, except we don't have alignments yet, and so can
    //  get by with a shallow copy.
    //
    if (_numExons >= maxExons) {
      maxExons *= 2;
      sim4polishExon *nnn = new sim4polishExon [maxExons];
      memcpy(nnn, _exons, sizeof(sim4polishExon) * _numExons);
      delete [] _exons;
      _exons = nnn;
    }

    _exons[_numExons] = exon;

    _exons[_numExons]._intronOrientation = SIM4_INTRON_NONE;

    if ((lines[cl][lengths[cl]-2] == '-') && (lines[cl][lengths[cl]-1] == '>'))
      _exons[_numExons]._intronOrientation = SIM4_INTRON_POSITIVE;
    if ((lines[cl][lengths[cl]-2] == '<') && (lines[cl][lengths[cl]-1] == '-'))
      _exons[_numExons]._intronOrientation = SIM4_INTRON_NEGATIVE;
    if ((lines[cl][lengths[cl]-2] == '-') && (lines[cl][lengths[cl]-1] == '-'))
      _exons[_numExons]._intronOrientation = SIM4_INTRON_AMBIGUOUS;
    if ((lines[cl][lengths[cl]-2] == '=') && (lines[cl][lengths[cl]-1] == '='))
      _exons[_numExons]._intronOrientation = SIM4_INTRON_GAP;

    _exons[_numExons]._estAlignment = 0L;
    _exons[_numExons]._genAlignment = 0L;

    _numCovered += _exons[_numExons]._estTo - _exons[_numExons]._estFrom + 1;

    _numExons++;

    cl++;
  }

  if (_numExons == 0) {
    fprintf(stderr, "sim4polish::s4p_linesToPolish()--  line "u32bitFMT": '%s'\n", lineNumber, lines[cl]);
    fprintf(stderr, "sim4polish::s4p_linesToPolish()--  WARNING: found ZERO exons?\n");
  }

  _querySeqIdentity = s4p_percentCoverageApprox();

  //  Now, if we are not at 'sim4end', assume that there are alignment lines for each exon.
  //
  //  We used to check that we didn't hit 'sim4end' before reading all the alignment lines, and if
  //  we did, we'd compain about it and remove the alignment lines.  Too much work.
  //
  if (strcmp(lines[cl], "sim4end") != 0) {
    for (u32bit el=0; el<_numExons; el++) {
      _exons[el]._estAlignment = new char [lengths[cl] + 1];
      strcpy(_exons[el]._estAlignment, lines[cl]);
      cl++;

      _exons[el]._genAlignment = new char [lengths[cl] + 1];
      strcpy(_exons[el]._genAlignment, lines[cl]);
      cl++;
    }
  }
}



//  Convert string str to an array of lines, call function to parse those lines.
//
//  NOTE that this is DESTRUCTIVE to the string!
void
sim4polish::s4p_stringToPolish(char *s) {

  //  Clear this polish.

  _numExons = 0;

  delete [] _comment;     _comment = 0L;
  delete [] _estDefLine;  _estDefLine = 0L;
  delete [] _genDefLine;  _genDefLine = 0L;
  delete [] _exons;       _exons = 0L;

  //  Decide the type of record we're reading.

  //  Read it.

  if ((s == 0L) || (s[0] == 0))
    return;

  //  Count the number of lines in the string

  u32bit  numLines = 0;
  u32bit  curLine  = 0;

  for (u32bit r=0; s[r] != 0; r++)
    if (s[r] == '\n')
      numLines++;

  if (numLines == 0)
    return;

  //  Allocate space for the lines

  char   **lines   = new char * [numLines + 1];
  u32bit  *lengths = new u32bit [numLines + 1];

  //  Set pointers to the character after each \n, change \n's to 0's.  Even though we set it, the
  //  last \n doesn't contribute a line.  The code is simpler if we treat it as a line though.

  lines[curLine++] = s;

  for (u32bit r=0; s[r] != 0; r++) {
    if (s[r] == '\n') {
      s[r] = 0;
      lines[curLine++] = s + r + 1;
    }
  }

  assert(curLine == numLines + 1);

  for (u32bit r=0; r<numLines; r++)
    lengths[r] = strlen(lines[r]);

  s4p_linesToPolish(0, numLines, lines, lengths);

  delete [] lines;
  delete [] lengths;
}
