#include "sim4polish.H"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>

#include <assert.h>

void
sim4polish::s4p_linesToPolishS4DB(u32bit startPosition,
                                  u32bit maxLines,
                                  char **lines,
                                  u32bit *lengths) {
  char             mOri[65];
  char             sOri[65];

  assert(_comment    == 0L);
  assert(_estDefLine == 0L);
  assert(_genDefLine == 0L);
  assert(_exons      == 0L);
  assert(_numExons   == 0);

  if (strcmp(lines[0], "sim4begin")) {
    fprintf(stderr, "sim4polish::s4p_linesToPolishS4DB()-- Invalid sim4db format, got '%s' instead of sim4begin.  Cannot convert.\n",
            lines[0]);
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
    fprintf(stderr, "sim4polish::s4p_linesToPolishS4DB()--  byte "u32bitFMT": '%s'\n", startPosition, lines[cl]);
    fprintf(stderr, "sim4polish::s4p_linesToPolishS4DB()--  Expecting description line, found %d tokens instead of 12.\n", r);
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
      fprintf(stderr, "sim4polish::s4p_linesToPolishS4DB()--  byte "u32bitFMT": '%s'\n", startPosition, lines[cl]);
      fprintf(stderr, "sim4polish::s4p_linesToPolishS4DB()--  unknown match orientation\n");
      break;
  }

  switch (sOri[2]) {
    case 'r': _strandOrientation = SIM4_STRAND_POSITIVE; break;
    case 'v': _strandOrientation = SIM4_STRAND_NEGATIVE; break;
    case 'k': _strandOrientation = SIM4_STRAND_UNKNOWN; break;
    case 't': _strandOrientation = SIM4_STRAND_INTRACTABLE; break;
    case 'i': _strandOrientation = SIM4_STRAND_FAILED; break;
    default:
      fprintf(stderr, "sim4polish::s4p_linesToPolishS4DB()--  byte "u32bitFMT": '%s'\n", startPosition, lines[cl]);
      fprintf(stderr, "sim4polish::s4p_linesToPolishS4DB()--  unknown strand orientation\n");
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
    fprintf(stderr, "sim4polish::s4p_linesToPolishS4DB()--  byte "u32bitFMT": '%s'\n", startPosition, lines[cl]);
    fprintf(stderr, "sim4polish::s4p_linesToPolishS4DB()--  WARNING: found ZERO exons?\n");
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


//  NOTE: This alters the lines array, with strtok()
void
sim4polish::s4p_linesToPolishGFF3(u32bit startPosition,
                                  u32bit maxLines,
                                  char **lines,
                                  u32bit *lengths) {
  char             mOri;
  char             sOri;
  char            *clptr;
  int              matchID;
  char            *tok, *crttok;
  int              dummy1, dummy2;
  char             dummybuf[1000];

  u32bit           r;
  bool             ok = true;

  assert(_comment    == 0L);
  assert(_estDefLine == 0L);
  assert(_genDefLine == 0L);
  assert(_exons      == 0L);
  assert(_numExons   == 0);

  //  Don't need to store matchID; re-assigned when file changes

  u32bit cl = 0;
  for (cl=0; lines[cl] && (lines[cl][0]=='#'); cl++);
  if (lines[cl] == NULL) {
    fprintf(stderr, "sim4polish::s4p_linesToPolishGFF3()-- Empty record.  Cannot convert (%s).\n", lines[0]);
    return;
  }

  if (!strcmp(lines[0], "\tsim4db\tmRNA")) {
    fprintf(stderr, "sim4polish::s4p_linesToPolishGFF3()-- Invalid GFF3 format, got '%s' instead of GFF3 mRNA line.  Cannot convert.\n",
            lines[0]);
    return;
  }

  cl = 0; while (lines[cl] && (lines[cl][0] == '#')) cl++;
  if (lines[cl] == NULL) {
    fprintf(stderr, "sim4polish::s4p_linesToPolishGFF3()--  ERROR: Critical error when reading GFF3 record. Skipping.\n");
    return;
  }


  //  Scan mRNA line

  _genDefLine = new char [lengths[cl]];

  r = sscanf(lines[cl], ""u32bitFMT":%s\tsim4db\tmRNA\t"u32bitFMT"\t"u32bitFMT"\t"u32bitFMT"\t%c\t.\t",
                        &_genID, _genDefLine, &dummy1, &dummy2, &_percentIdentity, &sOri);
  if (r != 6) {
    fprintf(stderr, "sim4polish::s4p_linesToPolishGFF3()--  byte "u32bitFMT": '%s'\n", startPosition, lines[cl]);
    fprintf(stderr, "sim4polish::s4p_linesToPolishGFF3()--  Expecting description line, found %d tokens instead of 6.\n", r);
  }

  switch (sOri) {
    case '+' : _strandOrientation = SIM4_STRAND_POSITIVE; break;
    case '-' : _strandOrientation = SIM4_STRAND_NEGATIVE; break;
    case '.' : _strandOrientation = SIM4_STRAND_UNKNOWN;  break;
    default  :  ok = false;
  } 
  

  if (ok == true) {
    //  skip over the first eight columns in the GFF3 format
    clptr = lines[cl];
    while (*clptr!='\t') clptr++; clptr++;
    while (*clptr!='\t') clptr++; clptr++;
    while (*clptr!='\t') clptr++; clptr++;
    while (*clptr!='\t') clptr++; clptr++;
    while (*clptr!='\t') clptr++; clptr++;
    while (*clptr!='\t') clptr++; clptr++;
    while (*clptr!='\t') clptr++; clptr++;
    while (*clptr!='\t') clptr++; clptr++;

    tok = strtok(clptr, "\n");

    crttok = strtok(tok, ";");
    while (crttok) {
      if (!strncmp(crttok, "ID=sim4db", 9)) {
        r = sscanf(crttok, "ID=sim4db"u32bitFMT"", &matchID);
        if (r != 1)  ok = false;
      } else if (!strncmp(crttok, "Name", 4)) {
        if (_estDefLine == 0L)
          _estDefLine = new char [lengths[cl]];
        r = sscanf(crttok, "Name="u32bitFMT":%s", &_estID, _estDefLine);
        if (r != 2)  ok = false;
      } else if (!strncmp(crttok, "Target", 6)) {
        if (_estDefLine == 0L)
          _estDefLine = new char [lengths[cl]];
        r = sscanf(crttok, "Target="u32bitFMT":%s "u32bitFMT" "u32bitFMT" %c", &_estID, _estDefLine, &dummy1, &dummy2, &mOri);
        if (r != 5)  ok = false;
        if (mOri == '+') _matchOrientation = SIM4_MATCH_FORWARD;
        else
          if (mOri == '-') _matchOrientation = SIM4_MATCH_COMPLEMENT;
          else
            ok = false;
      } else if (!strncmp(crttok, "targetLen", 9)) { 
        r = sscanf(crttok, "targetLen="u32bitFMT"", &_estLen);
        if (r != 1)  ok = false;
      } else if (!strncmp(crttok, "pA", 2)) {
        r = sscanf(crttok, "pA="u32bitFMT"", &_estPolyA);
        if (r != 1)  ok = false;
      } else if (!strncmp(crttok, "pT", 2)) {
        r = sscanf(crttok, "pT="u32bitFMT"", &_estPolyT);
        if (r != 1)  ok = false; 
      } else if (!strncmp(crttok, "genRegion", 9)) {
        r = sscanf(crttok, "genRegion="u32bitFMT"-"u32bitFMT"", &_genRegionOffset, &dummy1);
        if (r != 2)  ok = false;
        else 
          _genRegionLength = dummy1 - _genRegionOffset + 1;
      }

      crttok = strtok(NULL, ";");
    }

    //  Check that we read what we should have read so far 
    if ((ok == false) || !_estDefLine  || !_genDefLine || !_estLen || !_matchOrientation || !_strandOrientation) {
      fprintf(stderr, "sim4polish::s4p_linesToPolishGFF3()--  byte "u32bitFMT": '%s'\n", startPosition, lines[cl]);
      fprintf(stderr, "sim4polish::s4p_linesToPolishGFF3()--  Expecting mRNA description line, %s.\n", (ok==false) ? "failed":"incomplete");
    }
  }
 

  //
  //  While we get exons, make exons.
  //

  sim4polishExon    exon;
  u32bit            maxExons = 1024;

  _numExons = 0;
  _exons    = new sim4polishExon [maxExons];

  _numCovered = 0;

  cl++; while (lines[cl] && (lines[cl][0] == '#')) cl++;

  while (lines[cl] && strstr(lines[cl], "\tsim4db\texon\t")) {

    ok = true;

    exon._intronOrientation = SIM4_INTRON_NONE;

    r = sscanf(lines[cl], ""u32bitFMT":%s\tsim4db\texon\t"u32bitFMT"\t"u32bitFMT"\t"u32bitFMT"\t%c\t.\t",
                        &dummy1, dummybuf, &exon._genFrom, &exon._genTo, &exon._percentIdentity, &sOri);
    if (r != 6) {
      fprintf(stderr, "sim4polish::s4p_linesToPolishGFF3()--  byte "u32bitFMT": '%s'\n", startPosition, lines[cl]);
      fprintf(stderr, "sim4polish::s4p_linesToPolishGFF3()--  Expecting exon description line, found %d tokens instead of 6.\n", r);
    }

    if ((dummy1 != _genID) || strcmp(dummybuf, _genDefLine) ||
        (sOri != '+') && (sOri != '-') && (sOri != '.'))
      ok = false;

    if (ok) {
      clptr = lines[cl];
      while (*clptr!='\t') clptr++; clptr++;
      while (*clptr!='\t') clptr++; clptr++;
      while (*clptr!='\t') clptr++; clptr++;
      while (*clptr!='\t') clptr++; clptr++;
      while (*clptr!='\t') clptr++; clptr++;
      while (*clptr!='\t') clptr++; clptr++;
      while (*clptr!='\t') clptr++; clptr++;
      while (*clptr!='\t') clptr++; clptr++;

      tok = strtok(clptr, "\n");

      crttok = strtok(tok, ";");
      while (crttok) {
        if (!strncmp(crttok, "Parent=sim4db", 13)) {
          r = sscanf(crttok, "Parent=sim4db"u32bitFMT"", &dummy1);
          if ((r != 1) || (dummy1 != matchID))  ok = false;

        } else if (!strncmp(crttok, "Target=", 7)) {
          r = sscanf(crttok, "Target=%s "u32bitFMT" "u32bitFMT" %c", &dummybuf, &exon._estFrom, &exon._estTo, &mOri);
          if ((r != 4) ||
              ((mOri == '+') && (_matchOrientation == SIM4_MATCH_COMPLEMENT)) ||
              ((mOri == '-') && (_matchOrientation == SIM4_MATCH_FORWARD))) 
            ok = false;

        } else if (!strncmp(crttok, "nMatches=", 9)) {
          r = sscanf(crttok, "nMatches="u32bitFMT"", &exon._numMatches);
          if (r != 1)  ok = false;
        } else if (!strncmp(crttok, "Gap=", 4)) {
           ; // Handle this later or, better yet, just skip alignment

        } else if (!strncmp(crttok, "intron=", 7)) {
          r = sscanf(crttok, "intron=%s", &dummybuf);
          if (r != 1)  ok = false;
          if (!strcmp(dummybuf, "->"))
            exon._intronOrientation = SIM4_INTRON_POSITIVE;
          else if (!strcmp(dummybuf, "<-"))
            exon._intronOrientation = SIM4_INTRON_NEGATIVE;
          else if (!strcmp(dummybuf, "--"))
            exon._intronOrientation = SIM4_INTRON_AMBIGUOUS;
          else if (!strcmp(dummybuf, "=="))
            exon._intronOrientation = SIM4_INTRON_GAP;
          else 
            ok = false;
        }

        crttok = strtok(NULL, ";");
      }
    }

    //  Check that we read what we should have read so far
    if ((ok == false) || !exon._estFrom  || !exon._estTo || !exon._numMatches) {
      fprintf(stderr, "sim4polish::s4p_linesToPolishGFF3()--  byte "u32bitFMT": '%s'\n", startPosition, lines[cl]);
      fprintf(stderr, "sim4polish::s4p_linesToPolishGFF3()--  Expecting exon description line, %s.\n", (ok==false) ? "failed":"incomplete");
    }

    //  Now load everything into the real exons array:
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

    _exons[_numExons]._numMatchesN = 0;    // Most likely!

    _exons[_numExons]._estAlignment = 0L;
    _exons[_numExons]._genAlignment = 0L;

    _numCovered  += _exons[_numExons]._estTo - _exons[_numExons]._estFrom + 1;
    _numMatches  += _exons[_numExons]._numMatches;
    _numMatchesN += _exons[_numExons]._numMatchesN;

    _numExons++;

    cl++;
    while (lines[cl] && (lines[cl][0] == '#')) cl++;
  }

  if (_numExons == 0) {
    fprintf(stderr, "sim4polish::s4p_linesToPolishGFF3()--  byte "u32bitFMT": '%s'\n", startPosition, lines[cl]);
    fprintf(stderr, "sim4polish::s4p_linesToPolishGFF3()--  WARNING: found ZERO exons?\n");
  }

  _querySeqIdentity = s4p_percentCoverageApprox();

}

