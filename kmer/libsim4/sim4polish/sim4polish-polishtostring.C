#include "sim4polish.H"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <assert.h>

//#define DEBUG_CIGAR

const char *mOriFWD = "forward";
const char *mOriCMP = "complement";
const char *mOriERR = "error";
const char *mOriDEF = "UNKNOWN";

const char *sOriFWD = "forward";
const char *sOriREV = "reverse";
const char *sOriUNK = "unknown";
const char *sOriINT = "intractable";
const char *sOriABT = "aborted";
const char *sOriERR = "error";
const char *sOriDEF = "UNKNOWN";

const char *iOriPOS = " ->";
const char *iOriNEG = " <-";
const char *iOriAMB = " --";
const char *iOriGAP = " ==";
const char *iOriERR = " ??";
const char *iOriNOO = "";

sim4polishStyle sim4polishStyleDefault = sim4polishS4DB;
u32bit          sim4polishPolishID     = 0;




char *
encodeGap(char *ref, char *tgt) {

  if ((ref == 0L) || (tgt == 0L))
    return(0L);

  u32bit lenref = strlen(ref);
  u32bit lentgt = strlen(tgt);

  assert(lenref == lentgt);

  char   *gap    = new char [3 * lenref];
  char   *gpp    = gap;

  char    gaptyp = 0;
  u32bit  gapcnt = 0;

  for (u32bit i=0; i<lenref; i++) {
    if        ((ref[i] == '-') && (tgt[i] != '-')) {
      if (gaptyp != 'I') {
        if (gaptyp != 0) {
          sprintf(gpp, "%c"u32bitFMT" ", gaptyp, gapcnt);
          while (*gpp) gpp++;
        }
        gaptyp = 'I';
        gapcnt = 0;
      }
      gapcnt++;
    } else if ((ref[i] != '-') && (tgt[i] == '-')) {
      if (gaptyp != 'D') {
        if (gaptyp != 0) {
          sprintf(gpp, "%c"u32bitFMT" ", gaptyp, gapcnt);
          while (*gpp) gpp++;
        }
        gaptyp = 'D';
        gapcnt = 0;
      }
      gapcnt++;
    } else if ((ref[i] == '-') && (tgt[i] == '-')) {
      assert(0);
    } else {
      if (gaptyp != 'M') {
        if (gaptyp != 0) {
          sprintf(gpp, "%c"u32bitFMT" ", gaptyp, gapcnt);
          while (*gpp) gpp++;
        }
        gaptyp = 'M';
        gapcnt = 0;
      }
      gapcnt++;
    }
  }

  if (gaptyp != 0) {
    sprintf(gpp, "%c"u32bitFMT"", gaptyp, gapcnt);
    while (*gpp) gpp++;
  }

#ifdef DEBUG_CIGAR
  fprintf(stderr, "REF=%s\n", ref);
  fprintf(stderr, "TGT=%s\n", tgt);
  fprintf(stderr, "GAP=%s\n", gap);
  fprintf(stderr, "---\n");
#endif

  return(gap);
}





char *
sim4polish::s4p_polishToString(sim4polishStyle style) {
  char *ret = NULL;

  if (_numExons == 0)
    return(ret);

  switch (sim4polishStyleDefault) {
    case sim4polishS4DB:
      ret = s4p_polishToStringS4DB();
      break;
    case sim4polishGFF3:
      ret = s4p_polishToStringGFF3();
      break;
    case sim4polishATAC:
      ret = s4p_polishToStringATAC();
      break;
    default:
      fprintf(stderr, "s4p_polishToString()-- unknown style default='%d' style='%d'\n",
              sim4polishStyleDefault, style);
      exit(1);
  }

  return(ret);
}

char *
sim4polish::s4p_polishToStringS4DB(void) {
  const char   *mOri = mOriDEF;
  const char   *sOri = sOriDEF;
  const char   *iOri = iOriERR;

  //  Make a decent estimate of how much space we'll need to store the string
  //
  u32bit spaceNeeded = (1024 + 128 * _numExons +
                        ((_comment)    ? strlen(_comment)    : 0) +
                        ((_estDefLine) ? strlen(_estDefLine) : 0) +
                        ((_genDefLine) ? strlen(_genDefLine) : 0));

  for (u32bit i=0; i<_numExons; i++)
    if (_exons[i]._estAlignment)
      spaceNeeded += 2 * strlen(_exons[i]._estAlignment);

  char *outs = new char [spaceNeeded];
  char *outc = outs;

  switch (_matchOrientation) {
    case SIM4_MATCH_FORWARD:     mOri = mOriFWD;  break;
    case SIM4_MATCH_COMPLEMENT:  mOri = mOriCMP;  break;
    case SIM4_MATCH_ERROR:       mOri = mOriERR;  break;
    default:
      fprintf(stderr, "sim4reader: Unknown matchOrientation '"u32bitFMT"' in printPolish()\n", _matchOrientation);
      mOri = mOriDEF;
      break;
  }

  switch (_strandOrientation) {
    case SIM4_STRAND_POSITIVE:    sOri = sOriFWD;  break;
    case SIM4_STRAND_NEGATIVE:    sOri = sOriREV;  break;
    case SIM4_STRAND_UNKNOWN:     sOri = sOriUNK;  break;
    case SIM4_STRAND_INTRACTABLE: sOri = sOriINT;  break;
    case SIM4_STRAND_FAILED:      sOri = sOriABT;  break;
    case SIM4_STRAND_ERROR:       sOri = sOriERR;  break;
    default:
      fprintf(stderr, "sim4reader: Unknown strandOrientation '"u32bitFMT"' in printPolish()\n", _matchOrientation);
      sOri = sOriDEF;
      break;
  }

  sprintf(outc, "sim4begin\n"u32bitFMT"["u32bitFMT"-"u32bitFMT"-"u32bitFMT"] "u32bitFMT"["u32bitFMT"-"u32bitFMT"] <"u32bitFMT"-"u32bitFMT"-"u32bitFMT"-%s-%s>\n",
          _estID, _estLen, _estPolyA, _estPolyT,
          _genID, _genRegionOffset, _genRegionLength,
          _numMatches, _numMatchesN, _percentIdentity, mOri, sOri);
  while (*outc)  outc++;

  if (_comment) {
    sprintf(outc, "comment=%s\n", _comment);
    while (*outc)  outc++;
  }

  if (_estDefLine) {
    sprintf(outc, "edef=%s\n", _estDefLine);
    while (*outc)  outc++;
  }

  if (_genDefLine) {
    sprintf(outc, "ddef=%s\n", _genDefLine);
    while (*outc)  outc++;
  }

  for (u32bit i=0; i<_numExons; i++) {
    switch (_exons[i]._intronOrientation) {
      case SIM4_INTRON_POSITIVE:    iOri = iOriPOS;  break;
      case SIM4_INTRON_NEGATIVE:    iOri = iOriNEG;  break;
      case SIM4_INTRON_AMBIGUOUS:   iOri = iOriAMB;  break;
      case SIM4_INTRON_GAP:         iOri = iOriGAP;  break;
      case SIM4_INTRON_ERROR:       iOri = iOriERR;  break;
      default:                      iOri = iOriNOO;  break;
    }

    sprintf(outc, ""u32bitFMT"-"u32bitFMT" ("u32bitFMT"-"u32bitFMT") <"u32bitFMT"-"u32bitFMT"-"u32bitFMT">%s\n",
            _exons[i]._estFrom, _exons[i]._estTo,
            _exons[i]._genFrom, _exons[i]._genTo,
            _exons[i]._numMatches, _exons[i]._numMatchesN, _exons[i]._percentIdentity, iOri);

    while (*outc)  outc++;
  }

  for (u32bit i=0; i<_numExons; i++) {
    if (_exons[i]._estAlignment) {
      strcpy(outc, _exons[i]._estAlignment);
      while (*outc)  outc++;
      *outc++ = '\n';
    }
    if (_exons[i]._genAlignment) {
      strcpy(outc, _exons[i]._genAlignment);
      while (*outc)  outc++;
      *outc++ = '\n';
    }
  }

  strcpy(outc, "sim4end\n");

  return(outs);
}

char *
sim4polish::s4p_polishToStringGFF3(void) {

  //  9 columns, tab separated
  //  tab, newline, cr and control MUST be escaped
  //  reserved letters:  ; = % & ,
  //  spaces ARE ALLOWED in fields
  //  undefined values should use '.'
  //
  //  1 seqid, genome name (a-zA-Z0-9.:^*$@!+_?-|), no whitespace (??) and not begin with >
  //  2 source ("sim4cc")
  //  3 type ("mRNA" or "exon")
  //  4 begin, 1-based
  //  5 end, zero-length start=end, to the right of this base
  //  6 score, ill-defined
  //  7 strand
  //  8 phase
  //  9 attributes
  //      ID        (unique within scope of file)
  //      Name      (display name)
  //      Parent    ()
  //      Target
  //      Gap
  //      Derives_from
  //      Note
  //      Dbxref
  //      Ontology_term
  //      Is_circular
  //


  //  Make a decent estimate of how much space we'll need to store the string
  //
  u32bit spaceNeeded = (1024 + 128 * _numExons +
                        ((_comment)    ? strlen(_comment)    : 0) +
                        ((_estDefLine) ? strlen(_estDefLine) : 0) +
                        ((_genDefLine) ? strlen(_genDefLine) : 0));

  for (u32bit i=0; i<_numExons; i++)
    if (_exons[i]._estAlignment)
      spaceNeeded += 2 * strlen(_exons[i]._estAlignment);

  char   *outs = new char [spaceNeeded];
  char   *outc = outs;

  //  Find extents of this match.
  u32bit  estbgn = _exons[0]._estFrom;
  u32bit  estend = _exons[_numExons-1]._estTo;
  u32bit  genbgn = _exons[0]._genFrom;
  u32bit  genend = _exons[_numExons-1]._genTo;

  for (u32bit i=0; i<_numExons; i++) {
    if (_exons[i]._genFrom < genbgn)  genbgn = _exons[i]._genFrom;
    if (_exons[i]._genTo   < genbgn)  genbgn = _exons[i]._genTo;
    if (genend < _exons[i]._genFrom)  genend = _exons[i]._genFrom;
    if (genend < _exons[i]._genTo)    genend = _exons[i]._genTo;

    if (_exons[i]._genFrom < estbgn)  estbgn = _exons[i]._genFrom;
    if (_exons[i]._genTo   < estbgn)  estbgn = _exons[i]._genTo;
    if (estend < _exons[i]._genFrom)  estend = _exons[i]._genFrom;
    if (estend < _exons[i]._genTo)    estend = _exons[i]._genTo;
  }

  //  Find the orientation
  char    mOri = '?';

  if (_matchOrientation == SIM4_MATCH_FORWARD)     mOri = '+';
  if (_matchOrientation == SIM4_MATCH_COMPLEMENT)  mOri = '-';

  //  Get rid of spaces in the names (and do it non-destructively).

  u32bit  estDefSpace = 0;
  u32bit  genDefSpace = 0;

  while ((_estDefLine[estDefSpace]) && (isspace(_estDefLine[estDefSpace]) == 0))
    estDefSpace++;
  while ((_genDefLine[genDefSpace]) && (isspace(_genDefLine[genDefSpace]) == 0))
    genDefSpace++;

  char estDefChar = _estDefLine[estDefSpace];
  char genDefChar = _genDefLine[genDefSpace];

  _estDefLine[estDefSpace] = 0;
  _genDefLine[genDefSpace] = 0;

  //  The main mRNA match line.

  sprintf(outc, "%s\tsim4cc\tmRNA\t"u32bitFMT"\t"u32bitFMT"\t.\t%c\t.\t",
          _genDefLine, genbgn, genend, mOri);
  while (*outc)  outc++;

  sprintf(outc, "ID=sim4cc"u32bitFMT";Name=%s;Target=%s "u32bitFMT" "u32bitFMT" %c\n",
          sim4polishPolishID, _estDefLine, _estDefLine, estbgn, estend, mOri);
  while (*outc)  outc++;

  //  Exons.

  for (u32bit i=0; i<_numExons; i++) {
    char *gap = encodeGap(_exons[i]._genAlignment, _exons[i]._estAlignment);

    sprintf(outc, "%s\tsim4cc\texon\t"u32bitFMT"\t"u32bitFMT"\t.\t%c\t.\t",
            _genDefLine, _exons[i]._genFrom, _exons[i]._genTo, mOri);
    while (*outc)  outc++;

    if (gap)
      sprintf(outc, "Parent=sim4cc"u32bitFMT";Target=%s "u32bitFMT" "u32bitFMT" %c;Gap=%s\n",
              sim4polishPolishID, _estDefLine, _exons[i]._estFrom, _exons[i]._estTo, '-', gap);
    else
      sprintf(outc, "Parent=sim4cc"u32bitFMT";Target=%s "u32bitFMT" "u32bitFMT" %c\n",
              sim4polishPolishID, _estDefLine, _exons[i]._estFrom, _exons[i]._estTo, '-');
    while (*outc)  outc++;

    delete [] gap;
  }

  sim4polishPolishID++;

  _estDefLine[estDefSpace] = estDefChar;
  _genDefLine[genDefSpace] = genDefChar;

  return(outs);
}

char *
sim4polish::s4p_polishToStringATAC(void) {
  return(0L);
}



