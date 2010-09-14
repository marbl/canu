#include "sim4polish.H"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

char*
sim4polish::s4p_polishToString(void) {
  u32bit        spaceNeeded = 0;
  char         *outc = NULL;
  char         *outs = NULL;
  const char   *mOri = "";
  const char   *sOri = "";
  const char   *iOri = "";

  if (_numExons == 0)
    return(0L);

  //  Make a decent estimate of how much space we'll need to
  //  store the string
  //
  spaceNeeded = 1024 + 128 * _numExons;

  if (_comment)     spaceNeeded += strlen(_comment);
  if (_estDefLine)  spaceNeeded += strlen(_estDefLine);
  if (_genDefLine)  spaceNeeded += strlen(_genDefLine);

  for (u32bit i=0; i<_numExons; i++) {
    if (_exons[i]._estAlignment)
      spaceNeeded += 2 * strlen(_exons[i]._estAlignment);
  }

  outc = outs = new char [spaceNeeded];

  switch (_matchOrientation) {
    case SIM4_MATCH_FORWARD:
      mOri = "forward";
      break;
    case SIM4_MATCH_COMPLEMENT:
      mOri = "complement";
      break;
    case SIM4_MATCH_ERROR:
      mOri = "error";
      break;
    default:
      fprintf(stderr, "sim4reader: Unknown matchOrientation '"u32bitFMT"' in printPolish()\n", _matchOrientation);
      mOri = "UNKNOWN";
      break;
  }

  switch (_strandOrientation) {
    case SIM4_STRAND_POSITIVE:
      sOri = "forward";
      break;
    case SIM4_STRAND_NEGATIVE:
      sOri = "reverse";
      break;
    case SIM4_STRAND_UNKNOWN:
      sOri = "unknown";
      break;
    case SIM4_STRAND_INTRACTABLE:
      sOri = "intractable";
      break;
    case SIM4_STRAND_FAILED:
      sOri = "aborted";
      break;
    case SIM4_STRAND_ERROR:
      sOri = "error";
      break;
    default:
      fprintf(stderr, "sim4reader: Unknown strandOrientation '"u32bitFMT"' in printPolish()\n", _matchOrientation);
      sOri = "UNKNOWN";
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
      case SIM4_INTRON_POSITIVE:
        iOri = " ->";
        break;
      case SIM4_INTRON_NEGATIVE:
        iOri = " <-";
        break;
      case SIM4_INTRON_AMBIGUOUS:
        iOri = " --";
        break;
      case SIM4_INTRON_GAP:
        iOri = " ==";
        break;
      case SIM4_INTRON_ERROR:
        iOri = " ??";
        break;
      default:
        iOri = "";
        break;
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
      *outc = '\n';
      outc++;
    }
    if (_exons[i]._genAlignment) {
      strcpy(outc, _exons[i]._genAlignment);
      while (*outc)  outc++;
      *outc = '\n';
      outc++;
    }
  }

  strcpy(outc, "sim4end\n");

  return(outs);
}

