#include "sim4polish.H"

void
sim4polish::s4p_printPolish(FILE *O, u32bit flags) {
  char const  *mOri = "";
  char const  *sOri = "";

  if (_numExons == 0)
    return;

  fprintf(O, "sim4begin\n");

  switch (_matchOrientation) {
    case SIM4_MATCH_FORWARD:
      mOri = "forward";
      break;
    case SIM4_MATCH_COMPLEMENT:
      mOri = "complement";
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
    default:
      fprintf(stderr, "sim4reader: Unknown strandOrientation '"u32bitFMT"' in printPolish()\n", _matchOrientation);
      sOri = "UNKNOWN";
      break;
  }

  fprintf(O, ""u32bitFMT"["u32bitFMT"-"u32bitFMT"-"u32bitFMT"] "u32bitFMT"["u32bitFMT"-"u32bitFMT"] <"u32bitFMT"-"u32bitFMT"-"u32bitFMT"-%s-%s>\n",
          _estID, _estLen, _estPolyA, _estPolyT,
          _genID, _genRegionOffset, _genRegionLength,
          _numMatches, _numMatchesN, _percentIdentity, mOri, sOri);

  if (_comment)
    fprintf(O, "comment=%s\n", _comment);

  if ((_estDefLine) && ((flags & S4P_PRINTPOLISH_NODEFS) == 0))
    fprintf(O, "edef=%s\n", _estDefLine);

  if ((_genDefLine) && ((flags & S4P_PRINTPOLISH_NODEFS) == 0))
    fprintf(O, "ddef=%s\n", _genDefLine);

  for (u32bit i=0; i<_numExons; i++) {
    fprintf(O, ""u32bitFMT"-"u32bitFMT" ("u32bitFMT"-"u32bitFMT") <"u32bitFMT"-"u32bitFMT"-"u32bitFMT">",
            _exons[i]._estFrom, _exons[i]._estTo,
            _exons[i]._genFrom, _exons[i]._genTo,
            _exons[i]._numMatches, _exons[i]._numMatchesN, _exons[i]._percentIdentity);

    switch (_exons[i]._intronOrientation) {
      case SIM4_INTRON_POSITIVE:
        fprintf(O, " ->\n");
        break;
      case SIM4_INTRON_NEGATIVE:
        fprintf(O, " <-\n");
        break;
      case SIM4_INTRON_AMBIGUOUS:
        fprintf(O, " --\n");
        break;
      case SIM4_INTRON_GAP:
        fprintf(O, " ==\n");
        break;
      case SIM4_INTRON_NONE:
        fprintf(O, "\n");
        break;
      default:
        fprintf(O, "\n");
        break;
    }
  }

  if ((flags & S4P_PRINTPOLISH_NOALIGNS) == 0)
    for (u32bit i=0; i<_numExons; i++) {
      if (_exons[i]._estAlignment)
        fprintf(O, "%s\n", _exons[i]._estAlignment);
      if (_exons[i]._genAlignment)
        fprintf(O, "%s\n", _exons[i]._genAlignment);
    }

  fprintf(O, "sim4end\n");
}
