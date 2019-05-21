
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' (http://kmer.sourceforge.net)
 *  both originally distributed by Applera Corporation under the GNU General
 *  Public License, version 2.
 *
 *  Canu branched from Celera Assembler at its revision 4587.
 *  Canu branched from the kmer project at its revision 1994.
 *
 *  This file is derived from:
 *
 *    src/stores/gkLibrary.C
 *    src/stores/gkStore.C
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2014-NOV-26 to 2015-AUG-10
 *      are Copyright 2014-2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2015-OCT-09
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *    Sergey Koren beginning on 2015-DEC-09
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "sqStore.H"

//
//  sqLibrary is lightweight, except for three functions that need to parse strings
//

void
sqLibrary::sqLibrary_parsePreset(char *p) {

  if (strcasecmp(p, "contig") == 0) {
    _readCorrection   = SQ_CORRECTION_NONE;
    _readType         = SQ_READTYPE_CONTIG;
    return;
  }

  if (strcasecmp(p, "pacbio-raw") == 0) {
    _readCorrection   = SQ_CORRECTION_CONSENSUS;
    _readType         = SQ_READTYPE_PACBIO_RAW;
    _checkForSubReads = true;
    return;
  }

  if (strcasecmp(p, "pacbio-corrected") == 0) {
    _readCorrection   = SQ_CORRECTION_NONE;
    _readType         = SQ_READTYPE_PACBIO_CORRECTED;
    return;
  }

  if (strcasecmp(p, "nanopore-raw") == 0) {
    _readCorrection   = SQ_CORRECTION_CONSENSUS;
    _readType         = SQ_READTYPE_NANOPORE_RAW;
    _checkForSubReads = true;
    return;
  }

  if (strcasecmp(p, "nanopore-corrected") == 0) {
    _readCorrection   = SQ_CORRECTION_NONE;
    _readType         = SQ_READTYPE_NANOPORE_CORRECTED;
    return;
  }

  if (strcasecmp(p, "pacbio-hifi") == 0) {
    _readCorrection   = SQ_CORRECTION_NONE;
    _readType         = SQ_READTYPE_PACBIO_HIFI;
    return;
  }

  fprintf(stderr, "sqLibrary::sqLibrary_parsePreset()--  ERROR: unknown preset '%s'\n", p);
  exit(1);
}


void
sqLibrary::sqLibrary_setReadType(char *t) {

  if      (strcasecmp(t, "generic") == 0)
    _readType = SQ_READTYPE_GENERIC;

  else if (strcasecmp(t, "contig") == 0)
    _readType = SQ_READTYPE_CONTIG;

  else if (strcasecmp(t, "pacbio_raw") == 0)
    _readType = SQ_READTYPE_PACBIO_RAW;

  else if (strcasecmp(t, "pacbio_corrected") == 0)
    _readType = SQ_READTYPE_PACBIO_CORRECTED;

  else if (strcasecmp(t, "nanopore_raw") == 0)
    _readType = SQ_READTYPE_NANOPORE_RAW;

  else if (strcasecmp(t, "nanopore_corrected") == 0)
    _readType = SQ_READTYPE_NANOPORE_CORRECTED;

  else if (strcasecmp(t, "pacbio_hifi") == 0)
    _readType = SQ_READTYPE_PACBIO_HIFI;

  else
    fprintf(stderr, "sqLibrary::sqLibrary_setReadType()--  ERROR: unknown read type '%s'\n",
            t), exit(1);
}


char const *
sqLibrary::sqLibrary_readTypeString(void) {
  switch (_readType) {
    case SQ_READTYPE_GENERIC:             return("generic");             break;
    case SQ_READTYPE_CONTIG:              return("contig");              break;
    case SQ_READTYPE_PACBIO_RAW:          return("pacbio-raw");          break;
    case SQ_READTYPE_PACBIO_CORRECTED:    return("pacbio-corrected");    break;
    case SQ_READTYPE_NANOPORE_RAW:        return("nanopore-raw");        break;
    case SQ_READTYPE_NANOPORE_CORRECTED:  return("nanopore-corrected");  break;
    case SQ_READTYPE_PACBIO_HIFI:         return("pacbio-hifi");         break;
    default:                              return("invalid");             break;
  }
}


void
sqLibrary::sqLibrary_setReadCorrection(char *t) {

  if      (strcasecmp(t, "none") == 0)
    _finalTrim = SQ_CORRECTION_NONE;

  else if (strcasecmp(t, "consensus") == 0)
    _finalTrim = SQ_CORRECTION_CONSENSUS;

  else if (strcasecmp(t, "mer") == 0)
    _finalTrim = SQ_CORRECTION_MER;

  else
    fprintf(stderr, "sqLibrary::sqLibrary_setReadCorrection()--  ERROR: unknown read correction '%s'\n",
            t), exit(1);
}


char const *
sqLibrary::sqLibrary_readCorrectionString(void) {
  switch (_readCorrection) {
    case SQ_CORRECTION_NONE:              return("none");                break;
    case SQ_CORRECTION_CONSENSUS:         return("consensus");           break;
    case SQ_CORRECTION_MER:               return("mer");                 break;
    default:                              return("invalid");             break;
  }
}


void
sqLibrary::sqLibrary_setFinalTrim(char *t) {

  if      (strcasecmp(t, "none") == 0)
    _finalTrim = SQ_FINALTRIM_NONE;

  else if (strcasecmp(t, "largest") == 0)
    _finalTrim = SQ_FINALTRIM_LARGEST_COVERED;

  else if (strcasecmp(t, "bestedge") == 0)
    _finalTrim = SQ_FINALTRIM_BEST_EDGE;

  else
    fprintf(stderr, "sqLibrary::sqLibrary_setFinalTrim()--  ERROR: unknown final trim '%s'\n",
            t), exit(1);
}


char const *
sqLibrary::sqLibrary_finalTrimString(void) {
  switch (_finalTrim) {
    case SQ_FINALTRIM_NONE:               return("none");                break;
    case SQ_FINALTRIM_LARGEST_COVERED:    return("largestCovered");      break;
    case SQ_FINALTRIM_BEST_EDGE:          return("bestEdge");            break;
    default:                              return("invalid");             break;
  }
}
