
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 1999-2004, Applera Corporation. All rights reserved.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received (LICENSE.txt) a copy of the GNU General Public
 * License along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *************************************************************************/

static char const *rcsid = "$Id: AS_GKP_errors.c,v 1.19 2012-02-03 21:47:58 brianwalenz Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

#include "AS_global.h"
#include "AS_GKP_include.h"

#include <vector>

using namespace std;

//  A rather unwieldy error handling system.  We want to count the
//  number of times we see each error, and try to make it somewhat
//  easy to add a new error.

#define AS_GKP_NUM_ERRORS  128

class AS_GKP_ePL {
 public:
  AS_GKP_ePL() {
    memset(errorCs, 0, sizeof(uint32) * AS_GKP_NUM_ERRORS);
  };
  ~AS_GKP_ePL() {
  };

  uint32      errorCs[AS_GKP_NUM_ERRORS];
};

static char         const *errorMs[AS_GKP_NUM_ERRORS] = {0};
static char         const *errorSs[AS_GKP_NUM_ERRORS] = {0};
static uint32              errorCs[AS_GKP_NUM_ERRORS] = {0};
static vector<AS_GKP_ePL>  libError;

void
AS_GKP_reportError(int error, uint32 libIID, ...) {
  va_list ap;

  if (errorMs[0] == 0) {

    //
    //  The error message reported in the log file
    //

    errorMs[AS_GKP_BAT_ZERO_UID           ] = "# BAT Error: Batch has zero or no UID; can't add it.\n";
    errorMs[AS_GKP_BAT_EXISTS             ] = "# BAT Error: Batch %s exists, can't add it again.\n";

    errorMs[AS_GKP_FRG_INVALID_CHAR_SEQ   ] = "# FRG Error: Fragment %s invalid char %c at position "F_SIZE_T" in sequence.\n";
    errorMs[AS_GKP_FRG_INVALID_CHAR_QLT   ] = "# FRG Error: Fragment %s invalid char %c at position " F_SIZE_T " in quality.\n";
    errorMs[AS_GKP_FRG_INVALID_LENGTH     ] = "# FRG Error: Fragment %s sequence length %d != %d quality length.\n";
    errorMs[AS_GKP_FRG_ZERO_UID           ] = "# FRG Error: Fragment has UID of zero; can't add it.\n";
    errorMs[AS_GKP_FRG_EXISTS             ] = "# FRG Error: Fragment %s exists; can't add it again.\n";
    errorMs[AS_GKP_FRG_SEQ_TOO_LONG       ] = "# FRG Error: Fragment %s sequence length %d > %d max allowed sequence length.\n";
    errorMs[AS_GKP_FRG_SEQ_TOO_SHORT      ] = "# FRG Error: Fragment %s sequence length %d < %d min allowed sequence length.\n";
    errorMs[AS_GKP_FRG_CLR_BGN            ] = "# FRG Alert: Fragment %s invalid clear range (%d,%d) valid range is (0,%d) -- fixing clear_rng.bgn.\n";
    errorMs[AS_GKP_FRG_CLR_END            ] = "# FRG Alert: Fragment %s invalid clear range (%d,%d) valid range is (0,%d) -- fixing clear_rng.end.\n";
    errorMs[AS_GKP_FRG_CLR_TOO_SHORT      ] = "# FRG Error: Fragment %s clear range length %d < %d min allowed length.\n";
    errorMs[AS_GKP_FRG_CLR_INVALID        ] = "# FRG Error: Fragment %s clear range clr is invalid; begin = %d > end = %d.\n";
    errorMs[AS_GKP_FRG_UNKNOWN_LIB        ] = "# FRG Error: Fragment %s references unknown library %s.\n";
    errorMs[AS_GKP_FRG_LOADED_DELETED     ] = "# FRG Alert: Fragment %s loaded, but marked as deleted due to errors previously reported.\n";
    errorMs[AS_GKP_FRG_DOESNT_EXIST       ] = "# FRG Error: Fragment %s does not exist, can't delete it.\n";
    errorMs[AS_GKP_FRG_HAS_MATE           ] = "# FRG Error: Fragment %s has mate pair relationship, can't delete it.\n";
    errorMs[AS_GKP_FRG_UNKNOWN_ACTION     ] = "# FRG Error: invalid action %c.\n";

    errorMs[AS_GKP_ILL_NOT_SEQ_START_LINE ] = "# ILL Error: File '%s': seq name '%s' is not a sequence start line.\n";
    errorMs[AS_GKP_ILL_NOT_QLT_START_LINE ] = "# ILL Error: File '%s': qlt name '%s' is not a quality start line.\n";
    errorMs[AS_GKP_ILL_SEQ_QLT_NAME_DIFFER] = "# ILL Error: File '%s': seq/qlt names differ; seq='%s' qlt='%s'\n";
    errorMs[AS_GKP_ILL_SEQ_QLT_LEN_DIFFER ] = "# ILL Error: File '%s': seq/qlt lengths differ for read '%s'; seq=%d qlt=%d\n";
    errorMs[AS_GKP_ILL_SEQ_TOO_LONG       ] = "# ILL Alert: Fragment %s is longer than gkpShortReadLength bases, truncating from "F_U32" to "F_U32" bases.\n";
    errorMs[AS_GKP_ILL_CANT_OPEN_INPUT    ] = "# ILL Error: Couldn't open Illumina file '%s' for reading: %s\n";
    errorMs[AS_GKP_ILL_BAD_QV             ] = "# ILL Alert: Fragment %s has "F_U32" invalid QVs.\n";

    errorMs[AS_GKP_LIB_ILLEGAL_MEAN_STDDEV] = "# LIB Alert: Library %s has lllegal mean (%g) and standard deviation (%g); reset to mean 3000, stddev 300.\n";
    errorMs[AS_GKP_LIB_INVALID_MEAN       ] = "# LIB Alert: Library %s has invalid mean (%g); reset mean to 10 * stddev = %g.\n";
    errorMs[AS_GKP_LIB_INVALID_STDDEV     ] = "# LIB Alert: Library %s has invalid stddev (%g); reset stddev to 0.1 * mean = %g.\n";
    errorMs[AS_GKP_LIB_STDDEV_TOO_BIG     ] = "# LIB Alert: Library %s stddev (%g) too big for mean (%g); reset stddev to 0.1 * mean = %g.\n";
    errorMs[AS_GKP_LIB_STDDEV_TOO_SMALL   ] = "# LIB Alert: Library %s has suspicious mean (%g) and standard deviation (%g); reset stddev to 0.10 * mean = %g.\n";
    errorMs[AS_GKP_LIB_EXISTS             ] = "# LIB Error: Library %s,"F_IID" already exists; can't add it again.\n";
    errorMs[AS_GKP_LIB_ZERO_UID           ] = "# LIB Error: Library has zero or no UID; can't add it.\n";
    errorMs[AS_GKP_LIB_DOESNT_EXIST_UPDATE] = "# LIB Error: Library %s does not exist, can't update it.\n";
    errorMs[AS_GKP_LIB_UNKNOWN_ACTION     ] = "# LIB Error: invalid action %c.\n";

    errorMs[AS_GKP_LKG_SELF_LINK          ] = "# LKG Error: Can't make a link from fragment %s to itself.\n";
    errorMs[AS_GKP_LKG_UNSUPPORTED_TYPE   ] = "# LKG Error: Unsupported LKG type '%c' for frags %s,%s in library %s.\n";
    errorMs[AS_GKP_LKG_FRG_DOESNT_EXIST   ] = "# LKG Error: Fragment %s not previously defined.\n";
    errorMs[AS_GKP_LKG_FRG_DELETED        ] = "# LKG Error: Fragment %s is marked as deleted.\n";
    errorMs[AS_GKP_LKG_ALREADY_MATED      ] = "# LKG Error: Fragment %s,"F_IID" already has mate of iid="F_IID"; wanted to set to %s,"F_IID".\n";
    errorMs[AS_GKP_LKG_LIB_DOESNT_EXIST   ] = "# LKG Error: Library %s not previously defined.\n";
    errorMs[AS_GKP_LKG_DIFFERENT_LIB      ] = "# LKG Error: Fragment "F_IID" in lib "F_IID", different from fragment "F_IID" in lib "F_IID".\n";
    errorMs[AS_GKP_LKG_UNMATED_LIB        ] = "# LKG Error: Fragments "F_IID" and "F_IID" cannot be added to unmated library "F_IID".\n";
    errorMs[AS_GKP_LKG_DIFFERENT_ORIENT   ] = "# LKG Error: Fragments "F_IID" (mate type %c) and "F_IID" (mate type %c) not allowed in library "F_IID" (mate type %c).\n";
    errorMs[AS_GKP_LKG_UNKNOWN_ACTION     ] = "# LKG Error: Unknown action %c.\n";

    errorMs[AS_GKP_SFF_UID_ERROR          ] = "# SFF Error: 454 Universal Accession Number '%s' out of range.\n";
    errorMs[AS_GKP_SFF_ALREADY_EXISTS     ] = "# SFF Error: Fragment %s exists, can't add it again.\n";
    errorMs[AS_GKP_SFF_TOO_SHORT          ] = "# SFF Alert: Fragment %s too short (assembler minimum length is %d bases).\n";
    errorMs[AS_GKP_SFF_TOO_LONG           ] = "# SFF Alert: Fragment %s trimmed from %d bases to (assembler maximum length) %d bases.\n";
    errorMs[AS_GKP_SFF_N                  ] = "# SFF Alert: Fragment %s contains an N; deleted.\n";

    errorMs[AS_GKP_PLC_SAME_CONSTRAINT    ] = "# PLC Error: Can't constrain fragment %s by the same fragment %s on both sides.\n";
    errorMs[AS_GKP_PLC_SELF_CONSTRAINT    ] = "# PLC Error: Can't constrain fragment %s using itself.\n";
    errorMs[AS_GKP_PLC_FRG_DOESNT_EXIST   ] = "# PLC Error: Fragment %s not previously defined.\n";
    errorMs[AS_GKP_PLC_FRG_DELETED        ] = "# PLC Error: Fragment %s is marked as deleted.\n";
    errorMs[AS_GKP_PLC_ALREADY_CONSTRAINED] = "# PLC Error: Fragment %s,"F_IID" already has constraint enforced on it by %s,"F_IID" %s,"F_IID"; wanted to set to %s,"F_IID" %s"F_IID".\n";
    errorMs[AS_GKP_PLC_UNKNOWN_ACTION     ] = "# PLC Error: Unknown action %c.\n";

    errorMs[AS_GKP_UNKNOWN_MESSAGE        ] = "# GKP Error: Unknown message with type %s.\n";

    //
    //  A short description of the error, for summarizing the errors
    //  at the end of a run
    //

    errorSs[AS_GKP_BAT_ZERO_UID           ] = "# BAT Error: Batch has zero or no UID.\n";
    errorSs[AS_GKP_BAT_EXISTS             ] = "# BAT Error: Batch already exists.\n";

    errorSs[AS_GKP_FRG_INVALID_CHAR_SEQ   ] = "# FRG Error: invalid char in sequence.\n";
    errorSs[AS_GKP_FRG_INVALID_CHAR_QLT   ] = "# FRG Error: invalid char in quality.\n";
    errorSs[AS_GKP_FRG_INVALID_LENGTH     ] = "# FRG Error: sequence length != quality length.\n";
    errorSs[AS_GKP_FRG_ZERO_UID           ] = "# FRG Error: Fragment has UID of zero; can't add it.\n";
    errorSs[AS_GKP_FRG_EXISTS             ] = "# FRG Error: fragment exists; can't add it again.\n";
    errorSs[AS_GKP_FRG_SEQ_TOO_LONG       ] = "# FRG Error: sequence length > max allowed sequence length.\n";
    errorSs[AS_GKP_FRG_SEQ_TOO_SHORT      ] = "# FRG Error: sequence length < min allowed sequence length.\n";
    errorSs[AS_GKP_FRG_CLR_BGN            ] = "# FRG Alert: invalid clear range -- fixing clear_rng.bgn.\n";
    errorSs[AS_GKP_FRG_CLR_END            ] = "# FRG Alert: invalid clear range -- fixing clear_rng.end.\n";
    errorSs[AS_GKP_FRG_CLR_TOO_SHORT      ] = "# FRG Error: clear range length < min allowed length.\n";
    errorSs[AS_GKP_FRG_CLR_INVALID        ] = "# FRG Error: clear range clr is invalid.\n";
    errorSs[AS_GKP_FRG_UNKNOWN_LIB        ] = "# FRG Error: references unknown library.\n";
    errorSs[AS_GKP_FRG_LOADED_DELETED     ] = "# FRG Alert: loaded, but marked as deleted due to errors previously reported.\n";
    errorSs[AS_GKP_FRG_DOESNT_EXIST       ] = "# FRG Error: does not exist, can't delete.\n";
    errorSs[AS_GKP_FRG_HAS_MATE           ] = "# FRG Error: has mate pair relationship, can't delete.\n";
    errorSs[AS_GKP_FRG_UNKNOWN_ACTION     ] = "# FRG Error: invalid action.\n";

    errorSs[AS_GKP_ILL_NOT_SEQ_START_LINE ] = "# ILL Error: not a sequence start line.\n";
    errorSs[AS_GKP_ILL_NOT_QLT_START_LINE ] = "# ILL Error: not a quality start line.\n";
    errorSs[AS_GKP_ILL_SEQ_QLT_NAME_DIFFER] = "# ILL Error: seq/qlt names differ.\n";
    errorSs[AS_GKP_ILL_SEQ_QLT_LEN_DIFFER ] = "# ILL Error: seq/qlt lengths differ.\n";
    errorSs[AS_GKP_ILL_SEQ_TOO_LONG       ] = "# ILL Error: seq longer than longer than gkpShortReadLength bases, truncating.\n";
    errorSs[AS_GKP_ILL_CANT_OPEN_INPUT    ] = "# ILL Error: Couldn't open Illumina file for reading.\n";
    errorSs[AS_GKP_ILL_BAD_QV             ] = "# ILL Alert: invalid QV in read.\n";

    errorSs[AS_GKP_LIB_ILLEGAL_MEAN_STDDEV] = "# LIB Alert: lllegal mean and standard deviation; reset to mean 3000, stddev 300.\n";
    errorSs[AS_GKP_LIB_INVALID_MEAN       ] = "# LIB Alert: invalid mean; reset mean to 10 * stddev.\n";
    errorSs[AS_GKP_LIB_INVALID_STDDEV     ] = "# LIB Alert: invalid stddev; reset stddev to 0.1 * mean.\n";
    errorSs[AS_GKP_LIB_STDDEV_TOO_BIG     ] = "# LIB Alert: stddev too big for mean; reset stddev to 0.1 * mean.\n";
    errorSs[AS_GKP_LIB_STDDEV_TOO_SMALL   ] = "# LIB Alert: suspicious mean and standard deviation; reset stddev to 0.10 * mean.\n";
    errorSs[AS_GKP_LIB_EXISTS             ] = "# LIB Alert: already exists; can't add it again.\n";
    errorSs[AS_GKP_LIB_ZERO_UID           ] = "# LIB Error: zero or no UID; can't add it.\n";
    errorSs[AS_GKP_LIB_DOESNT_EXIST_UPDATE] = "# LIB Error: does not exist, can't update it.\n";
    errorSs[AS_GKP_LIB_UNKNOWN_ACTION     ] = "# LIB Error: invalid action.\n";

    errorSs[AS_GKP_LKG_SELF_LINK          ] = "# LKG Error: link from fragment to itself.\n";
    errorSs[AS_GKP_LKG_UNSUPPORTED_TYPE   ] = "# LKG Error: Unsupported link type.\n";
    errorSs[AS_GKP_LKG_FRG_DOESNT_EXIST   ] = "# LKG Error: Fragment not previously defined.\n";
    errorSs[AS_GKP_LKG_FRG_DELETED        ] = "# LKG Error: Fragment is marked as deleted.\n";

    errorSs[AS_GKP_LKG_ALREADY_MATED      ] = "# LKG Error: Fragment already mated.\n";

    errorSs[AS_GKP_LKG_LIB_DOESNT_EXIST   ] = "# LKG Error: Library not previously defined.\n";
    errorSs[AS_GKP_LKG_DIFFERENT_LIB      ] = "# LKG Error: Fragments to mate are in different libraries.\n";
    errorSs[AS_GKP_LKG_UNMATED_LIB        ] = "# LKG Error: Fragments to mate are in an unmated library.\n";
    errorSs[AS_GKP_LKG_DIFFERENT_ORIENT   ] = "# LKG Error: Fragments to mate are not supported in library of different mate type.\n";
    errorSs[AS_GKP_LKG_UNKNOWN_ACTION     ] = "# LKG Error: Unknown action.\n";

    errorSs[AS_GKP_SFF_UID_ERROR          ] = "# SFF Error: 454 Universal Accession Number out of range.\n";
    errorSs[AS_GKP_SFF_ALREADY_EXISTS     ] = "# SFF Error: Fragment exists, can't add it again.\n";

    errorSs[AS_GKP_SFF_TOO_SHORT          ] = "# SFF Alert: Fragment too short.\n";
    errorSs[AS_GKP_SFF_TOO_LONG           ] = "# SFF Alert: Fragment trimmed to assembler maximum length.\n";
    errorSs[AS_GKP_SFF_N                  ] = "# SFF Alert: Fragment contains an N; deleted.\n";

    errorSs[AS_GKP_PLC_SAME_CONSTRAINT    ] = "# PLC Error: Fragment constrained by the same fragment on both sides.\n";
    errorSs[AS_GKP_PLC_SELF_CONSTRAINT    ] = "# PLC Error: Fragment constrained by itself.\n";
    errorSs[AS_GKP_PLC_FRG_DOESNT_EXIST   ] = "# PLC Error: Fragment not previously defined.\n";
    errorSs[AS_GKP_PLC_FRG_DELETED        ] = "# PLC Error: Fragment is marked as deleted.\n";
    errorSs[AS_GKP_PLC_ALREADY_CONSTRAINED] = "# PLC Error: Fragment is already constrained.\n";
    errorSs[AS_GKP_PLC_UNKNOWN_ACTION     ] = "# PLC Error: invalid action.\n";

    errorSs[AS_GKP_UNKNOWN_MESSAGE        ] = "# GKP Error: Unknown message.\n";
  }

  va_start(ap, libIID);

  vfprintf(errorFP, errorMs[error], ap);

  errorCs[error]++;

  if (libError.size() < libIID)
    libError.resize(libIID + 1);

  libError[libIID].errorCs[error]++;

  va_end(ap);
}


bool
AS_GKP_summarizeErrors(char *gkpStoreName) {
  gkStore      *gkp    = new gkStore(gkpStoreName, FALSE, FALSE, TRUE);
  gkStoreStats *stats  = new gkStoreStats(gkp);

  //  Report the stats -- this is the same as 'gatekeeper -dumpinfo'.

  char  N[FILENAME_MAX];

  sprintf(N, "%s.info", gkpStoreName);

  FILE *F = fopen(N, "w");
  if (errno)
    fprintf(stderr, "WARNING: failed to open '%s' for write: %s\n", N, strerror(errno));

  if (F) {
    fprintf(F, "libIID\tbgnIID\tendIID\tactive\tdeleted\tmated\ttotLen\tclrLen\tlibName\n");

    fprintf(F, "0\t%d\t%d\t"F_U32"\t"F_U32"\t"F_U32"\t"F_U64"\t"F_U64"\tGLOBAL\n",
            1, stats->numActiveFrag + stats->numDeletedFrag, stats->numActiveFrag, stats->numDeletedFrag, stats->numMatedFrag, stats->readLength, stats->clearLength);

    for (uint32 j=0; j<gkp->gkStore_getNumLibraries() + 1; j++)
      fprintf(F, F_U32"\t"F_U32"\t"F_U32"\t"F_U32"\t"F_U32"\t"F_U32"\t"F_U64"\t"F_U64"\t%s\n",
              j,
              stats->lowestIID[j],
              stats->highestIID[j],
              stats->numActivePerLib[j],
              stats->numDeletedPerLib[j],
              stats->numMatedPerLib[j],
              stats->readLengthPerLib[j],
              stats->clearLengthPerLib[j],
              (j == 0 ? "LegacyUnmatedReads" : gkp->gkStore_getLibrary(j)->libraryName));

    fclose(F);
  }

  //  Deal with errors.

  uint32 nerrs = 0;

  for (uint32 i=0; i<AS_GKP_NUM_ERRORS; i++)
    nerrs += errorCs[i];

  if (nerrs == 0) {
    fprintf(stderr, "GKP finished with no alerts or errors.\n\n");
    return(0);
  }

  fprintf(stderr, "GKP finished with %d alerts or errors:\n", nerrs);

  for (uint32 i=0; i<AS_GKP_NUM_ERRORS; i++)
    if (errorCs[i])
      fprintf(stderr, "%d\t%s", errorCs[i], errorSs[i]);

  fprintf(stderr, "\n");

  //  Gatekeeper never dies.  We're always successful.  Maybe we should die only on very serious
  //  errors.

  bool   fatal = false;

  //  This particular error indicates either a library record didn't load properly, or your input
  //  frags might have a serious problem.
  //
  if ((errorCs[AS_GKP_FRG_UNKNOWN_LIB]))
    fatal = true;

  //  If any of the per-library errors are set, fail if they are a significant fraction of the
  //  library.  Currently, only FASTQ libraries set these.  Currently, we do not know the number of
  //  reads per library without scanning the entire store.
  //
  for (uint32 ll=1; ll < libError.size(); ll++) {
    uint32  nf = stats->numActivePerLib[ll] + stats->numDeletedPerLib[ll];
    uint32  ne = 0;

    for (uint32 i=0; i < AS_GKP_NUM_ERRORS; i++)
      ne += libError[ll].errorCs[i];

    if (ne >= nf * 0.05) {
      fprintf(stderr, "ERROR: library IID %u '%s' has %.2f%% errors or warnings.\n", ll, gkp->gkStore_getLibrary(ll)->libraryName, 100.0 * ne / nf);
      fatal = true;
    } else if (ne > 0) {
      fprintf(stderr, "WARNING: library IID %u '%s' has %.2f%% errors or warnings.\n", ll, gkp->gkStore_getLibrary(ll)->libraryName, 100.0 * ne / nf);
    }
  }

  delete stats;
  delete gkp;

  return(fatal);
}
