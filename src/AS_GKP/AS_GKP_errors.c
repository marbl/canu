
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

static char const *rcsid = "$Id: AS_GKP_errors.c,v 1.3 2007-11-08 12:38:12 brianwalenz Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

#include "AS_global.h"
#include "AS_GKP_include.h"

//  A rather unwieldy error handling system.  We want to count the
//  number of times we see each error, and try to make it somewhat
//  easy to add a new error.

#define AS_GKP_NUM_ERRORS  128

static char   *errorMs[AS_GKP_NUM_ERRORS] = {0};
static char   *errorSs[AS_GKP_NUM_ERRORS] = {0};
static uint32  errorCs[AS_GKP_NUM_ERRORS] = {0};


void
AS_GKP_reportError(int error, ...) {
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
    errorMs[AS_GKP_FRG_UNKNOWN_LIB        ] = "# FRG Error: Fragment %s references unknown library %s.\n";
    errorMs[AS_GKP_FRG_LOADED_DELETED     ] = "# FRG Alert: Fragment %s loaded, but marked as deleted due to errors previously reported.\n";
    errorMs[AS_GKP_FRG_DOESNT_EXIST       ] = "# FRG Error: Fragment %s does not exist, can't delete it.\n";
    errorMs[AS_GKP_FRG_HAS_MATE           ] = "# FRG Error: Fragment %s has mate pair relationship, can't delete it.\n";
    errorMs[AS_GKP_FRG_UNKNOWN_ACTION     ] = "# FRG Error: invalid action %c.\n";

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
    errorMs[AS_GKP_LKG_UNKNOWN_ACTION     ] = "# LKG Error: Unknown action %c.\n";

    errorMs[AS_GKP_SFF_UID_ERROR          ] = "# SFF Error: 454 Universal Accession Number '%s' out of range.\n";
    errorMs[AS_GKP_SFF_ALREADY_EXISTS     ] = "# SFF Error: Fragment %s exists, can't add it again.\n";

    errorMs[AS_GKP_UNKNOWN_MESSAGE        ] = "# GKP Error: Unknown message with type %s.\n";

    //
    //  A short description of the error, for summarizing the errors
    //  at the end of a run

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
    errorSs[AS_GKP_FRG_UNKNOWN_LIB        ] = "# FRG Error: references unknown library.\n";
    errorSs[AS_GKP_FRG_LOADED_DELETED     ] = "# FRG Alert: loaded, but marked as deleted due to errors previously reported.\n";
    errorSs[AS_GKP_FRG_DOESNT_EXIST       ] = "# FRG Error: does not exist, can't delete.\n";
    errorSs[AS_GKP_FRG_HAS_MATE           ] = "# FRG Error: has mate pair relationship, can't delete.\n";
    errorSs[AS_GKP_FRG_UNKNOWN_ACTION     ] = "# FRG Error: invalid action.\n";

    errorSs[AS_GKP_LIB_ILLEGAL_MEAN_STDDEV] = "# LIB Alert: lllegal mean and standard deviation; reset to mean 3000, stddev 300.\n";
    errorSs[AS_GKP_LIB_INVALID_MEAN       ] = "# LIB Alert: invalid mean; reset mean to 10 * stddev.\n";
    errorSs[AS_GKP_LIB_INVALID_STDDEV     ] = "# LIB Alert: invalid stddev; reset stddev to 0.1 * mean.\n";
    errorSs[AS_GKP_LIB_STDDEV_TOO_BIG     ] = "# LIB Alert: stddev too big for mean; reset stddev to 0.1 * mean.\n";
    errorSs[AS_GKP_LIB_STDDEV_TOO_SMALL   ] = "# LIB Alert: suspicious mean and standard deviation; reset stddev to 0.10 * mean.\n";
    errorSs[AS_GKP_LIB_EXISTS             ] = "# LIB Error: already exists; can't add it again.\n";
    errorSs[AS_GKP_LIB_ZERO_UID           ] = "# LIB Error: zero or no UID; can't add it.\n";
    errorSs[AS_GKP_LIB_DOESNT_EXIST_UPDATE] = "# LIB Error: does not exist, can't update it.\n";
    errorSs[AS_GKP_LIB_UNKNOWN_ACTION     ] = "# LIB Error: invalid action.\n";

    errorSs[AS_GKP_LKG_SELF_LINK          ] = "# LKG Error: link from fragment to itself.\n";
    errorSs[AS_GKP_LKG_UNSUPPORTED_TYPE   ] = "# LKG Error: Unsupported link type.\n";
    errorSs[AS_GKP_LKG_FRG_DOESNT_EXIST   ] = "# LKG Error: Fragment not previously defined.\n";
    errorMs[AS_GKP_LKG_FRG_DELETED        ] = "# LKG Error: Fragment is marked as deleted.\n";

    errorSs[AS_GKP_LKG_ALREADY_MATED      ] = "# LKG Error: Fragment already mated.\n";

    errorSs[AS_GKP_LKG_LIB_DOESNT_EXIST   ] = "# LKG Error: Library not previously defined.\n";
    errorSs[AS_GKP_LKG_DIFFERENT_LIB      ] = "# LKG Error: Fragments to mate are in different libraries.\n";
    errorSs[AS_GKP_LKG_UNKNOWN_ACTION     ] = "# LKG Error: Unknown action.\n";

    errorSs[AS_GKP_SFF_UID_ERROR          ] = "# SFF Error: 454 Universal Accession Number out of range.\n";
    errorSs[AS_GKP_SFF_ALREADY_EXISTS     ] = "# SFF Error: Fragment exists, can't add it again.\n";


    errorSs[AS_GKP_UNKNOWN_MESSAGE        ] = "# GKP Error: Unknown message.\n";
  }

  va_start(ap, error);

  vfprintf(errorFP, errorMs[error], ap);
  errorCs[error]++;

  va_end(ap);
}


int
AS_GKP_summarizeErrors(void) {
  int nerrs = 0;
  int i;

  for (i=0; i<AS_GKP_NUM_ERRORS; i++)
    nerrs += errorCs[i];

  if (nerrs) {
    fprintf(stderr, "GKP finished with %d errors:\n", nerrs);
    for (i=0; i<AS_GKP_NUM_ERRORS; i++)
      if (errorCs[i])
        fprintf(stderr, "%d\t%s", errorCs[i], errorSs[i]);
  } else {
    fprintf(stderr, "GKP finished with no errors.\n");
  }

  //  Gatekeeper never dies.  We're always successful.  Maybe we
  //  should die only on very serious errors.
  //
  //  This particular error indicates either a library record didn't
  //  load properly, or your input frags might have a serious problem.
  //
  int   fatal = 0;

  if ((errorCs[AS_GKP_FRG_UNKNOWN_LIB]))
    fatal = 1;

  return(fatal);
}
