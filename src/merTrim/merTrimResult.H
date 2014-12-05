
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 2010, J. Craig Venter Institute
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

static const char *rcsid_MERTRIMRESULT_H = "$Id$";

#include <stdio.h>
#include <stdlib.h>

#include "AS_global.H"

#ifndef MERTRIMRESULT_H
#define MERTRIMRESULT_H


class mertrimResult {
public:
  mertrimResult() {
    readIID = 0;
    deleted = 0;
    clrBgn  = 0;
    clrEnd  = 0;
    chimer  = 0;
    chmBgn  = 0;
    chmEnd  = 0;
  };

  void   print(FILE *F) {
    if (F == NULL)
      return;
    if (chimer)
      fprintf(F, F_IID"\t"F_U32"\t"F_U32"\tchimer\t"F_U32"\t"F_U32"%s\n",
              readIID, clrBgn, clrEnd, chmBgn, chmEnd,
              (deleted) ? "\tdeleted" : "");
    else
      fprintf(F, F_IID"\t"F_U32"\t"F_U32"%s\n",
              readIID, clrBgn, clrEnd,
              (deleted) ? "\tdeleted" : "");
  };

  void   writeResult(FILE *W) {
    if (W == NULL)
      return;
    AS_UTL_safeWrite(W, this, "merTrimResult", sizeof(mertrimResult), 1);
  };

  bool   readResult(FILE *R) {
    if (R == NULL)
      return(false);
    if (!feof(R))
      AS_UTL_safeRead(R, this, "merTrimResult", sizeof(mertrimResult), 1);
    return(feof(R) == false);
  };

  AS_IID readIID;
  uint32 deleted;
  uint32 clrBgn;
  uint32 clrEnd;
  uint32 chimer;
  uint32 chmBgn;
  uint32 chmEnd;
};


#endif  //  MERTRIMRESULT_H
