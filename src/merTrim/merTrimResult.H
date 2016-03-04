
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
 *    src/AS_MER/merTrimResult.H
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2011-AUG-22 to 2013-AUG-01
 *      are Copyright 2011,2013 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz on 2014-DEC-05
 *      are Copyright 2014 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2016-JAN-11
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

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
