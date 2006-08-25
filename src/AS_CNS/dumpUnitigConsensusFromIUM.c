
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

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "AS_global.h"
#include "AS_MSG_pmesg.h"

int
main(int argc, char **argv) {

  if (argc != 1) {
    fprintf(stderr, "usage: %s < messages.cgi > unitigs.fasta\n");
    exit(1);
  }

  MesgReader   ReadMesg_AS = (MesgReader)InputFileType_AS(stdin);
  GenericMesg *pmesg;

  while (EOF != ReadMesg_AS(stdin, &pmesg)) {
    if (pmesg->t == MESG_IUM) {
      IntUnitigMesg  *ium_mesg  = (IntUnitigMesg *)pmesg->m;

      int   ungapped_unitig_length = 0;
      int   inew = 0;
      int   iold = 0;

      assert(ium_mesg->length == strlen(ium_mesg->consensus));
      assert(ium_mesg->consensus[ium_mesg->length] == '\0');

      /* De-gap the sequence in-place. */

      for(iold = 0; iold < ium_mesg->length; iold++) {
        char ch = ium_mesg->consensus[iold];
        if (ch != '-') {
          assert(ch == 'A' || ch == 'C' || ch == 'G' || ch == 'T' || ch == '\0');
          ium_mesg->consensus[inew] = ium_mesg->consensus[iold];
          inew++;
        }
      }
      ungapped_unitig_length = inew;
      ium_mesg->consensus[ungapped_unitig_length] = '\0';

      //  You could extend this to include the coverage stat
      //
      if ((ium_mesg->coverage_stat >= -1000000.0) &&
          (ium_mesg->num_frags > 1))
        fprintf(stdout,">unitig"F_IID" length=%d num_frags="F_IID" Astat=%.2f\n%s\n",
                ium_mesg->iaccession,
                ungapped_unitig_length,
                ium_mesg->num_frags,
                ium_mesg->coverage_stat,
                ium_mesg->consensus);
    }
  }

  return(0);
}
