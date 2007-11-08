
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
#include <string.h>

#include "AS_global.h"
#include "AS_MSG_pmesg.h"

void
usage(char *name) {
  fprintf(stderr, "usage: %s [options] < messages.cgi > unitigs.fasta\n", name);
  fprintf(stderr, "  -n <numfrags>   print only unitigs with at least numfrags reads\n");
  fprintf(stderr, "                  in them.  Default is 0 (dump all unitigs).\n");
}


int
main(int argc, char **argv) {
  int  arg      = 1;
  int  error    = 0;
  int  numfrags = 0;

  while (arg < argc) {
    if        (strcmp(argv[arg], "-n") == 0) {
      numfrags = atoi(argv[++arg]);
    } else if (strcmp(argv[arg], "-h") == 0) {
      error = 1;
    } else {
      fprintf(stderr, "Invalid option '%s'\n", argv[arg]);
      error = 1;
    }
    arg++;
  }

  if (error) {
    usage(argv[0]);
    exit(1);
  }

  GenericMesg *pmesg;

  while (EOF != ReadProtoMesg_AS(stdin, &pmesg)) {
    if (pmesg->t == MESG_UTG) {
      SnapUnitigMesg  *utg_mesg  = (SnapUnitigMesg *)pmesg->m;

      int   ungapped_unitig_length = 0;
      int   inew = 0;
      int   iold = 0;

      assert(utg_mesg->length == strlen(utg_mesg->consensus));
      assert(utg_mesg->consensus[utg_mesg->length] == '\0');

      /* De-gap the sequence in-place. */

      for(iold = 0; iold < utg_mesg->length; iold++) {
        char ch = utg_mesg->consensus[iold];
        if (ch != '-') {
          assert(ch == 'A' || ch == 'C' || ch == 'G' || ch == 'T' || ch == '\0');
          utg_mesg->consensus[inew] = utg_mesg->consensus[iold];
          inew++;
        }
      }
      ungapped_unitig_length = inew;
      utg_mesg->consensus[ungapped_unitig_length] = '\0';

      if (utg_mesg->num_frags > numfrags)
        fprintf(stdout,">unitig%s length=%d num_frags="F_IID" Astat=%.2f\n%s\n",
                AS_UID_toString(utg_mesg->eaccession),
                ungapped_unitig_length,
                utg_mesg->num_frags,
                utg_mesg->coverage_stat,
                utg_mesg->consensus);
    }

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

      if (ium_mesg->num_frags > numfrags)
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
