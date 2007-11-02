
/**************************************************************************
 * This file is based on 
 * ... the Celera Assembler, a software program that 
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 1999-2004, Applera Corporation. All rights reserved.
 * 
 * The Celera Assembler is free software; you can redistribute it and/or modify
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
#include <unistd.h>

#include "AS_global.h"
#include "AS_MSG_pmesg.h"

int
main(int argc, char **argv) {
  GenericMesg      *pmesg;
  SnapConConMesg   *contig;
  char              status[265];
  int               inclDegenerate = 0;
  int               onlyDegenerate = 0;

  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-d") == 0) {
      inclDegenerate = 1;
    } else if (strcmp(argv[arg], "-D") == 0) {
      inclDegenerate = 1;
      onlyDegenerate = 1;
    } else {
      fprintf(stderr, "unknown options '%s'\n", argv[arg]);
      err = 1;
    }
    arg++;
  }

  if (isatty(fileno(stdin)) || (err > 0)) {
    fprintf(stderr,"Usage: %s [-d | -D] < asmfile > contigs_fasta_file\n", argv[0]);
    fprintf(stderr, "  -d (-D)    also (only) print degenerate contigs\n");
    exit(1);
  }
  
  while (ReadProtoMesg_AS(stdin, &pmesg) != EOF) {
    if (pmesg->t ==MESG_CCO)  {
      contig = pmesg->m;

      //  By "definition", a degenerate contig has one unitig and is
      //  unplaced.
      //
      int isDeg = 0;
      if ((contig->num_unitigs == 1) &&
          (contig->placed      == AS_UNPLACED))
        isDeg = 1;

      if (((isDeg == 1) && (inclDegenerate == 0)) ||
          ((isDeg == 0) && (onlyDegenerate == 1)))
        continue;

      assert(strlen(contig->consensus) == contig->length);

      int src = 0;
      int dst = 0;
      while (src < contig->length) {
        if (contig->consensus[src] != '-') {
          if (src != dst)
            contig->consensus[dst] = contig->consensus[src];
          dst++;
        }
        src++;
      }

      contig->consensus[dst] = 0;

      if (contig->placed  == AS_PLACED)
        strcpy(status, "placed=true");
      else if (contig->placed  == AS_UNPLACED)
        strcpy(status, "placed=false");
      else
        strcpy(status, "placed=????");

      if (contig->num_unitigs == 1)
        strcat(status, " degenerate=true");
      else
        strcat(status, " degenerate=false");

      AS_UTL_writeFastA(stdout,
                        contig->consensus, strlen(contig->consensus),
                        ">contig"F_UID","F_IID" %s\n",
                        contig->eaccession, contig->iaccession, status);
    }
  }

  return(0);
}
