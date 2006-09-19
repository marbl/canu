
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

#include "AS_global.h"

int
main(int argc, char **argv) {
  MesgReader        reader;
  GenericMesg      *pmesg;
  SnapConConMesg   *contig;
  char              status[265];

  if (isatty(fileno(stdin)) || (argc > 1)) {
    fprintf(stderr,"Usage: %s < asmfile > contigs_fasta_file\n", argv[0]);
    exit(1);
  }
  
  reader = (MesgReader)InputFileType_AS(stdin);

  while (reader(stdin, &pmesg) != EOF) {
    if (pmesg->t ==MESG_CCO)  {
      contig = pmesg->m;
      int src, dst;

      assert(strlen(contig->consensus) == contig->length);

      src = 0;
      dst = 0;
      while (src < contig->length) {
        if (contig->consensus[src] != '-') {
          if (src != dst)
            contig->consensus[dst] = contig->consensus[src];
          dst++;
        }
        src++;
      }

      contig->consensus[dst] = 0;

      strcpy(status, "???");
      if (contig->placed  == AS_PLACED)
        strcpy(status, "true");
      if (contig->placed  == AS_UNPLACED)
        strcpy(status, "false");

      printf(">contig"F_UID","F_IID" placed=%s\n%s\n",
             contig->eaccession,
             contig->iaccession,
             status,
             contig->consensus);
    }
  }

  return(0);
}
