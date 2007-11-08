
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
/* $Id: remove_fragment.c,v 1.14 2007-11-08 12:38:13 brianwalenz Exp $ */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <ctype.h>

#include "AS_global.h"
#include "AS_MSG_pmesg.h"
#include "AS_UTL_Hash.h"

int
main(int argc, char ** argv) {

  char            line[1024] = {0};

  HashTable_AS    *frag_uids_desired = CreateScalarHashTable_AS(1024);
  HashTable_AS    *frag_uids_found   = CreateScalarHashTable_AS(1024);

  GenericMesg    *gen = NULL;
  FragMesg       *frg = NULL;
  LinkMesg       *lkg = NULL;

  int arg = 1;
  int err = 0;
  while (arg < argc) {
    if (strcmp(argv[arg], "-f") == 0) {
      errno = 0;
      FILE *F = fopen(argv[++arg], "r");
      if (errno) {
      }
      while (fgets(line, 1024, F)) {
        AS_UID  uid = AS_UID_lookup(line, NULL);
        InsertInHashTable_AS(frag_uids_desired, AS_UID_toInteger(uid), 0, 0, 0);
      }
      fclose(F);
    } else {
      AS_UID  uid = AS_UID_lookup(argv[arg], NULL);
      InsertInHashTable_AS(frag_uids_desired, AS_UID_toInteger(uid), 0, 0, 0);
    }
    arg++;
  }
  if (err) {
    fprintf(stderr, "usage: %s [-f UIDfile] UID UID .... < in > out\n", argv[0]);
    fprintf(stderr, "  -f UIDfile    one uid per line\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Strips messages (FRG, LKG only) from the input file.\n");
    exit(1);
  }

  
  while (ReadProtoMesg_AS(stdin, &gen) != EOF) {
    if (gen->t == MESG_FRG) {
      frg = (FragMesg *)gen->m;
      if (ExistsInHashTable_AS(frag_uids_desired, AS_UID_toInteger(frg->eaccession), 0)) {
        InsertInHashTable_AS(frag_uids_found, AS_UID_toInteger(frg->eaccession), 0, 0, 0);
        WriteProtoMesg_AS(stdout, gen);
      } else {
        //WriteProtoMesg_AS(stdout, gen);
      }
    } else if (gen->t == MESG_LKG) {
      lkg = (LinkMesg *)gen->m;
      if (ExistsInHashTable_AS(frag_uids_found, AS_UID_toInteger(lkg->frag1), 0) &&
          ExistsInHashTable_AS(frag_uids_found, AS_UID_toInteger(lkg->frag2), 0)) {
        WriteProtoMesg_AS(stdout, gen);
      } else {
        //WriteProtoMesg_AS(stdout, gen);
      }
    } else {
      WriteProtoMesg_AS(stdout, gen);
    }
  }

  return(0);
}
