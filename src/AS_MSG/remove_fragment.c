
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

const char *mainid = "$Id: remove_fragment.c,v 1.19 2009-11-24 21:42:11 brianwalenz Exp $";

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
  char           *inputName    = NULL;
  char           *strippedName = NULL;
  char           *removedName  = NULL;
  char           *UIDname      = NULL;

  FILE           *I = NULL;
  FILE           *S = NULL;
  FILE           *R = NULL;

  char            line[1024] = {0};

  HashTable_AS    *frag_uids_desired = CreateScalarHashTable_AS();
  HashTable_AS    *frag_uids_found   = CreateScalarHashTable_AS();

  GenericMesg    *gen = NULL;
  FragMesg       *frg = NULL;
  LinkMesg       *lkg = NULL;

  uint32          frgRemoved = 0;
  uint32          lkgRemoved = 0;

  argc = AS_configure(argc, argv);

  int arg = 1;
  int err = 0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-f") == 0) {
      UIDname = argv[++arg];

    } else if (strcmp(argv[arg], "-i") == 0) {
      inputName = argv[++arg];

    } else if (strcmp(argv[arg], "-o") == 0) {
      strippedName = argv[++arg];

    } else if (strcmp(argv[arg], "-O") == 0) {
      removedName = argv[++arg];

    } else {
      fprintf(stderr, "%s: unknown option '%s'\n", argv[0], argv[arg]);
      err++;
    }

    arg++;
  }
  if ((err) ||
      (inputName == NULL) ||
      ((strippedName == NULL) && (removedName == NULL))) {
    fprintf(stderr, "usage: %s -f UIDfile -i original.frg [-o stripped.frg] [-O removed.frg]\n", argv[0]);
    fprintf(stderr, "  -f UIDfile       one uid per line\n");
    fprintf(stderr, "  -i original.frg  fragments input\n");
    fprintf(stderr, "  -o stripped.frg  fragments NOT listed in UIDfile are saved here\n");
    fprintf(stderr, "  -O removed.frg   fragments     listed in UIDfile are saved here\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Removes fragments and links of supplied UIDs from the input fragment file.\n");
    exit(1);
  }

  errno = 0;
  FILE *F = fopen(UIDname, "r");
  if (errno) {
    fprintf(stderr, "%s: failed to open '%s' for reading UIDs: %s\n", argv[0], UIDname, strerror(errno));
    exit(1);
  }
  while (fgets(line, 1024, F)) {
    AS_UID  uid = AS_UID_load(line);
    InsertInHashTable_AS(frag_uids_desired, AS_UID_toInteger(uid), 0, 0, 0);
  }
  fclose(F);

  if (inputName) {
    I = fopen(inputName, "r");
    if (errno) {
      fprintf(stderr, "%s: failed to open input '%s' for reading: %s\n", argv[0], inputName, strerror(errno));
      exit(1);
    }
  }

  if (strippedName) {
    S = fopen(strippedName, "w");
    if (errno) {
      fprintf(stderr, "%s: failed to open output '%s' for writing: %s\n", argv[0], strippedName, strerror(errno));
      exit(1);
    }
  }

  if (removedName) {
    R = fopen(removedName, "w");
    if (errno) {
      fprintf(stderr, "%s: failed to open output '%s' for writing: %s\n", argv[0], removedName, strerror(errno));
      exit(1);
    }
  }


  while (ReadProtoMesg_AS(I, &gen) != EOF) {
    if (gen->t == MESG_FRG) {
      frg = (FragMesg *)gen->m;
      if (ExistsInHashTable_AS(frag_uids_desired, AS_UID_toInteger(frg->eaccession), 0)) {
        frgRemoved++;
        InsertInHashTable_AS(frag_uids_found, AS_UID_toInteger(frg->eaccession), 0, 0, 0);
        if (R)
          WriteProtoMesg_AS(R, gen);
      } else {
        if (S)
          WriteProtoMesg_AS(S, gen);
      }
    } else if (gen->t == MESG_LKG) {
      lkg = (LinkMesg *)gen->m;
      if (ExistsInHashTable_AS(frag_uids_found, AS_UID_toInteger(lkg->frag1), 0) &&
          ExistsInHashTable_AS(frag_uids_found, AS_UID_toInteger(lkg->frag2), 0)) {
        lkgRemoved++;
        if (R)
          WriteProtoMesg_AS(R, gen);
      } else {
        if (S)
          WriteProtoMesg_AS(S, gen);
      }
    } else {
      if (R)
        WriteProtoMesg_AS(R, gen);
      if (S)
        WriteProtoMesg_AS(S, gen);
    }
  }

  fclose(I);
  if (R)
    fclose(R);
  if (S)
    fclose(S);

  fprintf(stderr, "Removed %u fragments and %u links.\n", frgRemoved, lkgRemoved);

  return(0);
}
