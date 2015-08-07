
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

const char *mainid = "$Id: tgStoreDump.C 6974 2015-07-01 18:27:05Z walenzb $";

#include "AS_global.H"

#include "gkStore.H"
#include "tgStore.H"



void
operationBuild(char   *buildName,
               char   *tigName,
               uint32  tigVers) {

  errno = 0;
  FILE *F = fopen(buildName, "r");
  if (errno)
    fprintf(stderr, "Failed to open '%s' for reading: %s\n", buildName, strerror(errno)), exit(1);

  if (AS_UTL_fileExists(tigName, TRUE, TRUE)) {
    fprintf(stderr, "ERROR: '%s' exists, and I will not clobber an existing store.\n", tigName);
    exit(1);
  }

  tgStore *tigStore = new tgStore(tigName);
  tgTig    *tig      = new tgTig();

  for (int32 v=1; v<tigVers; v++)
    tigStore->nextVersion();

  while (tig->loadLayout(F) == true) {
    if (tig->numberOfChildren() == 0)
      continue;

    //  The log isn't correct.  For new tigs (all of these are) we don't know the
    //  id until after it is added.  Further, if these come with id's already set,
    //  they can't be added to a new store -- they don't exist.

#if 0
    fprintf(stderr, "INSERTING tig %d (%d children) (originally ID %d)\n",
            tig->tigID(), tig->numberOfChildren(), oID);
#endif

    tigStore->insertTig(tig, false);
  }

  fclose(F);

  delete tig;
  delete tigStore;
}





int
main (int argc, char **argv) {
  char         *gkpName        = NULL;
  char         *tigName        = NULL;
  int           tigVers        = -1;

  argc = AS_configure(argc, argv);

  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-G") == 0) {
      gkpName = argv[++arg];

    } else if (strcmp(argv[arg], "-T") == 0) {
      tigName = argv[++arg];
      tigVers = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-r") == 0) {
      AS_UTL_decodeRange(argv[++arg], tigIDbgn, tigIDend);
      tigIDset    = true;



    } else if (strcmp(argv[arg], "-B") == 0) {
      opType = OPERATION_BUILD;
      buildName = argv[++arg];


    } else if (strcmp(argv[arg], "-R") == 0) {
      opType = OPERATION_REPLACE;
      replaceName = argv[++arg];

    } else if (strcmp(argv[arg], "-N") == 0) {
      sameVersion = false;



    } else {
      fprintf(stderr, "%s: Unknown option '%s'\n", argv[0], argv[arg]);
      err++;
    }

    arg++;
  }
  if ((err) || (gkpName == NULL) || (tigName == NULL)) {
    fprintf(stderr, "usage: %s -G <gkpStore> -T <tigStore> <v> [opts]\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "  -G <gkpStore>         Path to the gatekeeper store\n");
    fprintf(stderr, "  -T <tigStore> <v>     Path to the tigStore, version, to use\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -r <file>             Load the tig(s) in 'file', replacing whatever is in the store\n");

    fprintf(stderr, "  -B <layout-file>      Construct a new store with unitigs in 'layout-file'.  Store versions\n");
    fprintf(stderr, "                        before that specified on the '-t' option are created but are empty.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -R <layout>           Replace a multialign with this one (type and id are from the layout)\n");
    fprintf(stderr, "                        The multialign is replaced in version <v> from -t.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -N                    Replace a multialign in the next version of the store.  This option is\n");
    fprintf(stderr, "                        needed if the version of the store to add a multialign does not exist.\n");
    fprintf(stderr, "                        The multialign is replaced in version <v>+1 from -t.\n");
    fprintf(stderr, "\n");

    exit(1);
  }

  //  To add a new multialign: In the layout, assign an id of -1 to the multialign (e.g., "unitig
  //  -1" or "contig -1").  Use -R to 'replace' this unitig in the store.  The store will assign the
  //  next unitig/contig ID to this new multialign.  WARNING!  The new multialign MUST be added to
  //  the latest version.
  //
  //  To delete a multialign: Remove ALL FRG and UTG lines, and set data.num_frags and
  //  data.num_unitigs to zero.  Use -R to 'replae' this unitig in the store.
  //  EXCEPT the code below will ignore treat these as EOF.

  if ((opType == OPERATION_BUILD) && (buildName != NULL)) {
    operationBuild(buildName, tigName, tigVers);
    exit(0);
  }



  gkStore *gkpStore = new gkStore(gkpName);
  tgStore *tigStore = new tgStore(tigName, tigVers);


  if ((opType == OPERATION_REPLACE) && (replaceName != NULL)) {
    if (tigIDset) {
      fprintf(stderr, "ERROR:  -R is incompatible with -c, -u, -C and -U.  Did you mean -cp or -up instead?\n");
      exit(1);
    }

    errno = 0;
    FILE         *F = fopen(replaceName, "r");
    if (errno)
      fprintf(stderr, "Failed to open '%s': %s\n", replaceName, strerror(errno)), exit(1);

    fprintf(stderr, "Reading layouts from '%s'.\n", replaceName);

    delete tigStore;

    if (sameVersion)
      tigStore = new tgStore(tigName, tigVers, true, true, false);  //  default
    else
      tigStore = new tgStore(tigName, tigVers, true, false, append);

    tgTig  *tig = new tgTig;

    while (tig->loadLayout(F) == true) {
      if (tig->numberOfChildren() == 0) {
        if (tigStore->isDeleted(tig->tigID()) == true) {
          fprintf(stderr, "DELETING tig %d -- ALREADY DELETED\n", tig->tigID());
        } else {
          fprintf(stderr, "DELETING tig %d\n", tig->tigID());
          tigStore->deleteTig(tig->tigID());
        }
      } else {
        tigStore->insertTig(tig, false);
        fprintf(stderr, "INSERTING tig %d\n", tig->tigID());
      }
    }

    fprintf(stderr, "Reading layouts from '%s' completed.\n", replaceName);

    delete tig;

    fclose(F);
  }

  delete gkpStore;
  delete tigStore;

  exit(0);
}
