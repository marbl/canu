
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

const char *mainid = "$Id: tigStore.C,v 1.2 2009-10-07 08:23:50 brianwalenz Exp $";

#include "AS_global.h"
#include "MultiAlign.h"
#include "MultiAlignStore.h"
#include "MultiAlignment_CNS.h"
#include "MultiAlignment_CNS_private.h"

#define DUMP_PROPERTIES  0x01
#define DUMP_FRAGS       0x02
#define DUMP_UNITIGS     0x04
#define DUMP_CONSENSUS   0x08
#define DUMP_LAYOUT      0x10
#define DUMP_MULTIALIGN  0x20

#define DUMP_UNITIGLIST  1
#define DUMP_CONTIGLIST  2
#define DUMP_EDIT        3

void
changeProperties(MultiAlignStore *tigStore,
                 char            *editName) {
  char  editLine[1024];

  errno = 0;
  FILE *editFile = fopen(editName, "r");
  if (errno)
    fprintf(stderr, "Failed to open '%s' for reading: %s\n", editName, strerror(errno)), exit(1);

  fgets(editLine, 1024, editFile);
  while (!feof(editFile)) {
    chomp(editLine);

    //  Decode the string into three pieces, the operation, the tig id and the value.  Note that
    //  these are just pointers to the start of the piece, and that the pieces are NOT
    //  zero-terminated (so we can print an error of the whole line if needed).
    //
    char   *op = editLine;
    char   *tp = editLine;
    char   *vp = editLine;

    //  Skip whitespace before the operation
    while (isspace(*op) == 1)
      op++;

    //  Skip the operation, then skip whitespace before the tig.
    tp = op;
    while (isspace(*tp) == 0)
      tp++;
    while (isspace(*tp) == 1)
      tp++;

    //  Skip the tig, then skip whitespace before the value.
    vp = tp;
    while (isspace(*vp) == 0)
      vp++;
    while (isspace(*vp) == 1)
      vp++;

    int32   tid = atoi(tp);

    if        (strncmp(op, "unitig_coverage_stat", 20) == 0) {
      tigStore->setUnitigCoverageStat(tid, atof(vp));

    } else if (strncmp(op, "unitig_microhet_prob", 20) == 0) {
      tigStore->setUnitigMicroHetProb(tid, atof(vp));

    } else if (strncmp(op, "unitig_status", 13) == 0) {
      UnitigStatus  st = tigStore->getUnitigStatus(tid);
      switch (*vp) {
        case AS_UNIQUE:      st = AS_UNIQUE;      break;
        case AS_NOTREZ:      st = AS_NOTREZ;      break;
        case AS_SEP:         st = AS_SEP;         break;
        case AS_UNASSIGNED:  st = AS_UNASSIGNED;  break;
        default:
          fprintf(stderr, "unknown unitig_status in '%s'\n", editLine);
          break;
      }
      tigStore->setUnitigStatus(tid, st);

    } else if (strncmp(op, "unitig_unique_rept", 18) == 0) {
      UnitigFUR ur = tigStore->getUnitigFUR(tid);
      switch (*vp) {
        case AS_FORCED_NONE:    ur = AS_FORCED_NONE;    break;
        case AS_FORCED_UNIQUE:  ur = AS_FORCED_UNIQUE;  break;
        case AS_FORCED_REPEAT:  ur = AS_FORCED_REPEAT;  break;
        default:
          fprintf(stderr, "unknown unitig_unique_rept in '%s'\n", editLine);
          break;
      }
      tigStore->setUnitigFUR(tid, ur);

    } else if (strncmp(op, "contig_status", 13) == 0) {
      ContigStatus  st = tigStore->getContigStatus(tid);
      switch (*vp) {
        case AS_PLACED:      st = AS_PLACED;    break;
        case AS_UNPLACED:    st = AS_UNPLACED;  break;
        default:
          fprintf(stderr, "unknown contig_status in '%s'\n", editLine);
          break;
      }
      tigStore->setContigStatus(tid, st);

    } else {
      fprintf(stderr, "unknown edit '%s'\n", editLine);
    }

    fgets(editLine, 1024, editFile);
  }
}


int
main (int argc, char **argv) {
  char          tmpName[FILENAME_MAX] = {0};
  char         *gkpName     = NULL;
  char         *tigName     = NULL;
  int           tigVers     = -1;
  int           tigPartU    = 0;
  int           tigPartC    = 0;
  int32         tigID       = 0;
  int32         tigIsUnitig = 0;
  uint32        dumpType    = 0;
  uint32        dumpFlags   = 0;
  char         *editName    = NULL;
  char         *replaceName = NULL;
  MultiAlignT  *ma          = NULL;

  argc = AS_configure(argc, argv);

  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-g") == 0) {
      gkpName = argv[++arg];

    } else if (strcmp(argv[arg], "-t") == 0) {
      tigName = argv[++arg];
      tigVers = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-up") == 0) {
      tigPartU = atoi(argv[++arg]);
    } else if (strcmp(argv[arg], "-cp") == 0) {
      tigPartC = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-u") == 0) {
      tigID       = atoi(argv[++arg]);
      tigIsUnitig = TRUE;

    } else if (strcmp(argv[arg], "-c") == 0) {
      tigID       = atoi(argv[++arg]);
      tigIsUnitig = FALSE;

    } else if (strcmp(argv[arg], "-d") == 0) {
      arg++;

      if      (strncmp(argv[arg], "properties", 1) == 0)
        dumpFlags |= DUMP_PROPERTIES;

      else if (strncmp(argv[arg], "frags", 1) == 0)
        dumpFlags |= DUMP_FRAGS;

      else if (strncmp(argv[arg], "unitigs", 1) == 0)
        dumpFlags |= DUMP_UNITIGS;

      else if (strncmp(argv[arg], "consensus", 1) == 0)
        dumpFlags |= DUMP_CONSENSUS;

      else if (strncmp(argv[arg], "layout", 1) == 0)
        dumpFlags |= DUMP_LAYOUT;

      else if (strncmp(argv[arg], "multialign", 1) == 0)
        dumpFlags |= DUMP_MULTIALIGN;

      else
        fprintf(stderr, "%s: Unknown dump option '-d %s'\n", argv[0], argv[arg]);

    } else if (strcmp(argv[arg], "-D") == 0) {
      arg++;

      if      (strncmp(argv[arg], "unitiglist", 1) == 0)
        dumpType = DUMP_UNITIGLIST;

      else if (strncmp(argv[arg], "contiglist", 1) == 0)
        dumpType = DUMP_CONTIGLIST;

      else if (strncmp(argv[arg], "properties", 1) == 0)
        dumpType = DUMP_EDIT;

      else
        fprintf(stderr, "%s: Unknown dump option '-D %s'\n", argv[0], argv[arg]);

    } else if (strcmp(argv[arg], "-E") == 0) {
      editName = argv[++arg];

    } else if (strcmp(argv[arg], "-R") == 0) {
      replaceName = argv[++arg];

    } else {
      fprintf(stderr, "%s: Unknown option '%s'\n", argv[0], argv[arg]);
      err++;
    }

    arg++;
  }
  if ((err) || (gkpName == NULL) || (tigName == NULL)) {
    fprintf(stderr, "usage: %s -g <gkpStore> -t <tigStore> <v> [opts]\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "  -g <gkpStore>           Path to the gatekeeper store\n");
    fprintf(stderr, "  -t <tigStore> <v>       Path to the tigStore, and version to use\n");
    fprintf(stderr, "  -up <p>                 ...limit to unitigs in partition <p>\n");
    fprintf(stderr, "  -cp <p>                 ...limit to contigs in partition <p>\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -D <operation>        Dump something about the store\n");
    fprintf(stderr, "     unitiglist         ...a list of the unitigs in the store\n");
    fprintf(stderr, "     contiglist         ...a list of the contigs in the store\n");
    fprintf(stderr, "     properties         ...a list of properties for ALL multialigns in the store (for -E)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -u id                 Unitig to dump (for -d option)\n");
    fprintf(stderr, "  -c id                 Contig to dump (for -d option)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -d <operation>        Dump something about a multialign (-c or -u) in the store\n");
    fprintf(stderr, "     properties         ...properties\n");
    fprintf(stderr, "     frags              ...a list of fragments\n");
    fprintf(stderr, "     unitigs            ...a list of unitigs\n");
    fprintf(stderr, "     consensus          ...the consensus sequence\n");
    fprintf(stderr, "     layout             ...the layout\n");
    fprintf(stderr, "     multialign         ...the full multialignment\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -E <editFile>         Change properties of multialigns\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -R <layout>           Replace a multialign with this one (type and id are from the layout)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\n");
    exit(1);
  }

  gkpStore = new gkStore(gkpName, FALSE, FALSE);
  tigStore = new MultiAlignStore(tigName, tigVers, tigPartU, tigPartC, FALSE, FALSE);

  //
  //
  //

  if (editName) {
    delete tigStore;
    tigStore = new MultiAlignStore(tigName, tigVers, tigPartU, tigPartC, TRUE, TRUE);
    changeProperties(tigStore, editName);
    goto finish;
  }

  if (replaceName) {
    errno = 0;
    FILE         *F = fopen(replaceName, "r");
    if (errno)
      fprintf(stderr, "Failed to open '%s': %s\n", replaceName, strerror(errno)), exit(1);

    MultiAlignT  *ma       = CreateEmptyMultiAlignT();
    bool          isUnitig = LoadMultiAlignFromHuman(ma, F);

    fclose(F);

    delete tigStore;
    tigStore = new MultiAlignStore(tigName, tigVers, tigPartU, tigPartC, TRUE, TRUE);

    tigStore->insertMultiAlign(ma, isUnitig, TRUE);

    goto finish;
  }

  //
  //
  //

  if (dumpType == DUMP_UNITIGLIST) {
    tigStore->dumpMultiAlignRTable(true);
    goto finish;
  }

  if (dumpType == DUMP_CONTIGLIST) {
    tigStore->dumpMultiAlignRTable(false);
    goto finish;
  }

  if (dumpType == DUMP_EDIT) {
    for (uint32 i=0; i<tigStore->numUnitigs(); i++) {
      if (tigStore->isDeleted(i, TRUE) == false) {
        fprintf(stdout, "unitig_coverage_stat %8u %f\n", i, tigStore->getUnitigCoverageStat(i));
        fprintf(stdout, "unitig_microhet_prob %8u %f\n", i, tigStore->getUnitigMicroHetProb(i));
        fprintf(stdout, "unitig_status        %8u %c\n", i, tigStore->getUnitigStatus(i));
        fprintf(stdout, "unitig_unique_rept   %8u %c\n", i, tigStore->getUnitigFUR(i));
      }
    }
    for (uint32 i=0; i<tigStore->numContigs(); i++) {
      if (tigStore->isDeleted(i, FALSE) == false) {
        fprintf(stdout, "contig_status        %8u %c\n", i, tigStore->getContigStatus(i));
      }
    }
    goto finish;
  }

  //
  //
  //

  ma = tigStore->loadMultiAlign(tigID, tigIsUnitig);

  if (ma == NULL) {
    fprintf(stderr, "MultiAlign not loaded.\n");
    exit(0);
  }

  if (dumpFlags & DUMP_PROPERTIES) {
    fprintf(stderr, "maID                "F_S32"\n", ma->maID);
    fprintf(stderr, "unitigCoverageStat  %f\n",      ma->data.unitig_coverage_stat);
    fprintf(stderr, "unitigMicrohetProb  %f\n",      ma->data.unitig_microhet_prob);
    fprintf(stderr, "unitigStatus        %c/%d\n",   ma->data.unitig_status, ma->data.unitig_status);
    fprintf(stderr, "unitigFUR           %c/%d\n",   ma->data.unitig_unique_rept, ma->data.unitig_unique_rept);
    fprintf(stderr, "contigStatus        %c/%d\n",   ma->data.contig_status, ma->data.contig_status);
    fprintf(stderr, "numFrags            "F_U32" (vs "F_U32")\n", ma->data.num_frags, GetNumIntMultiPoss(ma->f_list));
    fprintf(stderr, "numUnitigs          "F_U32" (vs "F_U32")\n", ma->data.num_unitigs, GetNumIntUnitigPoss(ma->u_list));

    tigStore->dumpMultiAlignR(tigID, tigIsUnitig);
  }

  if (dumpFlags & DUMP_FRAGS) {
    for (int32 i=0; i<GetNumIntMultiPoss(ma->f_list); i++) {
      IntMultiPos *imp = GetIntMultiPos(ma->f_list, i);

      fprintf(stdout, "FRG %7d %5d,%5d\n",
              imp->ident, imp->position.bgn, imp->position.end);
    }
  }

  if (dumpFlags & DUMP_UNITIGS) {
    for (int32 i=0; i<GetNumIntUnitigPoss(ma->u_list); i++) {
      IntUnitigPos *iup = GetIntUnitigPos(ma->u_list, i);

      fprintf(stdout, "UTG %7d %5d,%5d\n",
              iup->ident, iup->position.bgn, iup->position.end);
    }
  }

  if (dumpFlags & DUMP_CONSENSUS) {
    if (ma->consensus) {
      fprintf(stderr, "cns=%d (incl null)\n", GetNumchars(ma->consensus));
      fprintf(stderr, "cns=%s\n", Getchar(ma->consensus, 0));
    }
  }

  if (dumpFlags & DUMP_LAYOUT) {
    DumpMultiAlignForHuman(stdout, ma, tigIsUnitig);
  }

  if (dumpFlags & DUMP_MULTIALIGN) {
    int showQV   = 0;
    int showDots = 1;

    PrintMultiAlignT(stdout, ma, gkpStore, showQV, showDots, AS_READ_CLEAR_OBTCHIMERA);
  }

 finish:
  delete gkpStore;
  delete tigStore;

  exit(0);
}
