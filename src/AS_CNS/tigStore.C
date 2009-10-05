
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

const char *mainid = "$Id: tigStore.C,v 1.1 2009-10-05 22:49:42 brianwalenz Exp $";

#include "AS_global.h"
#include "MultiAlign.h"
#include "MultiAlignStore.h"
#include "MultiAlignment_CNS.h"
#include "MultiAlignment_CNS_private.h"

#define DUMP_PROPERTIES  0x01
#define DUMP_FRAGS       0x02
#define DUMP_UNITIGS     0x04
#define DUMP_CONSENSUS   0x08
#define DUMP_MULTIALIGN  0x16

#define DUMP_UNITIGLIST  1
#define DUMP_CONTIGLIST  2

int
main (int argc, char **argv) {
  char   tmpName[FILENAME_MAX] = {0};

  char  *gkpName = NULL;

  char  *tigName = NULL;
  int    tigVers = -1;
  int    tigPart = -1;

  int32  tigID       = 0;
  int32  tigIsUnitig = 0;

  uint32 dumpType    = 0;
  uint32 dumpFlags   = 0;

  argc = AS_configure(argc, argv);

  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-g") == 0) {
      gkpName = argv[++arg];

    } else if (strcmp(argv[arg], "-t") == 0) {
      tigName = argv[++arg];
      tigVers = atoi(argv[++arg]);
      tigPart = atoi(argv[++arg]);

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

      else if (strncmp(argv[arg], "multialign", 1) == 0)
        dumpFlags |= DUMP_MULTIALIGN;

      else
        fprintf(stderr, "%s: Unknown dump option '-d %s'\n", argv[0], argv[arg]);

    } else if (strcmp(argv[arg], "-D") == 0) {
      arg++;

      if      (strncmp(argv[arg], "unitiglist", 1) == 0)
        dumpType |= DUMP_UNITIGLIST;

      else if (strncmp(argv[arg], "contiglist", 1) == 0)
        dumpType |= DUMP_CONTIGLIST;

      else
        fprintf(stderr, "%s: Unknown dump option '-D %s'\n", argv[0], argv[arg]);

    } else {
      fprintf(stderr, "%s: Unknown option '%s'\n", argv[0], argv[arg]);
      err++;
    }

    arg++;
  }
  if ((err) || (gkpName == NULL) || (tigName == NULL)) {
    fprintf(stderr, "usage: %s -g gkpStore -t tigStore version partition [opts]\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "  -D <operation>      Dump something about the store\n");
    fprintf(stderr, "     unitiglist       ...a list of the unitigs in the store\n");
    fprintf(stderr, "     contiglist       ...a list of the contigs in the store\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -d <operation>      Dump something about a multialign in the store\n");
    fprintf(stderr, "     properties       ...properties\n");
    fprintf(stderr, "     frags            ...a list of fragments\n");
    fprintf(stderr, "     unitigs          ...a list of unitigs\n");
    fprintf(stderr, "     consensus        ...the consensus sequence\n");
    fprintf(stderr, "     multialign       ...the full multialignment\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -u id               Unitig to dump\n");
    fprintf(stderr, "  -c id               Contig to dump\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\n");
    exit(1);
  }

  gkpStore = new gkStore(gkpName, FALSE, FALSE);
  tigStore = new MultiAlignStore(tigName, tigVers, tigPart, 0, TRUE);

  fprintf(stderr, "tigVers=%d tigPart=%d\n", tigVers, tigPart);


  //
  //
  //

  if (dumpType == DUMP_UNITIGLIST) {
    tigStore->dumpMultiAlignRTable(true);
    exit(0);
  }

  if (dumpType == DUMP_CONTIGLIST) {
    tigStore->dumpMultiAlignRTable(false);
    exit(0);
  }

  //
  //
  //

  MultiAlignT *ma = tigStore->loadMultiAlign(tigID, tigIsUnitig);

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

  if (dumpFlags & DUMP_MULTIALIGN) {
    int showQV   = 0;
    int showDots = 1;

    PrintMultiAlignT(stdout, ma, gkpStore, showQV, showDots, AS_READ_CLEAR_OBTCHIMERA);
  }

  exit(0);
}
