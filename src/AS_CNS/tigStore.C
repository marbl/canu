
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

const char *mainid = "$Id: tigStore.C,v 1.26 2013-01-31 16:38:21 brianwalenz Exp $";

#include "AS_global.h"
#include "AS_UTL_decodeRange.H"
#include "MultiAlign.h"
#include "MultiAlignStore.h"
#include "MultiAlignMatePairAnalysis.H"
#include "MultiAlignSizeAnalysis.H"
#include "MultiAlignment_CNS.h"
#include "MultiAlignment_CNS_private.h"

#define DUMP_PROPERTIES       1
#define DUMP_FRAGS            2
#define DUMP_UNITIGS          3
#define DUMP_CONSENSUS        4
#define DUMP_CONSENSUSGAPPED  5
#define DUMP_LAYOUT           6
#define DUMP_MULTIALIGN       7
#define DUMP_MATEPAIR         8
#define DUMP_SIZES            9

#define OPERATION_UNITIGLIST  1
#define OPERATION_CONTIGLIST  2
#define OPERATION_PROPERTIES  3
#define OPERATION_TIG         4
#define OPERATION_EDIT        5
#define OPERATION_REPLACE     6
#define OPERATION_BUILD       7
#define OPERATION_COMPRESS    8
#define OPERATION_FMAP        9

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
    while (isspace(*op))
      op++;

    //  Skip the operation, then skip whitespace before the tig.
    tp = op;
    while (!isspace(*tp))
      tp++;
    while (isspace(*tp))
      tp++;

    //  Skip the tig, then skip whitespace before the value.
    vp = tp;
    while (!isspace(*vp))
      vp++;
    while (isspace(*vp))
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



void
dumpProperties(MultiAlignStore *tigStore,
               int32 tigID,
               int32 tigIsUnitig,
               MultiAlignT *ma) {

  fprintf(stdout, "maID                "F_S32"\n", ma->maID);
  fprintf(stdout, "unitigCoverageStat  %f\n",      ma->data.unitig_coverage_stat);
  fprintf(stdout, "unitigMicrohetProb  %f\n",      ma->data.unitig_microhet_prob);
  fprintf(stdout, "unitigStatus        %c/%d\n",   ma->data.unitig_status, ma->data.unitig_status);
  fprintf(stdout, "unitigFUR           %c/%d\n",   ma->data.unitig_unique_rept, ma->data.unitig_unique_rept);
  fprintf(stdout, "contigStatus        %c/%d\n",   ma->data.contig_status, ma->data.contig_status);

  #if GCCONTENT
  float gcContent = 0.0;
  int ulen = 0;
  int glen = 0;

  if (ma->consensus)
  {
    char *cns = Getchar(ma->consensus, 0);

    if (cns && *cns)
    {
      int gcs = 0;

      while (*cns)
      {
        glen++;

        if (*cns != '-')
        {
          if (*cns == 'G' || *cns == 'C') { gcs++; }
          ulen++;
        }

        cns++;
      }

      if (ulen) { gcContent = (float) gcs / (float) ulen; }
    }
  }

  fprintf(stdout, "gcContent           %f\n",      gcContent);
  fprintf(stdout, "uLen                %d\n",      ulen);
  fprintf(stdout, "gLen                %d\n",      glen);
  #endif

  fprintf(stdout, "numFrags            "F_U32" (vs "F_U64")\n", ma->data.num_frags, (uint64)GetNumIntMultiPoss(ma->f_list));
  fprintf(stdout, "numUnitigs          "F_U32" (vs "F_U64")\n", ma->data.num_unitigs, (uint64)GetNumIntUnitigPoss(ma->u_list));

  tigStore->dumpMultiAlignR(tigID, tigIsUnitig);
}


void
dumpFrags(MultiAlignStore *tigStore,
          int32 tigID,
          int32 tigIsUnitig,
          MultiAlignT *ma) {

  for (uint32 i=0; i<GetNumIntMultiPoss(ma->f_list); i++) {
    IntMultiPos *imp = GetIntMultiPos(ma->f_list, i);

    fprintf(stdout, "FRG %7d %5d,%5d\n",
            imp->ident, imp->position.bgn, imp->position.end);
  }
}


void
dumpUnitigs(MultiAlignStore *tigStore,
            int32 tigID,
            int32 tigIsUnitig,
            MultiAlignT *ma) {

  for (uint32 i=0; i<GetNumIntUnitigPoss(ma->u_list); i++) {
    IntUnitigPos *iup = GetIntUnitigPos(ma->u_list, i);

    fprintf(stdout, "UTG %7d %5d,%5d\n",
            iup->ident, iup->position.bgn, iup->position.end);
  }
}


void
dumpConsensus(MultiAlignStore *tigStore,
              int32 tigID,
              int32 tigIsUnitig,
              MultiAlignT *ma,
              bool withGaps) {

  if (ma->consensus == NULL)
    return;

  char *cns = Getchar(ma->consensus, 0);

  if ((cns == NULL) || (cns[0] == 0))
    return;

  if (withGaps == false) {
    char *o = cns;
    char *n = cns;

    while (*n) {
      if (*n != '-')
        *o++ = *n;
      n++;
    }

    *o = 0;
  }

  if (tigIsUnitig)
    fprintf(stdout, ">utg%d len="F_U64" reads="F_U32" status=%c microHet=%.2f covStat=%.2f\n%s\n",
            ma->maID, GetNumchars(ma->consensus) - 1, ma->data.num_frags,
            ma->data.unitig_status,
            ma->data.unitig_microhet_prob,
            ma->data.unitig_coverage_stat,
            cns);
  else
    fprintf(stdout, ">ctg%d len="F_U64" reads="F_U32" unitigs="F_U32" status=%c\n%s\n",
            ma->maID, GetNumchars(ma->consensus) - 1, ma->data.num_frags,
            ma->data.num_unitigs,
            ma->data.contig_status,
            cns);
}



void
operationBuild(char *buildName, char *tigName,  int tigVers) {
  uint32  utgID = 0;
  uint32  ctgID = 0;
  uint32  orgID = 0;

  MultiAlignStore  *tigStore = NULL;
  MultiAlignT       *ma       = CreateEmptyMultiAlignT();
  bool               isUnitig = false;

  errno = 0;
  FILE         *F = fopen(buildName, "r");
  if (errno)
    fprintf(stderr, "Failed to open '%s': %s\n", buildName, strerror(errno)), exit(1);

  if (AS_UTL_fileExists(tigName, TRUE, TRUE)) {
    fprintf(stderr, "ERROR: '%s' exists, and I will not clobber an existing store.\n", tigName);
    exit(1);
    tigStore = new MultiAlignStore(tigName, tigVers, 0, 0, TRUE, TRUE, FALSE);
  } else {
    tigStore = new MultiAlignStore(tigName);

    for (int32 v=1; v<tigVers; v++)
      tigStore->nextVersion();
  }

  while (LoadMultiAlignFromHuman(ma, isUnitig, F) == true) {
    if (ma->data.num_frags + ma->data.num_unitigs == 0)
      continue;

    orgID = ma->maID;

    if (isUnitig)
      ma->maID = utgID;
    else
      ma->maID = ctgID;

#if 0
    fprintf(stderr, "INSERTING %s %d (%d frags %d unitigs) (originally ID %d)\n",
            (isUnitig) ? "unitig" : "contig",
            ma->maID,
            ma->data.num_frags, ma->data.num_unitigs,
            orgID);
#endif

    tigStore->insertMultiAlign(ma, isUnitig, FALSE);

    if (isUnitig)
      utgID++;
    else
      ctgID++;
  }

  fclose(F);

  delete tigStore;
}



void
operationCompress(char *tigName, int tigVers) {
  MultiAlignStore  *tigStore = new MultiAlignStore(tigName, tigVers, 0, 0, FALSE, FALSE, FALSE);
  bool              isUnitig = TRUE;

  uint32            nUtgErrors = 0;
  uint32            nCtgErrors = 0;

  uint32            nUtgCompress = 0;
  uint32            nCtgCompress = 0;

  //  Pass 0:  Fail if this isn't the latest version.  If we try to compress something that isn't the
  //  latest version, versions after this still point to the uncompressed tigs.

  //if (tigStore->

  //  Pass 1:  Check that we aren't going to pull a tig out of the future and place it in the past.

  for (uint32 ti=0; ti<tigStore->numUnitigs(); ti++) {
    if (tigStore->isDeleted(ti, isUnitig))
      continue;

    if (tigStore->getUnitigVersion(ti) > tigVers) {
      fprintf(stderr, "WARNING:  Attempt to move future unitig "F_U32" from version "F_U32" to previous version %d.\n",
              ti, tigStore->getUnitigVersion(ti), tigVers);
      nUtgErrors++;
    } else if (tigStore->getUnitigVersion(ti) < tigVers) {
      nUtgCompress++;
    }
  }

  isUnitig = FALSE;

  for (uint32 ti=0; ti<tigStore->numContigs(); ti++) {
    if (tigStore->isDeleted(ti, isUnitig))
      continue;

    if (tigStore->getContigVersion(ti) > tigVers) {
      fprintf(stderr, "WARNING:  Attempt to move future contig "F_U32" from version "F_U32" to previous version %d.\n",
              ti, tigStore->getContigVersion(ti), tigVers);
      nCtgErrors++;
    } else if (tigStore->getContigVersion(ti) < tigVers) {
      nCtgCompress++;
    }
  }

  if (nUtgErrors + nCtgErrors > 0) {
    fprintf(stderr, "Store can't be compressed; probably trying to compress to something that isn't the latest version.\n");
    fprintf(stderr, "  "F_U32" unitigs failed; "F_U32" compressable\n", nUtgErrors, nUtgCompress);
    fprintf(stderr, "  "F_U32" contigs failed; "F_U32" compressable\n", nCtgErrors, nCtgCompress);
    delete tigStore;
    exit(1);
  }

  //  Pass 2:  Actually do the moves

  if (nUtgCompress + nCtgCompress > 0) {
    delete tigStore;
    tigStore = new MultiAlignStore(tigName, tigVers, 0, 0, TRUE, TRUE, FALSE);
  }

  if (nUtgCompress > 0) {
    isUnitig = TRUE;

    fprintf(stderr, "Compressing "F_U32" unitigs into version %d\n", nUtgCompress, tigVers);

    for (uint32 ti=0; ti<tigStore->numUnitigs(); ti++) {
      if ((ti % 1000000) == 0)
        fprintf(stderr, "tig %d\n", ti);

      if (tigStore->isDeleted(ti, isUnitig)) {
        continue;
      }

      if (tigStore->getUnitigVersion(ti) == tigVers)
        continue;

      MultiAlignT *ma = tigStore->loadMultiAlign(ti, isUnitig);

      if (ma == NULL) {
        //fprintf(stderr, "WARNING: unitig "F_U32" is NULL.\n", ti);
        continue;
      }

      tigStore->insertMultiAlign(ma, isUnitig, TRUE);
      tigStore->unloadMultiAlign(ti, isUnitig);
    }
  }

  if (nCtgCompress > 0) {
    isUnitig = FALSE;

    fprintf(stderr, "Compressing "F_U32" contigs into version %d\n", nCtgCompress, tigVers);

    for (uint32 ti=0; ti<tigStore->numContigs(); ti++) {
      if ((ti % 1000000) == 0)
        fprintf(stderr, "tig %d\n", ti);

      if (tigStore->isDeleted(ti, isUnitig))
        continue;

      if (tigStore->getContigVersion(ti) == tigVers)
        continue;

      MultiAlignT *ma = tigStore->loadMultiAlign(ti, isUnitig);

      if (ma == NULL) {
        //fprintf(stderr, "WARNING: contig "F_U32" is NULL.\n", ti);
        continue;
      }

      tigStore->insertMultiAlign(ma, isUnitig, TRUE);
      tigStore->unloadMultiAlign(ti, isUnitig);
    }
  }

  //  Now clean up the older files

  if (nUtgCompress + nCtgCompress > 0) {
    for (uint32 version=1; version<tigVers; version++) {
      fprintf(stderr, "Purge version "F_U32".\n", version);
      tigStore->purgeVersion(version);
    }
  }

  //  And the newer files

  delete tigStore;
}


void
operationFmap(char   *tigName,
              int32   tigVers,
              uint32  minFmap,
              uint32  maxFmap,
              bool    tigIsUnitig,
              uint32  tigID,
              bool    dumpAll) {
  MultiAlignStore  *tigStore = new MultiAlignStore(tigName, tigVers, 0, 0, FALSE, FALSE, FALSE);

  uint32  bgn = tigID;
  uint32  end = tigID + 1;

  if (dumpAll == TRUE) {
    bgn = 0;
    end = (tigIsUnitig) ? tigStore->numUnitigs() : tigStore->numContigs();
  }

  for (int32 ti=bgn; ti<end; ti++) {
    if (tigStore->isDeleted(ti, tigIsUnitig))
      continue;

    if ((tigStore->getNumFrags(ti, tigIsUnitig) < minFmap) ||
        (maxFmap < tigStore->getNumFrags(ti, tigIsUnitig)))
      continue;

    MultiAlignT  *ma    = tigStore->loadMultiAlign(ti, tigIsUnitig);
    uint32        fiMax = GetNumIntMultiPoss(ma->f_list);

    for (uint32 fi=0; fi<fiMax; fi++) {
      IntMultiPos  *imp = GetIntMultiPos(ma->f_list, fi);

      fprintf(stdout, F_U32"\t"F_U32"\t"F_S32"\t"F_S32"\n",
              imp->ident, ti, imp->position.bgn, imp->position.end);
    }

    tigStore->unloadMultiAlign(ti, tigIsUnitig);
  }
}




int
main (int argc, char **argv) {
  char          tmpName[FILENAME_MAX] = {0};
  char         *gkpName        = NULL;
  char         *tigName        = NULL;
  int           tigVers        = -1;
  int           tigPartU       = 0;
  int           tigPartC       = 0;
  uint32        tigID          = 0;
  uint32        tigIDbgn       = 0;
  uint32        tigIDend       = UINT32_MAX;
  int32         tigIsUnitig    = TRUE;
  uint32        opType         = 0;
  uint32        dumpFlags      = 0;
  uint32        dumpAll        = 0;
  char         *editName       = NULL;
  char         *replaceName    = NULL;
  bool          replaceInPlace = TRUE;
  char         *buildName      = NULL;

  uint32        minFmap        = 0;
  uint32        maxFmap        = UINT32_MAX;

  MultiAlignT  *ma             = NULL;
  int           showQV         = 0;
  int           showDots       = 1;

  matePairAnalysis  *mpa       = NULL;
  char              *mpaPrefix = NULL;

  sizeAnalysis      *siz       = NULL;
  char              *sizPrefix = NULL;
  uint64             sizSize   = 0;

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
      AS_UTL_decodeRange(argv[++arg], tigIDbgn, tigIDend);
      tigID       = tigIDbgn;
      tigIsUnitig = TRUE;

    } else if (strcmp(argv[arg], "-U") == 0) {
      dumpAll     = TRUE;
      tigIsUnitig = TRUE;

    } else if (strcmp(argv[arg], "-c") == 0) {
      AS_UTL_decodeRange(argv[++arg], tigIDbgn, tigIDend);
      tigID       = tigIDbgn;
      tigIsUnitig = FALSE;

    } else if (strcmp(argv[arg], "-C") == 0) {
      dumpAll     = TRUE;
      tigIsUnitig = FALSE;

    } else if (strcmp(argv[arg], "-d") == 0) {
      arg++;

      opType = OPERATION_TIG;

      if      (strncmp(argv[arg], "properties", 1) == 0)
        dumpFlags = DUMP_PROPERTIES;

      else if (strncmp(argv[arg], "frags", 1) == 0)
        dumpFlags = DUMP_FRAGS;

      else if (strncmp(argv[arg], "unitigs", 1) == 0)
        dumpFlags = DUMP_UNITIGS;

      else if (strncmp(argv[arg], "consensus", 1) == 0)
        dumpFlags = DUMP_CONSENSUS;

      else if (strncmp(argv[arg], "consensusgapped", 1) == 0)
        dumpFlags = DUMP_CONSENSUSGAPPED;

      else if (strncmp(argv[arg], "layout", 1) == 0)
        dumpFlags = DUMP_LAYOUT;

      else if (strncmp(argv[arg], "multialign", 2) == 0)
        dumpFlags = DUMP_MULTIALIGN;

      else if (strncmp(argv[arg], "matepair", 2) == 0)
        dumpFlags = DUMP_MATEPAIR;

      else if (strncmp(argv[arg], "sizes", 2) == 0)
        dumpFlags = DUMP_SIZES;
      
      else
        fprintf(stderr, "%s: Unknown dump option '-d %s'\n", argv[0], argv[arg]);

    } else if (strcmp(argv[arg], "-D") == 0) {
      arg++;

      if      (strncmp(argv[arg], "unitiglist", 1) == 0)
        opType = OPERATION_UNITIGLIST;

      else if (strncmp(argv[arg], "contiglist", 1) == 0)
        opType = OPERATION_CONTIGLIST;

      else if (strncmp(argv[arg], "properties", 1) == 0)
        opType = OPERATION_PROPERTIES;

      else
        fprintf(stderr, "%s: Unknown dump option '-D %s'\n", argv[0], argv[arg]);

    } else if (strcmp(argv[arg], "-E") == 0) {
      opType = OPERATION_EDIT;
      editName = argv[++arg];

    } else if (strcmp(argv[arg], "-B") == 0) {
      opType = OPERATION_BUILD;
      buildName = argv[++arg];


    } else if (strcmp(argv[arg], "-R") == 0) {
      opType = OPERATION_REPLACE;
      replaceName = argv[++arg];

    } else if (strcmp(argv[arg], "-N") == 0) {
      replaceInPlace = FALSE;

    } else if (strcmp(argv[arg], "-compress") == 0) {
      opType = OPERATION_COMPRESS;

    } else if (strcmp(argv[arg], "-fmap") == 0) {
      opType  = OPERATION_FMAP;
      minFmap = atoi(argv[++arg]);
      maxFmap = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-w") == 0) {
      MULTIALIGN_PRINT_WIDTH = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-s") == 0) {
      MULTIALIGN_PRINT_SPACING = atoi(argv[++arg]);
      sizSize                  = atoll(argv[arg]);

    } else if (strcmp(argv[arg], "-o") == 0) {
      mpaPrefix = argv[++arg];

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
    fprintf(stderr, "  -u id[-id]            Unitig to dump (for -d option); if A-B, dump tigs from id A to id B, inclusive\n");
    fprintf(stderr, "  -c id[-id]            Contig to dump (for -d option); if A-B, dump tigs from id A to id B, inclusive\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -U                    Dump ALL unitigs (for -d option)\n");
    fprintf(stderr, "  -C                    Dump ALL contigs (for -d option)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -d <operation>        Dump something about a multialign (-c or -u) in the store\n");
    fprintf(stderr, "     properties         ...properties\n");
    fprintf(stderr, "     frags              ...a list of fragments\n");
    fprintf(stderr, "     unitigs            ...a list of unitigs\n");
    fprintf(stderr, "     consensus          ...the consensus sequence\n");
    fprintf(stderr, "     consensusgapped    ...the consensus sequence, with gaps as indicated in the multialignment\n");
    fprintf(stderr, "     layout             ...the layout\n");
    fprintf(stderr, "     multialign         ...the full multialignment\n");
    fprintf(stderr, "     matepair           ...an analysis of the mate pairs\n");
    fprintf(stderr, "     sizes              ...an analysis of sizes of the tigs\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -E <editFile>         Change properties of multialigns\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -B <layout-file>      Construct a new store with unitigs in 'layout-file'.  Store versions\n");
    fprintf(stderr, "                        before that specified on the '-t' option are created but are empty.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -R <layout>           Replace a multialign with this one (type and id are from the layout)\n");
    fprintf(stderr, "  -N                    Replace a multialign in the next version of the store.  This option is\n");
    fprintf(stderr, "                        rarely useful, but is needed if the version of the store to add a multialign\n");
    fprintf(stderr, "                        does not exist.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -compress             Move tigs from earlier versions into the specified version.  This removes\n");
    fprintf(stderr, "                        historical versions of unitigs/contigs, and can save tremendous storage space,\n");
    fprintf(stderr, "                        but makes it impossible to back up the assembly past the specified versions\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -fmap min max         Output a map of fragment id to unitig id, for unitigs with min <= x <= max fragments.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  For '-d multialign':\n");
    fprintf(stderr, "  -w width              Width of the page.\n");
    fprintf(stderr, "  -s spacing            Spacing between reads on the same line.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  For '-d matepair':\n");
    fprintf(stderr, "  -o prefix             Output files will be written to 'prefix.*' in the current directory.\n");
    fprintf(stderr, "                        (defaults to 'tigStore' (the -t option) if not set.)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  For '-d sizes':\n");
    fprintf(stderr, "  -s genomesize         Denominator to use for n50 computation\n");
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
  //
  //  One can change partitioning by deleting a multialign from one partition and adding it to
  //  another partition.  Doing so WILL cause consensus to fail, as consensus is expecting a
  //  specific set of fragments in each partition.
  //
  //  It is not possible to add a new partition:
  //  MultiAlignStore::MultiAlignStore()-- ERROR, didn't find any unitigs or contigs in the store.  Correct version?


  if ((opType == OPERATION_BUILD) && (buildName != NULL)) {
    operationBuild(buildName, tigName, tigVers);
    exit(0);
  }


  if (opType == OPERATION_COMPRESS) {
    operationCompress(tigName, tigVers);
    exit(0);
  }


  if (opType == OPERATION_FMAP) {
    operationFmap(tigName, tigVers, minFmap, maxFmap, tigIsUnitig, tigID, dumpAll);
    exit(0);
  }


  gkpStore = new gkStore(gkpName, FALSE, FALSE);
  tigStore = new MultiAlignStore(tigName, tigVers, tigPartU, tigPartC, FALSE, FALSE, FALSE);


  if ((opType == OPERATION_EDIT) && (editName != NULL)) {
    delete tigStore;
    tigStore = new MultiAlignStore(tigName, tigVers, tigPartU, tigPartC, TRUE, TRUE, FALSE);
    changeProperties(tigStore, editName);
  }


  if ((opType == OPERATION_REPLACE) && (replaceName != NULL)) {
    if (tigID) {
      fprintf(stderr, "ERROR:  -R is incompatible with -c and -u.  Did you mean -cp or -up instead?\n");
      exit(1);
    }

    errno = 0;
    FILE         *F = fopen(replaceName, "r");
    if (errno)
      fprintf(stderr, "Failed to open '%s': %s\n", replaceName, strerror(errno)), exit(1);

    fprintf(stderr, "Reading layouts from '%s'.\n", replaceName);

    delete tigStore;
    tigStore = new MultiAlignStore(tigName, tigVers, tigPartU, tigPartC, TRUE, replaceInPlace, !replaceInPlace);

    MultiAlignT  *ma       = CreateEmptyMultiAlignT();
    bool          isUnitig = false;

    while (LoadMultiAlignFromHuman(ma, isUnitig, F) == true) {
      if (ma->data.num_frags + ma->data.num_unitigs == 0) {
        if (tigStore->isDeleted(ma->maID, isUnitig) == true) {
          fprintf(stderr, "DELETING %s %d -- ALREADY DELETED\n", (isUnitig) ? "unitig" : "contig", ma->maID);
        } else {
          fprintf(stderr, "DELETING %s %d\n", (isUnitig) ? "unitig" : "contig", ma->maID);
          tigStore->deleteMultiAlign(ma->maID, isUnitig);
        }
      } else {
        tigStore->insertMultiAlign(ma, isUnitig, FALSE);
        fprintf(stderr, "INSERTING %s %d\n", (isUnitig) ? "unitig" : "contig", ma->maID);
      }
    }

    fprintf(stderr, "Reading layouts from '%s' completed.\n", replaceName);

    fclose(F);
  }


  if (opType == OPERATION_UNITIGLIST) {
    tigStore->dumpMultiAlignRTable(true);
  }


  if (opType == OPERATION_CONTIGLIST) {
    tigStore->dumpMultiAlignRTable(false);
  }


  if (opType == OPERATION_PROPERTIES) {
    for (uint32 i=0; i<tigStore->numUnitigs(); i++) {
      if (tigStore->isDeleted(i, TRUE) == false) {
        fprintf(stdout, "unitig_coverage_stat %8u %d\n", i, tigStore->getUnitigCoverageStat(i));
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
  }


  if (opType == OPERATION_TIG) {
    int32  bgn = tigID;
    int32  end = tigID + 1;

    if (tigIDbgn != tigIDend) {
      bgn = tigIDbgn;
      end = tigIDend + 1;
    }

    if (dumpAll == TRUE) {
      bgn = 0;
      end = (tigIsUnitig) ? tigStore->numUnitigs() : tigStore->numContigs();
    }

    if ((tigIsUnitig == true) && (tigStore->numUnitigs() <= tigID)) {
      fprintf(stderr, "ERROR: only "F_U32" unitigs in the store; can't dump requested unitig "F_S32"\n",
              tigStore->numUnitigs(), tigID);
      bgn = 0;
      end = 0;
    }

    if ((tigIsUnitig == false) && (tigStore->numContigs() <= tigID)) {
      fprintf(stderr, "ERROR: only "F_U32" contigs in the store; can't dump requested contig "F_S32"\n",
              tigStore->numContigs(), tigID);
      bgn = 0;
      end = 0;
    }

    if (dumpFlags == DUMP_MATEPAIR) {
      mpa = new matePairAnalysis(gkpName);
      if (mpaPrefix == NULL)
        mpaPrefix = tigName;
    }

    if (dumpFlags == DUMP_SIZES) {
      siz = new sizeAnalysis(sizSize);
      if (sizPrefix == NULL)
        sizPrefix = tigName;
    }

    for (tigID=bgn; tigID<end; tigID++) {
      ma = tigStore->loadMultiAlign(tigID, tigIsUnitig);

      if (ma == NULL)
        continue;

      if (dumpFlags == DUMP_PROPERTIES)
        dumpProperties(tigStore, tigID, tigIsUnitig, ma);

      if (dumpFlags == DUMP_FRAGS)
        dumpFrags(tigStore, tigID, tigIsUnitig, ma);

      if (dumpFlags == DUMP_UNITIGS)
        dumpUnitigs(tigStore, tigID, tigIsUnitig, ma);

      if (dumpFlags == DUMP_CONSENSUS)
        dumpConsensus(tigStore, tigID, tigIsUnitig, ma, false);

      if (dumpFlags == DUMP_CONSENSUSGAPPED)
        dumpConsensus(tigStore, tigID, tigIsUnitig, ma, true);

      if (dumpFlags == DUMP_LAYOUT)
        DumpMultiAlignForHuman(stdout, ma, tigIsUnitig);

      if (dumpFlags == DUMP_MULTIALIGN)
        PrintMultiAlignT(stdout, ma, gkpStore, showQV, showDots, (tigIsUnitig) ? AS_READ_CLEAR_OBTCHIMERA : AS_READ_CLEAR_LATEST);

      if (dumpFlags == DUMP_MATEPAIR)
        mpa->evaluateTig(ma);

      if (dumpFlags == DUMP_SIZES)
        siz->evaluateTig(ma, tigIsUnitig);

      tigStore->unloadMultiAlign(tigID, tigIsUnitig);
    }
  }

  if (mpa) {
    mpa->finalize();
    mpa->printSummary(stdout);
    mpa->writeUpdate(mpaPrefix);
    mpa->drawPlots(mpaPrefix);
    delete mpa;
  }

  if (siz) {
    siz->finalize();
    siz->printSummary(stdout);
    delete siz;
  }

  delete gkpStore;
  delete tigStore;

  exit(0);
}
