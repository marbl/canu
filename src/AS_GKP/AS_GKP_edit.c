
/**************************************************************************
 * This file is part of Celera Assembler, a software program that 
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 2007, J. Craig Venter Institute
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

static char CM_ID[] = "$Id: AS_GKP_edit.c,v 1.1 2007-08-10 06:50:40 brianwalenz Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <ctype.h>

#include "AS_global.h"
#include "AS_UTL_fileIO.h"
#include "AS_GKP_include.h"
#include "AS_PER_gkpStore.h"



//  perl's chomp is pretty nice
#define chomp(S) { char *t=(S); while (*t) t++; t--; while (isspace(*t)) *t--=0; }
#define munch(S) { while (*(S) && isspace(*(S))) (S)++; }



static
void
setClear(GateKeeperFragmentRecord *gkfr, char *E, uint32 which) {
  gkfr->clearBeg[which] = strtoul(E, &E, 10);
  munch(E);
  gkfr->clearEnd[which] = strtoul(E, &E, 10);
}



void
updateVectorClear(char *vectorClearFile, char *gkpStoreName) {
  char          line[256];
  int           nlines  = 0;
  int           nupdate = 0;
  FILE         *F = NULL;
  char          N[FILENAME_MAX] = {0};

  errno = 0;
  FILE   *v = fopen(vectorClearFile, "r");
  if (errno) {
    fprintf(stderr, "couldn't open '%s' to read vector clear ranges: %s\n",
            vectorClearFile, strerror(errno));
    exit(1);
  }

  GateKeeperStore *gkpStore = openGateKeeperStore(gkpStoreName, TRUE);

  //  Suck in the whole frg info file, rather than doing random
  //  access all over the place.  Don't forget that there is a header
  //  on said file.  XXX -- really, someone should improve the
  //  GenericStore interface so that one could loadPartialStore()
  //  and dump it back.
#warning Violating genericStore encapsulation.

  uint32                     nfrg = getNumGateKeeperFragments(gkpStore);
  uint32                     nfrr = 0;
  GateKeeperFragmentRecord  *frg  = (GateKeeperFragmentRecord *)safe_malloc(nfrg * sizeof(GateKeeperFragmentRecord));

  sprintf(N, "%s/frg", gkpStoreName);

  errno = 0;
  F = fopen(N, "r");
  if (errno) {
    fprintf(stderr, "couldn't open '%s' to read fragment info: %s\n",
            vectorClearFile, strerror(errno));
    exit(1);
  }
  fprintf(stderr, "Loading "F_U32" fragment info's from '%s'\n", nfrg, N);
  AS_UTL_fseek(F, sizeof(StoreStat), SEEK_SET);
  nfrr = AS_UTL_safeRead(F, frg, "GKPvectorLoaderHack", sizeof(GateKeeperFragmentRecord), nfrg);
  if (nfrr != nfrg)
    fprintf(stderr, "ERROR:  Read only "F_U32" records instead of "F_U32".\n", nfrr, nfrg);
  assert(nfrr == nfrg);
  fclose(F);

  fgets(line, 256, v);
  while (!feof(v)) {
    char          *pine = line;

    CDS_UID_t uid = STR_TO_UID(pine, &pine, 10);
    int       l   = strtol(pine, &pine, 10);
    int       r   = strtol(pine, &pine, 10);

    if (uid == 0) {
      fprintf(stderr, "unexpected line: %s", line);
    } else {
      CDS_IID_t  iid = getGatekeeperUIDtoIID(gkpStore, uid, NULL);
      if (iid) {
        iid--;  // because an iid is not an index into an array
        frg[iid].hasVectorClear = 1;
        if (l < r) {
          frg[iid].clearBeg[AS_READ_CLEAR_VEC] = l - 1;  //  Assume they are base-based.
          frg[iid].clearEnd[AS_READ_CLEAR_VEC] = r;
        } else {
          frg[iid].clearBeg[AS_READ_CLEAR_VEC] = r - 1;
          frg[iid].clearEnd[AS_READ_CLEAR_VEC] = l;
        }
        nupdate++;
      }
      nlines++;
    }
    fgets(line, 256, v);
  }

  fclose(v);

  closeGateKeeperStore(gkpStore);

  fprintf(stderr, "in %d lines, updated %d fragments.\n", nlines, nupdate);

  errno = 0;
  F = fopen(N, "r+");
  if (errno) {
    fprintf(stderr, "couldn't open '%s' to write fragment info: %s\n",
            vectorClearFile, strerror(errno));
    exit(1);
  }
  fprintf(stderr, "Writing "F_U32" fragment info's from '%s'\n", nfrg, N);
  AS_UTL_fseek(F, sizeof(StoreStat), SEEK_SET);
  AS_UTL_safeWrite(F, frg, "GKPvectorLoaderHack", sizeof(GateKeeperFragmentRecord), nfrg);
  fclose(F);

  exit(0);
}




void
editStore(char *editsFileName, char *gkpStoreName) {

  if (testOpenGateKeeperStore(gkpStoreName, FALSE) == 0) {
    fprintf(stderr, "failed to open store '%s', exit.\n", gkpStoreName);
    exit(1);
  }

  GateKeeperStore *gkp      = NULL;
  fragRecord       fr       = {0};

  gkp = openGateKeeperStore(gkpStoreName, FALSE);
  if (gkp == NULL) {
    fprintf(stderr, "Failed to open %s\n", gkpStoreName);
    exit(1);
  }

  CDS_IID_t        lastElem = getLastElemFragStore(gkp) + 1;
  FILE            *F        = NULL;
  char             L[1024];
  char            *E;

  //  "frg uid UID THING DATA"
  //  "frg iid IID THING DATA"
  //  "lib uid UID THING DATA"
  //  "lib iid IID THING DATA"
  //
  //  E.g., "frg uid 1047118958955 lib 1099982711595"
  //        "lib iid 1 mean 4300.0"
  //        "lib iid 1 stddev 400.0"
  //        "lib iid 1 distance 4300.0 400.0"
  //        "lib iid 1 comment sample=ak42 label=ak42 name=ak42-G-01-1p4-2p0"

  errno = 0;
  if (strcmp(editsFileName, "-") == 0)
    F = stdin;
  else
    F = fopen(editsFileName, "r");
  if (errno) {
    fprintf(stderr, "couldn't open edits file '%s': %s\n", editsFileName, strerror(errno));
    exit(1);
  }
  fgets(L, 1024, F);
  while (!feof(F)) {
    int        isFRG = 0;
    int        isLIB = 0;
    CDS_UID_t  UID   = 0;
    CDS_IID_t  IID   = 0;
    char      *ACT   = NULL;

    chomp(L);
    E = L;

    if        (strncmp("frg", E, 3) == 0) {
      E += 3;
      munch(E);
      isFRG = 1;
    } else if (strncmp("lib", E, 3) == 0) {
      E += 3;
      munch(E);
      isLIB = 1;
    } else {
      fprintf(stderr, "unknwon edit line format: '%s'\n", L);
      goto nextline;
    }

    if        (strncmp("uid", E, 3) == 0) {
      E += 3;
      munch(E);
      UID = STR_TO_UID(E, &E, 10);
      IID = getGatekeeperUIDtoIID(gkp, UID, NULL);
    } else if (strncmp("iid", E, 3) == 0) {
      E += 3;
      munch(E);
      UID = 0;
      IID = STR_TO_IID(E, &E, 10);
    } else {
      fprintf(stderr, "unknwon edit line format: '%s'\n", L);
      goto nextline;
    }

    if (IID == 0) {
      fprintf(stderr, "invalid id (UID="F_UID", IID="F_IID") in edit line: '%s'\n", UID, IID, L);
      goto nextline;
    }

    munch(E);
    ACT = E;

    while (*E && !isspace(*E))
      E++;
    *E++ = 0;
    munch(E);

    if (isFRG) {
      GateKeeperFragmentRecord gkfr;

      getIndexStore(gkpStore->frg, IID, &gkfr);

      if        (strcmp(ACT, "orig") == 0) {
        setClear(&gkfr, E, AS_READ_CLEAR_ORIG);
      } else if (strcmp(ACT, "qlt") == 0) {
        setClear(&gkfr, E, AS_READ_CLEAR_QLT);
      } else if (strcmp(ACT, "vec") == 0) {
        setClear(&gkfr, E, AS_READ_CLEAR_VEC);
      } else if (strcmp(ACT, "obtini") == 0) {
        setClear(&gkfr, E, AS_READ_CLEAR_OBTINI);
      } else if (strcmp(ACT, "obt") == 0) {
        setClear(&gkfr, E, AS_READ_CLEAR_OBT);
      } else if (strcmp(ACT, "utg") == 0) {
        setClear(&gkfr, E, AS_READ_CLEAR_UTG);
      } else if (strcmp(ACT, "ecr1") == 0) {
        setClear(&gkfr, E, AS_READ_CLEAR_ECR1);
      } else if (strcmp(ACT, "ecr2") == 0) {
        setClear(&gkfr, E, AS_READ_CLEAR_ECR2);
      } else if (strcmp(ACT, "mate") == 0) {
        gkfr.mateIID = getGatekeeperUIDtoIID(gkp, STR_TO_UID(E, &E, 10), NULL);
      } else if (strcmp(ACT, "uid") == 0) {
        gkfr.readUID = STR_TO_UID(E, &E, 10);  //  I _really_ hope you know what you're doing
      } else if (strcmp(ACT, "lib") == 0) {
        gkfr.libraryIID = getGatekeeperUIDtoIID(gkp, STR_TO_UID(E, &E, 10), NULL);
      } else if (strcmp(ACT, "plate") == 0) {
        gkfr.plateUID = STR_TO_UID(E, &E, 10);
      } else if (strcmp(ACT, "plateloc") == 0) {
        gkfr.plateLocation = strtoul(E, &E, 10);
      } else if (strcmp(ACT, "isnonrandom") == 0) {
        gkfr.nonrandom = (strtoul(E, &E, 10) == 0) ? 0 : 1;
      } else if (strcmp(ACT, "isdeleted") == 0) {
        gkfr.deleted = (strtoul(E, &E, 10) == 0) ? 0 : 1;
      } else if (strcmp(ACT, "status") == 0) {
        gkfr.status = strtoul(E, &E, 10);
      } else if (strcmp(ACT, "orient") == 0) {
        gkfr.orientation = STR_TO_UID(E, &E, 10);
      } else {
        fprintf(stderr, "invalid frg action in edit line: '%s'\n", L);
        goto nextline;
      }

      setIndexStore(gkpStore->frg, IID, &gkfr);
    }

    if (isLIB) {
      GateKeeperLibraryRecord  gklr;

      getIndexStore(gkpStore->lib, IID, &gklr);

      if        (strcmp(ACT, "mean") == 0) {
        gklr.mean   = atof(E);
      } else if (strcmp(ACT, "stddev") == 0) {
        gklr.stddev = atof(E);
      } else if (strcmp(ACT, "distance") == 0) {
        gklr.mean   = atof(E);
        while (*E && !isspace(*E))
          E++;
        munch(E);
        gklr.stddev = atof(E);
      } else if (strcmp(ACT, "comment") == 0) {
        memset(gklr.comment, 0, AS_PER_COMMENT_LEN);
        strncpy(gklr.comment, E, AS_PER_COMMENT_LEN);

      } else if (strcmp(ACT, "hpsisflowgram") == 0) {
        gklr.hpsIsFlowGram = (strtoul(E, &E, 10) == 0) ? 0 : 1;
      } else if (strcmp(ACT, "hpsispeakspacing") == 0) {
        gklr.hpsIsPeakSpacing = (strtoul(E, &E, 10) == 0) ? 0 : 1;
      } else if (strcmp(ACT, "donottrusthomopolymerruns") == 0) {
        gklr.doNotTrustHomopolymerRuns = (strtoul(E, &E, 10) == 0) ? 0 : 1;
      } else if (strcmp(ACT, "donotoverlaptrim") == 0) {
        gklr.doNotOverlapTrim = (strtoul(E, &E, 10) == 0) ? 0 : 1;

      } else if (strcmp(ACT, "isnotrandom") == 0) {
        gklr.isNotRandom = (strtoul(E, &E, 10) == 0) ? 0 : 1;
      } else if (strcmp(ACT, "orientation") == 0) {
        gklr.orientation = (strtoul(E, &E, 10) == 0) ? 0 : 1;

      } else if (strcmp(ACT, "allfragsdeleted") == 0) {
      } else if (strcmp(ACT, "allfragsnonrandom") == 0) {
      } else if (strcmp(ACT, "allfragsunmated") == 0) {

      } else {
        fprintf(stderr, "invalid lib action in edit line: '%s'\n", L);
        goto nextline;
      }

      setIndexStore(gkpStore->lib, IID, &gklr);
    }


  nextline:
    fgets(L, 1024, F);
  }

  exit(0);
}
