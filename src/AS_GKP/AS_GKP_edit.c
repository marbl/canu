
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

static char CM_ID[] = "$Id: AS_GKP_edit.c,v 1.2 2007-08-13 05:48:47 brianwalenz Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <ctype.h>

#include "AS_global.h"
#include "AS_UTL_fileIO.h"
#include "AS_GKP_include.h"
#include "AS_PER_gkpStore.h"



//  perl's chomp is pretty nice
#define chomp(S)  { char *t=(S); while (*t) t++; t--; while (isspace(*t)) *t--=0; }
#define munch(S)  { while (*(S) &&  isspace(*(S))) (S)++; }
#define crunch(S) { while (*(S) && !isspace(*(S))) (S)++; }



static
void
setClear(GateKeeperFragmentRecord *gkfr, char *E, uint32 which, int update) {
  int b = gkfr->clearBeg[which];
  int e = gkfr->clearEnd[which];
  gkfr->clearBeg[which] = strtoul(E, &E, 10);
  munch(E);
  gkfr->clearEnd[which] = strtoul(E, &E, 10);
  if (update)
    fprintf(stdout, "frg uid "F_UID" %s %d %d -> %d %d\n",
            gkfr->readUID, AS_READ_CLEAR_NAMES[which],
            b, e,
            gkfr->clearBeg[which], gkfr->clearEnd[which]);
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
editStore(char *editsFileName, char *gkpStoreName, int update) {
  FILE            *F        = NULL;
  char             L[1024]  = {0};
  char            *E        = NULL;

  int              errors   = 0;

  if (testOpenGateKeeperStore(gkpStoreName, FALSE) == 0) {
    fprintf(stderr, "failed to open store '%s', exit.\n", gkpStoreName);
    exit(1);
  }

  gkpStore = openGateKeeperStore(gkpStoreName, update);
  if (gkpStore == NULL) {
    fprintf(stderr, "Failed to open %s\n", gkpStoreName);
    exit(1);
  }

  CDS_IID_t        lastElem = getLastElemFragStore(gkpStore) + 1;


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
    int        isFRG     = 0;
    int        isLIB     = 0;
    CDS_UID_t  UID       = 0;
    CDS_IID_t  IID       = 0;
    char       ACT[1024] = {0};

    chomp(L);
    E = L;
    munch(E);

    if ((E[0] == '#') || (E[0] == ';') || (E[0] == 0)) {
      if (update)
        fprintf(stdout, "%s\n", L);
      goto nextline;
    }

    if        (strncasecmp("frg", E, 3) == 0) {
      E += 3;
      munch(E);
      isFRG = 1;
    } else if (strncasecmp("lib", E, 3) == 0) {
      E += 3;
      munch(E);
      isLIB = 1;
    } else {
      fprintf(stderr, "unknwon edit line format: '%s'\n", L);
      goto nextline;
    }

    if        (strncasecmp("uid", E, 3) == 0) {
      E += 3;
      munch(E);
      UID = STR_TO_UID(E, &E, 10);
      IID = getGatekeeperUIDtoIID(gkpStore, UID, NULL);
    } else if (strncasecmp("iid", E, 3) == 0) {
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
      errors++;
      goto nextline;
    }

    //  Stripped out thing-type and id-type.
    munch(E);

    //  Copy the action to ACT, then terminate after the first word.
    //
    strcpy(ACT, E);
    {
      char *t = ACT;
      crunch(t);
      *t = 0;
    }

    //  Advance E past the action.
    //
    crunch(E);
    munch(E);

    //fprintf(stderr, "ACT='%s' E='%s'\n", ACT, E);

    if (isFRG) {
      GateKeeperFragmentRecord gkfr = {0};

      getIndexStore(gkpStore->frg, IID, &gkfr);

      if        (strcasecmp(ACT, AS_READ_CLEAR_NAMES[AS_READ_CLEAR_ORIG]) == 0) {
        setClear(&gkfr, E, AS_READ_CLEAR_ORIG, update);
      } else if (strcasecmp(ACT, AS_READ_CLEAR_NAMES[AS_READ_CLEAR_QLT]) == 0) {
        setClear(&gkfr, E, AS_READ_CLEAR_QLT, update);
      } else if (strcasecmp(ACT, AS_READ_CLEAR_NAMES[AS_READ_CLEAR_VEC]) == 0) {
        setClear(&gkfr, E, AS_READ_CLEAR_VEC, update);
      } else if (strcasecmp(ACT, AS_READ_CLEAR_NAMES[AS_READ_CLEAR_OBTINI]) == 0) {
        setClear(&gkfr, E, AS_READ_CLEAR_OBTINI, update);
      } else if (strcasecmp(ACT, AS_READ_CLEAR_NAMES[AS_READ_CLEAR_OBT]) == 0) {
        setClear(&gkfr, E, AS_READ_CLEAR_OBT, update);
      } else if (strcasecmp(ACT, AS_READ_CLEAR_NAMES[AS_READ_CLEAR_UTG]) == 0) {
        setClear(&gkfr, E, AS_READ_CLEAR_UTG, update);
      } else if (strcasecmp(ACT, AS_READ_CLEAR_NAMES[AS_READ_CLEAR_ECR1]) == 0) {
        setClear(&gkfr, E, AS_READ_CLEAR_ECR1, update);
      } else if (strcasecmp(ACT, AS_READ_CLEAR_NAMES[AS_READ_CLEAR_ECR2]) == 0) {
        setClear(&gkfr, E, AS_READ_CLEAR_ECR2, update);
      } else if (strcasecmp(ACT, "mateiid") == 0) {
        CDS_IID_t o = gkfr.mateIID;
        gkfr.mateIID = STR_TO_IID(E, &E, 10);
        if (update)
          fprintf(stdout, "frg uid "F_UID" mateiid "F_IID" -> mateiid "F_IID"\n",
                  gkfr.readUID, o, gkfr.mateIID);
      } else if (strcasecmp(ACT, "mateuid") == 0) {
        CDS_IID_t o = gkfr.mateIID;
        CDS_UID_t n = STR_TO_UID(E, &E, 10);
        gkfr.mateIID = getGatekeeperUIDtoIID(gkpStore, n, NULL);
        if (update)
          fprintf(stdout, "frg uid "F_UID" mateiid "F_IID" -> mateiid "F_IID" mateuid "F_UID"\n",
                  gkfr.readUID, o, gkfr.mateIID, n);
      } else if (strcasecmp(ACT, "readuid") == 0) {
        CDS_UID_t o = gkfr.readUID;
        gkfr.readUID = STR_TO_UID(E, &E, 10);  //  I _really_ hope you know what you're doing
        if (update)
          fprintf(stdout, "frg iid "F_IID" readuid "F_UID" -> "F_UID"\n",
                  gkfr.readIID, o, gkfr.readUID);
      } else if (strcasecmp(ACT, "libiid") == 0) {
        CDS_IID_t o = gkfr.libraryIID;
        gkfr.libraryIID = STR_TO_UID(E, &E, 10);
        if (update)
          fprintf(stdout, "frg uid "F_UID" libiid "F_IID" -> libiid "F_IID"\n",
                  gkfr.readUID, o, gkfr.libraryIID);
      } else if (strcasecmp(ACT, "libuid") == 0) {
        CDS_IID_t o = gkfr.libraryIID;
        CDS_UID_t n = STR_TO_UID(E, &E, 10);
        gkfr.libraryIID = getGatekeeperUIDtoIID(gkpStore, n, NULL);
        if (update)
          fprintf(stdout, "frg uid "F_UID" libiid "F_IID" -> libiid "F_IID" libuid "F_UID"\n",
                  gkfr.readUID, o, gkfr.libraryIID, n);
      } else if (strcasecmp(ACT, "plate") == 0) {
        CDS_UID_t o = gkfr.plateUID;
        gkfr.plateUID = STR_TO_UID(E, &E, 10);
        if (update)
          fprintf(stdout, "frg uid "F_UID" plate "F_UID" -> "F_UID"\n",
                  gkfr.readUID, o, gkfr.plateUID);
      } else if (strcasecmp(ACT, "platelocation") == 0) {
        uint32 o = gkfr.plateLocation;
        gkfr.plateLocation = strtoul(E, &E, 10);
        if (update)
          fprintf(stdout, "frg uid "F_UID" platelocation "F_U32" -> "F_U32"\n",
                  gkfr.readUID, o, gkfr.plateLocation);
      } else if (strcasecmp(ACT, "isnonrandom") == 0) {
        uint32 o = gkfr.nonrandom;
        if      ((E[0] == '1') || (E[0] == 't') || (E[0] == 'T'))
          gkfr.nonrandom = 1;
        else if ((E[0] == '0') || (E[0] == 'f') || (E[0] == 'F'))
          gkfr.nonrandom = 0;
        else {
          fprintf(stderr, "invalid frg isnonrandom flag in edit line: '%s'\n", L);
          errors++;
          goto nextline;
        }
        if (update)
          fprintf(stdout, "frg uid "F_UID" isnonrandom "F_U32" -> "F_U32"\n",
                  gkfr.readUID, o, gkfr.nonrandom);
      } else if (strcasecmp(ACT, "isdeleted") == 0) {
        uint32 o = gkfr.deleted;
        if      ((E[0] == '1') || (E[0] == 't') || (E[0] == 'T'))
          gkfr.deleted = 1;
        else if ((E[0] == '0') || (E[0] == 'f') || (E[0] == 'F'))
          gkfr.deleted = 0;
        else {
          fprintf(stderr, "invalid frg isdeleted flag in edit line: '%s'\n", L);
          errors++;
          goto nextline;
        }
        if (update)
          fprintf(stdout, "frg uid "F_UID" isdeleted "F_U32" -> "F_U32"\n",
                  gkfr.readUID, o, gkfr.deleted);
      } else if (strcasecmp(ACT, "status") == 0) {
        uint32 o = gkfr.status;
        uint32 i;
        for (i=0; i<9; i++) {
          if (tolower(AS_READ_STATUS_NAMES[i][0]) == tolower(E[0])) {
            gkfr.status = i;
            break;
          }
        }
        if (i == 9) {
          fprintf(stderr, "invalid frg status in edit line: '%s'\n", L);
          errors++;
          goto nextline;
        }
        if (update)
          fprintf(stdout, "frg uid "F_UID" status %s -> %s\n",
                  gkfr.readUID, AS_READ_STATUS_NAMES[o], AS_READ_STATUS_NAMES[gkfr.status]);
      } else if (strcasecmp(ACT, "orientation") == 0) {
        uint32 o = gkfr.orientation;
        uint32 i;
        for (i=0; i<5; i++)
          if (tolower(AS_READ_ORIENT_NAMES[i][0]) == tolower(E[0])) {
            gkfr.orientation = i;
            break;
          }
        if (i == 5) {
          fprintf(stderr, "invalid frg orientation in edit line: '%s'\n", L);
          errors++;
          goto nextline;
        }
        if (update)
          fprintf(stdout, "frg uid "F_UID" orientation %s -> %s\n",
                  gkfr.readUID, AS_READ_ORIENT_NAMES[o], AS_READ_ORIENT_NAMES[gkfr.orientation]);
      } else {
        fprintf(stderr, "invalid frg action in edit line: '%s'\n", L);
        errors++;
        goto nextline;
      }

      if (update)
        setIndexStore(gkpStore->frg, IID, &gkfr);
    }

    if (isLIB) {
      GateKeeperLibraryRecord  gklr = {0};

      getIndexStore(gkpStore->lib, IID, &gklr);

      if        (strcasecmp(ACT, "mean") == 0) {
        double m = gklr.mean;
        gklr.mean   = atof(E);
        fprintf(stdout, "lib uid "F_UID" mean %f -> %f\n",
                gklr.libraryUID, m, gklr.mean);
      } else if (strcasecmp(ACT, "stddev") == 0) {
        double s = gklr.stddev;
        gklr.stddev = atof(E);
        fprintf(stdout, "lib uid "F_UID" mean %f -> %f\n",
                gklr.libraryUID, s, gklr.stddev);
      } else if (strcasecmp(ACT, "distance") == 0) {
        double m = gklr.mean;
        double s = gklr.stddev;
        gklr.mean   = atof(E);
        crunch(E);
        munch(E);
        gklr.stddev = atof(E);
        fprintf(stdout, "lib uid "F_UID" distance %f %f -> %f %f\n",
                gklr.libraryUID, m, s, gklr.mean, gklr.stddev);
      } else if (strcasecmp(ACT, "comment") == 0) {
        fprintf(stdout, "lib uid "F_UID" comment \"%s\" -> \"%s\"\n",
                gklr.libraryUID, gklr.comment, E);
        memset(gklr.comment, 0, AS_PER_COMMENT_LEN);
        strncpy(gklr.comment, E, AS_PER_COMMENT_LEN);

      //  Lots of boilerplate here for T/F flags.  Sigh.

      } else if (strcasecmp(ACT, "hpsisflowgram") == 0) {
        uint32 o = gklr.hpsIsFlowGram;
        if      ((E[0] == '1') || (E[0] == 't') || (E[0] == 'T'))
          gklr.hpsIsFlowGram = 1;
        else if ((E[0] == '0') || (E[0] == 'f') || (E[0] == 'F'))
          gklr.hpsIsFlowGram = 0;
        else {
          fprintf(stderr, "invalid lib hpsisflowgram flag in edit line: '%s'\n", L);
          errors++;
          goto nextline;
        }
        fprintf(stdout, "lib uid "F_UID" hpsisflowgram %c -> %c\n",
                gklr.libraryUID, (o) ? 'T' : 'F', (gklr.hpsIsFlowGram) ? 'T' : 'F');
      } else if (strcasecmp(ACT, "hpsispeakspacing") == 0) {
        uint32 o = gklr.hpsIsPeakSpacing;
        if      ((E[0] == '1') || (E[0] == 't') || (E[0] == 'T'))
          gklr.hpsIsPeakSpacing = 1;
        else if ((E[0] == '0') || (E[0] == 'f') || (E[0] == 'F'))
          gklr.hpsIsPeakSpacing = 0;
        else {
          fprintf(stderr, "invalid lib hpsispeakspacing flag in edit line: '%s'\n", L);
          errors++;
          goto nextline;
        }
        fprintf(stdout, "lib uid "F_UID" hpsispeakspacing %c -> %c\n",
                gklr.libraryUID, (o) ? 'T' : 'F', (gklr.hpsIsPeakSpacing) ? 'T' : 'F');
      } else if (strcasecmp(ACT, "donottrusthomopolymerruns") == 0) {
        uint32 o = gklr.doNotTrustHomopolymerRuns;
        if      ((E[0] == '1') || (E[0] == 't') || (E[0] == 'T'))
          gklr.doNotTrustHomopolymerRuns = 1;
        else if ((E[0] == '0') || (E[0] == 'f') || (E[0] == 'F'))
          gklr.doNotTrustHomopolymerRuns = 0;
        else {
          fprintf(stderr, "invalid lib donottrushhomopolymerruns flag in edit line: '%s'\n", L);
          errors++;
          goto nextline;
        }
        fprintf(stdout, "lib uid "F_UID" donottrushhomopolymerruns %c -> %c\n",
                gklr.libraryUID, (o) ? 'T' : 'F', (gklr.doNotTrustHomopolymerRuns) ? 'T' : 'F');
      } else if (strcasecmp(ACT, "donotoverlaptrim") == 0) {
        uint32 o = gklr.doNotOverlapTrim;
        if      ((E[0] == '1') || (E[0] == 't') || (E[0] == 'T'))
          gklr.doNotOverlapTrim = 1;
        else if ((E[0] == '0') || (E[0] == 'f') || (E[0] == 'F'))
          gklr.doNotOverlapTrim = 0;
        else {
          fprintf(stderr, "invalid lib donotoverlaptrim flag in edit line: '%s'\n", L);
          errors++;
          goto nextline;
        }
        fprintf(stdout, "lib uid "F_UID" donotoverlaptrim %c -> %c\n",
                gklr.libraryUID, (o) ? 'T' : 'F', (gklr.doNotOverlapTrim) ? 'T' : 'F');

      } else if (strcasecmp(ACT, "isnotrandom") == 0) {
        uint32 o = gklr.isNotRandom;
        if      ((E[0] == '1') || (E[0] == 't') || (E[0] == 'T'))
          gklr.isNotRandom = 1;
        else if ((E[0] == '0') || (E[0] == 'f') || (E[0] == 'F'))
          gklr.isNotRandom = 0;
        else {
          fprintf(stderr, "invalid lib isnotrandom flag in edit line: '%s'\n", L);
          errors++;
          goto nextline;
        }
        fprintf(stdout, "lib uid "F_UID" isnotrandom %c -> %c\n",
                gklr.libraryUID, (o) ? 'T' : 'F', (gklr.isNotRandom) ? 'T' : 'F');
      } else if (strcasecmp(ACT, "orientation") == 0) {
        uint32 o = gklr.orientation;
        uint32 i;
        for (i=0; i<5; i++)
          if (tolower(AS_READ_ORIENT_NAMES[i][0]) == tolower(E[0])) {
            gklr.orientation = i;
            break;
          }
        if (i == 5) {
          fprintf(stderr, "invalid lkg orientation in edit line: '%s'\n", L);
          errors++;
          goto nextline;
        }
        if (update)
          fprintf(stdout, "lkg uid "F_UID" orientation %s -> %s\n",
                  gklr.libraryUID, AS_READ_ORIENT_NAMES[o], AS_READ_ORIENT_NAMES[gklr.orientation]);
      } else if (strcasecmp(ACT, "allfragsdeleted") == 0) {
        assert(0);
      } else if (strcasecmp(ACT, "allfragsnonrandom") == 0) {
        assert(0);
      } else if (strcasecmp(ACT, "allfragsunmated") == 0) {
        assert(0);

      } else {
        fprintf(stderr, "invalid lib action in edit line: '%s'\n", L);
        errors++;
        goto nextline;
      }

      if (update)
        setIndexStore(gkpStore->lib, IID, &gklr);
    }


  nextline:
    fgets(L, 1024, F);
  }

  if (errors)
    fprintf(stderr, "%d errors detected.\n", errors);

  exit(errors != 0);
}
