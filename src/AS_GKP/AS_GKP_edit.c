
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

static char const *rcsid = "$Id: AS_GKP_edit.c,v 1.11 2008-06-27 06:29:16 brianwalenz Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <ctype.h>

#include "AS_global.h"
#include "AS_UTL_fileIO.h"
#include "AS_GKP_include.h"
#include "AS_PER_gkpStore.h"


#define EDIT_ALL_CLR "ALL"


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

  GateKeeperStore          *gkpStore = openGateKeeperStore(gkpStoreName, TRUE);
  fragRecord                fr;

  gkpStore->frg = convertStoreToMemoryStore(gkpStore->frg);

  fgets(line, 256, v);
  while (!feof(v)) {
    char          *pine = line;

    AS_UID    uid = AS_UID_lookup(pine, &pine);
    int       l   = strtol(pine, &pine, 10);
    int       r   = strtol(pine, &pine, 10);
    int       ll;
    int       rr;

    if (AS_UID_isDefined(uid) == FALSE) {
      fprintf(stderr, "unexpected line: %s", line);
    } else {
      AS_IID     iid = getGatekeeperUIDtoIID(gkpStore, uid, NULL);

      if (iid) {
        getFrag(gkpStore, iid, &fr, FRAG_S_INF);

        fr.gkfr.hasVectorClear = 1;

        //  Silently flip coordinates if needed.
        if (l < r) {
          ll = l - 1;
          rr = r;
        } else {
          ll = r - 1;
          rr = ll;
        }

        //  Not silently fix invalid coords.
        if ((ll < 0) ||
            (rr < 0) ||
            (ll >= AS_READ_MAX_LEN) ||
            (rr >= AS_READ_MAX_LEN)) {
          if (ll <  0)                  ll = 0;
          if (rr >= AS_READ_MAX_LEN)    ll = AS_READ_MAX_LEN-1;
          if (ll <  0)                  rr = 0;
          if (rr >= AS_READ_MAX_LEN)    rr = AS_READ_MAX_LEN-1;
          fprintf(stderr, "WARNING:  Fixing vector clear range for %s from (%d,%d) to (%d,%d).\n", l, r, ll, rr);
        }

        fr.gkfr.clearBeg[AS_READ_CLEAR_VEC] = ll;
        fr.gkfr.clearEnd[AS_READ_CLEAR_VEC] = rr;

        setFrag(gkpStore, iid, &fr);
        nupdate++;
      }
      nlines++;
    }
    fgets(line, 256, v);
  }

  fclose(v);

  closeGateKeeperStore(gkpStore);

  fprintf(stderr, "in %d lines, updated %d fragments.\n", nlines, nupdate);

  exit(0);
}





static
void
setClear(GateKeeperFragmentRecord *gkfr, char *E, uint32 which, int update) {
  int b = gkfr->clearBeg[which];
  int e = gkfr->clearEnd[which];
  int start = 0;
  int end   = 0;
  int updateAll = 0;


  /* possible formats are CLR start end or CLR ALL start end
   * here, we grab the next word and see if it is equal to all or not
   */
  char *next = E;
  crunch(E);
  // put terminator character for next
  {
    char *t = E;
    munch(E);
    *t = '\0';
  }

  if (strcasecmp(next, EDIT_ALL_CLR) == 0) {
     // update all clear ranges
     start = strtoul(E, &E, 10);
     munch(E);
     end = strtoul(E, &E, 10);
     updateAll = 1;
  } else {
    gkfr->clearBeg[which] = strtoul(next, NULL, 10);
    munch(E);
    gkfr->clearEnd[which] = strtoul(E, &E, 10);
  }

  if (which == AS_READ_CLEAR_VEC) {
    gkfr->hasVectorClear = 1;
  }
  if (which == AS_READ_CLEAR_QLT)
    gkfr->hasQualityClear = 1;

  if (updateAll) {
     for (; which <= AS_READ_CLEAR_LATEST; which++) {
       gkfr->clearBeg[which] = start;
       gkfr->clearEnd[which] = end;
   }
   which = AS_READ_CLEAR_LATEST;
  }

  if (update)
    fprintf(stdout, "frg uid %s %s %d %d -> %d %d\n",
            AS_UID_toString(gkfr->readUID),
            AS_READ_CLEAR_NAMES[which],
            b, e,
            gkfr->clearBeg[which], gkfr->clearEnd[which]);
}




static
void
allFrags(GateKeeperStore *gkpStore,
         AS_IID           IID,
         char             action,
         int              flag,
         int              update) {
  GateKeeperFragmentRecord gkfr = {0};
  uint32                   i;

  int64  firstElem = getFirstElemStore(gkpStore->frg);
  int64  lastElem  = getLastElemStore(gkpStore->frg);

  if (update)
    fprintf(stderr, "delete all frags in lib "F_IID" (%d,%d)\n",
            IID, firstElem, lastElem);

  for (i=firstElem; i<lastElem; i++) {
    getIndexStore(gkpStore->frg, i, &gkfr);
    if (gkfr.libraryIID == IID) {
      if        (action == 'd') {
        gkfr.deleted = flag;
      } else if (action == 'r') {
        gkfr.nonrandom = flag;
      } else if (action == 'u') {
        gkfr.mateIID = 0;
      } else if (action == 's') {
        gkfr.status = flag;
      } else if (action == 'o') {
        gkfr.orientation = flag;
      } else {
        assert(0);
      }
      if (update)
        setIndexStore(gkpStore->frg, i, &gkfr);
    }
  }
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

  AS_IID           lastElem = getLastElemFragStore(gkpStore) + 1;


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
    AS_UID     UID       = AS_UID_undefined();
    AS_IID     IID       = 0;
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
      UID = AS_UID_lookup(E, &E);
      IID = getGatekeeperUIDtoIID(gkpStore, UID, NULL);
    } else if (strncasecmp("iid", E, 3) == 0) {
      E += 3;
      munch(E);
      UID = AS_UID_undefined();
      IID = AS_IID_fromString(E, &E);
    } else {
      fprintf(stderr, "unknwon edit line format: '%s'\n", L);
      goto nextline;
    }

    if (IID == 0) {
      fprintf(stderr, "invalid id (UID=%s, IID="F_IID") in edit line: '%s'\n", AS_UID_toString(UID), IID, L);
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

      if (IID > getNumGateKeeperFragments(gkpStore)) {
        fprintf(stderr, "invalid frg iid "F_IID" in edit line: '%s'\n", IID, L);
        errors++;
        goto nextline;
      }

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
        AS_IID    o = gkfr.mateIID;
        gkfr.mateIID = AS_IID_fromString(E, &E);
        if (update)
          fprintf(stdout, "frg uid %s mateiid "F_IID" -> mateiid "F_IID"\n",
                  AS_UID_toString(gkfr.readUID), o, gkfr.mateIID);
      } else if (strcasecmp(ACT, "mateuid") == 0) {
        AS_IID    o = gkfr.mateIID;
        AS_UID    n = AS_UID_lookup(E, &E);
        gkfr.mateIID = getGatekeeperUIDtoIID(gkpStore, n, NULL);
        if (update)
          fprintf(stdout, "frg uid %s mateiid "F_IID" -> mateiid "F_IID" mateuid %s\n",
                  AS_UID_toString(gkfr.readUID), o, gkfr.mateIID, AS_UID_toString(n));
      } else if (strcasecmp(ACT, "readuid") == 0) {
        AS_UID    o = gkfr.readUID;
        gkfr.readUID = AS_UID_lookup(E, &E);  //  I _really_ hope you know what you're doing
        if (update)
          fprintf(stdout, "frg iid "F_IID" readuid %s -> %s\n",
                  gkfr.readIID, AS_UID_toString(o), AS_UID_toString(gkfr.readUID));
      } else if (strcasecmp(ACT, "libiid") == 0) {
        AS_IID    o = gkfr.libraryIID;
        gkfr.libraryIID = AS_IID_fromString(E, &E);
        if (update)
          fprintf(stdout, "frg uid %s libiid "F_IID" -> libiid "F_IID"\n",
                  AS_UID_toString(gkfr.readUID), o, gkfr.libraryIID);
      } else if (strcasecmp(ACT, "libuid") == 0) {
        AS_IID    o = gkfr.libraryIID;
        AS_UID    n = AS_UID_lookup(E, &E);
        gkfr.libraryIID = getGatekeeperUIDtoIID(gkpStore, n, NULL);
        if (update)
          fprintf(stdout, "frg uid %s libiid "F_IID" -> libiid "F_IID" libuid %s\n",
                  AS_UID_toString(gkfr.readUID), o, gkfr.libraryIID, AS_UID_toString(n));
      } else if (strcasecmp(ACT, "plate") == 0) {
        AS_UID    o = gkfr.plateUID;
        gkfr.plateUID = AS_UID_lookup(E, &E);
        if (update)
          fprintf(stdout, "frg uid %s plate %s -> %s\n",
                  AS_UID_toString(gkfr.readUID), AS_UID_toString(o), AS_UID_toString(gkfr.plateUID));
      } else if (strcasecmp(ACT, "platelocation") == 0) {
        uint32 o = gkfr.plateLocation;
        gkfr.plateLocation = strtoul(E, &E, 10);
        if (update)
          fprintf(stdout, "frg uid %s platelocation "F_U32" -> "F_U32"\n",
                  AS_UID_toString(gkfr.readUID), o, gkfr.plateLocation);
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
          fprintf(stdout, "frg uid %s isnonrandom "F_U32" -> "F_U32"\n",
                  AS_UID_toString(gkfr.readUID), o, gkfr.nonrandom);
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
          fprintf(stdout, "frg uid %s isdeleted "F_U32" -> "F_U32"\n",
                  AS_UID_toString(gkfr.readUID), o, gkfr.deleted);
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
          fprintf(stdout, "frg uid %s status %s -> %s\n",
                  AS_UID_toString(gkfr.readUID), AS_READ_STATUS_NAMES[o], AS_READ_STATUS_NAMES[gkfr.status]);
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
          fprintf(stdout, "frg uid %s orientation %s -> %s\n",
                  AS_UID_toString(gkfr.readUID), AS_READ_ORIENT_NAMES[o], AS_READ_ORIENT_NAMES[gkfr.orientation]);
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

      if (IID > getNumGateKeeperLibraries(gkpStore)) {
        fprintf(stderr, "invalid lib iid "F_IID" in edit line: '%s'\n", IID, L);
        errors++;
        goto nextline;
      }

      getIndexStore(gkpStore->lib, IID, &gklr);

      if        (strcasecmp(ACT, "mean") == 0) {
        double m = gklr.mean;
        gklr.mean   = atof(E);
        if (update)
          fprintf(stdout, "lib uid %s mean %f -> %f\n",
                  AS_UID_toString(gklr.libraryUID), m, gklr.mean);
      } else if (strcasecmp(ACT, "stddev") == 0) {
        double s = gklr.stddev;
        gklr.stddev = atof(E);
        if (update)
          fprintf(stdout, "lib uid %s mean %f -> %f\n",
                  AS_UID_toString(gklr.libraryUID), s, gklr.stddev);
      } else if (strcasecmp(ACT, "distance") == 0) {
        double m = gklr.mean;
        double s = gklr.stddev;
        gklr.mean   = atof(E);
        crunch(E);
        munch(E);
        gklr.stddev = atof(E);
        if (update)
          fprintf(stdout, "lib uid %s distance %f %f -> %f %f\n",
                  AS_UID_toString(gklr.libraryUID), m, s, gklr.mean, gklr.stddev);
      } else if (strcasecmp(ACT, "comment") == 0) {
        if (update)
          fprintf(stdout, "lib uid %s comment \"%s\" -> \"%s\"\n",
                  AS_UID_toString(gklr.libraryUID), gklr.comment, E);
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
        if (update)
          fprintf(stdout, "lib uid %s hpsisflowgram %c -> %c\n",
                  AS_UID_toString(gklr.libraryUID), (o) ? 'T' : 'F', (gklr.hpsIsFlowGram) ? 'T' : 'F');
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
        if (update)
          fprintf(stdout, "lib uid %s hpsispeakspacing %c -> %c\n",
                  AS_UID_toString(gklr.libraryUID), (o) ? 'T' : 'F', (gklr.hpsIsPeakSpacing) ? 'T' : 'F');
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
        if (update)
          fprintf(stdout, "lib uid %s donottrushhomopolymerruns %c -> %c\n",
                  AS_UID_toString(gklr.libraryUID), (o) ? 'T' : 'F', (gklr.doNotTrustHomopolymerRuns) ? 'T' : 'F');
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
        if (update)
          fprintf(stdout, "lib uid %s donotoverlaptrim %c -> %c\n",
                  AS_UID_toString(gklr.libraryUID), (o) ? 'T' : 'F', (gklr.doNotOverlapTrim) ? 'T' : 'F');
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
        if (update)
          fprintf(stdout, "lib uid %s isnotrandom %c -> %c\n",
                  AS_UID_toString(gklr.libraryUID), (o) ? 'T' : 'F', (gklr.isNotRandom) ? 'T' : 'F');
        allFrags(gkpStore, IID, 'r', gklr.isNotRandom, update);
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
          fprintf(stdout, "lib uid %s orientation %s -> %s\n",
                  AS_UID_toString(gklr.libraryUID), AS_READ_ORIENT_NAMES[o], AS_READ_ORIENT_NAMES[gklr.orientation]);
        allFrags(gkpStore, IID, 'o', gklr.orientation, update);
      } else if (strcasecmp(ACT, "allfragsdeleted") == 0) {
        uint32 maketrue = 0;
        if      ((E[0] == '1') || (E[0] == 't') || (E[0] == 'T'))
          maketrue = 1;
        else if ((E[0] == '0') || (E[0] == 'f') || (E[0] == 'F'))
          maketrue = 0;
        else {
          fprintf(stderr, "invalid lib allfragsdeleted flag in edit line: '%s'\n", L);
          errors++;
          goto nextline;
        }
        allFrags(gkpStore, IID, 'd', maketrue, update);
        IID = 0;
      } else if (strcasecmp(ACT, "allfragsnonrandom") == 0) {
        uint32 maketrue = 0;
        if      ((E[0] == '1') || (E[0] == 't') || (E[0] == 'T'))
          maketrue = 1;
        else if ((E[0] == '0') || (E[0] == 'f') || (E[0] == 'F'))
          maketrue = 0;
        else {
          fprintf(stderr, "invalid lib allfragsnonrandom flag in edit line: '%s'\n", L);
          errors++;
          goto nextline;
        }
        allFrags(gkpStore, IID, 'r', maketrue, update);
        IID = 0;
      } else if (strcasecmp(ACT, "allfragsunmated") == 0) {
        uint32 maketrue = 0;
        if      ((E[0] == '1') || (E[0] == 't') || (E[0] == 'T')) {
          allFrags(gkpStore, IID, 'u', 1, update);
        } else {
          fprintf(stderr, "invalid lib allfragsunmated flag in edit line: '%s'\n", L);
          errors++;
          goto nextline;
        }
        IID = 0;
      } else {
        fprintf(stderr, "invalid lib action in edit line: '%s'\n", L);
        errors++;
        goto nextline;
      }

      //  IID == 0 if we did a "allfrags" op.  The lib doesn't need to
      //  be (and shouldn't get) updated.

      if ((IID > 0) && (update))
        setIndexStore(gkpStore->lib, IID, &gklr);
    }


  nextline:
    fgets(L, 1024, F);
  }

  closeGateKeeperStore(gkpStore);

  if (errors) {
    fprintf(stderr, "%d errors detected.\n", errors);
    exit(1);
  }

  fprintf(stderr, "Success!\n", errors);

  exit(0);
}
