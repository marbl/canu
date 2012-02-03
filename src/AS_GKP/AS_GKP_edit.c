
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

static char const *rcsid = "$Id: AS_GKP_edit.c,v 1.31 2012-02-03 21:47:58 brianwalenz Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <ctype.h>

#include "AS_global.h"
#include "AS_UTL_fileIO.h"
#include "AS_GKP_include.h"
#include "AS_PER_gkpStore.h"


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

  gkStore          *gkpStore = new gkStore(gkpStoreName, FALSE, TRUE);
  gkFragment        fr;
  
  // enable the clear range
  gkpStore->gkStore_enableClearRange(AS_READ_CLEAR_VEC);

  fgets(line, 256, v);
  while (!feof(v)) {
    char          *pine = line;
    char          *lstr;
    char          *rstr;

    chomp(line);

    AS_UID    uid = AS_UID_lookup(pine, &pine);
    munch(pine);
    int       l   = strtol(lstr = pine, &pine, 10);
    munch(pine);
    int       r   = strtol(rstr = pine, &pine, 10);
    int       ll;
    int       rr;
    int       mm;

    if (!isdigit(*lstr) || !isdigit(*rstr) || (AS_UID_isDefined(uid) == FALSE)) {
      fprintf(stderr, "unexpected line: %s\n", line);
    } else {
      AS_IID     iid = gkpStore->gkStore_getUIDtoIID(uid, NULL);

      if (iid) {
        gkpStore->gkStore_getFragment(iid, &fr, GKFRAGMENT_INF);

        mm = fr.gkFragment_getSequenceLength();

        //  Silently flip coordinates if needed.
        if (l < r) {
          ll = l - 1;
          rr = r;
        } else {
          ll = r - 1;
          rr = ll;
        }

        //  Now silently fix invalid coords.
        if ((ll < 0) ||
            (rr < 0) ||
            (ll > mm) ||
            (rr > mm)) {
          if (ll <  0)    ll = 0;
          if (rr > mm)    ll = mm;
          if (ll <  0)    rr = 0;
          if (rr > mm)    rr = mm;
          fprintf(stderr, "WARNING:  Fixing vector clear range for '%s' to (%d,%d).\n", line, ll, rr);
        }

        fr.gkFragment_setClearRegion(ll, rr, AS_READ_CLEAR_VEC);
        gkpStore->gkStore_setFragment(&fr);

        nupdate++;
      }
      nlines++;
    }
    fgets(line, 256, v);
  }

  fclose(v);

  delete gkpStore;

  fprintf(stderr, "in %d lines, updated %d fragments.\n", nlines, nupdate);

  exit(0);
}



void
revertClearRange(char *clearRegionName, char *gkpStoreName) {
  gkStore    *gkpStore = new gkStore(gkpStoreName, FALSE, TRUE);
  gkFragment  fr;
  uint32      br, er;  //  Begin, End, of the range to revert to
  uint32      bl, el;  //  Begin, End, of the latest
  uint32      which = gkStore_decodeClearRegionLabel(clearRegionName);

  if (which == AS_READ_CLEAR_ERROR)
    fprintf(stderr, "invalid clear region label %s\n", clearRegionName), exit(1);

  fr.gkFragment_enableGatekeeperMode(gkpStore);

  for (int32 i=1; i<gkpStore->gkStore_getNumFragments(); i++) {
    gkpStore->gkStore_getFragment(i, &fr, GKFRAGMENT_INF);

    fr.gkFragment_getClearRegion(br, er, which);
    fr.gkFragment_getClearRegion(bl, el, AS_READ_CLEAR_LATEST);

    //  If the latest is different, reset the clear region.  It looks like nonsense (why set the
    //  clear range to the value it already has?!) but it also updates the latest clear range to
    //  this value, which is what we want to do.
    //
    if (((br != bl) || (er != el)) &&
        (br <= er)) {
      fr.gkFragment_setClearRegion(br, er, which);
      gkpStore->gkStore_setFragment(&fr);
    }
  }

  for (which++; which < AS_READ_CLEAR_NUM; which++)
    gkpStore->gkStore_purgeClearRange(which);

  delete gkpStore;

  exit(0);
}




static
void
setClear(gkFragment *fr, char *E, uint32 which, int update, int verbose) {
  uint32 bo, eo, bn, en;

  //  Code used to exist to update "ALL" clear ranges, but that makes
  //  no sense (I hope) any more.  This was used by toggleAndRun, but
  //  now we can just nuke the appropriate clear range files in the
  //  store.

  bn = strtoul(E, &E, 10);
  munch(E);
  en = strtoul(E, &E, 10);

  fr->gkFragment_getClearRegion(bo, eo, which);
  fr->gkFragment_setClearRegion(bn, en, which);

  if (verbose)
    fprintf(stdout, "frg uid %s %s %d %d -> %d %d\n",
            AS_UID_toString(fr->gkFragment_getReadUID()),
            AS_READ_CLEAR_NAMES[which],
            bo, eo, bn, en);
}




static
void
allFrags(gkStore *gkpStore,
         AS_IID           IID,
         char             action,
         int              flag,
         int              update,
         int              verbose) {
  gkFragment  fr;
  uint32      i;

  int64  firstElem = 1;
  int64  lastElem  = gkpStore->gkStore_getNumFragments();

  if (verbose) {
    fprintf(stderr, "update all frags in lib "F_IID" ("F_S64","F_S64")\n",
            IID, firstElem, lastElem);    
  }

  fr.gkFragment_enableGatekeeperMode(gkpStore);
  for (i=firstElem; i<lastElem; i++) {
    gkpStore->gkStore_getFragment(i, &fr, GKFRAGMENT_INF);

    if (fr.gkFragment_getLibraryIID() == IID) {
      if        (action == 'd') {
        fr.gkFragment_setIsDeleted(flag);

      } else if (action == 'r') {
        fr.gkFragment_setIsNonRandom(flag);

      } else if (action == 'u') {
        fr.gkFragment_setMateIID(0);

      } else if (action == 'o') {
        fr.gkFragment_setOrientation(flag);

      } else {
        assert(0);
      }

      if (update)
        gkpStore->gkStore_setFragment(&fr);
    }
  }
}



void
editStore(char *editsFileName, char *gkpStoreName, int update) {
  FILE            *F        = NULL;
  char             L[1024]  = {0};
  char            *E        = NULL;

  int              errors   = 0;
  int              verbose  = 1;

  gkpStore = new gkStore(gkpStoreName, FALSE, update);

  AS_IID           lastElem = gkpStore->gkStore_getNumFragments() + 1;


  //  "frg uid UID THING DATA"
  //  "frg iid IID THING DATA"
  //  "lib uid UID THING DATA"
  //  "lib iid IID THING DATA"
  //
  //  E.g., "frg uid 1047118958955 lib 1099982711595"
  //        "lib iid 1 mean 4300.0"
  //        "lib iid 1 stddev 400.0"
  //        "lib iid 1 distance 4300.0 400.0"

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
      if (verbose)
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
      IID = gkpStore->gkStore_getUIDtoIID(UID, NULL);
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
      gkFragment fr;
      fr.gkFragment_enableGatekeeperMode(gkpStore);

      if (IID > gkpStore->gkStore_getNumFragments()) {
        fprintf(stderr, "invalid frg iid "F_IID" in edit line: '%s'\n", IID, L);
        errors++;
        goto nextline;
      }

      gkpStore->gkStore_getFragment(IID, &fr, GKFRAGMENT_INF);

      uint32  clr = gkStore_decodeClearRegionLabel(ACT);

      if        (clr != AS_READ_CLEAR_ERROR) {
        setClear(&fr, E, clr, update, verbose);
      } else if (strcasecmp(ACT, "mateiid") == 0) {
        AS_IID    o = fr.gkFragment_getMateIID();
        fr.gkFragment_setMateIID(AS_IID_fromString(E, &E));
        if (verbose)
          fprintf(stdout, "frg uid %s mateiid "F_IID" -> mateiid "F_IID"\n",
                  AS_UID_toString(fr.gkFragment_getReadUID()), o, fr.gkFragment_getMateIID());
      } else if (strcasecmp(ACT, "mateuid") == 0) {
        AS_IID    o = fr.gkFragment_getMateIID();
        AS_UID    n = AS_UID_lookup(E, &E);
        fr.gkFragment_setMateIID(gkpStore->gkStore_getUIDtoIID(n, NULL));
        if (verbose)
          fprintf(stdout, "frg uid %s mateiid "F_IID" -> mateiid "F_IID" mateuid %s\n",
                  AS_UID_toString(fr.gkFragment_getReadUID()), o, fr.gkFragment_getMateIID(), AS_UID_toString(n));
      } else if (strcasecmp(ACT, "readuid") == 0) {
        AS_UID    o = fr.gkFragment_getReadUID();
        fr.gkFragment_setReadUID(AS_UID_lookup(E, &E));  //  I _really_ hope you know what you're doing
        if (verbose)
          fprintf(stdout, "frg iid "F_IID" readuid %s -> %s\n",
                  fr.gkFragment_getReadIID(), AS_UID_toString(o), AS_UID_toString(fr.gkFragment_getReadUID()));
      } else if (strcasecmp(ACT, "libiid") == 0) {
        AS_IID    o = fr.gkFragment_getLibraryIID();
        fr.gkFragment_setLibraryIID(AS_IID_fromString(E, &E));
        if (verbose)
          fprintf(stdout, "frg uid %s libiid "F_IID" -> libiid "F_IID"\n",
                  AS_UID_toString(fr.gkFragment_getReadUID()), o, fr.gkFragment_getLibraryIID());
      } else if (strcasecmp(ACT, "libuid") == 0) {
        AS_IID    o = fr.gkFragment_getLibraryIID();
        AS_UID    n = AS_UID_lookup(E, &E);
        fr.gkFragment_setLibraryIID(gkpStore->gkStore_getUIDtoIID(n, NULL));
        if (verbose)
          fprintf(stdout, "frg uid %s libiid "F_IID" -> libiid "F_IID" libuid %s\n",
                  AS_UID_toString(fr.gkFragment_getReadUID()), o, fr.gkFragment_getLibraryIID(), AS_UID_toString(n));
      } else if (strcasecmp(ACT, "isnonrandom") == 0) {
        uint32 o = fr.gkFragment_getIsNonRandom();
        if      ((E[0] == '1') || (E[0] == 't') || (E[0] == 'T')) {
          fr.gkFragment_setIsNonRandom(1);
        } else if ((E[0] == '0') || (E[0] == 'f') || (E[0] == 'F')) {
          fr.gkFragment_setIsNonRandom(0);
        } else {
          fprintf(stderr, "invalid frg isnonrandom flag in edit line: '%s'\n", L);
          errors++;
          goto nextline;
        }
        if (verbose)
          fprintf(stdout, "frg uid %s isnonrandom "F_U32" -> "F_U32"\n",
                  AS_UID_toString(fr.gkFragment_getReadUID()), o, fr.gkFragment_getIsNonRandom());
      } else if (strcasecmp(ACT, "isdeleted") == 0) {
        uint32 o = fr.gkFragment_getIsDeleted();
        if      ((E[0] == '1') || (E[0] == 't') || (E[0] == 'T')) {
          fr.gkFragment_setIsDeleted(1);
        } else if ((E[0] == '0') || (E[0] == 'f') || (E[0] == 'F')) {
          fr.gkFragment_setIsDeleted(0);
        } else {
          fprintf(stderr, "invalid frg isdeleted flag in edit line: '%s'\n", L);
          errors++;
          goto nextline;
        }
        if (verbose)
          fprintf(stdout, "frg uid %s isdeleted "F_U32" -> "F_U32"\n",
                  AS_UID_toString(fr.gkFragment_getReadUID()), o, fr.gkFragment_getIsDeleted());
      } else if (strcasecmp(ACT, "orientation") == 0) {
        uint32 o = fr.gkFragment_getOrientation();
        uint32 i;
        for (i=0; i<5; i++)
          if (tolower(AS_READ_ORIENT_NAMES[i][0]) == tolower(E[0])) {
            fr.gkFragment_setOrientation(i);
            break;
          }
        if (i == 5) {
          fprintf(stderr, "invalid frg orientation in edit line: '%s'\n", L);
          errors++;
          goto nextline;
        }
        if (verbose)
          fprintf(stdout, "frg uid %s orientation %s -> %s\n",
                  AS_UID_toString(fr.gkFragment_getReadUID()), AS_READ_ORIENT_NAMES[o], AS_READ_ORIENT_NAMES[fr.gkFragment_getOrientation()]);
      } else {
        fprintf(stderr, "invalid frg action in edit line: '%s'\n", L);
        errors++;
        goto nextline;
      }

      if (update)
        gkpStore->gkStore_setFragment(&fr);
    }

    if (isLIB) {
      gkLibrary  gklr;

      if (IID > gkpStore->gkStore_getNumLibraries()) {
        fprintf(stderr, "invalid lib iid "F_IID" in edit line: '%s'\n", IID, L);
        errors++;
        goto nextline;
      }

      gkpStore->gkStore_getLibrary(IID, &gklr);

      if        (strcasecmp(ACT, "mean") == 0) {
        double m = gklr.mean;
        gklr.mean   = atof(E);
        if (verbose)
          fprintf(stdout, "lib uid %s mean %f -> %f\n",
                  AS_UID_toString(gklr.libraryUID), m, gklr.mean);
      } else if (strcasecmp(ACT, "stddev") == 0) {
        double s = gklr.stddev;
        gklr.stddev = atof(E);
        if (verbose)
          fprintf(stdout, "lib uid %s mean %f -> %f\n",
                  AS_UID_toString(gklr.libraryUID), s, gklr.stddev);
      } else if (strcasecmp(ACT, "distance") == 0) {
        double m = gklr.mean;
        double s = gklr.stddev;
        gklr.mean   = atof(E);
        crunch(E);
        munch(E);
        gklr.stddev = atof(E);
        if (verbose)
          fprintf(stdout, "lib uid %s distance %f %f -> %f %f\n",
                  AS_UID_toString(gklr.libraryUID), m, s, gklr.mean, gklr.stddev);

      //  Lots of boilerplate here for T/F flags.  Sigh.

#define setBoolean(XX)                                                              \
      } else if (strcasecmp(ACT, #XX) == 0) {                                       \
        uint32 o = gklr.XX;                                                         \
        if      ((E[0] == '1') || (E[0] == 't') || (E[0] == 'T'))                   \
          gklr.XX = 1;                                                              \
        else if ((E[0] == '0') || (E[0] == 'f') || (E[0] == 'F'))                   \
          gklr.XX = 0;                                                              \
        else {                                                                      \
          fprintf(stderr, "invalid lib XX flag in edit line: '%s'\n", L);           \
          errors++;                                                                 \
          goto nextline;                                                            \
        }                                                                           \
        if (verbose)                                                                \
          fprintf(stdout, "lib uid %s XX %c -> %c\n",                               \
                  AS_UID_toString(gklr.libraryUID),                                 \
                  (o) ? 'T' : 'F',                                                  \
                  (gklr.XX) ? 'T' : 'F');                                           \

      setBoolean(forceBOGunitigger)
      setBoolean(isNotRandom)
      setBoolean(doNotTrustHomopolymerRuns)

      setBoolean(doTrim_initialNone)
      setBoolean(doTrim_initialMerBased)
      setBoolean(doTrim_initialFlowBased)
      setBoolean(doTrim_initialQualityBased)

      setBoolean(doRemoveDuplicateReads)

      setBoolean(doTrim_finalLargestCovered)
      setBoolean(doTrim_finalEvidenceBased)

      setBoolean(doRemoveSpurReads)
      setBoolean(doRemoveChimericReads)
 
      setBoolean(doConsensusCorrection)

      setBoolean(forceShortReadFormat)

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
        if (verbose)
          fprintf(stdout, "lib uid %s orientation %s -> %s\n",
                  AS_UID_toString(gklr.libraryUID), AS_READ_ORIENT_NAMES[o], AS_READ_ORIENT_NAMES[gklr.orientation]);
        allFrags(gkpStore, IID, 'o', gklr.orientation, update, verbose);
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
        allFrags(gkpStore, IID, 'd', maketrue, update, verbose);
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
        allFrags(gkpStore, IID, 'r', maketrue, update, verbose);
        IID = 0;
      } else if (strcasecmp(ACT, "allfragsunmated") == 0) {
        uint32 maketrue = 0;
        if      ((E[0] == '1') || (E[0] == 't') || (E[0] == 'T')) {
          allFrags(gkpStore, IID, 'u', 1, update, verbose);
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
        gkpStore->gkStore_setLibrary(IID, &gklr);
    }


  nextline:
    fgets(L, 1024, F);
  }

  delete gkpStore;

  if (errors) {
    fprintf(stderr, "%d errors detected.\n", errors);
    exit(1);
  }

  fprintf(stderr, "Success!\n");

  exit(0);
}
