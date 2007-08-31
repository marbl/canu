
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

static char const *rcsid = "$Id: AS_GKP_rearrange.c,v 1.6 2007-08-31 21:06:16 brianwalenz Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "AS_global.h"
#include "AS_GKP_include.h"
#include "AS_PER_gkpStore.h"
#include "AS_PER_encodeSequenceQuality.h"

static char   encodedsequence[AS_FRAG_MAX_LEN+1] = {0};

void
addFragToStore(GateKeeperStore           *gkp,
               GateKeeperFragmentRecord  *gkf,
               char                      *seq,
               char                      *qlt,
               char                      *hps,
               char                      *src) {
  StoreStat   stats;

  gkf->seqLen = strlen(seq);
  gkf->hpsLen = strlen(hps);
  gkf->srcLen = strlen(src);

  statsStore(gkp->seq, &stats);
  gkf->seqOffset = stats.lastElem;

  statsStore(gkp->qlt, &stats);
  gkf->qltOffset = stats.lastElem;

  statsStore(gkp->hps, &stats);
  gkf->hpsOffset = stats.lastElem;

  statsStore(gkp->src, &stats);
  gkf->srcOffset = stats.lastElem;

  setGatekeeperUIDtoIID(gkp, gkf->readUID, gkf->readIID, AS_IID_FRG);
  appendIndexStore(gkp->frg, gkf);

  appendVLRecordStore(gkp->seq, seq, gkf->seqLen);

  encodeSequenceQuality(encodedsequence, seq, qlt);
  appendVLRecordStore(gkp->qlt, encodedsequence, gkf->seqLen);

  appendVLRecordStore(gkp->hps, hps, gkf->hpsLen);
  appendVLRecordStore(gkp->src, src, gkf->srcLen);
}



void
rearrangeStore(char *uidFileName, char *gkpStoreName, char *newStoreName) {

  if (testOpenGateKeeperStore(gkpStoreName, FALSE) == 0) {
    fprintf(stderr, "failed to open old store '%s', exit.\n", gkpStoreName);
    exit(1);
  }

  if (testOpenGateKeeperStore(newStoreName, FALSE)) {
    fprintf(stderr, "new store '%s' exists, exit.\n", newStoreName);
    exit(1);
  }

  GateKeeperStore *gkp      = NULL;
  GateKeeperStore *nst      = NULL;
  fragRecord       fr       = {0};
  uint32           i;

  gkp = openGateKeeperStore(gkpStoreName, FALSE);
  if (gkp == NULL) {
    fprintf(stderr, "Failed to open %s\n", gkpStoreName);
    exit(1);
  }

  CDS_IID_t        lastElem = getLastElemFragStore(gkp) + 1;
  FILE            *F        = NULL;
  char             L[1024];

  //  Wasteful...  Sigh.  order[] keeps a mapping of new iid -> old
  //  iid.  If order[i]==0 then missing[i] will contain the uid (used
  //  for UIDs that aren't in the original store, but are in the
  //  list).

  uint32           orderLen   = 0;
  CDS_IID_t       *newIIDorder      = (CDS_IID_t *)safe_calloc(lastElem, sizeof(CDS_IID_t));
  CDS_UID_t       *newIIDordertoUID = (CDS_UID_t *)safe_calloc(lastElem, sizeof(CDS_UID_t));
  CDS_IID_t       *oldIIDtonewIID   = (CDS_IID_t *)safe_calloc(lastElem, sizeof(CDS_IID_t));

  //  Read the list of UIDs, and generate a list of IIDs that we want to dump.

  errno = 0;
  if (strcmp(uidFileName, "-") == 0)
    F = stdin;
  else
    F = fopen(uidFileName, "r");
  if (errno) {
    fprintf(stderr, "couldn't open uid file '%s': %s\n", uidFileName, strerror(errno));
    exit(1);
  }
  fgets(L, 1024, F);
  while (!feof(F)) {
    CDS_UID_t      uid = STR_TO_UID(L, 0L, 10);
    CDS_IID_t      iid = getGatekeeperUIDtoIID(gkp, uid, NULL);

    newIIDorder[orderLen]       = iid;
    oldIIDtonewIID[iid]         = orderLen + 1;
    newIIDordertoUID[orderLen]  = uid;

    if (iid == 0)
      fprintf(stderr, "UID "F_UID" doesn't exist in original store, adding empty (deleted) fragment in it's place (new IID="F_IID").\n", uid, orderLen);

    orderLen++;

    fgets(L, 1024, F);
  }
  if (F != stdin)
    fclose(F);

  //  VERY important.  Old IID 0 is new IID 0, which we probably set
  //  to something else above.
  //
  oldIIDtonewIID[0] = 0;


  //
  //  Make a new store, stuff in all the old reads.  Simple, right?
  //

  nst = createGateKeeperStore(newStoreName);
  if (nst == NULL) {
    fprintf(stderr, "Failed to create %s\n", newStoreName);
    exit(1);
  }

  //  Stuff in all the libraries.  We can just copy from the old to
  //  the new.
  //
  for (i=1; getGateKeeperLibrary(gkp, i) != NULL; i++) {
    appendIndexStore(nst->lib, getGateKeeperLibrary(gkp, i));
    setGatekeeperUIDtoIID(nst, getGateKeeperLibrary(gkp, i)->libraryUID, getLastElemStore(nst->lib), AS_IID_LIB);
  }


  //  Stuff in all the reads.
  //
  for (i=0; i<orderLen; i++) {
    if (newIIDorder[i] > 0) {
      //  A valid fragment iid pointer
      getFrag(gkp, newIIDorder[i], &fr, FRAG_S_ALL);

      //fprintf(stderr, "i+1=%d readIID=%d oldIIDtonewIID=%d\n", i+1, fr.gkfr.readIID, oldIIDtonewIID[fr.gkfr.readIID]);
      assert(i+1             == oldIIDtonewIID[fr.gkfr.readIID]);

      //fprintf(stderr, "readUID "F_UID"  newIIDorderUID "F_UID"\n", fr.gkfr.readUID, newIIDordertoUID[i]);
      assert(fr.gkfr.readUID == newIIDordertoUID[i]);
    } else {
      //  Nope, fragment in our UID list isn't in the source store,
      //  add a deleted bogus fragment.
      //
      clr_fragRecord(&fr);
      fr.gkfr.readUID = newIIDordertoUID[i];
      fr.gkfr.deleted = 1;
    }

    //  Everybody gets a new IID.  Which is obvious once one remembers
    //  that we're adding new frags (deleted, yes, but still new) to
    //  our new store.
    //
    fr.gkfr.readIID = i + 1;
    fr.gkfr.mateIID = oldIIDtonewIID[fr.gkfr.mateIID];

    addFragToStore(nst, &fr.gkfr, fr.seq, fr.qlt, fr.hps, fr.src);
  }

  //  Done!

  closeGateKeeperStore(gkp);
  closeGateKeeperStore(nst);
}
