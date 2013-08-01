
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

static char *rcsid = "$Id: AS_PER_gkStore_clearRange.C,v 1.16 2011-08-22 04:49:11 brianwalenz Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>
#include <sys/stat.h>

#include "AS_global.H"
#include "AS_PER_genericStore.H"
#include "AS_PER_gkpStore.H"
#include "AS_PER_encodeSequenceQuality.H"
#include "AS_UTL_fileIO.H"


uint32
gkStore_decodeClearRegionLabel(const char *label) {
  uint32 i;
  for (i=0; i<AS_READ_CLEAR_NUM; i++)
    if (strcasecmp(label, AS_READ_CLEAR_NAMES[i]) == 0)
      return(i);
  return(AS_READ_CLEAR_ERROR);
}


static
char *
gkClearRange_makeName(gkStore *gkp, uint32 readType, uint32 clearType) {
  static char filePath[FILENAME_MAX];

  sprintf(filePath, "%s/clr-%s-%02d-%s",
          gkp->gkStore_path(),
          AS_READ_TYPE_NAMES[readType],
          clearType,
          AS_READ_CLEAR_NAMES[clearType]);

  return(filePath);
}


gkClearRange::gkClearRange(gkStore *gkp_, uint32 clearType_, uint32 create_) {

  gkp          = gkp_;

  clearType    = clearType_;
  create       = create_;

  pkconfigured = 0;
  pkdirty      = 0;
  pkmaxiid     = 0;
  pk           = NULL;

  nmconfigured = 0;
  nmdirty      = 0;
  nmmaxiid     = 0;
  nm           = NULL;

  sbconfigured = 0;
  sbdirty      = 0;
  sbmaxiid     = 0;
  sb           = NULL;
}


gkClearRange::~gkClearRange() {

  if ((pkdirty) && (pk != NULL)) {
    char *filePath = gkClearRange_makeName(gkp, GKFRAGMENT_PACKED, clearType);

    errno = 0;
    FILE *F = fopen(filePath, "w");
    if (errno)
      fprintf(stderr, "gkClearRange::~gkClearRange()-- failed to write clear range file '%s': %s\n", filePath, strerror(errno)), exit(1);

    AS_UTL_safeWrite(F, pk, "gkClearRange::~glClearRange()--sm", sizeof(uint8), gkp->inf.numPacked * 2 + 2);

    fclose(F);
  }

  if ((nmdirty) && (nm != NULL)) {
    char *filePath = gkClearRange_makeName(gkp, GKFRAGMENT_NORMAL, clearType);

    errno = 0;
    FILE *F = fopen(filePath, "w");
    if (errno)
      fprintf(stderr, "gkClearRange::~gkClearRange()-- failed to write clear range file '%s': %s\n", filePath, strerror(errno)), exit(1);

    AS_UTL_safeWrite(F, nm, "gkClearRange::~glClearRange()--nm", sizeof(uint16), gkp->inf.numNormal * 2 + 2);

    fclose(F);
  }

  if ((sbdirty) && (sb != NULL)) {
    char *filePath = gkClearRange_makeName(gkp, GKFRAGMENT_STROBE, clearType);

    errno = 0;
    FILE *F = fopen(filePath, "w");
    if (errno)
      fprintf(stderr, "gkClearRange::~gkClearRange()-- failed to write clear range file '%s': %s\n", filePath, strerror(errno)), exit(1);

    AS_UTL_safeWrite(F, sb, "gkClearRange::~glClearRange()--sb", sizeof(uint32), gkp->inf.numStrobe * 2 + 2);

    fclose(F);
  }

  delete [] pk;
  delete [] nm;
  delete [] sb;
}



void
gkClearRange::gkClearRange_purge(void) {
  char *filePath = NULL;

  filePath = gkClearRange_makeName(gkp, GKFRAGMENT_PACKED, clearType);
  if (AS_UTL_fileExists(filePath, FALSE, FALSE)) {
    fprintf(stderr, "gkStore: purging clear region '%s'\n", filePath);
    unlink(filePath);
  }

  filePath = gkClearRange_makeName(gkp, GKFRAGMENT_NORMAL, clearType);
  if (AS_UTL_fileExists(filePath, FALSE, FALSE)) {
    fprintf(stderr, "gkStore: purging clear region '%s'\n", filePath);
    unlink(filePath);
  }

  filePath = gkClearRange_makeName(gkp, GKFRAGMENT_STROBE, clearType);
  if (AS_UTL_fileExists(filePath, FALSE, FALSE)) {
    fprintf(stderr, "gkStore: purging clear region '%s'\n", filePath);
    unlink(filePath);
  }

  delete [] pk;
  delete [] nm;
  delete [] sb;

  pkconfigured = 0;
  pkdirty      = 0;
  pkmaxiid     = 0;
  pk           = NULL;

  nmconfigured = 0;
  nmdirty      = 0;
  nmmaxiid     = 0;
  nm           = NULL;

  sbconfigured = 0;
  sbdirty      = 0;
  sbmaxiid     = 0;
  sb           = NULL;
}


void
gkClearRange::gkClearRange_getClearRegion(gkFragment *fr, uint32& begin, uint32& end) {

  //  We can be configured (pk is valid) but still ask for a fragment
  //  outside the range.  If the first fragment (any early fragment)
  //  is added to the store with a vector clear, but no other
  //  fragments have a vector clear, the clear range is configured,
  //  but too small (we don't add undefined clear ranges).

  if (fr->type == GKFRAGMENT_PACKED) {
    if (!pkconfigured)
      gkClearRange_configurePacked();
    if ((pk) && (fr->tiid <= pkmaxiid)) {
      begin = pk[2*fr->tiid+0];
      end   = pk[2*fr->tiid+1];
    } else {
      begin = 1;
      end   = 0;
    }
  }
  if (fr->type == GKFRAGMENT_NORMAL) {
    if (!nmconfigured)
      gkClearRange_configureNormal();
    if ((nm) && (fr->tiid <= nmmaxiid)) {
      begin = nm[2*fr->tiid+0];
      end   = nm[2*fr->tiid+1];
    } else {
      begin = 1;
      end   = 0;
    }
  }
  if (fr->type == GKFRAGMENT_STROBE) {
    if (!sbconfigured)
      gkClearRange_configureStrobe();
    if ((sb) && (fr->tiid <= sbmaxiid)) {
      begin = sb[2*fr->tiid+0];
      end   = sb[2*fr->tiid+1];
    } else {
      begin = 1;
      end   = 0;
    }
  }
}



void
gkClearRange::gkClearRange_setClearRegion(gkFragment *fr, uint32  begin, uint32  end) {

  if (begin > end)
    fprintf(stderr, "ERROR: fragment %s,%d clear begin (%d) > end (%d).\n",
            AS_UID_toString(fr->gkFragment_getReadUID()), fr->gkFragment_getReadIID(), begin, end);

  if (begin > fr->gkFragment_getSequenceLength())
    fprintf(stderr, "ERROR: fragment %s,%d clear begin (%d) > sequence length (%d).\n",
            AS_UID_toString(fr->gkFragment_getReadUID()), fr->gkFragment_getReadIID(), begin, fr->gkFragment_getSequenceLength());

  if (end   > fr->gkFragment_getSequenceLength())
    fprintf(stderr, "ERROR: fragment %s,%d clear end (%d) > sequence length (%d).\n",
            AS_UID_toString(fr->gkFragment_getReadUID()), fr->gkFragment_getReadIID(), end, fr->gkFragment_getSequenceLength());

  assert(begin <= end);
  assert(begin <= fr->gkFragment_getSequenceLength());
  assert(end   <= fr->gkFragment_getSequenceLength());

  if (fr->type == GKFRAGMENT_PACKED) {
    if (!pkconfigured)
      gkClearRange_configurePacked();
    assert(pk != NULL);
    pkdirty = 1;
    pk[2*fr->tiid+0] = begin;
    pk[2*fr->tiid+1] = end;
  }
  if (fr->type == GKFRAGMENT_NORMAL) {
    if (nm == NULL)
      gkClearRange_configureNormal();
    assert(nm != NULL);
    nmdirty = 1;
    nm[2*fr->tiid+0] = begin;
    nm[2*fr->tiid+1] = end;
  }
  if (fr->type == GKFRAGMENT_STROBE) {
    if (sb == NULL)
      gkClearRange_configureStrobe();
    assert(sb != NULL);
    sbdirty = 1;
    sb[2*fr->tiid+0] = begin;
    sb[2*fr->tiid+1] = end;
  }
}





void
gkClearRange::gkClearRange_configurePacked(void) {
  char *filePath = gkClearRange_makeName(gkp, GKFRAGMENT_PACKED, clearType);

  if (AS_UTL_fileExists(filePath, FALSE, FALSE)) {
    pkdirty  = 0;
    pkmaxiid = gkp->inf.numPacked;
    pk       = NULL;

    errno = 0;
    FILE *F = fopen(filePath, "r");
    if (errno)
      fprintf(stderr, "gkClearRange::glClearRange()-- failed to load clear range file '%s': %s\n", filePath, strerror(errno)), exit(1);

    pk = new uint8 [pkmaxiid * 2 + 2];
    AS_UTL_safeRead(F, pk, "gkClearRange::glClearRange()--pk", sizeof(uint8), pkmaxiid * 2 + 2);

    fclose(F);
  } else if (create) {
    pkmaxiid = (gkp->inf.numPacked > 0) ? gkp->inf.numPacked : 1048576;
    pkdirty  = 1;
    pk       = new uint8 [pkmaxiid * 2 + 2];

    gkFragment  fr;

    fr.type = GKFRAGMENT_PACKED;

    //  Normal clear ranges copy the latest clear range on create.  The 'taint' region is invalid on create.

    for (uint32 iid=1; iid<=gkp->inf.numPacked; iid++) {
      getIndexStore(gkp->fpk, iid, &fr.fr.packed);
      pk[2*iid+0] = (clearType == AS_READ_CLEAR_TNT) ? 1 : fr.gkFragment_getClearRegionBegin();
      pk[2*iid+1] = (clearType == AS_READ_CLEAR_TNT) ? 0 : fr.gkFragment_getClearRegionEnd();
    }
  } else {
    //  Nothing to configure.  Not told to create, and file doesn't exist.  We'll return the
    //  invalid clear range (1,0) on any accesses.
  }

  pkconfigured = 1;
}


void
gkClearRange::gkClearRange_configureNormal(void) {
  char *filePath = gkClearRange_makeName(gkp, GKFRAGMENT_NORMAL, clearType);

  assert(nm == NULL);
  assert(nmconfigured == 0);

  if (AS_UTL_fileExists(filePath, FALSE, FALSE)) {
    nmdirty   = 0;
    nmmaxiid  = gkp->inf.numNormal;
    nm        = NULL;

    errno = 0;
    FILE *F = fopen(filePath, "r");
    if (errno)
      fprintf(stderr, "gkClearRange::glClearRange()-- failed to load clear range file '%s': %s\n", filePath, strerror(errno)), exit(1);

    nm = new uint16 [nmmaxiid * 2 + 2];
    AS_UTL_safeRead(F, nm, "gkClearRange::glClearRange()--nm", sizeof(uint16), nmmaxiid * 2 + 2);

    fclose(F);
  } else if (create) {
    nmmaxiid = (gkp->inf.numNormal > 0) ? gkp->inf.numNormal : 1048576;
    nmdirty  = 1;
    nm       = new uint16 [nmmaxiid * 2 + 2];

    gkFragment  fr;

    fr.type = GKFRAGMENT_NORMAL;

    for (uint32 iid=1; iid<=gkp->inf.numNormal; iid++) {
      getIndexStore(gkp->fnm, iid, &fr.fr.normal);
      nm[2*iid+0] = (clearType == AS_READ_CLEAR_TNT) ? 1 : fr.gkFragment_getClearRegionBegin();
      nm[2*iid+1] = (clearType == AS_READ_CLEAR_TNT) ? 0 : fr.gkFragment_getClearRegionEnd();
    }
  } else {
    //  Nothing to configure.  Not told to create, and file doesn't exist.  We'll return the
    //  invalid clear range (1,0) on any accesses.
  }

  nmconfigured = 1;
}


void
gkClearRange::gkClearRange_configureStrobe(void) {
  char *filePath = gkClearRange_makeName(gkp, GKFRAGMENT_STROBE, clearType);

  if (AS_UTL_fileExists(filePath, FALSE, FALSE)) {
    sbdirty   = 0;
    sbmaxiid  = gkp->inf.numStrobe;
    sb        = NULL;

    errno = 0;
    FILE *F = fopen(filePath, "r");
    if (errno)
      fprintf(stderr, "gkClearRange::glClearRange()-- failed to load clear range file '%s': %s\n", filePath, strerror(errno)), exit(1);

    sb = new uint32 [sbmaxiid * 2 + 2];
    AS_UTL_safeRead(F, sb, "gkClearRange::glClearRange()--sb", sizeof(uint32), sbmaxiid * 2 + 2);

    fclose(F);
  } else if (create) {
    sbmaxiid = (gkp->inf.numStrobe > 0) ? gkp->inf.numStrobe : 1048576;
    sbdirty  = 1;
    sb       = new uint32 [sbmaxiid * 2 + 2];

    gkFragment  fr;

    fr.type = GKFRAGMENT_STROBE;

    for (uint32 iid=1; iid<=gkp->inf.numStrobe; iid++) {
      getIndexStore(gkp->fsb, iid, &fr.fr.strobe);
      sb[2*iid+0] = (clearType == AS_READ_CLEAR_TNT) ? 1 : fr.gkFragment_getClearRegionBegin();
      sb[2*iid+1] = (clearType == AS_READ_CLEAR_TNT) ? 0 : fr.gkFragment_getClearRegionEnd();
    }
  } else {
    //  Nothing to configure.  Not told to create, and file doesn't exist.  We'll return the
    //  invalid clear range (1,0) on any accesses.
  }

  sbconfigured = 1;
}




void
gkClearRange::gkClearRange_makeSpacePacked(AS_IID tiid, uint32 bgn, uint32 end) {

  if (tiid <= pkmaxiid)
    //  Already got space.
    return;

  if ((bgn > end) && (pk == NULL))
    //  Nobody has used this clear range yet, and we don't want it either.
    return;

  uint64  newpkmaxiid = (pkmaxiid == 0) ? (1048576) : (pkmaxiid * 2);
  uint8  *newpk       = NULL;

  while (newpkmaxiid < tiid)
    newpkmaxiid *= 2;

  newpk = new uint8 [newpkmaxiid * 2 + 2];

  if (pkmaxiid > 0)
    memcpy(newpk, pk, sizeof(uint8) * pkmaxiid * 2 + 2);

  for (uint64 i=pkmaxiid; i <= newpkmaxiid; i++) {
    newpk[2*i+0] = 1;
    newpk[2*i+1] = 0;
  }

  delete [] pk;
  pk       = newpk;
  pkmaxiid = newpkmaxiid;

  pkconfigured = 1;
}


void
gkClearRange::gkClearRange_makeSpaceNormal(AS_IID tiid, uint32 bgn, uint32 end) {

  if (tiid < nmmaxiid)
    return;

  if ((bgn > end) && (nm == NULL))
    return;

  uint64  newnmmaxiid = (nmmaxiid == 0) ? (1048576) : (nmmaxiid * 2);
  uint16  *newnm      = NULL;

  while (newnmmaxiid < tiid)
    newnmmaxiid *= 2;

  newnm = new uint16 [newnmmaxiid * 2 + 2];

  if (nmmaxiid > 0)
    memcpy(newnm, nm, sizeof(uint16) * nmmaxiid * 2 + 2);

  for (uint64 i=nmmaxiid; i <= newnmmaxiid; i++) {
    newnm[2*i+0] = 1;
    newnm[2*i+1] = 0;
  }

  delete [] nm;
  nm       = newnm;
  nmmaxiid = newnmmaxiid;

  nmconfigured = 1;
}


void
gkClearRange::gkClearRange_makeSpaceStrobe(AS_IID tiid, uint32 bgn, uint32 end) {

  if (tiid < sbmaxiid)
    return;

  if ((bgn > end) && (sb == NULL))
    return;

  uint64  newsbmaxiid = (sbmaxiid == 0) ? (1048576) : (sbmaxiid * 2);
  uint32  *newsb      = NULL;

  while (newsbmaxiid < tiid)
    newsbmaxiid *= 2;

  newsb = new uint32 [newsbmaxiid * 2 + 2];

  if (sbmaxiid > 0)
    memcpy(newsb, sb, sizeof(uint32) * sbmaxiid * 2 + 2);

  for (uint64 i=sbmaxiid; i <= newsbmaxiid; i++) {
    newsb[2*i+0] = 1;
    newsb[2*i+1] = 0;
  }

  delete [] sb;
  sb       = newsb;
  sbmaxiid = newsbmaxiid;

  sbconfigured = 1;
}






void
gkStore::gkStore_enableClearRange(uint32 which) {
  assert(partmap == NULL);
  clearRange[which]->gkClearRange_enableCreate();
}


void
gkStore::gkStore_purgeClearRange(uint32 which) {
  assert(partmap == NULL);
  clearRange[which]->gkClearRange_purge();
}
