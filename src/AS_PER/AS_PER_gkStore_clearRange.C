
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

static char *rcsid = "$Id: AS_PER_gkStore_clearRange.C,v 1.6 2009-10-26 13:20:26 brianwalenz Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>
#include <sys/stat.h>

#include "AS_global.h"
#include "AS_PER_genericStore.h"
#include "AS_PER_gkpStore.h"
#include "AS_PER_encodeSequenceQuality.h"
#include "AS_UTL_fileIO.h"


uint32
gkStore_decodeClearRegionLabel(const char *label) {
  uint32 i;
  for (i=0; i<AS_READ_CLEAR_NUM; i++)
    if (strcasecmp(label, AS_READ_CLEAR_NAMES[i]) == 0)
      return(i);
  return(AS_READ_CLEAR_ERROR);
};


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

  smconfigured = 0;
  smdirty      = 0;
  smmaxiid     = 0;
  sm           = NULL;

  mdconfigured = 0;
  mddirty      = 0;
  mdmaxiid     = 0;
  md           = NULL;

  lgconfigured = 0;
  lgdirty      = 0;
  lgmaxiid     = 0;
  lg           = NULL;
}


gkClearRange::~gkClearRange() {

  if ((smdirty) && (sm != NULL)) {
    char *filePath = gkClearRange_makeName(gkp, GKFRAGMENT_PACKED, clearType);

    errno = 0;
    FILE *F = fopen(filePath, "w");
    if (errno)
      fprintf(stderr, "gkClearRange::~gkClearRange()-- failed to write clear range file '%s': %s\n", filePath, strerror(errno)), exit(1);

    AS_UTL_safeWrite(F, sm, "gkClearRange::~glClearRange()--sm", sizeof(uint8), gkp->inf.numShort * 2 + 2);

    fclose(F);
  }

  if ((mddirty) && (md != NULL)) {
    char *filePath = gkClearRange_makeName(gkp, GKFRAGMENT_NORMAL, clearType);

    errno = 0;
    FILE *F = fopen(filePath, "w");
    if (errno)
      fprintf(stderr, "gkClearRange::~gkClearRange()-- failed to write clear range file '%s': %s\n", filePath, strerror(errno)), exit(1);

    AS_UTL_safeWrite(F, md, "gkClearRange::~glClearRange()--md", sizeof(uint16), gkp->inf.numMedium * 2 + 2);

    fclose(F);
  }

  if ((lgdirty) && (lg != NULL)) {
    char *filePath = gkClearRange_makeName(gkp, GKFRAGMENT_STROBE, clearType);

    errno = 0;
    FILE *F = fopen(filePath, "w");
    if (errno)
      fprintf(stderr, "gkClearRange::~gkClearRange()-- failed to write clear range file '%s': %s\n", filePath, strerror(errno)), exit(1);

    AS_UTL_safeWrite(F, lg, "gkClearRange::~glClearRange()--lg", sizeof(uint32), gkp->inf.numLong * 2 + 2);

    fclose(F);
  }

  delete [] sm;
  delete [] md;
  delete [] lg;
};



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

  delete [] sm;
  delete [] md;
  delete [] lg;

  smconfigured = 0;
  smdirty      = 0;
  smmaxiid     = 0;
  sm           = NULL;

  mdconfigured = 0;
  mddirty      = 0;
  mdmaxiid     = 0;
  md           = NULL;

  lgconfigured = 0;
  lgdirty      = 0;
  lgmaxiid     = 0;
  lg           = NULL;
}


void
gkClearRange::gkClearRange_getClearRegion(gkFragment *fr, uint32& begin, uint32& end) {

  //  We can be configured (sm is valid) but still ask for a fragment
  //  outside the range.  If the first fragment (any early fragment)
  //  is added to the store with a vector clear, but no other
  //  fragments have a vector clear, the clear range is configured,
  //  but too small (we don't add undefined clear ranges).

  if (fr->type == GKFRAGMENT_PACKED) {
    if (!smconfigured)
      gkClearRange_configureShort();
    if ((sm) && (fr->tiid <= smmaxiid)) {
      begin = sm[2*fr->tiid+0];
      end   = sm[2*fr->tiid+1];
    } else {
      begin = 1;
      end   = 0;
    }
  }
  if (fr->type == GKFRAGMENT_NORMAL) {
    if (!mdconfigured)
      gkClearRange_configureMedium();
    if ((md) && (fr->tiid <= mdmaxiid)) {
      begin = md[2*fr->tiid+0];
      end   = md[2*fr->tiid+1];
    } else {
      begin = 1;
      end   = 0;
    }
  }
  if (fr->type == GKFRAGMENT_STROBE) {
    if (!lgconfigured)
      gkClearRange_configureLong();
    if ((lg) && (fr->tiid <= lgmaxiid)) {
      begin = lg[2*fr->tiid+0];
      end   = lg[2*fr->tiid+1];
    } else {
      begin = 1;
      end   = 0;
    }
  }
}



void
gkClearRange::gkClearRange_setClearRegion(gkFragment *fr, uint32  begin, uint32  end) {
  if (fr->type == GKFRAGMENT_PACKED) {
    if (!smconfigured)
      gkClearRange_configureShort();
    assert(sm != NULL);
    smdirty = 1;
    sm[2*fr->tiid+0] = begin;
    sm[2*fr->tiid+1] = end;
  }
  if (fr->type == GKFRAGMENT_NORMAL) {
    if (md == NULL)
      gkClearRange_configureMedium();
    assert(md != NULL);
    mddirty = 1;
    md[2*fr->tiid+0] = begin;
    md[2*fr->tiid+1] = end;
  }
  if (fr->type == GKFRAGMENT_STROBE) {
    if (lg == NULL)
      gkClearRange_configureLong();
    assert(lg != NULL);
    lgdirty = 1;
    lg[2*fr->tiid+0] = begin;
    lg[2*fr->tiid+1] = end;
  }
}





void
gkClearRange::gkClearRange_configureShort(void) {
  char *filePath = gkClearRange_makeName(gkp, GKFRAGMENT_PACKED, clearType);

  if (AS_UTL_fileExists(filePath, FALSE, TRUE)) {
    smdirty  = 0;
    smmaxiid = gkp->inf.numShort;
    sm       = NULL;

    errno = 0;
    FILE *F = fopen(filePath, "r");
    if (errno)
      fprintf(stderr, "gkClearRange::glClearRange()-- failed to load clear range file '%s': %s\n", filePath, strerror(errno)), exit(1);

    sm = new uint8 [smmaxiid * 2 + 2];
    AS_UTL_safeRead(F, sm, "gkClearRange::glClearRange()--sm", sizeof(uint8), smmaxiid * 2 + 2);

    fclose(F);
  } else if (create) {
    smmaxiid = (gkp->inf.numShort > 0) ? gkp->inf.numShort : 1048576;
    smdirty  = 1;
    sm       = new uint8 [smmaxiid * 2 + 2];

    gkFragment  fr;
    gkStream   *gs = new gkStream(gkp,
                                  1,
                                  gkp->inf.numShort,
                                  GKFRAGMENT_INF);
    while (gs->next(&fr)) {
      sm[2*fr.tiid+0] = fr.gkFragment_getClearRegionBegin();
      sm[2*fr.tiid+1] = fr.gkFragment_getClearRegionEnd();
    }
    delete gs;
  }

  smconfigured = 1;
}


void
gkClearRange::gkClearRange_configureMedium(void) {
  char *filePath = gkClearRange_makeName(gkp, GKFRAGMENT_NORMAL, clearType);

  assert(md == NULL);
  assert(mdconfigured == 0);

  if (AS_UTL_fileExists(filePath, FALSE, TRUE)) {
    mddirty   = 0;
    mdmaxiid  = gkp->inf.numMedium;
    md        = NULL;

    errno = 0;
    FILE *F = fopen(filePath, "r");
    if (errno)
      fprintf(stderr, "gkClearRange::glClearRange()-- failed to load clear range file '%s': %s\n", filePath, strerror(errno)), exit(1);

    md = new uint16 [mdmaxiid * 2 + 2];
    AS_UTL_safeRead(F, md, "gkClearRange::glClearRange()--md", sizeof(uint16), mdmaxiid * 2 + 2);

    fclose(F);
  } else if (create) {
    mdmaxiid = (gkp->inf.numMedium > 0) ? gkp->inf.numMedium : 1048576;
    mddirty  = 1;
    md       = new uint16 [mdmaxiid * 2 + 2];

    gkFragment  fr;
    gkStream   *gs = new gkStream(gkp,
                                  gkp->inf.numShort + 1,
                                  gkp->inf.numShort + gkp->inf.numMedium,
                                  GKFRAGMENT_INF);
    while (gs->next(&fr)) {
      md[2*fr.tiid+0] = fr.gkFragment_getClearRegionBegin();
      md[2*fr.tiid+1] = fr.gkFragment_getClearRegionEnd();
    }
    delete gs;
  }

  mdconfigured = 1;
}


void
gkClearRange::gkClearRange_configureLong(void) {
  char *filePath = gkClearRange_makeName(gkp, GKFRAGMENT_STROBE, clearType);

  if (AS_UTL_fileExists(filePath, FALSE, TRUE)) {
    lgdirty   = 0;
    lgmaxiid  = gkp->inf.numLong;
    lg        = NULL;

    errno = 0;
    FILE *F = fopen(filePath, "r");
    if (errno)
      fprintf(stderr, "gkClearRange::glClearRange()-- failed to load clear range file '%s': %s\n", filePath, strerror(errno)), exit(1);

    lg = new uint32 [lgmaxiid * 2 + 2];
    AS_UTL_safeRead(F, lg, "gkClearRange::glClearRange()--lg", sizeof(uint32), lgmaxiid * 2 + 2);

    fclose(F);
  } else if (create) {
    lgmaxiid = (gkp->inf.numLong > 0) ? gkp->inf.numLong : 1048576;
    lgdirty  = 1;
    lg       = new uint32 [lgmaxiid * 2 + 2];

    gkFragment  fr;
    gkStream   *gs = new gkStream(gkp,
                                  gkp->inf.numShort + gkp->inf.numMedium + 1,
                                  gkp->inf.numShort + gkp->inf.numMedium + gkp->inf.numLong,
                                  GKFRAGMENT_INF);
    while (gs->next(&fr)) {
      lg[2*fr.tiid+0] = fr.gkFragment_getClearRegionBegin();
      lg[2*fr.tiid+1] = fr.gkFragment_getClearRegionEnd();
    }
    delete gs;
  }

  lgconfigured = 1;
}




void
gkClearRange::gkClearRange_makeSpaceShort(AS_IID tiid, uint32 bgn, uint32 end) {

  if (tiid <= smmaxiid)
    //  Already got space.
    return;

  if ((bgn > end) && (sm == NULL))
    //  Nobody has used this clear range yet, and we don't want it either.
    return;

  uint32  newsmmaxiid = (smmaxiid == 0) ? (1048576) : (smmaxiid * 2);
  uint8  *newsm       = NULL;

  while (newsmmaxiid < tiid)
    newsmmaxiid *= 2;

  newsm = new uint8 [newsmmaxiid * 2 + 2];

  if (smmaxiid > 0)
    memcpy(newsm, sm, sizeof(uint8) * smmaxiid * 2 + 2);

  for (uint32 i=smmaxiid; i <= newsmmaxiid; i++) {
    newsm[2*i+0] = 1;
    newsm[2*i+1] = 0;
  }

  delete [] sm;
  sm       = newsm;
  smmaxiid = newsmmaxiid;

  smconfigured = 1;
}


void
gkClearRange::gkClearRange_makeSpaceMedium(AS_IID tiid, uint32 bgn, uint32 end) {

  if (tiid < mdmaxiid)
    return;

  if ((bgn > end) && (md == NULL))
    return;

  uint32  newmdmaxiid = (mdmaxiid == 0) ? (1048576) : (mdmaxiid * 2);
  uint16  *newmd      = NULL;

  while (newmdmaxiid < tiid)
    newmdmaxiid *= 2;

  newmd = new uint16 [newmdmaxiid * 2 + 2];

  if (mdmaxiid > 0)
    memcpy(newmd, md, sizeof(uint16) * mdmaxiid * 2 + 2);

  for (uint32 i=mdmaxiid; i <= newmdmaxiid; i++) {
    newmd[2*i+0] = 1;
    newmd[2*i+1] = 0;
  }

  delete [] md;
  md       = newmd;
  mdmaxiid = newmdmaxiid;

  mdconfigured = 1;
}


void
gkClearRange::gkClearRange_makeSpaceLong(AS_IID tiid, uint32 bgn, uint32 end) {

  if (tiid < lgmaxiid)
    return;

  if ((bgn > end) && (lg == NULL))
    return;

  uint32  newlgmaxiid = (lgmaxiid == 0) ? (1048576) : (lgmaxiid * 2);
  uint32  *newlg      = NULL;

  while (newlgmaxiid < tiid)
    newlgmaxiid *= 2;

  newlg = new uint32 [newlgmaxiid * 2 + 2];

  if (lgmaxiid > 0)
    memcpy(newlg, lg, sizeof(uint32) * lgmaxiid * 2 + 2);

  for (uint32 i=lgmaxiid; i <= newlgmaxiid; i++) {
    newlg[2*i+0] = 1;
    newlg[2*i+1] = 0;
  }

  delete [] lg;
  lg       = newlg;
  lgmaxiid = newlgmaxiid;

  lgconfigured = 1;
}
