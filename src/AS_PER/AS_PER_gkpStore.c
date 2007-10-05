
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

static char CM_ID[] = "$Id: AS_PER_gkpStore.c,v 1.41 2007-10-05 06:25:50 brianwalenz Exp $";

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

static
int
fileExists(const char *path,
           int directory,
           int readwrite) {
  struct stat  s;
  int          r;

  errno = 0;
  r = stat(path, &s);
  if (errno) {
    //fprintf(stderr, "Failed to stat() file '%s'; assumed to not exist.\n", path);
    return(0);
  }

  if (directory == 1) {
    if ((readwrite == 0) &&
        (s.st_mode & S_IFDIR) &&
        (s.st_mode & (S_IRUSR | S_IRGRP | S_IROTH)) &&
        (s.st_mode & (S_IXUSR | S_IXGRP | S_IXOTH)))
      return(1);
    if ((readwrite == 1) &&
        (s.st_mode & S_IFDIR) &&
        (s.st_mode & (S_IRUSR | S_IRGRP | S_IROTH)) &&
        (s.st_mode & (S_IWUSR | S_IWGRP | S_IWOTH)) &&
        (s.st_mode & (S_IXUSR | S_IXGRP | S_IXOTH)))
      return(1);
  }

  if (directory == 0) {
    if ((readwrite == 0) &&
        (s.st_mode & (S_IRUSR | S_IRGRP | S_IROTH)))
      return(1);
    if ((readwrite == 1) &&
        (s.st_mode & (S_IRUSR | S_IRGRP | S_IROTH)) &&
        (s.st_mode & (S_IWUSR | S_IWGRP | S_IWOTH)))
      return(1);
  }

  return(0);
}


int
testOpenGateKeeperStore(const char *path,
                        int   writable) {

  char name[FILENAME_MAX];
  int  fileCount = 0;

  if (fileExists(path, 1, writable)) {
    sprintf(name,"%s/gkp", path);
    fileCount += fileExists(name, 0, writable);

    sprintf(name,"%s/bat", path);
    fileCount += fileExists(name, 0, writable);

    sprintf(name,"%s/frg", path);
    fileCount += fileExists(name, 0, writable);

    sprintf(name,"%s/lib", path);
    fileCount += fileExists(name, 0, writable);

    sprintf(name,"%s/seq", path);
    fileCount += fileExists(name, 0, writable);

    sprintf(name,"%s/qlt", path);
    fileCount += fileExists(name, 0, writable);

    sprintf(name,"%s/hps", path);
    fileCount += fileExists(name, 0, writable);

    sprintf(name,"%s/src", path);
    fileCount += fileExists(name, 0, writable);

    sprintf(name,"%s/map", path);
    fileCount += fileExists(name, 0, writable);

  }

  return(fileCount == NUM_GKP_FILES);
}


GateKeeperStore *
openGateKeeperStore(const char *path,
                    int   writable) {

  char              name[FILENAME_MAX];
  FILE             *gkpinfo;

  GateKeeperStore  *gkpStore = (GateKeeperStore *)safe_calloc(1, sizeof(GateKeeperStore));

  strcpy(gkpStore->storePath, path);

  gkpStore->bat = NULL;
  gkpStore->frg = NULL;
  gkpStore->lib = NULL;

  gkpStore->seq = NULL;
  gkpStore->qlt = NULL;
  gkpStore->hps = NULL;
  gkpStore->src = NULL;

  gkpStore->UIDtoIID = NULL;
  
  gkpStore->partnum = -1;
  gkpStore->partfrg = NULL;
  gkpStore->partqlt = NULL;
  gkpStore->parthps = NULL;
  gkpStore->partsrc = NULL;
  gkpStore->partmap = NULL;


  sprintf(name,"%s/gkp", gkpStore->storePath);
  errno = 0;
  gkpinfo = fopen(name, "r");
  if (errno) {
    fprintf(stderr, "failed to open gatekeeper store '%s': %s\n", name, strerror(errno));
    exit(1);
  }

  if (1 != AS_UTL_safeRead(gkpinfo, &gkpStore->gkp, "openGateKeeperStore:header", sizeof(GateKeeperStoreInfo), 1)) {
    fprintf(stderr, "failed to open gatekeeper store '%s': couldn't read the header (%s)\n", name, strerror(errno));
    exit(1);
  }
  fclose(gkpinfo);

  if (gkpStore->gkp.gkpMagic != 1) {
    fprintf(stderr, "invalid magic!\n");
  }
  if (gkpStore->gkp.gkpVersion != 1) {
    fprintf(stderr, "invalid version!\n");
  }
  if (gkpStore->gkp.gkpBatchRecordSize != sizeof(GateKeeperBatchRecord)) {
    fprintf(stderr, "ERROR!  Store built unsing GateKeeperBatchRecord of size %d bytes.\n", gkpStore->gkp.gkpBatchRecordSize);
    fprintf(stderr, "        Code compiled with GateKeeperBatchRecord of size %d bytes.\n", (int)sizeof(GateKeeperBatchRecord));
    assert(gkpStore->gkp.gkpBatchRecordSize == sizeof(GateKeeperBatchRecord));
  }
  if (gkpStore->gkp.gkpLibraryRecordSize != sizeof(GateKeeperLibraryRecord)) {
    fprintf(stderr, "ERROR!  Store built unsing GateKeeperLibraryRecord of size %d bytes.\n", gkpStore->gkp.gkpLibraryRecordSize);
    fprintf(stderr, "        Code compiled with GateKeeperLibraryRecord of size %d bytes.\n", (int)sizeof(GateKeeperLibraryRecord));
    assert(gkpStore->gkp.gkpLibraryRecordSize == sizeof(GateKeeperLibraryRecord));
  }
  if (gkpStore->gkp.gkpFragmentRecordSize != sizeof(GateKeeperFragmentRecord)) {
    fprintf(stderr, "ERROR!  Store built unsing GateKeeperFragmentRecord of size %d bytes.\n", gkpStore->gkp.gkpFragmentRecordSize);
    fprintf(stderr, "        Code compiled with GateKeeperFragmentRecord of size %d bytes.\n", (int)sizeof(GateKeeperFragmentRecord));
    assert(gkpStore->gkp.gkpFragmentRecordSize == sizeof(GateKeeperFragmentRecord));
  }


  char  mode[4];
  if (writable)
    strcpy(mode, "r+");
  else
    strcpy(mode, "r");

  if (writable >= 0) {
    sprintf(name,"%s/bat", gkpStore->storePath);
    gkpStore->bat   = openStore(name, mode);

    sprintf(name,"%s/frg", gkpStore->storePath);
    gkpStore->frg   = openStore(name, mode);

    sprintf(name,"%s/lib", gkpStore->storePath);
    gkpStore->lib   = openStore(name, mode);

    sprintf(name,"%s/seq", gkpStore->storePath);
    gkpStore->seq = openStore(name, mode);

    sprintf(name,"%s/qlt", gkpStore->storePath);
    gkpStore->qlt = openStore(name, mode);

    sprintf(name,"%s/hps", gkpStore->storePath);
    gkpStore->hps = openStore(name, mode);

    sprintf(name,"%s/src", gkpStore->storePath);
    gkpStore->src = openStore(name, mode);

    if ((NULL == gkpStore->bat) ||
        (NULL == gkpStore->frg) ||
        (NULL == gkpStore->lib) ||
        (NULL == gkpStore->seq) ||
        (NULL == gkpStore->qlt) ||
        (NULL == gkpStore->hps) ||
        (NULL == gkpStore->src)) {
      fprintf(stderr,"**** Failure to open Gatekeeper Store ...\n");
      assert(0);
    }

    if (writable) {
      char  name[FILENAME_MAX];
      sprintf(name,"%s/map", gkpStore->storePath);
      gkpStore->UIDtoIID = LoadUIDtoIIDHashTable_AS(name);
    }
  }

  return(gkpStore);
}


GateKeeperStore *
createGateKeeperStore(const char *path) {
  char   name[FILENAME_MAX];
  FILE  *gkpinfo;

  GateKeeperStore  *gkpStore = (GateKeeperStore *)safe_calloc(1, sizeof(GateKeeperStore));

  strcpy(gkpStore->storePath, path);

  gkpStore->bat = NULL;
  gkpStore->frg = NULL;
  gkpStore->lib = NULL;

  gkpStore->seq = NULL;
  gkpStore->qlt = NULL;
  gkpStore->hps = NULL;
  gkpStore->src = NULL;

  gkpStore->UIDtoIID = NULL;

  gkpStore->partnum = -1;
  gkpStore->partfrg = NULL;
  gkpStore->partqlt = NULL;
  gkpStore->parthps = NULL;
  gkpStore->partsrc = NULL;
  gkpStore->partmap = NULL;

  AS_UTL_mkdir(path);

  gkpStore->gkp.gkpMagic              = 1;
  gkpStore->gkp.gkpVersion            = 1;
  gkpStore->gkp.gkpBatchRecordSize    = sizeof(GateKeeperBatchRecord);
  gkpStore->gkp.gkpLibraryRecordSize  = sizeof(GateKeeperLibraryRecord);
  gkpStore->gkp.gkpFragmentRecordSize = sizeof(GateKeeperFragmentRecord);

  sprintf(name,"%s/gkp", path);
  errno = 0;
  gkpinfo = fopen(name, "w");
  if (errno) {
    fprintf(stderr, "failed to create gatekeeper store '%s': %s\n", name, strerror(errno));
    exit(1);
  }

  AS_UTL_safeWrite(gkpinfo, &gkpStore->gkp, "createGateKeeperStore:header", sizeof(GateKeeperStoreInfo), 1);
  if (fclose(gkpinfo)) {
    fprintf(stderr, "failed to create gatekeeper store '%s': %s\n", name, strerror(errno));
    exit(1);
  }

  sprintf(name,"%s/bat", path);
  gkpStore->bat = createIndexStore(name, "bat", sizeof(GateKeeperBatchRecord), 1);

  sprintf(name,"%s/frg", path);
  gkpStore->frg = createIndexStore(name, "frg", sizeof(GateKeeperFragmentRecord), 1);

  sprintf(name,"%s/lib", path);
  gkpStore->lib = createIndexStore(name, "lib", sizeof(GateKeeperLibraryRecord), 1);

  sprintf(name,"%s/seq", path);
  gkpStore->seq = createVLRecordStore(name, "seq", MAX_SEQ_LENGTH);

  sprintf(name,"%s/qlt", path);
  gkpStore->qlt = createVLRecordStore(name, "qlt", MAX_SEQ_LENGTH);

  sprintf(name,"%s/hps", path);
  gkpStore->hps = createVLRecordStore(name, "hps", MAX_HPS_LENGTH);

  sprintf(name,"%s/src", path);
  gkpStore->src = createVLRecordStore(name, "src", MAX_SRC_LENGTH);

  sprintf(name,"%s/map", path);
  gkpStore->UIDtoIID = CreateScalarHashTable_AS(32 * 1024);
  SaveHashTable_AS(name, gkpStore->UIDtoIID);

  return(gkpStore);
}


void
closeGateKeeperStore(GateKeeperStore *gkpStore) {
  char  name[FILENAME_MAX];
  FILE *gkpinfo;

  if (gkpStore == NULL)
    return;

  sprintf(name,"%s/gkp", gkpStore->storePath);
  errno = 0;
  gkpinfo = fopen(name, "w");
  if (errno) {
    fprintf(stderr, "failed to write gatekeeper store into to '%s': %s\n", name, strerror(errno));
    exit(1);
  }

  AS_UTL_safeWrite(gkpinfo, &gkpStore->gkp, "closeGateKeeperStore:header", sizeof(GateKeeperStoreInfo), 1);
  if (fclose(gkpinfo)) {
    fprintf(stderr, "failed to close gatekeeper store '%s': %s\n", name, strerror(errno));
    exit(1);
  }

  if(gkpStore->bat != NULL)
    closeStore(gkpStore->bat);

  if(gkpStore->frg != NULL)
    closeStore(gkpStore->frg);

  if(gkpStore->lib != NULL)
    closeStore(gkpStore->lib);

  if(gkpStore->seq != NULL)
    closeStore(gkpStore->seq);

  if(gkpStore->qlt != NULL)
    closeStore(gkpStore->qlt);

  if(gkpStore->hps != NULL)
    closeStore(gkpStore->hps);

  if(gkpStore->src != NULL)
    closeStore(gkpStore->src);

  if(gkpStore->UIDtoIID != NULL)
    DeleteHashTable_AS(gkpStore->UIDtoIID);

  if(gkpStore->partfrg != NULL)
    closeStore(gkpStore->partfrg);

  if(gkpStore->partqlt != NULL)
    closeStore(gkpStore->partqlt);

  if(gkpStore->parthps != NULL)
    closeStore(gkpStore->parthps);

  if(gkpStore->partsrc != NULL)
    closeStore(gkpStore->partsrc);

  if (gkpStore->partmap != NULL)
    DeleteHashTable_AS(gkpStore->partmap);


  safe_free(gkpStore);
}



////////////////////////////////////////////////////////////////////////////////




GateKeeperStore *createGateKeeperPartition(const char *path, uint32 partnum) {
  char       name[FILENAME_MAX];

  GateKeeperStore *gkp = openGateKeeperStore(path, -1);

  gkp->partnum = partnum;

  sprintf(name,"%s/frg.%03d", gkp->storePath, partnum);
  gkp->partfrg = createIndexStore(name, "partfrg", sizeof(GateKeeperFragmentRecord), 1);

  sprintf(name,"%s/qlt.%03d", gkp->storePath, partnum);
  gkp->partqlt = createVLRecordStore(name, "partqlt", MAX_SEQ_LENGTH);

  sprintf(name,"%s/hps.%03d", gkp->storePath, partnum);
  gkp->parthps = createVLRecordStore(name, "parthps", MAX_HPS_LENGTH);

  sprintf(name,"%s/src.%03d", gkp->storePath, partnum);
  gkp->partsrc = createVLRecordStore(name, "partsrc", MAX_SRC_LENGTH);

  return(gkp);
}




void       loadGateKeeperPartition(GateKeeperStore *gkp, uint32 partnum) {
  char       name[FILENAME_MAX];
  StoreStat  stats;
  int        i;

  if (gkp->partnum != -1) {
    fprintf(stderr, "WARNING:  Throwing out partition %d to load %d\n",
            gkp->partnum, partnum);

    if(gkp->partfrg != NULL)
      closeStore(gkp->partfrg);

    if(gkp->partqlt != NULL)
      closeStore(gkp->partqlt);

    if(gkp->parthps != NULL)
      closeStore(gkp->parthps);

    if(gkp->partsrc != NULL)
      closeStore(gkp->partsrc);

    if (gkp->partmap != NULL)
      DeleteHashTable_AS(gkp->partmap);
  }

  sprintf(name,"%s/frg.%03d", gkp->storePath, partnum);
  if (fileExists(name, 0, FALSE) == 0) {
    fprintf(stderr, "loadGateKeeperPartition()--  Partition %d doesn't exist; normal store used instead.\n", partnum);
    return;
  }

  gkp->partnum = partnum;

  //  load all our data

  sprintf(name,"%s/frg.%03d", gkp->storePath, partnum);
  gkp->partfrg = loadStorePartial(name, 0, 0);

  sprintf(name,"%s/qlt.%03d", gkp->storePath, partnum);
  gkp->partqlt = loadStorePartial(name, 0, 0);

  sprintf(name,"%s/hps.%03d", gkp->storePath, partnum);
  gkp->parthps = loadStorePartial(name, 0, 0);

  sprintf(name,"%s/src.%03d", gkp->storePath, partnum);
  gkp->partsrc = loadStorePartial(name, 0, 0);

  //  zip through the frg and build a map from iid to the frg record

  statsStore(gkp->partfrg, &stats);

  gkp->partmap = CreateScalarHashTable_AS(stats.lastElem + 1);

  for(i = stats.firstElem; i <= stats.lastElem; i++) {
    GateKeeperFragmentRecord *p = getIndexStorePtr(gkp->partfrg, i);

    if (InsertInHashTable_AS(gkp->partmap,
                             (uint64)p->readIID, 0,
                             (uint64)p, 0) != HASH_SUCCESS)
      assert(0);
  }
}


////////////////////////////////////////////////////////////////////////////////


void clearGateKeeperBatchRecord(GateKeeperBatchRecord *g) {
  memset(g, 0, sizeof(GateKeeperBatchRecord));
}

void clearGateKeeperLibraryRecord(GateKeeperLibraryRecord *g) {
  memset(g, 0, sizeof(GateKeeperLibraryRecord));
}

void clearGateKeeperFragmentRecord(GateKeeperFragmentRecord *g) {
  memset(g, 0, sizeof(GateKeeperFragmentRecord));
}


////////////////////////////////////////////////////////////////////////////////


static
int
AS_PER_decodeLibraryFeaturesBoolean(char *feature, char *value) {
  int  ret = 0;

  //  Decodes a string with 0/1, false/true, no/yes into an integer flag.

  switch (value[0]) {
    case '0':
    case 'f':
    case 'F':
    case 'n':
    case 'N':
      ret = 0;
      break;
    case '1':
    case 't':
    case 'T':
    case 'y':
    case 'Y':
      ret = 1;
      break;
    default:
      fprintf(stderr, "AS_PER_decodeLibraryFeatures()-- Found feature '%s' but has unknown boolean value '%s'\n",
              feature, value);
      break;
  }

  return(ret);
}



void
AS_PER_decodeLibraryFeatures(GateKeeperLibraryRecord *gkpl,
                             LibraryMesg             *lmesg) {
  int f;
  for (f=0; f<lmesg->num_features; f++) {
    char *fea = lmesg->features[f];
    char *val = lmesg->values[f];


    //  isNotRandom --
    if        (strcasecmp(fea, "isNotRandom") == 0) {
      gkpl->isNotRandom = AS_PER_decodeLibraryFeaturesBoolean("isNotRandom", val);
    }

    //  doNotOverlapTrim -- 
    else if (strcasecmp(fea, "doNotOverlapTrim") == 0) {
      gkpl->doNotOverlapTrim = AS_PER_decodeLibraryFeaturesBoolean("doNotOverlapTrim", val);
    }

    //  doNotTrustHomopolymerRuns -- 
    else if (strcasecmp(fea, "doNotTrustHomopolymerRuns") == 0) {
      gkpl->doNotTrustHomopolymerRuns = AS_PER_decodeLibraryFeaturesBoolean("doNotTrustHomopolymerRuns", val);
    }

    //  hpsIsPeakSpacing -- 
    else if (strcasecmp(fea, "hpsIsPeakSpacing") == 0) {
      gkpl->hpsIsPeakSpacing = AS_PER_decodeLibraryFeaturesBoolean("hpsIsPeakSpacing", val);
    }

    //  hpsIsFlowGram -- 
    else if (strcasecmp(fea, "hpsIsFlowGram") == 0) {
      gkpl->hpsIsFlowGram = AS_PER_decodeLibraryFeaturesBoolean("hpsIsFlowGram", val);
    }

    else {
      fprintf(stderr, "AS_PER_decodeLibraryFeatures()-- Found feature '%s' but don't understand it.\n",
              fea);
    }
  }
}


void
AS_PER_encodeLibraryFeaturesCleanup(LibraryMesg *lmesg) {
  while (lmesg->num_features > 0) {
    lmesg->num_features--;
    safe_free(lmesg->features[lmesg->num_features]);
    safe_free(lmesg->values  [lmesg->num_features]);
  }
  safe_free(lmesg->features);
  safe_free(lmesg->values);
}


void
AS_PER_encodeLibraryFeatures(GateKeeperLibraryRecord *gkpl,
                             LibraryMesg             *lmesg) {

  //  Examine the gkpl, allocate space to encode the features into
  //  features/values, return the number of features encoded.
  //
  //  Be sure to call AS_PER_encodeLibraryFeaturesCleanup to properly
  //  cleanup the LibraryMesg after it is written!

  //  We can hardcode the maximum number of features we expect to be
  //  writing.  Otherwise, we should count the number of features we
  //  want to encode, allocate....but what a pain.
  //
  lmesg->num_features = 0;
  lmesg->features     = (char **)safe_malloc(5 * sizeof(char*));
  lmesg->values       = (char **)safe_malloc(5 * sizeof(char*));

  int    nf  = 0;
  char **fea = lmesg->features;
  char **val = lmesg->values;

  //  Mostly for debugging, but just might be generally a
  //  GoodThing(tm) to always specify optional features.
  int    alwaysEncode = 1;

  if (gkpl->isNotRandom || alwaysEncode) {
    fea[nf] = (char *)safe_malloc(32 * sizeof(char));
    val[nf] = (char *)safe_malloc(32 * sizeof(char));
    sprintf(fea[nf], "isNotRandom");
    sprintf(val[nf], "%d", gkpl->isNotRandom);
    nf++;
  }

  if (gkpl->doNotOverlapTrim || alwaysEncode) {
    fea[nf] = (char *)safe_malloc(32 * sizeof(char));
    val[nf] = (char *)safe_malloc(32 * sizeof(char));
    sprintf(fea[nf], "doNotOverlapTrim");
    sprintf(val[nf], "%d", gkpl->doNotOverlapTrim);
    nf++;
  }

  if (gkpl->doNotTrustHomopolymerRuns || alwaysEncode) {
    fea[nf] = (char *)safe_malloc(32 * sizeof(char));
    val[nf] = (char *)safe_malloc(32 * sizeof(char));
    sprintf(fea[nf], "doNotTrustHomopolymerRuns");
    sprintf(val[nf], "%d", gkpl->doNotTrustHomopolymerRuns);
    nf++;
  }

  if (gkpl->hpsIsPeakSpacing || alwaysEncode) {
    fea[nf] = (char *)safe_malloc(32 * sizeof(char));
    val[nf] = (char *)safe_malloc(32 * sizeof(char));
    sprintf(fea[nf], "hpsIsPeakSpacing");
    sprintf(val[nf], "%d", gkpl->hpsIsPeakSpacing);
    nf++;
  }

  if (gkpl->hpsIsFlowGram || alwaysEncode) {
    fea[nf] = (char *)safe_malloc(32 * sizeof(char));
    val[nf] = (char *)safe_malloc(32 * sizeof(char));
    sprintf(fea[nf], "hpsIsFlowGram");
    sprintf(val[nf], "%d", gkpl->hpsIsFlowGram);
    nf++;
  }

  lmesg->num_features = nf;
}


////////////////////////////////////////////////////////////////////////////////


fragRecord *new_fragRecord(void) {
  fragRecord *fr = (fragRecord *)safe_malloc(sizeof(fragRecord));
  clr_fragRecord(fr);
  return(fr);
}

void        del_fragRecord(fragRecord *fr) {
  safe_free(fr);
}

void        clr_fragRecord(fragRecord *fr) {
  clearGateKeeperFragmentRecord(&fr->gkfr);
  fr->hasSEQ = 0;
  fr->hasQLT = 0;
  fr->hasHPS = 0;
  fr->hasSRC = 0;
  fr->seq[0] = 0;
  fr->qlt[0] = 0;
  fr->hps[0] = 0;
  fr->src[0] = 0;
}

//  Set clear region 'which' and all later clear regions.  You are
//  explicitly not allowed to set the original clear range.
//
void        setFragRecordClearRegion(fragRecord *fr,
                                     uint32 start,
                                     uint32 end,
                                     uint32 which) {
  assert(which >  AS_READ_CLEAR_VEC);
  assert(which <= AS_READ_CLEAR_LATEST);
  for (; which <= AS_READ_CLEAR_LATEST; which++) {
    fr->gkfr.clearBeg[which] = start;
    fr->gkfr.clearEnd[which] = end;
  }
}


void        getFragRecordClearRegion(fragRecord *fr, uint32 *start, uint32 *end, uint32 which) {
  if (which == AS_READ_CLEAR_UNTRIM) {
    *start = 0;
    *end   = fr->gkfr.seqLen;
  } else {
    assert(which <= AS_READ_CLEAR_LATEST);
    *start = fr->gkfr.clearBeg[which];
    *end   = fr->gkfr.clearEnd[which];
  }
}


uint32      getFragRecordClearRegionBegin(fragRecord *fr, uint32 which) {
  if (which == AS_READ_CLEAR_UNTRIM) {
    return(0);
  } else {
    assert(which <= AS_READ_CLEAR_LATEST);
    return(fr->gkfr.clearBeg[which]);
  }
}


uint32      getFragRecordClearRegionEnd  (fragRecord *fr, uint32 which) {
  if (which == AS_READ_CLEAR_UNTRIM) {
    return(fr->gkfr.seqLen);
  } else {
    assert(which <= AS_READ_CLEAR_LATEST);
    return(fr->gkfr.clearEnd[which]);
  }
}


////////////////////////////////////////////////////////////////////////////////


static
void
getFragData(GateKeeperStore *gkp, fragRecord *fr, int streamFlags) {
  VLSTRING_SIZE_T   actualLength = 0;
  StoreStruct      *store        = NULL;

  fr->hasSEQ = 0;
  fr->hasQLT = 0;
  fr->hasHPS = 0;
  fr->hasSRC = 0;

  fr->seq[0] = 0;
  fr->qlt[0] = 0;
  fr->hps[0] = 0;
  fr->src[0] = 0;

  if ((streamFlags & FRAG_S_SEQ) &&
      !(streamFlags & FRAG_S_QLT)) {
    fr->hasSEQ = 1;
    if (fr->gkfr.seqLen > 0) {
      assert(gkp->partmap == NULL);
      getVLRecordStore(gkp->seq,
                       fr->gkfr.seqOffset,
                       fr->seq,
                       VLSTRING_MAX_SIZE,
                       &actualLength);
      fr->seq[actualLength] = 0;
    }
  }
  if (streamFlags & FRAG_S_QLT) {
    assert(fr->hasSEQ == 0);
    fr->hasSEQ = 1;
    fr->hasQLT = 1;
    if (fr->gkfr.seqLen > 0) {
      store = gkp->qlt;
      if (gkp->partmap)
        store = gkp->partqlt;
      getVLRecordStore(store,
                       fr->gkfr.qltOffset,
                       fr->qlt,
                       VLSTRING_MAX_SIZE,
                       &actualLength);
      fr->qlt[actualLength] = 0;
      decodeSequenceQuality(fr->qlt, fr->seq, fr->qlt);
    }
  }
  if (streamFlags & FRAG_S_HPS) {
    fr->hasHPS = 1;
    if (fr->gkfr.hpsLen > 0) {
      store = gkp->hps;
      if (gkp->partmap)
        store = gkp->parthps;
      getVLRecordStore(store,
                       fr->gkfr.hpsOffset,
                       fr->hps,
                       VLSTRING_MAX_SIZE,
                       &actualLength);
      fr->hps[actualLength] = 0;
    }
  }
  if (streamFlags & FRAG_S_SRC) {
    fr->hasSRC = 1;
    if (fr->gkfr.srcLen > 0) {
      store = gkp->src;
      if (gkp->partmap)
        store = gkp->partsrc;
      getVLRecordStore(store,
                       fr->gkfr.srcOffset,
                       fr->src,
                       VLSTRING_MAX_SIZE,
                       &actualLength);
      fr->src[actualLength] = 0;
    }
  }
}


void    getFrag(GateKeeperStore *gkp, CDS_IID_t iid, fragRecord *fr, int32 flags) {
  if (gkp->partmap == NULL) {
    getIndexStore(gkp->frg, iid, &fr->gkfr);
  } else {
    GateKeeperFragmentRecord *gkfr;

    gkfr = (GateKeeperFragmentRecord *)LookupValueInHashTable_AS(gkp->partmap, iid, 0);
    if (gkfr == NULL) {
      fprintf(stderr, "getFrag()-- ERROR!  IID "F_IID" not in partition!\n", iid);
      assert(0);
    }

    memcpy(&fr->gkfr, gkfr, sizeof(GateKeeperFragmentRecord));
  }

  getFragData(gkp, fr, flags);
}


void    setFrag(GateKeeperStore *gkp, CDS_IID_t iid, fragRecord *fr) {
  assert(gkp->partmap == NULL);
  setIndexStore(gkp->frg, iid, &fr->gkfr);
}

void    delFrag(GateKeeperStore *gkp, CDS_IID_t iid) {
  GateKeeperFragmentRecord   gkfr;
  CDS_IID_t                  miid;

  assert(gkp->partmap == NULL);

  //  Delete fragment with iid from the store.  If the fragment has a
  //  mate, remove the mate relationship from both fragmentss.

  getIndexStore(gkp->frg, iid, &gkfr);
  miid = gkfr.mateIID;
  gkfr.deleted = 1;
  gkfr.mateIID = 0;
  setIndexStore(gkp->frg, iid, &gkfr);

  if (miid > 0) {
    getIndexStore(gkp->frg, miid, &gkfr);
    gkfr.mateIID = 0;
    setIndexStore(gkp->frg, miid, &gkfr);
  }
}



////////////////////////////////////////////////////////////////////////////////


FragStream      *openFragStream(GateKeeperStore *gkp, int flags) {
  FragStream  *fs = (FragStream *)safe_malloc(sizeof(FragStream));
  fs->gkp   = gkp;
  fs->frg   = NULL;
  fs->seq   = NULL;
  fs->qlt   = NULL;
  fs->hps   = NULL;
  fs->src   = NULL;
  fs->flags = flags;

#if 0
  int bufferSize = 1048576;
  fs->frgBuffer = (char *)safe_malloc(sizeof(char) * bufferSize);
  fs->seqBuffer = (char *)safe_malloc(sizeof(char) * bufferSize);
  fs->qltBuffer = (char *)safe_malloc(sizeof(char) * bufferSize);
  fs->hpsBuffer = (char *)safe_malloc(sizeof(char) * bufferSize);
  fs->seqBuffer = (char *)safe_malloc(sizeof(char) * bufferSize);
#else
  int bufferSize = 0;
  fs->frgBuffer = NULL;
  fs->seqBuffer = NULL;
  fs->qltBuffer = NULL;
  fs->hpsBuffer = NULL;
  fs->seqBuffer = NULL;
#endif

  fs->frg   = openStream(fs->gkp->frg, fs->frgBuffer, bufferSize);

  if ((fs->flags & FRAG_S_SEQ) && !(fs->flags & FRAG_S_QLT))
    fs->seq   = openStream(fs->gkp->seq, fs->seqBuffer, bufferSize);
  if (fs->flags & FRAG_S_QLT)
    fs->qlt   = openStream(fs->gkp->qlt, fs->qltBuffer, bufferSize);
  if (fs->flags & FRAG_S_HPS)
    fs->hps   = openStream(fs->gkp->hps, fs->hpsBuffer, bufferSize);
  if (fs->flags & FRAG_S_SRC)
    fs->src   = openStream(fs->gkp->src, fs->srcBuffer, bufferSize);

  return(fs);
}


void             resetFragStream(FragStream *fs, int64 startIndex, int64 endIndex) {
  int64  seqOffset, qltOffset, hpsOffset, srcOffset;

  if (startIndex == STREAM_FROMSTART) {
    seqOffset = STREAM_FROMSTART;
    qltOffset = STREAM_FROMSTART;
    hpsOffset = STREAM_FROMSTART;
    srcOffset = STREAM_FROMSTART;
  } else {

    //  Read the frag the hard way so we can find the offsets the other
    //  stores need to use.

    GateKeeperFragmentRecord  gkpf;

    getIndexStore(fs->gkp->frg, startIndex, &gkpf);

    seqOffset = gkpf.seqOffset;
    qltOffset = gkpf.qltOffset;
    hpsOffset = gkpf.hpsOffset;
    srcOffset = gkpf.srcOffset;
  }

  resetStream(fs->frg, startIndex, endIndex);

  if ((fs->flags & FRAG_S_SEQ) && !(fs->flags & FRAG_S_QLT))
    resetStream(fs->seq, seqOffset, STREAM_UNTILEND);
  if (fs->flags & FRAG_S_QLT)
    resetStream(fs->qlt, qltOffset, STREAM_UNTILEND);
  if (fs->flags & FRAG_S_HPS)
    resetStream(fs->hps, hpsOffset, STREAM_UNTILEND);
  if (fs->flags & FRAG_S_SRC)
    resetStream(fs->src, srcOffset, STREAM_UNTILEND);
}


void             closeFragStream(FragStream *fs) {
  closeStream(fs->frg);
  if ((fs->flags & FRAG_S_SEQ) && !(fs->flags & FRAG_S_QLT))
    closeStream(fs->seq);
  if (fs->flags & FRAG_S_QLT)
    closeStream(fs->qlt);
  if (fs->flags & FRAG_S_HPS)
    closeStream(fs->hps);
  if (fs->flags & FRAG_S_SRC)
    closeStream(fs->src);

  safe_free(fs->frgBuffer);
  safe_free(fs->seqBuffer);
  safe_free(fs->qltBuffer);
  safe_free(fs->hpsBuffer);
  safe_free(fs->seqBuffer);
}



int              nextFragStream(FragStream *fs, fragRecord *fr) {
  VLSTRING_SIZE_T  actualLength = 0;

  if (nextStream(fs->frg, &fr->gkfr) == 0)
    return(0);

  //  So we can use the stream, we can't use getFragData() here.  We
  //  need to duplicate it.

  fr->seq[0] = 0;
  fr->qlt[0] = 0;
  fr->hps[0] = 0;
  fr->src[0] = 0;

  fr->hasSEQ = 0;
  fr->hasQLT = 0;
  fr->hasHPS = 0;
  fr->hasSRC = 0;

  if ((fs->flags & FRAG_S_SEQ) && !(fs->flags & FRAG_S_QLT)) {
    fr->hasSEQ = 1;
    nextVLRecordStream(fs->seq, fr->seq, MAX_SEQ_LENGTH, &actualLength);
    fr->seq[actualLength] = 0;
  }

  if (fs->flags & FRAG_S_QLT) {
    assert(fr->hasSEQ == 0);
    fr->hasSEQ = 1;
    fr->hasQLT = 1;
    nextVLRecordStream(fs->qlt, fr->qlt, MAX_SEQ_LENGTH, &actualLength);
    fr->qlt[actualLength] = 0;
    decodeSequenceQuality(fr->qlt, fr->seq, fr->qlt);
  }

  if (fs->flags & FRAG_S_HPS) {
    fr->hasHPS = 1;
    nextVLRecordStream(fs->hps, fr->hps, MAX_HPS_LENGTH, &actualLength);
    fr->hps[actualLength] = 0;
  }

  if (fs->flags & FRAG_S_SRC) {
    fr->hasSRC = 1;
    nextVLRecordStream(fs->src, fr->src, MAX_SRC_LENGTH, &actualLength);
    fr->src[actualLength] = 0;
  }

  return(1);
}



void
loadGateKeeperStorePartial(GateKeeperStore *gkp,
                           int64            firstElem,
                           int64            lastElem,
                           int              flags) {
  GateKeeperFragmentRecord   gkf;

  int64  seqFirst = 0, seqLast = 0;
  int64  qltFirst = 0, qltLast = 0;
  int64  hpsFirst = 0, hpsLast = 0;
  int64  srcFirst = 0, srcLast = 0;

  if (firstElem != 0) {
    getIndexStore(gkp->frg, firstElem, &gkf);

    seqFirst = gkf.seqOffset;
    qltFirst = gkf.qltOffset;
    hpsFirst = gkf.hpsOffset;
    srcFirst = gkf.srcOffset;
  }

  if ((lastElem != 0) &&
      (lastElem != getLastElemFragStore(gkp))) {
    getIndexStore(gkp->frg, lastElem + 1, &gkf);

    seqLast = gkf.seqOffset;
    qltLast = gkf.qltOffset;
    hpsLast = gkf.hpsOffset;
    srcLast = gkf.srcOffset;
  }

  if ((flags & FRAG_S_SEQ) && !(flags & FRAG_S_QLT))
    gkp->seq = convertStoreToPartialMemoryStore(gkp->seq, seqFirst, seqLast);

  if (flags & FRAG_S_QLT)
    gkp->qlt = convertStoreToPartialMemoryStore(gkp->qlt, qltFirst, qltLast);

  if (flags & FRAG_S_HPS)
    gkp->hps = convertStoreToPartialMemoryStore(gkp->hps, hpsFirst, hpsLast);

  if (flags & FRAG_S_SRC)
    gkp->src = convertStoreToPartialMemoryStore(gkp->src, srcFirst, srcLast);

  //  Making the frag info a memory store is easy!  We do this last,
  //  since we need to get a frag that isn't included here -- being
  //  lastElem+1.  The cost of doing this last is two fragment loads
  //  that otherwise we could do from the memory store.
  //
  gkp->frg = convertStoreToPartialMemoryStore(gkp->frg, firstElem, lastElem);
}
