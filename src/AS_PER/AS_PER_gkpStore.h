
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

/* 	$Id: AS_PER_gkpStore.h,v 1.16 2007-02-13 17:59:34 brianwalenz Exp $	 */

#ifndef AS_PER_GKPFRGSTORE_H
#define AS_PER_GKPFRGSTORE_H

#include <sys/types.h>
#include <time.h>

#include "AS_MSG_pmesg.h"
#include "AS_PER_genericStore.h"
#include "AS_UTL_PHash.h"

#define NULL_LINK 0

#define AS_PER_NAME_LEN      256
#define AS_PER_COMMENT_LEN   256

//  The following counts represent the number of records of each type
//  PRIOR to processing this batch.  To get the range of batch i, take
//  the difference between batch i+1 and batch i.
//
typedef struct {
  CDS_UID_t      UID;
  char           name[AS_PER_NAME_LEN];
  char           comment[AS_PER_COMMENT_LEN];
  uint64         created;

  unsigned int   deleted:1;
  unsigned int   spare:31;

  int32          numFragments;
  int32          numLibraries;
  int32          numLibraries_s;
} GateKeeperBatchRecord;

#define AS_GKP_ORIENT_UNKNOWN    0x00
#define AS_GKP_ORIENT_INNIE      0x01
#define AS_GKP_ORIENT_OUTTIE     0x02
#define AS_GKP_ORIENT_NORMAL     0x03
#define AS_GKP_ORIENT_ANTINORMAL 0x04

typedef struct {
  CDS_UID_t      UID;

  char           name[AS_PER_NAME_LEN];
  char           comment[AS_PER_COMMENT_LEN];
  uint64         created;

  unsigned int   deleted:1;
  unsigned int   redefined:1;
  unsigned int   orientation:3;
  unsigned int   spare:28;

  double         mean;
  double         stddev;
  
  unsigned int   numFeatures;

  CDS_IID_t      prevInstanceID;    // Previous definitions are linked by this reference
  CDS_IID_t      prevID;            // If redefined == TRUE, the original ID of this

  uint16         birthBatch;        // This entry is valid
  uint16         deathBatch;        // [birthBatch, deatchBatch)
} GateKeeperLibraryRecord;


#define AS_GKP_STATUS_G   0x00
#define AS_GKP_STATUS_B   0x01
#define AS_GKP_STATUS_U   0x02
#define AS_GKP_STATUS_W   0x03
#define AS_GKP_STATUS_X   0x04
#define AS_GKP_STATUS_V   0x05
#define AS_GKP_STATUS_E   0x06
#define AS_GKP_STATUS_I   0x07
#define AS_GKP_STATUS_R   0x08

typedef struct{
  CDS_UID_t        UID;
 
  CDS_IID_t        readIID;
  CDS_IID_t        mateIID;

  CDS_IID_t        libraryIID;
  CDS_UID_t        plateUID;
  unsigned int     plateLocation;

  unsigned int     deleted:1;
  unsigned int     nonrandom:1;
  unsigned int     status:4;
  unsigned int     hasQLT:1;
  unsigned int     hasHPS:1;
  unsigned int     hasOVLclr:1;
  unsigned int     hasCNSclr:1;
  unsigned int     hasCGWclr:1;

  unsigned int     orientation:3;  //  copied from library

  //  If someone ever adds "type" to a read, search for getReadType in
  //  the source.  You'll need to add type to there too.

  unsigned int     spare:20;

  VLSTRING_SIZE_T  clrSta;
  VLSTRING_SIZE_T  clrEnd;

  VLSTRING_SIZE_T  ovlSta;
  VLSTRING_SIZE_T  ovlEnd;

  VLSTRING_SIZE_T  cnsSta;
  VLSTRING_SIZE_T  cnsEnd;

  VLSTRING_SIZE_T  cgwSta;
  VLSTRING_SIZE_T  cgwEnd;

  uint64           seqOffset;
  uint64           qltOffset;
  uint64           hpsOffset;
  uint64           srcOffset;

  uint16           birthBatch;         /* This entry is valid */
  uint16           deathBatch;         /* [birthBatch, deatchBatch) */
} GateKeeperFragmentRecord;



#define MAX_SEQ_LENGTH (AS_READ_MAX_LEN + 1)
#define MAX_HPS_LENGTH (AS_READ_MAX_LEN + 1)
#define MAX_SRC_LENGTH (512 + sizeof(int32) + sizeof(int64))

#define FRAG_S_INF 0x01
#define FRAG_S_SEQ 0x02
#define FRAG_S_QLT 0x04
#define FRAG_S_HPS 0x08
#define FRAG_S_SRC 0x10
#define FRAG_S_ALL 0x1f

typedef struct {
  GateKeeperFragmentRecord   gkfr;
  uint                       flags;
  char                       seq[MAX_SEQ_LENGTH];
  char                       qlt[MAX_SEQ_LENGTH];
  char                       hps[MAX_HPS_LENGTH];
  char                       src[MAX_SRC_LENGTH];
} ReadStruct;

typedef ReadStruct* ReadStructp;




static char getLinkOrientation(GateKeeperFragmentRecord *gkpf){

  switch (gkpf->orientation){
    case AS_GKP_ORIENT_UNKNOWN:
      return '?';
    case AS_GKP_ORIENT_INNIE:
      return 'I';
    case AS_GKP_ORIENT_OUTTIE:
      return 'O';
    case AS_GKP_ORIENT_NORMAL:
      return 'N';
    case AS_GKP_ORIENT_ANTINORMAL:
      return 'A';
    default:
      return '-';
  }

  return '-';
}

// Stores


#define INDEXSTORE_DEF_EXTEND(type)\
static int deleteAndMark ## type ## Store(type ## Store fs, int index, int batchID){\
  type ## Record dr;\
  getIndexStore(fs,index,&dr); \
  dr.deleted = TRUE;\
  dr.deathBatch = batchID;\
  setIndexStore(fs,index,&dr);\
  return(0);\
}

#define INDEXSTORE_DEF(type)\
typedef StoreHandle type ## Store;\
static int commit ## type ## Store(type ## Store sh){\
  return commitStore(sh);\
}\
static type ## Store reset ## type ## Store(type ## Store sh, int firstID){\
  return resetIndexStore(sh, firstID);\
}\
static int close ## type ## Store(type ## Store sh){\
  return closeStore(sh);\
}\
static int delete ## type ## Store(type ## Store fs, int index){\
  type ## Record dr;\
  getIndexStore(fs,index,&dr); \
  dr.deleted = TRUE;\
  setIndexStore(fs,index,&dr);\
  return(0);\
}\
static int get ## type ## Store(type ## Store fs, int index, type ## Record *dr){\
  return getIndexStore(fs,index,dr); \
}\
static int set ## type ## Store(type ## Store fs, int index, type ## Record *dr){\
  return setIndexStore(fs,index,dr); \
}\
static type ## Store create ## type ## Store(char *StorePath, char *ext, int firstID){\
  type ## Store s = createIndexStore(StorePath,ext, sizeof(type ## Record), 1, firstID);\
  return s;\
}\
static type ## Store open ## type ## Store(char *StorePath, char *rw){\
  return openStore(StorePath, rw);\
}\
static int append ## type ## Store(type ## Store store, type ## Record *element){\
  return appendIndexStore(store,element);\
}\
static int32 getNum ## type ## s(type ## Store store){\
  StoreStat stat;\
  statsStore(store, &stat);\
  return(stat.lastElem);\
}

INDEXSTORE_DEF(GateKeeperBatch);
INDEXSTORE_DEF(GateKeeperFragment);
INDEXSTORE_DEF_EXTEND(GateKeeperFragment);
INDEXSTORE_DEF(GateKeeperLibrary);
INDEXSTORE_DEF_EXTEND(GateKeeperLibrary);

#define NUM_GKP_FILES 10

// 1  is gatekeeper store info
// 2  is batches
// 3  is fragments
// 4  is libraries
// 5  is shadow libraries
// 6  is sequence
// 7  is quality
// 8  is homopolymer spacing and etc
// 9  is source
// 10 is gkp.phash

typedef struct {
  uint64    gkpMagic;
  uint64    gkpVersion;
} GateKeeperStoreInfo;

#define UID_NAMESPACE_AS 'U'

typedef struct {
  char                     storePath[FILENAME_MAX];

  GateKeeperStoreInfo      gkp;

  GateKeeperBatchStore     bat;
  GateKeeperFragmentStore  frg;
  GateKeeperLibraryStore   lib;    
  GateKeeperLibraryStore   lis;

  StoreHandle              seq;
  StoreHandle              qlt;
  StoreHandle              hps;

  StoreHandle              src;

  PHashTable_AS           *phs;
} GateKeeperStore;

////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////

GateKeeperStore *createGateKeeperStore(const char *path);
int              testOpenGateKeeperStore(const char *path, int writable);
GateKeeperStore *openGateKeeperStore(const char *path, int writable);
void             closeGateKeeperStore(GateKeeperStore *gkpStore);


void    clearGateKeeperBatchRecord(GateKeeperBatchRecord *g);
void    clearGateKeeperLibraryRecord(GateKeeperLibraryRecord *g);
void    clearGateKeeperFragmentRecord(GateKeeperFragmentRecord *g);

int     getFrag(GateKeeperStore *gkp, int64 iid, ReadStruct *rs, int32 flags);
int     setFrag(GateKeeperStore *gkp, int64 iid, ReadStruct *rs);
int     delFrag(GateKeeperStore *gkp, int64 iid);

#define getFirstElemFragStore(GKP)  getFirstElemStore((GKP)->frg)
#define getLastElemFragStore(GKP)   getLastElemStore((GKP)->frg)

////////////////////////////////////////////////////////////////////////////////

static
GateKeeperStore *loadFragStore(const char *path) {
  return(openGateKeeperStore(path, FALSE));
}

static
GateKeeperStore *loadFragStorePartial(const char *path, int64 firstElem, int64 lastElem) {
  return(openGateKeeperStore(path, FALSE));
}

////////////////////////////////////////////////////////////////////////////////


typedef struct {
  GateKeeperStore   *gkp;
  StreamHandle       frg;
  StreamHandle       seq;
  StreamHandle       qlt;
  StreamHandle       hps;
  StreamHandle       src;
} FragStream;

FragStream      *openFragStream(GateKeeperStore *gkp);
void             resetFragStream(FragStream *fs, int64 startIndex, int64 endIndex);
void             closeFragStream(FragStream *fs);

int64            getStartIndexFragStream(FragStream *fs);

int              nextFragStream(FragStream *fs, ReadStruct *rs, int streamFlags);

////////////////////////////////////////////////////////////////////////////////

typedef GateKeeperStore  tFragStorePartition;

static
tFragStorePartition *openFragStorePartition(char *fragStorePath, int32 partition, int loadData) {
  return(openGateKeeperStore(fragStorePath, FALSE));
};

static
void                 closeFragStorePartition(tFragStorePartition *partition) {
  closeGateKeeperStore(partition);
};

static
int                  getFragStorePartition(tFragStorePartition *partition,
                                           int32 indx,
                                           int32 getFlags,
                                           ReadStruct *rs) {
  return(getFrag(partition, indx, rs, getFlags));
};

#endif
