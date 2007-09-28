
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

static char const *rcsid = "$Id: AS_GKP_sff.c,v 1.1 2007-09-28 07:31:22 brianwalenz Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <ctype.h>

#include "AS_global.h"
#include "AS_UTL_fileIO.h"
#include "AS_GKP_include.h"
#include "AS_PER_gkpStore.h"
#include "AS_PER_encodeSequenceQuality.h"

#define SFF_KEY_SEQUENCE_MAX         64

#define SFF_NAME_LENGTH_MAX         256
#define SFF_NUMBER_OF_FLOWS_MAX     512
#define SFF_NUMBER_OF_BASES_MAX    2048  //  The assembler itself cannot handle longer


typedef struct {
  //  The next block is read in one swoop from the sff file.  DO NOT MODIFY!
  uint32   magic_number;
  char     version[4];
  uint64   index_offset;
  uint32   index_length;
  uint32   number_of_reads;
  uint16   header_length;
  uint16   key_length;
  uint16   number_of_flows_per_read;
  uint8    flowgram_format_code;

  char     flow_chars[SFF_NUMBER_OF_FLOWS_MAX];     //  h->number_of_flows_per_read
  char     key_sequence[SFF_KEY_SEQUENCE_MAX];      //  h->key_length

  uint32   swap_endianess;
} sffHeader;

typedef struct {
  //  The next block is read in one swoop from the sff file.  DO NOT MODIFY!
  uint16   read_header_length;
  uint16   name_length;
  uint32   number_of_bases;
  uint16   clip_quality_left;
  uint16   clip_quality_right;
  uint16   clip_adapter_left;
  uint16   clip_adapter_right;

  char     name[SFF_NAME_LENGTH_MAX];                     //  r->name_length

  uint16   flowgram_values[SFF_NUMBER_OF_FLOWS_MAX];      //  h->number_of_flows_per_read
  uint8    flow_index_per_base[SFF_NUMBER_OF_BASES_MAX];  //  r->number_of_bases
  char     bases[SFF_NUMBER_OF_BASES_MAX];                //  r->number_of_bases
  uint8    quality_scores[SFF_NUMBER_OF_BASES_MAX];       //  r->number_of_bases

  char     quality[SFF_NUMBER_OF_BASES_MAX];              //  quality_scores converted to CA-format qv
} sffRead;


inline
uint64
uint64Swap(uint64 x) {
  x = ((x >>  8) & 0x00ff00ff00ff00ffLLU) | ((x <<  8) & 0xff00ff00ff00ff00LLU);
  x = ((x >> 16) & 0x0000ffff0000ffffLLU) | ((x << 16) & 0xffff0000ffff0000LLU);
  x = ((x >> 32) & 0x00000000ffffffffLLU) | ((x << 32) & 0xffffffff00000000LLU);
  return(x);
}

inline
uint32
uint32Swap(uint32 x) {
  x = ((x >>  8) & 0x00ff00ff) | ((x <<  8) & 0xff00ff00);
  x = ((x >> 16) & 0x0000ffff) | ((x << 16) & 0xffff0000);
  return(x);
}

inline
uint16
uint16Swap(uint16 x) {
  x = ((x >>  8) & 0x000000ff) | ((x <<  8) & 0x0000ff00);
  return(x);
}


static
void
readsff_header(FILE *sff, sffHeader *h) {

  memset(h, 0, sizeof(sffHeader));

  AS_UTL_safeRead(sff, h, "readsff_header_1", 31, 1);

  if (h->magic_number != 0x2e736666) {
    h->swap_endianess           = 1;
    h->magic_number             = uint32Swap(h->magic_number);
    h->index_offset             = uint64Swap(h->index_offset);
    h->index_length             = uint32Swap(h->index_length);
    h->number_of_reads          = uint32Swap(h->number_of_reads);
    h->header_length            = uint16Swap(h->header_length);
    h->key_length               = uint16Swap(h->key_length);
    h->number_of_flows_per_read = uint16Swap(h->number_of_flows_per_read);
  }

  assert(h->magic_number == 0x2e736666);
  assert(h->number_of_flows_per_read < SFF_NUMBER_OF_FLOWS_MAX);
  assert(h->key_length < SFF_KEY_SEQUENCE_MAX);
  
  AS_UTL_safeRead(sff,  h->flow_chars,   "readsff_header_2", sizeof(char), h->number_of_flows_per_read);
  AS_UTL_safeRead(sff,  h->key_sequence, "readsff_header_3", sizeof(char), h->key_length);

  //h->flow_chars  [h->number_of_flows_per_read] = 0;
  //h->key_sequence[h->key_length] = 0;

  uint64  padding_length = h->header_length - 31 - h->number_of_flows_per_read - h->key_length;
  if (padding_length > 0) {
    uint64  junk;
    AS_UTL_safeRead(sff, &junk, "readsff_header_4", sizeof(char), padding_length);
  }

  //  The spec says the index might be here, however, all files I've
  //  seen have the index at the end of the file.  Don't get all worked
  //  up about how wasteful this block seems.
  //
  //  We read just because if we are a popen()ed file, we cannot seek.
  //
  //  AS_UTL_fseek(sff, h->index_length, SEEK_CUR);
  //
  if ((h->index_length > 0) && (h->index_offset == h->header_length)) {
    char *junk = (char *)safe_malloc(sizeof(char) * h->index_length);
    AS_UTL_safeRead(sff, junk, "readsff_index", sizeof(char), h->index_length);
    safe_free(junk);
  }

#if 0
  fprintf(stderr, "header: magic_number %8u\n", h->magic_number);
  fprintf(stderr, "header: version %d%d%d%d\n", h->version[0], h->version[1], h->version[2], h->version[3]);
  fprintf(stderr, "header: index_offset %lu  index_length %u\n", h->index_offset, h->index_length);
  fprintf(stderr, "header: number_of_reads %u\n", h->number_of_reads);
  fprintf(stderr, "header: header_length %u\n", h->header_length);
  fprintf(stderr, "header: key_length %u\n", h->key_length);
  fprintf(stderr, "header: number_of_flows_per_read %u\n", h->number_of_flows_per_read);
  fprintf(stderr, "header: flowgram_format_code %u\n", h->flowgram_format_code);
#endif
}


static
void
readsff_read(FILE *sff, sffHeader *h, sffRead *r) {

  memset(r, 0, sizeof(sffRead));

  AS_UTL_safeRead(sff, r, "readsff_read_1", 16, 1);

  if (h->swap_endianess) {
    r->read_header_length = uint16Swap(r->read_header_length);
    r->name_length        = uint16Swap(r->name_length);
    r->number_of_bases    = uint32Swap(r->number_of_bases);
    r->clip_quality_left  = uint16Swap(r->clip_quality_left);
    r->clip_quality_right = uint16Swap(r->clip_quality_right);
    r->clip_adapter_left  = uint16Swap(r->clip_adapter_left);
    r->clip_adapter_right = uint16Swap(r->clip_adapter_right);
  }

  assert(r->read_header_length < SFF_NAME_LENGTH_MAX);
  assert(r->number_of_bases < SFF_NUMBER_OF_BASES_MAX);

  AS_UTL_safeRead(sff, r->name, "readsff_read_2", sizeof(char), r->name_length);
  r->name[r->name_length] = 0;

  uint64  padding_length = r->read_header_length - 16 - r->name_length;
  if (padding_length > 0) {
    uint64  junk;
    AS_UTL_safeRead(sff, &junk, "readsff_read_3", sizeof(char), padding_length);
  }

  AS_UTL_safeRead(sff, r->flowgram_values,     "readsff_read_4", sizeof(uint16), h->number_of_flows_per_read);
  AS_UTL_safeRead(sff, r->flow_index_per_base, "readsff_read_5", sizeof(uint8),  r->number_of_bases);
  AS_UTL_safeRead(sff, r->bases,               "readsff_read_6", sizeof(char),   r->number_of_bases);
  AS_UTL_safeRead(sff, r->quality_scores,      "readsff_read_7", sizeof(uint8),  r->number_of_bases);

  int i;
  for (i=0; i<r->number_of_bases; i++)
    r->quality[i] = r->quality_scores[i] + '0';

  r->bases[r->number_of_bases] = 0;
  r->quality[r->number_of_bases] = 0;

  //  The padding_length is the number of bytes to make the above four
  //  chunks of data be of size that is divisible by 8.  The
  //  padding_length we compute directly below is the number of bytes
  //  we read past the last multiple of 8, and if that is non-zero, we
  //  need to read 8-padding_length bytes.
  //
  padding_length = (h->number_of_flows_per_read * sizeof(uint16) +
                    r->number_of_bases * sizeof(uint8) +
                    r->number_of_bases * sizeof(char) +
                    r->number_of_bases * sizeof(uint8)) % 8;
  if (padding_length > 0) {
    uint64  junk;
    AS_UTL_safeRead(sff, &junk, "readsff_read_8", sizeof(char), 8 - padding_length);
  }
}


//  Massage the sff into a new gatekeeper entry.
//
static
CDS_UID_t
readsff_constructUIDFromName(char *name, int constructReadUID) {
  CDS_UID_t   uid;
  char        base36[16];
  uint64      timestamp;
  uint64      rigname;
  uint64      region;
  uint64      position, x, y;

  base36[0] = name[0];
  base36[1] = name[1];
  base36[2] = name[2];
  base36[3] = name[3];
  base36[4] = name[4];
  base36[5] = name[5];
  base36[6] = 0;
  timestamp = strtoull(base36, NULL, 36);

  base36[0] = name[6];
  base36[1] = 0;
  rigname = strtoull(base36, NULL, 36) & 0x0000000f;

  base36[0] = name[7];
  base36[1] = name[8];
  base36[2] = 0;
  region = strtoul(base36, NULL, 10);

  base36[0] = name[9];
  base36[1] = name[10];
  base36[2] = name[11];
  base36[3] = name[12];
  base36[4] = name[13];
  base36[5] = 0;
  position = strtoul(base36, NULL, 36);

  x = position / 4096;
  y = position % 4096;

  int err = 0;

  if (timestamp > 1 << 31)
    err |= 1;
  if ((region == 0) || (region > 4))
    err |= 2;
  if (x > 16384)
    err |= 4;
  if (y > 4096)
    err |= 8;

  if (err)
    fprintf(stdout, "%s -- err %d -- timestamp:0x%08lx region:0x%08lx position:0x%08lx (x=%lu y=%lu)\n", name, err, timestamp, region, position, x, y);

  assert(4 + 31 + 3 + 14 + 12 == 64);

  if (constructReadUID) {
    uid   = 0;
    uid  |= rigname;    //  4 bits
    uid <<= 31;
    uid  |= timestamp;  //  31 bits
    uid <<= 3;
    uid  |= region;     //  3 bits
    uid <<= 14;
    uid  |= x;          //  14 bits
    uid <<= 12;
    uid  |= y;          //  12 bits
  } else {
    uid   = 0;
    uid  |= rigname;    //  4 bits
    uid <<= 31;
    uid  |= timestamp;  //  31 bits
    uid <<= 3;
    uid  |= region;     //  3 bits
    uid <<= 26;
    uid  |= 0x03ffffff; // 26 bits
  }

  return(uid);
}


static
CDS_IID_t
readsff_constructLibraryIIDFromName(char *name) {
  CDS_UID_t  uid = readsff_constructUIDFromName(name, 0);
  CDS_IID_t  iid = getGatekeeperUIDtoIID(gkpStore, uid, NULL);

  if (iid == 0) {
    GateKeeperLibraryRecord  gkl  = {0};

    gkl.libraryUID = uid;
    gkl.comment[0] = 0;

    gkl.spare2 = 0;
    gkl.spare1 = 0;

    gkl.hpsIsFlowGram    = 1;
    gkl.hpsIsPeakSpacing = 0;

    gkl.doNotTrustHomopolymerRuns  = 1;
    gkl.doNotOverlapTrim           = 0;
    gkl.isNotRandom                = 0;

    gkl.orientation = AS_READ_ORIENT_UNKNOWN;
    gkl.ZZZdeleted  = 0;

    gkl.mean   = 0.0;
    gkl.stddev = 0.0;

    appendIndexStore(gkpStore->lib, &gkl);
    setGatekeeperUIDtoIID(gkpStore, gkl.libraryUID, getLastElemStore(gkpStore->lib), AS_IID_LIB);

    iid = getLastElemStore(gkpStore->lib);

#if 0
    fprintf(stderr, "added library "F_UID" at iid "F_IID"\n", gkl.libraryUID, iid);
#endif
  }

  return(iid);
}


int
Load_SFF(FILE *sff) {

  sffHeader *h  = (sffHeader *)safe_calloc(sizeof(sffHeader), 1);
  sffRead   *r  = (sffRead   *)safe_calloc(sizeof(sffRead),   1);
  int        rn = 0;

  char       encodedsequence[AS_FRAG_MAX_LEN+1] = {0};

  readsff_header(sff, h);

  //  Construct a gkpLibraryRecord for this sff file.  Well, this is
  //  where we'd LIKE to do it, but since the sff doesn't give us any
  //  reasonable way to make a UID from the header, we defer until we
  //  get the first read.  Then, we use the read timestamp, hash and
  //  region to make a library.


  for (rn=0; rn < h->number_of_reads; rn++) {
    GateKeeperFragmentRecord gkf = {0};
    clearGateKeeperFragmentRecord(&gkf);

    readsff_read(sff, h, r);

    gkf.readUID = readsff_constructUIDFromName(r->name, 1);
    gkf.readIID = 0;
    gkf.mateIID = 0;

    if (getGatekeeperUIDtoIID(gkpStore, gkf.readUID, NULL)) {
      fprintf(errorFP, "# SFF Error: Fragment "F_UID" exists, can't add it again.\n",
              gkf.readUID);
      continue;
      //return(GATEKEEPER_FAILURE);
    }

    gkf.libraryIID  = readsff_constructLibraryIIDFromName(r->name);
    gkf.orientation = AS_READ_ORIENT_UNKNOWN;

    //  Set clear ranges
    //
    //  These are base-based.  If either value is 0, that means the
    //  value was not computed.
    //
    //  We have a policy decision here.  If only one of the ranges is
    //  set, we can either ignore both, or set the unset one to the
    //  maximum.

    int  clq = r->clip_quality_left;
    int  crq = r->clip_quality_right;
    int  which;

    if (clq == 0)  clq = 1;
    if (crq == 0)  crq = r->number_of_bases;

    for (which=0; which <= AS_READ_CLEAR_LATEST; which++) {
      gkf.clearBeg[which] = clq - 1;
      gkf.clearEnd[which] = crq - 1;
    }

    gkf.hasQualityClear = 0;
    gkf.clearBeg[AS_READ_CLEAR_QLT] = 0;
    gkf.clearEnd[AS_READ_CLEAR_QLT] = 0;

    if ((r->clip_quality_left > 0) && (r->clip_quality_right > 0)) {
      gkf.hasQualityClear = 1;
      gkf.clearBeg[AS_READ_CLEAR_QLT] = r->clip_quality_left  - 1;
      gkf.clearEnd[AS_READ_CLEAR_QLT] = r->clip_quality_right - 1;
    }

    gkf.hasVectorClear = 0;
    gkf.clearBeg[AS_READ_CLEAR_VEC] = 0;
    gkf.clearEnd[AS_READ_CLEAR_VEC] = 0;

    if ((r->clip_adapter_left > 0) && (r->clip_adapter_right > 0)) {
      gkf.hasVectorClear  = 1;
      gkf.clearBeg[AS_READ_CLEAR_VEC] = r->clip_adapter_left  - 1;
      gkf.clearEnd[AS_READ_CLEAR_VEC] = r->clip_adapter_right - 1;
    }


    //  Now add the fragment to the store
    //
    gkf.readIID = getLastElemStore(gkpStore->frg) + 1;

    gkf.seqLen = strlen(r->bases);
    gkf.hpsLen = 0;
    gkf.srcLen = strlen(r->name);

    {
      StoreStat   stats;

      statsStore(gkpStore->seq, &stats);
      gkf.seqOffset = stats.lastElem;

      statsStore(gkpStore->qlt, &stats);
      gkf.qltOffset = stats.lastElem;

      statsStore(gkpStore->hps, &stats);
      gkf.hpsOffset = stats.lastElem;

      statsStore(gkpStore->src, &stats);
      gkf.srcOffset = stats.lastElem;
    }

    setGatekeeperUIDtoIID(gkpStore, gkf.readUID, gkf.readIID, AS_IID_FRG);
    appendIndexStore(gkpStore->frg, &gkf);

    appendVLRecordStore(gkpStore->seq, r->bases, gkf.seqLen);

    encodeSequenceQuality(encodedsequence, r->bases, r->quality);
    appendVLRecordStore(gkpStore->qlt, encodedsequence, gkf.seqLen);

    appendVLRecordStore(gkpStore->hps, NULL,    0);
    appendVLRecordStore(gkpStore->src, r->name, gkf.srcLen);

#if 0
    fprintf(stderr, "Added '%s' of length %d  clears %d %d %d %d -- %d %d %d %d %d %d %d %d\n",
            r->name,
            r->number_of_bases,
            r->clip_quality_left, r->clip_quality_right,
            r->clip_adapter_left, r->clip_adapter_right,
            gkf.clearBeg[AS_READ_CLEAR_ORIG], gkf.clearEnd[AS_READ_CLEAR_ORIG],
            gkf.clearBeg[AS_READ_CLEAR_QLT],  gkf.clearEnd[AS_READ_CLEAR_QLT],
            gkf.clearBeg[AS_READ_CLEAR_VEC],  gkf.clearEnd[AS_READ_CLEAR_VEC],
            gkf.clearBeg[AS_READ_CLEAR_LATEST], gkf.clearEnd[AS_READ_CLEAR_LATEST]);
#endif
  }

  fprintf(stderr, "Added %d 454 reads.\n", rn);

  return(GATEKEEPER_SUCCESS);
}

