
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

static char const *rcsid = "$Id: AS_GKP_sff.c,v 1.9 2008-02-29 12:22:10 brianwalenz Exp $";

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

  char    *flow_chars;        //  h->number_of_flows_per_read
  char    *key_sequence;      //  h->key_length

  void    *data_block;
  uint32   data_block_len;

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

  char    *name;                 //  r->name_length
  uint16  *flowgram_values;      //  h->number_of_flows_per_read
  uint8   *flow_index_per_base;  //  r->number_of_bases
  char    *bases;                //  r->number_of_bases
  uint8   *quality_scores;       //  r->number_of_bases
  char    *quality;              //  quality_scores converted to CA-format qv

  void    *data_block;
  uint32   data_block_len;
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

  uint32 newlen = h->number_of_flows_per_read + h->key_length + 2;
  if (h->data_block_len < newlen) {
    h->data_block_len = newlen;
    h->data_block     = safe_realloc(h->data_block, h->data_block_len);
  }

  memset(h->data_block, 0, h->data_block_len);

  h->flow_chars   = h->data_block;
  h->key_sequence = h->data_block + (h->number_of_flows_per_read + 1) * sizeof(char);

  AS_UTL_safeRead(sff,  h->flow_chars,   "readsff_header_2", sizeof(char), h->number_of_flows_per_read);
  AS_UTL_safeRead(sff,  h->key_sequence, "readsff_header_3", sizeof(char), h->key_length);

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

  //  Can you say UGLY?  Hey, it's a lot better than what I originally came up with.

  uint32 ss[6];
  ss[0] = (r->name_length + 1)          * sizeof(char);
  ss[1] = (h->number_of_flows_per_read) * sizeof(uint16) + ss[0];
  ss[2] = (r->number_of_bases)          * sizeof(uint8)  + ss[1];
  ss[3] = (r->number_of_bases + 1)      * sizeof(char)   + ss[2];
  ss[4] = (r->number_of_bases + 1)      * sizeof(char)   + ss[3];
  ss[5] = (r->number_of_bases + 1)      * sizeof(char)   + ss[4];

  if (r->data_block_len < ss[5]) {
    r->data_block_len = ss[5];
    r->data_block     = safe_realloc(r->data_block, r->data_block_len);
  }

  memset(r->data_block, 0, r->data_block_len);

  r->name                 = r->data_block;
  r->flowgram_values      = r->data_block + ss[0];
  r->flow_index_per_base  = r->data_block + ss[1];
  r->bases                = r->data_block + ss[2];
  r->quality_scores       = r->data_block + ss[3];
  r->quality              = r->data_block + ss[4];

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



static
AS_UID
readsff_constructUIDFromName(char *name, int constructReadUID) {
  char  libname[32] = {0};

  assert(strlen(name) == 14);

  //  The read UID is just the name given by 454.  Thanks, 454!  The
  //  library UID is derived from the read UID by yanking out the
  //  position of the read, and replacing it with something obnoxious.
  //  We could have just truncated it to the last 9 letters, but
  //  that's kind of short.  Not that 15 is better.

  if (constructReadUID) {
    libname[0]  = name[0];
    libname[1]  = name[1];
    libname[2]  = name[2];
    libname[3]  = name[3];
    libname[4]  = name[4];
    libname[5]  = name[5];
    libname[6]  = name[6];
    libname[7]  = name[7];
    libname[8]  = name[8];
    libname[9]  = name[9];
    libname[10] = name[10];
    libname[11] = name[11];
    libname[12] = name[12];
    libname[13] = name[13];
    libname[14] = name[14];
    libname[15] = 0;
  } else {
    libname[0] = 'L';
    libname[1] = 'I';
    libname[2] = 'B';
    libname[3] = 'S';
    libname[4] = 'F';
    libname[5] = 'F';
    libname[6]  = name[0];
    libname[7]  = name[1];
    libname[8]  = name[2];
    libname[9]  = name[3];
    libname[10] = name[4];
    libname[11] = name[5];
    libname[12] = name[6];
    libname[13] = name[7];
    libname[14] = name[8];
    libname[15] = name[9];
    libname[16] = 0;
  }

  return(AS_UID_load(libname));
}


static
AS_IID
readsff_constructLibraryIIDFromName(char *name) {
  AS_UID  uid = readsff_constructUIDFromName(name, 0);
  AS_IID  iid = getGatekeeperUIDtoIID(gkpStore, uid, NULL);

  if (iid == 0) {
    GateKeeperLibraryRecord  gkl;

    clearGateKeeperLibraryRecord(&gkl);

    gkl.libraryUID = uid;

    gkl.hpsIsFlowGram    = 1;

    gkl.deletePerfectPrefixes      = 1;
    gkl.doNotTrustHomopolymerRuns  = 1;

    gkl.orientation = AS_READ_ORIENT_UNKNOWN;

    gkl.mean   = 0.0;
    gkl.stddev = 0.0;

    appendIndexStore(gkpStore->lib, &gkl);
    setGatekeeperUIDtoIID(gkpStore, gkl.libraryUID, getLastElemStore(gkpStore->lib), AS_IID_LIB);

    iid = getLastElemStore(gkpStore->lib);

    gkpStore->gkp.sffLibCreated++;
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

  gkpStore->gkp.sffInput += h->number_of_reads;

  for (rn=0; rn < h->number_of_reads; rn++) {
    GateKeeperFragmentRecord gkf = {0};
    clearGateKeeperFragmentRecord(&gkf);

    readsff_read(sff, h, r);

    gkf.readUID = readsff_constructUIDFromName(r->name, 1);
    gkf.readIID = 0;
    gkf.mateIID = 0;

    if (r->number_of_bases - h->key_length < AS_FRAG_MIN_LEN) {
      //  This isn't _really_ an error, and we'll notice it if we
      //  compare the sffInput to sffLoaded counts.
      //
      //gkpStore->gkp.sffErrors++;
      continue;
    }

    if (getGatekeeperUIDtoIID(gkpStore, gkf.readUID, NULL)) {
      AS_GKP_reportError(AS_GKP_SFF_ALREADY_EXISTS,
                         AS_UID_toString(gkf.readUID));
      gkpStore->gkp.sffErrors++;
      continue;
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
    //  maximum.  We set it to the maximum.

    int  clq = r->clip_quality_left;
    int  crq = r->clip_quality_right;
    int  which;

    assert((r->clip_quality_left == 0) || (h->key_length <= r->clip_quality_left));
    assert((r->clip_adapter_left == 0) || (h->key_length <= r->clip_adapter_left));

    if (clq == 0)  clq = h->key_length + 1;
    if (crq == 0)  crq = r->number_of_bases;

    for (which=0; which <= AS_READ_CLEAR_LATEST; which++) {
      gkf.clearBeg[which] = clq - h->key_length - 1;
      gkf.clearEnd[which] = crq - h->key_length;
    }

    if ((r->clip_quality_left > 0) && (r->clip_quality_right > 0)) {
      gkf.hasQualityClear = 1;
      gkf.clearBeg[AS_READ_CLEAR_QLT] = r->clip_quality_left  - h->key_length - 1;
      gkf.clearEnd[AS_READ_CLEAR_QLT] = r->clip_quality_right - h->key_length;
    } else if (r->clip_quality_left > 0) {
      gkf.hasQualityClear = 1;
      gkf.clearBeg[AS_READ_CLEAR_QLT] = r->clip_quality_left  - h->key_length - 1;
      gkf.clearEnd[AS_READ_CLEAR_QLT] = r->number_of_bases - h->key_length;
    } else if (r->clip_quality_right > 0) {
      gkf.hasQualityClear = 1;
      gkf.clearBeg[AS_READ_CLEAR_QLT] = 0;
      gkf.clearEnd[AS_READ_CLEAR_QLT] = r->clip_quality_right - h->key_length;
    } else {
      gkf.hasQualityClear = 0;
      gkf.clearBeg[AS_READ_CLEAR_QLT] = 0;
      gkf.clearEnd[AS_READ_CLEAR_QLT] = 0;
    }

    if ((r->clip_adapter_left > 0) && (r->clip_adapter_right > 0)) {
      gkf.hasVectorClear  = 1;
      gkf.clearBeg[AS_READ_CLEAR_VEC] = r->clip_adapter_left  - h->key_length - 1;
      gkf.clearEnd[AS_READ_CLEAR_VEC] = r->clip_adapter_right - h->key_length;
    } else if (r->clip_adapter_left > 0) {
      gkf.hasVectorClear  = 1;
      gkf.clearBeg[AS_READ_CLEAR_VEC] = r->clip_adapter_left  - h->key_length - 1;
      gkf.clearEnd[AS_READ_CLEAR_VEC] = r->number_of_bases - h->key_length;
    } else if (r->clip_adapter_right > 0) {
      gkf.hasVectorClear  = 1;
      gkf.clearBeg[AS_READ_CLEAR_VEC] = 0;
      gkf.clearEnd[AS_READ_CLEAR_VEC] = r->clip_adapter_right - h->key_length;
    } else {
      gkf.hasVectorClear = 0;
      gkf.clearBeg[AS_READ_CLEAR_VEC] = 0;
      gkf.clearEnd[AS_READ_CLEAR_VEC] = 0;
    }

    if (r->number_of_bases > AS_READ_MAX_LEN) {
      AS_GKP_reportError(AS_GKP_SFF_TOO_LONG, r->name, r->number_of_bases, AS_READ_MAX_LEN);
      gkpStore->gkp.sffWarnings++;
      
      r->number_of_bases = AS_READ_MAX_LEN;

      r->bases  [AS_READ_MAX_LEN + h->key_length] = 0;
      r->quality[AS_READ_MAX_LEN + h->key_length] = 0;

      for (which=0; which <= AS_READ_CLEAR_LATEST; which++) {
        if (gkf.clearBeg[which] > AS_READ_MAX_LEN)
          gkf.clearBeg[which] = AS_READ_MAX_LEN;
        if (gkf.clearEnd[which] > AS_READ_MAX_LEN)
          gkf.clearEnd[which] = AS_READ_MAX_LEN;
      }
    }
     

    //  Now add the fragment to the store
    //
    gkf.readIID = getLastElemStore(gkpStore->frg) + 1;

    gkf.seqLen = strlen(r->bases + h->key_length);
    gkf.hpsLen = 0;
    gkf.srcLen = 0;

    gkf.seqOffset = getLastElemStore(gkpStore->seq) + 1;
    gkf.qltOffset = getLastElemStore(gkpStore->qlt) + 1;
    gkf.hpsOffset = getLastElemStore(gkpStore->hps) + 1;
    gkf.srcOffset = getLastElemStore(gkpStore->src) + 1;

    setGatekeeperUIDtoIID(gkpStore, gkf.readUID, gkf.readIID, AS_IID_FRG);
    appendIndexStore(gkpStore->frg, &gkf);

    appendStringStore(gkpStore->seq, r->bases + h->key_length, gkf.seqLen);

    encodeSequenceQuality(encodedsequence,
                          r->bases + h->key_length,
                          r->quality + h->key_length);
    appendStringStore(gkpStore->qlt, encodedsequence, gkf.seqLen);

    appendStringStore(gkpStore->hps, NULL, 0);
    appendStringStore(gkpStore->src, NULL, 0);

    gkpStore->gkp.sffLoaded++;
  }

  fprintf(stderr, "Added %d 454 reads.\n", rn);

  safe_free(h->data_block);
  safe_free(h);
  safe_free(r->data_block);
  safe_free(r);

  return(0);
}

